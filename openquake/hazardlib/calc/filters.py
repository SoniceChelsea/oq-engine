# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2012-2017 GEM Foundation
#
# OpenQuake is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenQuake is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake. If not, see <http://www.gnu.org/licenses/>.
"""
Module :mod:`~openquake.hazardlib.calc.filters` contain filter functions for
calculators.

Filters are functions (or other callable objects) that should take generators
and return generators. There are two different kinds of filter functions:

1. Source-site filters. Those functions take a generator of two-item tuples,
   each pair consists of seismic source object (that is, an instance of
   a subclass of :class:`~openquake.hazardlib.source.base.BaseSeismicSource`)
   and a site collection (instance of
   :class:`~openquake.hazardlib.site.SiteCollection`).
2. Rupture-site filters. Those also take a generator of pairs, but in this
   case the first item in the pair is a rupture object (instance of
   :class:`~openquake.hazardlib.source.rupture.Rupture`). The second element in
   generator items is still site collection.

The purpose of both kinds of filters is to limit the amount of calculation
to be done based on some criteria, like the distance between the source
and the site. So common design feature of all the filters is the loop over
pairs of the provided generator, filtering the sites collection, and if
there are no items left in it, skipping the pair and continuing to the next
one. If some sites need to be considered together with that source / rupture,
the pair gets generated out, with a (possibly) :meth:`limited
<openquake.hazardlib.site.SiteCollection.filter>` site collection.

Consistency of filters' input and output stream format allows several filters
(obviously, of the same kind) to be chained together.

Filter functions should not make assumptions about the ordering of items
in the original generator or draw more than one pair at once. Ideally, they
should also perform reasonably fast (filtering stage that takes longer than
the actual calculation on unfiltered collection only decreases performance).

Module :mod:`openquake.hazardlib.calc.filters` exports one distance-based
filter function (see :func:`filter_sites_by_distance_to_rupture`) as well as
a "no operation" filter (`source_site_noop_filter`). There is
a class `SourceFilter` to determine the sites
affected by a given source: the second one uses an R-tree index and it is
faster if there are a lot of sources, i.e. if the initial time to prepare
the index can be compensed. Finally, there is a function
`filter_sites_by_distance_to_rupture` based on the Joyner-Boore distance.
"""
import sys
import logging
import collections
from contextlib import contextmanager
import numpy
from scipy.interpolate import interp1d
try:
    import rtree
except ImportError:
    rtree = None
from openquake.baselib.python3compat import raise_
from openquake.hazardlib.site import FilteredSiteCollection
from openquake.hazardlib.geo.utils import fix_lons_idl


@contextmanager
def context(src):
    """
    Used to add the source_id to the error message. To be used as

    with context(src):
        operation_with(src)

    Typically the operation is filtering a source, that can fail for
    tricky geometries.
    """
    try:
        yield
    except:
        etype, err, tb = sys.exc_info()
        msg = 'An error occurred with source id=%s. Error: %s'
        msg %= (src.source_id, err)
        raise_(etype, msg, tb)


def filter_sites_by_distance_to_rupture(rupture, integration_distance, sites):
    """
    Filter out sites from the collection that are further from the rupture
    than some arbitrary threshold.

    :param rupture:
        Instance of :class:`~openquake.hazardlib.source.rupture.Rupture`
        that was generated by :meth:
        `openquake.hazardlib.source.base.BaseSeismicSource.iter_ruptures`
        of an instance of this class.
    :param integration_distance:
        Threshold distance in km.
    :param sites:
        Instance of :class:`openquake.hazardlib.site.SiteCollection`
        to filter.
    :returns:
        Filtered :class:`~openquake.hazardlib.site.SiteCollection`.

    This function is similar to :meth:`openquake.hazardlib.source.base.BaseSeismicSource.filter_sites_by_distance_to_source`.
    The same notes about filtering criteria apply. Site
    should not be filtered out if it is not further than the integration
    distance from the rupture's surface projection along the great
    circle arc (this is known as Joyner-Boore distance, :meth:`
    openquake.hazardlib.geo.surface.base.BaseQuadrilateralSurface.get_joyner_boore_distance`).
    """
    jb_dist = rupture.surface.get_joyner_boore_distance(sites.mesh)
    return sites.filter(jb_dist <= integration_distance)


def getdefault(dic_with_default, key):
    """
    :param dic_with_default: a dictionary with a 'default' key
    :param key: a key that may be present in the dictionary or not
    :returns: the value associated to the key, or to 'default'
    """
    try:
        return dic_with_default[key]
    except KeyError:
        return dic_with_default['default']


def get_distances(rupture, mesh, param):
    """
    :param rupture: a rupture
    :param mesh: a mesh of points
    :param param: the kind of distance to compute (default rjb)
    :returns: an array of distances from the given mesh
    """
    if param == 'rrup':
        dist = rupture.surface.get_min_distance(mesh)
    elif param == 'rx':
        dist = rupture.surface.get_rx_distance(mesh)
    elif param == 'ry0':
        dist = rupture.surface.get_ry0_distance(mesh)
    elif param == 'rjb':
        dist = rupture.surface.get_joyner_boore_distance(mesh)
    elif param == 'rhypo':
        dist = rupture.hypocenter.distance_to_mesh(mesh)
    elif param == 'repi':
        dist = rupture.hypocenter.distance_to_mesh(mesh, with_depths=False)
    elif param == 'rcdpp':
        dist = rupture.get_cdppvalue(mesh)
    elif param == 'azimuth':
        dist = rupture.surface.get_azimuth(mesh)
    else:
        raise ValueError('Unknown distance measure %r' % param)
    return dist


class FarAwayRupture(Exception):
    """Raised if the rupture is outside the maximum distance for all sites"""


class IntegrationDistance(collections.Mapping):
    """
    Pickleable object wrapping a dictionary of integration distances per
    tectonic region type. The integration distances can be scalars or
    list of pairs (magnitude, distance). Here is an example using 'default'
    as tectonic region type sot that the same values will be used for all
    tectonic region types:

    >>> maxdist = IntegrationDistance({'default': [
    ...          (1, 10), (2, 20), (3, 30), (4, 40), (5, 100), (6, 200),
    ...          (7, 400), (8, 800)]})
    >>> maxdist('Some TRT', mag=5.5)
    array(150.0)

    It has also a method `.get_closest(sites, rupture)` returning the closest
    sites to the rupture and their distances. The integration distance can be
    missing if the sites have been already filtered (empty dictionary): in
    that case the method returns all the sites and all the distances.
    """
    def __init__(self, dic):
        self.dic = dic  # TRT -> float or list of pairs
        self.magdist = {}  # TRT -> (magnitudes, distances)
        for trt, value in dic.items():
            if isinstance(value, list):  # assume a list of pairs (mag, dist)
                value.sort()  # make sure the list is sorted by magnitude
                self.magdist[trt] = zip(*value)
            else:
                self.dic[trt] = float(value)

    def __call__(self, trt, mag=None):
        value = getdefault(self.dic, trt)
        if isinstance(value, float):  # scalar maximum distance
            return value
        elif mag is None:  # get the maximum magnitude distance
            return value[-1][1]
        elif not hasattr(self, 'interp'):
            self.interp = {}  # function cache
        try:
            md = self.interp[trt]  # retrieve from the cache
        except KeyError:  # fill the cache
            magdist = getdefault(self.magdist, trt)
            md = self.interp[trt] = interp1d(
                *magdist, bounds_error=False, fill_value='extrapolate')
        return md(mag)

    def get_closest(self, sites, rupture, distance_type='rrup'):
        """
        :param sites: a (Filtered)SiteColletion
        :param rupture: a rupture
        :param distance_type: default 'rjb'
        :returns: (close sites, close distances)
        :raises: a FarAwayRupture exception if the rupture is far away
        """
        distances = get_distances(rupture, sites.mesh, distance_type)
        if not self.dic:  # for sites already filtered
            return sites, distances
        mask = distances <= self(rupture.tectonic_region_type, rupture.mag)
        if mask.any():
            return sites.filter(mask), distances[mask]
        else:
            raise FarAwayRupture

    def __getitem__(self, trt):
        return self(trt)

    def __iter__(self):
        return iter(self.dic)

    def __len__(self):
        return len(self.dic)

    def __repr__(self):
        return repr(self.dic)


class SourceFilter(object):
    """
    The SourceFilter uses the rtree library if available. The index is
    generated at instantiation time and kept in memory. The filter should be
    instantiated only once per calculation, after the site collection is
    known. It should be used as follows::

      ss_filter = SourceFilter(sitecol, integration_distance)
      for src, sites in ss_filter(sources):
         do_something(...)

    As a side effect, sets the `.nsites` attribute of the source, i.e. the
    number of sites within the integration distance. Notice that SourceFilter
    instances can be pickled, but when unpickled the `use_rtree` flag is set to
    false and the index is lost: the reason is that libspatialindex indices
    cannot be properly pickled (https://github.com/Toblerity/rtree/issues/65).

    :param sitecol:
        :class:`openquake.hazardlib.site.SiteCollection` instance (or None)
    :param integration_distance:
        Threshold distance in km, this value gets passed straight to
        :meth:`openquake.hazardlib.source.base.BaseSeismicSource.filter_sites_by_distance_to_source`
        which is what is actually used for filtering.
    :param use_rtree:
        by default True, i.e. try to use the rtree module if available
    """
    def __init__(self, sitecol, integration_distance, use_rtree=True):
        self.integration_distance = (
            IntegrationDistance(integration_distance)
            if isinstance(integration_distance, dict)
            else integration_distance)
        self.sitecol = sitecol
        self.use_rtree = use_rtree and rtree and (
            integration_distance and sitecol is not None and
            sitecol.at_sea_level())
        if self.use_rtree:
            fixed_lons, self.idl = fix_lons_idl(sitecol.lons)
            self.index = rtree.index.Index(
                (sid, (lon, lat, lon, lat), None)
                for sid, lon, lat in
                zip(sitecol.sids, fixed_lons, sitecol.lats))
        if rtree is None:
            logging.info('Using distance filtering [no rtree]')

    def get_affected_box(self, src):
        """
        Get the enlarged bounding box of a source.

        :param src: a source object
        :returns: a bounding box (min_lon, min_lat, max_lon, max_lat)
        """
        maxdist = self.integration_distance[src.tectonic_region_type]
        min_lon, min_lat, max_lon, max_lat = src.get_bounding_box(maxdist)
        if self.idl:  # apply IDL fix
            if min_lon < 0 and max_lon > 0:
                return max_lon, min_lat, min_lon + 360, max_lat
            elif min_lon < 0 and max_lon < 0:
                return min_lon + 360, min_lat, max_lon + 360, max_lat
            elif min_lon > 0 and max_lon > 0:
                return min_lon, min_lat, max_lon, max_lat
            elif min_lon > 0 and max_lon < 0:
                return max_lon + 360, min_lat, min_lon, max_lat
        else:
            return min_lon, min_lat, max_lon, max_lat

    def get_rectangle(self, src):
        """
        :param src: a source object
        :returns: ((min_lon, min_lat), width, height), useful for plotting
        """
        min_lon, min_lat, max_lon, max_lat = self.get_affected_box(src)
        return (min_lon, min_lat), max_lon - min_lon, max_lat - min_lat

    def get_close_sites(self, source):
        """
        Returns the sites within the integration distance from the source,
        or None.
        """
        source_sites = list(self([source]))
        if source_sites:
            return source_sites[0][1]

    def __call__(self, sources, sites=None):
        if sites is None:
            sites = self.sitecol
        if self.sitecol is None:  # do not filter
            for source in sources:
                yield source, sites
            return
        for src in sources:
            if self.use_rtree:  # Rtree filtering, used in the controller
                box = self.get_affected_box(src)
                sids = numpy.array(sorted(self.index.intersection(box)))
                if len(sids):
                    src.nsites = len(sids)
                    yield src, FilteredSiteCollection(sids, sites.complete)
            elif not self.integration_distance:
                yield src, sites
            else:  # normal filtering, used in the workers
                maxdist = self.integration_distance(src.tectonic_region_type)
                with context(src):
                    s_sites = src.filter_sites_by_distance_to_source(
                        maxdist, sites)
                if s_sites is not None:
                    src.nsites = len(s_sites)
                    yield src, s_sites

    def __getstate__(self):
        return dict(integration_distance=self.integration_distance,
                    sitecol=self.sitecol, use_rtree=False)

source_site_noop_filter = SourceFilter(None, {})
