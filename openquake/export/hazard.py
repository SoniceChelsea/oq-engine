# Copyright (c) 2010-2012, GEM Foundation.
#
# OpenQuake is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenQuake is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake.  If not, see <http://www.gnu.org/licenses/>.

"""
Functionality for exporting and serializing hazard curve calculation results.
"""


import os

from collections import OrderedDict

from nhlib.calc import disagg
from nrml import writers as nrml_writers

from openquake import logs
from openquake.db import models
from openquake.export import core


LOG = logs.LOG


def export(output_id, target_dir):
    """
    Export the given hazard calculation output from the database to the
    specified directory.

    :param int output_id:
        ID of a :class:`openquake.db.models.Output`.
    :param str target_dir:
        Directory where output artifacts should be written.
    :returns:
        List of file names (including the full directory path) containing the
        exported results.

        The quantity and type of the result files depends on
        the type of output, as well as calculation parameters. (See the
        `output_type` attribute of :class:`openquake.db.models.Output`.)
    """
    output = models.Output.objects.get(id=output_id)

    export_fn = _export_fn_map().get(
        output.output_type, core._export_fn_not_implemented)

    return export_fn(output, os.path.expanduser(target_dir))


def _export_fn_map():
    """
    Creates a mapping from output type to export function.

    Each export function should implement a common interface and accept two
    arguments: a :class:`~openquake.db.models.Output` object and a target
    dir (`str`).

    Each function should return a list of the file names created by the export
    action.
    """

    fn_map = {
        'hazard_curve': export_hazard_curves,
        'gmf': export_gmf,
        'gmf_scenario': export_gmf_scenario,
        'ses': export_ses,
        'complete_lt_ses': export_ses,
        'complete_lt_gmf': export_gmf,
        'hazard_map': export_hazard_map,
        'disagg_matrix': export_disagg,
    }
    return fn_map


HAZARD_CURVES_FILENAME_FMT = 'hazard-curves-%(hazard_curve_id)s.xml'
HAZARD_MAP_FILENAME_FMT = 'hazard-map-%(hazard_map_id)s.xml'
GMF_FILENAME_FMT = 'gmf-%(gmf_coll_id)s.xml'
SES_FILENAME_FMT = 'ses-%(ses_coll_id)s.xml'
COMPLETE_LT_SES_FILENAME_FMT = 'complete-lt-ses-%(ses_coll_id)s.xml'
COMPLETE_LT_GMF_FILENAME_FMT = 'complete-lt-gmf-%(gmf_coll_id)s.xml'
GMF_SCENARIO_FMT = 'gmf-%(output_id)s.xml'


@core.makedirs
def export_hazard_curves(output, target_dir):
    """
    Export the specified hazard curve ``output`` to the ``target_dir``.

    :param output:
        :class:`openquake.db.models.Output` with an `output_type` of
        `hazard_curve`.
    :param str target_dir:
        Destination directory location for exported files.

    :returns:
        A list of exported file names (including the absolute path to each
        file).
    """
    hc = models.HazardCurve.objects.get(output=output.id)
    hcd = models.HazardCurveData.objects.filter(hazard_curve=hc.id)

    filename = HAZARD_CURVES_FILENAME_FMT % dict(hazard_curve_id=hc.id)
    path = os.path.abspath(os.path.join(target_dir, filename))

    if hc.lt_realization is not None:
        # If the curves are for a specified logic tree realization,
        # get the tree paths
        lt_rlz = hc.lt_realization
        smlt_path = core.LT_PATH_JOIN_TOKEN.join(lt_rlz.sm_lt_path)
        gsimlt_path = core.LT_PATH_JOIN_TOKEN.join(lt_rlz.gsim_lt_path)
    else:
        # These curves must be statistical aggregates
        smlt_path = None
        gsimlt_path = None

    metadata = {
        'quantile_value': hc.quantile,
        'statistics': hc.statistics,
        'smlt_path': smlt_path,
        'gsimlt_path': gsimlt_path,
        'sa_period': hc.sa_period,
        'sa_damping': hc.sa_damping,
        'investigation_time': hc.investigation_time,
        'imt': hc.imt,
        'imls': hc.imls,
    }
    writer = nrml_writers.HazardCurveXMLWriter(path, **metadata)
    writer.serialize(hcd)

    return [path]


@core.makedirs
def export_gmf(output, target_dir):
    """
    Export the GMF Collection specified by ``output`` to the ``target_dir``.

    :param output:
        :class:`openquake.db.models.Output` with an `output_type` of `gmf`.
    :param str target_dir:
        Destination directory location for exported files.

    :returns:
        A list of exported file names (including the absolute path to each
        file).
    """
    gmf_coll = models.GmfCollection.objects.get(output=output.id)
    lt_rlz = gmf_coll.lt_realization

    if output.output_type == 'complete_lt_gmf':
        filename = COMPLETE_LT_GMF_FILENAME_FMT % dict(gmf_coll_id=gmf_coll.id)

        # For the `complete logic tree` GMF, the LT paths are not relevant.
        sm_lt_path = None
        gsim_lt_path = None
    else:
        # output type should be `gmf`
        filename = GMF_FILENAME_FMT % dict(gmf_coll_id=gmf_coll.id)

        sm_lt_path = core.LT_PATH_JOIN_TOKEN.join(lt_rlz.sm_lt_path)
        gsim_lt_path = core.LT_PATH_JOIN_TOKEN.join(lt_rlz.gsim_lt_path)

    path = os.path.abspath(os.path.join(target_dir, filename))

    writer = nrml_writers.EventBasedGMFXMLWriter(
        path, sm_lt_path, gsim_lt_path)
    writer.serialize(gmf_coll)

    return [path]


@core.makedirs
def export_gmf_scenario(output, target_dir):
    """
    Export the GMFs specified by ``output`` to the ``target_dir``.

    :param output:
        :class:`openquake.db.models.Output`
        with an `output_type` of `gmf_scenario`.
    :param str target_dir:
        Destination directory location for exported files.

    :returns:
        A list of exported file names (including the absolute path to each
        file).
    """
    gmfs = models.GmfScenario.objects.get(output=output.id)
    filename = GMF_SCENARIO_FMT % dict(output_id=output.id)
    path = os.path.abspath(os.path.join(target_dir, filename))
    writer = nrml_writers.ScenarioGMFXMLWriter(path)
    writer.serialize(gmfs)
    return [path]


@core.makedirs
def export_ses(output, target_dir):
    """
    Export the Stochastic Event Set Collection specified by ``output`` to the
    ``target_dir``.

    :param output:
        :class:`openquake.db.models.Output` with an `output_type` of `ses`.
    :param str target_dir:
        Destination directory location for exported files.

    :returns:
        A list of exported file names (including the absolute path to each
        file).
    """
    ses_coll = models.SESCollection.objects.get(output=output.id)
    # lt_rlz can be `None` in the case of a `complete logic tree` SES
    lt_rlz = ses_coll.lt_realization

    if output.output_type == 'complete_lt_ses':
        filename = COMPLETE_LT_SES_FILENAME_FMT % dict(ses_coll_id=ses_coll.id)

        # For the `complete logic tree` SES, the LT paths are not relevant.
        sm_lt_path = None
        gsim_lt_path = None
    else:
        # output_type should be `ses`
        filename = SES_FILENAME_FMT % dict(ses_coll_id=ses_coll.id)

        sm_lt_path = core.LT_PATH_JOIN_TOKEN.join(lt_rlz.sm_lt_path)
        gsim_lt_path = core.LT_PATH_JOIN_TOKEN.join(lt_rlz.gsim_lt_path)

    path = os.path.abspath(os.path.join(target_dir, filename))

    writer = nrml_writers.SESXMLWriter(path, sm_lt_path, gsim_lt_path)
    writer.serialize(ses_coll)

    return [path]


@core.makedirs
def export_hazard_map(output, target_dir):
    """
    Export the specified hazard map ``output`` to the ``target_dir``.

    :param output:
        :class:`openquake.db.models.Output` with an `output_type` of
        `hazard_map`.
    :param str target_dir:
        Destination directory location for exported files.

    :returns:
        A list of exported file name (including the absolute path to each
        file).
    """
    hazard_map = models.HazardMap.objects.get(output=output)

    filename = HAZARD_MAP_FILENAME_FMT % dict(hazard_map_id=hazard_map.id)
    path = os.path.abspath(os.path.join(target_dir, filename))

    if hazard_map.lt_realization is not None:
        # If the maps are for a specified logic tree realization,
        # get the tree paths
        lt_rlz = hazard_map.lt_realization
        smlt_path = core.LT_PATH_JOIN_TOKEN.join(lt_rlz.sm_lt_path)
        gsimlt_path = core.LT_PATH_JOIN_TOKEN.join(lt_rlz.gsim_lt_path)
    else:
        # These maps must be constructed from mean or quantile curves
        smlt_path = None
        gsimlt_path = None

    metadata = {
        'quantile_value': hazard_map.quantile,
        'statistics': hazard_map.statistics,
        'smlt_path': smlt_path,
        'gsimlt_path': gsimlt_path,
        'sa_period': hazard_map.sa_period,
        'sa_damping': hazard_map.sa_damping,
        'investigation_time': hazard_map.investigation_time,
        'imt': hazard_map.imt,
        'poe': hazard_map.poe,
    }

    writer = nrml_writers.HazardMapXMLWriter(path, **metadata)
    writer.serialize(zip(hazard_map.lons, hazard_map.lats, hazard_map.imls))
    return [path]


class _DisaggMatrix(object):
    """
    A simple data model into which disaggregation matrix information can be
    packed. The :class:`nrml.hazard.writers.DisaggXMLWriter` expects a sequence
    of objects which match this interface.

    :param matrix:
        A n-dimensional numpy array representing a probability mass function
        produced by the disaggregation calculator. The calculator produces a 6d
        matrix, but the final results which are saved to XML are "subsets" of
        matrix showing contributions to hazard from different combinations of
        magnitude, distance, longitude, latitude, epsilon, and tectonic region
        type.
    :param dim_labels:
        Expected values are (as tuples, lists, or similar) one of the
        following, depending on the result `matrix` type:

        * ['Mag']
        * ['Dist']
        * ['TRT']
        * ['Mag', 'Dist']
        * ['Mag', 'Dist', 'Eps']
        * ['Lon', 'Lat']
        * ['Mag', 'Lon', 'Lat']
        * ['Lon', 'Lat', 'TRT']
    :param float poe:
        Probability of exceedence (specified in the calculation config file).
    :param float iml:
        Interpolated intensity value, corresponding to the ``poe``, extracted
        from the aggregated hazard curve (which was used as input to compute
        the ``matrix``).
    """

    def __init__(self, matrix, dim_labels, poe, iml):
        self.matrix = matrix
        self.dim_labels = dim_labels
        self.poe = poe
        self.iml = iml


@core.makedirs
def export_disagg(output, target_dir):
    """
    Export disaggregation histograms to the ``target_dir``.

    :param output:
        :class:`openquake.db.models.Output` with an `output_type` of
        `disagg_matrix`.
    :param str target_dir:
        Destination directory location for exported files.

    :returns:
        A list of exported file name (including the absolute path to each
        file).
    """
    # We expect 1 result per `Output`
    [disagg_result] = models.DisaggResult.objects.filter(output=output)
    lt_rlz = disagg_result.lt_realization

    filename = '%s.xml' % output.display_name
    path = os.path.abspath(os.path.join(target_dir, filename))

    pmf_map = OrderedDict([
        (('Mag', ), disagg.mag_pmf),
        (('Dist', ), disagg.dist_pmf),
        (('TRT', ), disagg.trt_pmf),
        (('Mag', 'Dist'), disagg.mag_dist_pmf),
        (('Mag', 'Dist', 'Eps'), disagg.mag_dist_eps_pmf),
        (('Lon', 'Lat'), disagg.lon_lat_pmf),
        (('Mag', 'Lon', 'Lat'), disagg.mag_lon_lat_pmf),
        (('Lon', 'Lat', 'TRT'), disagg.lon_lat_trt_pmf),
    ])

    writer_kwargs = dict(
        investigation_time=disagg_result.investigation_time,
        imt=disagg_result.imt,
        lon=disagg_result.location.x,
        lat=disagg_result.location.y,
        sa_period=disagg_result.sa_period,
        sa_damping=disagg_result.sa_damping,
        mag_bin_edges=disagg_result.mag_bin_edges,
        dist_bin_edges=disagg_result.dist_bin_edges,
        lon_bin_edges=disagg_result.lon_bin_edges,
        lat_bin_edges=disagg_result.lat_bin_edges,
        eps_bin_edges=disagg_result.eps_bin_edges,
        tectonic_region_types=disagg_result.trts,
        smlt_path=core.LT_PATH_JOIN_TOKEN.join(lt_rlz.sm_lt_path),
        gsimlt_path=core.LT_PATH_JOIN_TOKEN.join(lt_rlz.gsim_lt_path),
    )

    writer = nrml_writers.DisaggXMLWriter(path, **writer_kwargs)

    data = (_DisaggMatrix(pmf_fn(disagg_result.matrix), dim_labels,
                          disagg_result.poe, disagg_result.iml)
            for dim_labels, pmf_fn in pmf_map.iteritems())

    writer.serialize(data)

    return [path]
