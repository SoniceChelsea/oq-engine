# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2020 GEM Foundation
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
Module exports class:`MorikawaFujiwara2013`
"""

import numpy as np

from openquake.hazardlib.gsim.base import GMPE, CoeffsTable
from openquake.hazardlib import const
from openquake.hazardlib.imt import SA


class Mori2020(GMPE):
    """
    Implements the GMM from Mori et.al published as "Ground motion prediction equation for the Kathmandu Valley, Nepal based
on strong motion records during the 2015 Gorkha Nepal
earthquake sequence",Soil Dynamics and Earthquake Engineering 135(2020) 106208.
    """

    #: Supported tectonic region type is active shallow crust
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST

    #: Supported intensity measure types are spectral acceleration
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([
        SA
    ])

    #: Supported intensity measure component is orientation-independent
    #: measure :attr:`~openquake.hazardlib.const.IMC.RotD50`
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    #: Supported standard deviation types are inter-event, intra-event
    #: and total, see equation 2, pag 106.
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = set([
        const.StdDev.TOTAL,
    ])

    #: Required site parameters are:
    #: - d3200 - Depth of bedrock [m]
    REQUIRES_SITES_PARAMETERS = {'d3200'}

    #: Required rupture parameters are magnitude/
    REQUIRES_RUPTURE_PARAMETERS = {'mag'}

    #: Required distance measure is Rrup [km]
    REQUIRES_DISTANCES = {'rrup'}

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        C = self.COEFFS[imt]
        mw = rup.mag

        mag_term_exceptc = self._get_magnitude_term_exceptc(C, dists.rrup, mw)
        c_val = self._cvalue(C,sites.d3200)
        mean = mag_term_exceptc + c_val

        stddevs = self._compute_stddevs(C, len(sites.d3200), stddev_types)
        mean = np.log(10**mean/980.665)
        return mean, stddevs

    def _compute_stddevs(self, C, num_sites, stddev_types):
        """ Return total standard deviation (converted to base e) """
        stddevs = []
        for _ in stddev_types:
            stddevs.append(np.zeros(num_sites) + C['sigma'] * np.log(10))
        return stddevs

    def _get_magnitude_term_exceptc(self, C, rrup, mw):

        return (C['a']*(mw) + C['b'] * rrup + C['c'] -
                np.log10(rrup + C['d'] * 10.**(self.CONSTS['e']*mw)))

    def _cvalue(self, C, d3200):
        tmpp = np.ones_like(d3200) * C['p']
        tmpq = np.ones_like(d3200) * C['q']
        return tmpp * d3200 + tmpq


    COEFFS = CoeffsTable(sa_damping=5, table="""\
     IMT       a         b         d       p       q   sigma
    0.05  0.5411 -0.005398  0.008437     NaN     NaN  0.4897
    0.06  0.5549 -0.005602  0.009935     NaN     NaN  0.4657
    0.07  0.5411 -0.005783  0.012141     NaN     NaN  0.6473
    0.08  0.5440 -0.005908  0.014089     NaN     NaN  0.7441
    0.09  0.5460 -0.005977  0.015398     NaN     NaN  0.7660
     0.1  0.5605 -0.005942  0.015938     NaN     NaN  0.6764
    0.11  0.5513 -0.005838  0.014729     NaN     NaN  0.8040
    0.12  0.5509 -0.005734  0.013673     NaN     NaN  0.9029
    0.13  0.5357 -0.005602  0.012742     NaN     NaN  0.9941
    0.15  0.5478 -0.005408  0.011173     NaN     NaN  0.8719
    0.17  0.5815 -0.005198  0.009898     NaN     NaN  0.6527
     0.2  0.5935 -0.004829  0.008375     NaN     NaN  0.4841
    0.22  0.5744 -0.004623  0.007550     NaN     NaN  0.5493
    0.25  0.5271 -0.004309  0.006522     NaN     NaN  0.8545
     0.3  0.5816 -0.003880  0.005205     NaN     NaN  0.4096
    0.35  0.5299 -0.003504  0.004226     NaN     NaN  0.7306
     0.4  0.5237 -0.003134  0.003474     NaN     NaN  0.7353
    0.45  0.5424 -0.002863  0.002883     NaN     NaN  0.4184
     0.5  0.5447 -0.002642  0.002411     NaN     NaN  0.3931
     0.6  0.5994 -0.002285  0.001717     NaN     NaN  0.0337
     0.7  0.6102 -0.001921  0.001246     NaN     NaN  0.1811
     0.8  0.6390 -0.001580  0.000923     NaN     NaN  0.4920
     0.9  0.6843 -0.001357  0.000700     NaN     NaN  0.8187
       1  0.7106 -0.001256  0.000550  0.0014 -1.4695  0.4091
     1.1  0.7194 -0.001092  0.000453  0.0016 -1.6260  0.4074
     1.2  0.7089 -0.001019  0.000395  0.0017 -1.6302  0.4061
     1.3  0.7060 -0.000956  0.000368  0.0017 -1.6463  0.4046
     1.5  0.6600 -0.000002  0.000801  0.0017 -1.4361  0.4035
     1.7  0.6950 -0.000002  0.000707  0.0019 -1.7647  0.4007
       2  0.6710 -0.000002  0.000655  0.0021 -1.7460  0.3927
     2.2  0.6767 -0.000002  0.000591  0.0022 -1.8259  0.3883
     2.5  0.6918 -0.000003  0.000562  0.0025 -2.0804  0.3831
       3  0.6965 -0.000003  0.000470  0.0029 -2.3284  0.3775
     3.5  0.7164 -0.000004  0.000469  0.0027 -2.5713  0.3713
       4  0.7456 -0.000004  0.000561  0.0024 -2.7991  0.3646
     4.5  0.8042 -0.000005  0.000552  0.0020 -3.2100  0.3603
       5  0.8543 -0.000005  0.000510  0.0018 -3.5390  0.3552
     5.5  0.8665 -0.000006  0.000573  0.0016 -3.6477  0.3494
       6  0.8938 -0.000006  0.000651  0.0014 -3.8600  0.3428
     6.5  0.9085 -0.000007  0.000736  0.0014 -4.0125  0.3366
       7  0.9167 -0.000007  0.000811  0.0014 -4.1215  0.3300
     7.5  0.9223 -0.000008  0.000840  0.0013 -4.2077  0.3242
       8  0.9190 -0.000008  0.000873  0.0013 -4.2334  0.3185
     8.5  0.9129 -0.000009  0.000911  0.0013 -4.2545  0.3130
       9  0.9024 -0.000009  0.000990  0.0012 -4.2472  0.3090
     9.5  0.8902 -0.000010  0.001083  0.0012 -4.2158  0.3047
      10  0.8774 -0.000010  0.001171  0.0012 -4.1854  0.3007
 """)


    CONSTS = {
        "D0": 300.,
        "e": 0.5}
