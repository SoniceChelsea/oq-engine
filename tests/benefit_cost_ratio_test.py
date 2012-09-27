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

import unittest

from risklib.curve import Curve
from risklib.benefit_cost_ratio import (
    _bcr, _mean_loss_ratio_curve, _mean_loss)


class RiskCommonTestCase(unittest.TestCase):

    def test_compute_bcr(self):
        # numbers are proven to be correct
        eal_orig = 0.00838
        eal_retrofitted = 0.00587
        retrofitting_cost = 0.1
        interest = 0.05
        life_expectancy = 40
        expected_result = 0.43405

        result = _bcr(eal_orig, eal_retrofitted, interest,
                      life_expectancy, retrofitting_cost)
        self.assertAlmostEqual(result, expected_result, delta=2e-5)

    def test_loss_ratio_pe_mid_curve_computation(self):
        expected_curve = Curve([(0.0600, 0.1445),
                (0.1200, 0.0400), (0.1800, 0.0190), (0.2400, 0.0155),
                (0.3225, 0.0095)])

        loss_ratio_curve = Curve([(0, 0.3460), (0.06, 0.12),
                (0.12, 0.057), (0.18, 0.04),
                (0.24, 0.019), (0.3, 0.009), (0.45, 0)])

        self.assertEqual(expected_curve,
                _mean_loss_ratio_curve(loss_ratio_curve))

    def test_mean_loss_ratio_computation(self):
        loss_ratio_curve = Curve([(0, 0.3460), (0.06, 0.12),
                (0.12, 0.057), (0.18, 0.04),
                (0.24, 0.019), (0.3, 0.009), (0.45, 0)])

        self.assertAlmostEqual(0.023305,
                               _mean_loss(loss_ratio_curve), 3)
