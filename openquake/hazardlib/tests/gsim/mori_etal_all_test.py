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
import os
import unittest
from openquake.hazardlib.gsim.mori_etal_all import (
        Mori2020,
        )
from openquake.hazardlib.tests.gsim.utils import BaseGSIMTestCase


class Mori2020Test(BaseGSIMTestCase):
    GSIM_CLASS = Mori2020

    def test_mean(self):
        self.check('Mori20/mean_mori.csv', max_discrep_percentage=0.1)

    def test_sigma(self):
        self.check('Mori20/total_std_mori.csv', max_discrep_percentage=0.1)

