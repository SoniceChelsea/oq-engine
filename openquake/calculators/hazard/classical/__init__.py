# -*- coding: utf-8 -*-
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
The Classical Probabilistic Seismic Hazard Analysis (cPSHA) approach
allows calculation of hazard curves and hazard maps following the
classical integration procedure (**Cornell [1968]**, **McGuire [1976]**)
as formulated by **Field et al. [2003]**.

Sources:

* | Cornell, C. A. (1968).
  | Engineering seismic risk analysis.
  | Bulletin of the Seismological Society of America, 58:1583–1606.
* | Field, E. H., Jordan, T. H., and Cornell, C. A. (2003).
  | OpenSHA - A developing Community-Modeling
  | Environment for Seismic Hazard Analysis. Seism. Res. Lett., 74:406–419.
* | McGuire, K. K. (1976).
  | Fortran computer program for seismic risk analysis. Open-File report 76-67,
  | United States Department of the Interior, Geological Survey. 102 pages.


*******
Outputs
*******

* Hazard Curves
* Hazard Maps

Hazard Curves
=============

Hazard Curves are discrete functions which describe probability of ground
motion exceedance in a given time frame. Hazard curves are composed of several
key elements:

* **Intensity Measure Levels (IMLs)** - IMLs define the x-axis values (or
  "ordinates") of the curve. IMLs are defined with an Intensity Measure Type
  (IMT) unit. IMLs are a `strictly monotonically increasing sequence`.
* **Probabilitites of Exceedance (PoEs)** - PoEs define the y-axis values, (or
  "abscissae") of the curve. For each node in the curve, the PoE denotes the
  probability that ground motion will exceedence a given level in a given time
  span.
* **Intensity Measure Type (IMT)** - The unit of measurement for the defined
  IMLs.
* **Investigation time** - The period of time (in years) for an earthquake
  hazard study. It is important to consider the investigation time when
  analyzing hazard curve results, because one can logically conclude that, the
  longer the time span, there is greater probability of ground motion exceeding
  the given values.
* **Spectral Acceleration (SA) Period** - Optional; used only if the IMT is SA.
* **Spectral Acceleration (SA) Damping** - Optional; used only if the IMT is
  SA.

For a given calculation, hazard curves are computed for each logic tree
realization, each IMT/IML definition, and each geographical point of interest.
(In other words: If a calculation specifies 4 logic tree samples, a geometry
with 10 points of interest, and 3 IMT/IML definitions, 120 curves will be
computed.)

Another way to put it is:

``T = R * P * I``

where

* ``T`` is the total number of curves
* ``R`` is the total number of logic tree realizations
* ``P`` is the number of geographical points of interest
* ``I`` is the number of IMT/IML definitions

Statistical Curves
------------------

The classical hazard calculator is also capable of producing `mean` and
`quantile` curves. These aggregates are computed from the curves for a given
point and IMT over all logic tree realizations.

Mean Curves
^^^^^^^^^^^

When computing a mean hazard curve for a given point/IMT, there are two
possible approaches:

1. mean, unweighted
2. mean, weighted

Technically, both approaches are "weighted". In the first approach, however,
the weights are `implicit` and are taken into account in the process of logic
tree sampling. This approach is used in the case of random Monte-Carlo logic
tree sampling. The total of number of logic tree samples is defined by the user
with the `number_of_logic_tree_samples` configuration parameter.

In the second approach, the weights are explicit in the
caluclation of the mean. This approach is used in the case of end-branch logic
tree enumeration, whereby each possible logic tree path is traversed. (Each
logic tree path in this case defines a weight.) The total number of logic tree
samples in this case is determined by the total number of possible tree paths.
(To perform end-branch enumeration, the user must specify
`number_of_logic_tree_samples = 0` in the job configuration.

Quantile Curves
^^^^^^^^^^^^^^^

Similar to mean curves, quantiles curves can be produced for a given
point/IMT/quantile level (in the range [0.0, 1.0]), there are two possible
approaches:

1. quantile, unweighted
2. quantile, weighted

As with mean curves, `unweighted quantiles` are calculated when Monte-Carlo
logic tree sampling is used and `weighted quantiles` are calculated when logic
tree end-branch enumeration is used.
"""
