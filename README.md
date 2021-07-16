# Summary
 This package provides substitute functions for `GeoDestination`, `GeoDistance` and `GeoDirection`, which are built-in functions of _Methematica_. This package can solve geodesic problems more accurately with arbitrary precision calculation by designating option `WorkingPrecision`. This package also contains interface to GeographicLib via _J/Link_.

 Although this package and [GeographicLib](http://geographiclib.sourceforge.net/) refer to the same paper by Charles F. F. Karney, they differ in the way to calculate. GeographicLib uses functional expansion to calculate elliptic integrals, but this package uses _Mathematica_ built-in functions such as `EllipticE`. Additionally, inverse function of elliptic integral is represented by `InverseFunction`, and `FindRoot` is used to solve the inverse problems numerically. This enables arbitrary precision calculation with _Mathematica_.

 Accuracy of distance calculation is 0.1 m order for `WorkingPrecision->$MachinePrecision` and 0.001m order for `WorkingPrecision->30`, because calculation accuracy of `FindRoot` is a half of `WorkingPrecision`. This is much better than `GeoDistance[Method->"Vintency75"]`(10000m order), but less than GeographicLib(GeographicLib has 15nm accuracy).

# Installation
 To install the package, 
+ download the zip file,
+ extract,
+ move "Geodesic" directory to the directory _Mathematica_ can find such as `$UserBaseDirectory` or another one of `$Path`.

# Usage
 In your _Mathematica_ session,
+ evaluate ``<<Geodesic` `` to load this package,
+ use like built-in function e.g. `GeoDistance2[{lat1, lon1}, {lat2, lon2}]`,
+ or to use GeographicLib, write `GeoDistanceGL[{lat1, lon1}, {lat2, lon2}]`.

# License
 This software is released under the MIT License.
* Mathematica files: Copyright (c) 2014 Kei Misawa, http://opensource.org/licenses/mit-license.php
* Java compiled files from GeographicLib: Copyright (c) 2013 Charles F.F. Karney.

# Resources
* Charles F. F. Karney, "Algorithms for geodesics," J. Geodesy 87, 43-55 (2013), http://dx.doi.org/10.1007/s00190-012-0578-z
* Charles F. F. Karney, "Geodesics on an ellipsoid of revolution," Feb. 2011, http://arxiv.org/abs/1102.1215v1
