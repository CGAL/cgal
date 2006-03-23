// Copyright (c) 2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://lrineau@scm.gforge.inria.fr/svn/cgal/trunk/Surface_mesher/examples/Surface_mesher/implicit_function.C $
// $Id: implicit_function.C 29089 2006-03-06 15:17:07Z lrineau $
//
//
// Author(s)     : Steve OUDOT, Laurent RINEAU


///////////////// Definitions of several famous surfaces /////////////////
double sphere_function (double, double, double);  // (c=(0,0,0), r=1)
double ellipsoid_function (double, double, double);  // (c=(0,0,0), r=1)
double torus_function (double, double, double);  // (c=(0,0,0), r=2)
double chair_function (double, double, double);  // (c=(0,0,0), r=6)
double tanglecube_function (double, double, double);  // (c=(0,0,0), r=4)
double octic_function (double, double, double);  // (c=(0,0,0), r=2)
double heart_function (double, double, double);  // (c=(0,0,0), r=2)
double klein_function (double, double, double);  // (c=(0,0,0), r=4)
double ring_function (double, double, double);  // (c=(0,0,0), r=?)
double false_knot_function (double, double, double);  // (c=(0,0,0), r=1)
double knot1_function (double, double, double);  // (c=(0,0,0), r=4)
double knot2_function (double, double, double);  // (c=(0,0,0), r=4)
double knot3_function (double, double, double);  // (c=(0,0,0), r=4)
