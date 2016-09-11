// Copyright (c) 2009  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Boolean_set_operations_2/demo/Boolean_set_operations_2/typedefs.h $
// $Id: typedefs.h 37003 2007-03-10 16:55:12Z spion $
//
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_TYPEDEFS_H
#define CGAL_TYPEDEFS_H

//
// Linear polygons
//
typedef CGAL::Simple_cartesian<double>            Linear_kernel ;
typedef CGAL::Polygon_2<Linear_kernel>            Linear_polygon;
typedef CGAL::Polygon_with_holes_2<Linear_kernel> Linear_polygon_with_holes;

typedef Linear_kernel::Point_2 Linear_point ;

//
// Circlular polygons
//

#ifdef CGAL_USE_GMP

  typedef CGAL::Gmpq                     Base_nt;

#else

  typedef CGAL::Quotient<CGAL::MP_Float> Base_nt;

#endif

typedef CGAL::Lazy_exact_nt<Base_nt> Coord_type;


struct Gps_circular_kernel : public CGAL::Cartesian<Coord_type> {};

typedef CGAL::Gps_circle_segment_traits_2<Gps_circular_kernel> Circular_traits;
typedef Circular_traits::Curve_2                               Circular_curve;
typedef Circular_traits::X_monotone_curve_2                    Circular_X_monotone_curve;
typedef Circular_traits::Point_2                               Circular_point ;
typedef Circular_traits::Polygon_2                             Circular_polygon;
typedef CGAL::General_polygon_with_holes_2<Circular_polygon>   Circular_polygon_with_holes;
typedef CGAL::General_polygon_set_2<Circular_traits>           Circular_polygon_set;

typedef CGAL::Qt::Circular_set_graphics_item<Circular_polygon_set> Circular_GI;

typedef std::vector<Circular_polygon_with_holes>  Circular_region_source_container ;


//
// Bezier curves typedefs
//
#ifdef CGAL_USE_CORE

typedef CGAL::CORE_algebraic_number_traits            Bezier_nt_traits;
typedef Bezier_nt_traits::Rational                    Bezier_rational;
typedef Bezier_nt_traits::Algebraic                   Bezier_algebraic;

struct Bezier_rat_kernel  : public CGAL::Cartesian<Bezier_rational>  {};
struct Bezier_alg_kernel  : public CGAL::Cartesian<Bezier_algebraic> {};

struct Bezier_traits : public CGAL::Arr_Bezier_curve_traits_2<Bezier_rat_kernel, Bezier_alg_kernel, Bezier_nt_traits> {};
  
typedef Bezier_rat_kernel::Point_2                      Bezier_rat_point;
typedef Bezier_traits::Curve_2                          Bezier_curve;
typedef Bezier_traits::X_monotone_curve_2               Bezier_X_monotone_curve;
typedef Bezier_traits::Point_2                          Bezier_point;
typedef CGAL::Gps_traits_2<Bezier_traits>               Bezier_gps_traits;
typedef Bezier_gps_traits::General_polygon_2            Bezier_polygon;
typedef std::vector<Bezier_polygon>                     Bezier_polygon_vector ;
typedef Bezier_gps_traits::General_polygon_with_holes_2 Bezier_polygon_with_holes;
typedef CGAL::General_polygon_set_2<Bezier_gps_traits>  Bezier_polygon_set ;

typedef CGAL::Qt::Bezier_set_graphics_item<Bezier_polygon_set> Bezier_GI;

typedef std::vector<Bezier_curve>                Bezier_boundary_source ;
typedef std::vector<Bezier_boundary_source>      Bezier_region_source ;
typedef std::vector<Bezier_region_source>        Bezier_region_source_container ;

#endif

#endif
