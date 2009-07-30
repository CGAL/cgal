// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//
#ifndef CGAL_TYPEDEFS_H
#define CGAL_TYPEDEFS_H

typedef CGAL::Simple_cartesian<double>         Dbl_kernel ;
typedef CGAL::Polygon_with_holes_2<Dbl_kernel> Dbl_polygon_with_holes ;

#ifdef CGAL_USE_GMP

  typedef CGAL::Gmpq                     Base_nt;

#else

  typedef CGAL::Quotient<CGAL::MP_Float> Base_nt;

#endif

typedef CGAL::Lazy_exact_nt<Base_nt> Coord_type;


struct Kernel : public CGAL::Cartesian<Coord_type> {};

typedef Kernel::Segment_2                             Segment;
typedef Kernel::Point_2                               Point_2;
typedef Kernel::Circle_2                              Circle;
typedef Kernel::Iso_rectangle_2                       Iso_rectangle;

//
// Linear polygons
//
typedef CGAL::Polygon_2<Kernel>            Linear_polygon;
typedef CGAL::Polygon_with_holes_2<Kernel> Linear_polygon_with_holes;
typedef CGAL::Polygon_set_2<Kernel>        Linear_polygon_set;

typedef CGAL::Qt::GeneralPolygonSetGraphicsItem<Linear_polygon_set> Linear_GI;


//
// Circle-segment polygons
//

typedef CGAL::Gps_circle_segment_traits_2<Kernel>            Circular_traits;
typedef Circular_traits::Curve_2                             Circular_curve;
typedef Circular_traits::X_monotone_curve_2                  Circular_X_monotone_curve;
typedef Circular_traits::Point_2                             Circular_point ;
typedef Circular_traits::Polygon_2                           Circular_polygon;
typedef CGAL::General_polygon_with_holes_2<Circular_polygon> Circular_polygon_with_holes;
typedef CGAL::General_polygon_set_2<Circular_traits>         Circular_polygon_set;

//typedef Polygon_with_holes::Hole_const_iterator       Hole_const_iterator;

typedef CGAL::Qt::PolygonWithHolesGraphicsItem<Circular_polygon_with_holes
                                              , CGAL::Circular_polygon_with_holes_sampler<Dbl_polygon_with_holes>
                                              > Circular_GI;



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
typedef CGAL::Gps_traits_2<Bezier_traits>               Bezier_gps_traits;
typedef Bezier_gps_traits::General_polygon_2            Bezier_polygon;
typedef std::vector<Bezier_polygon>                     Bezier_polygon_vector ;
typedef Bezier_gps_traits::General_polygon_with_holes_2 Bezier_polygon_with_holes;
typedef CGAL::General_polygon_set_2<Bezier_gps_traits>  Bezier_polygon_set ;

typedef CGAL::Qt::GeneralPolygonSetGraphicsItem<Bezier_polygon_set
                                               ,CGAL::Bezier_polygon_with_holes_sampler<Dbl_polygon_with_holes>
                                               > 
                                               Bezier_GI;


#endif

#endif
