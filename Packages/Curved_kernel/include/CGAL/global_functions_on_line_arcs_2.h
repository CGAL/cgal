// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
//                     National & Kapodistrian University of Athens (Greece).
// All rights reserved.
//
// Authors : Monique Teillaud     <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion         <Sylvain.Pion@sophia.inria.fr>
//           Athanasios Kakargias <grad0460@di.uoa.gr>
// 
// Partially supported by INRIA's project "CALAMATA", a
// bilateral collaboration with National Kapodistrian University of
// Athens.
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (CGAL - Effective Computational Geometry for Curves and Surfaces) 
// All rights reserved.

// file : include/CGAL/global_functions_on_circular_arcs_2.h

#ifndef CGAL_GLOBAL_FUNCTIONS_ON_LINE_ARCS_2_H
#define CGAL_GLOBAL_FUNCTIONS_ON_LINE_ARCS_2_H

// global functions

#include <CGAL/Circular_arc_2.h>
#include <CGAL/Circular_arc_endpoint_2.h>
#include <CGAL/Line_arc_2.h>

namespace CGAL {

// TODO : Add the other ones...

//Circles

//template< class CK >
//inline
//CGAL::Comparison_result 
//compare_x(const Line_arc_2<CK> &A1, const bool b1, 
//	  const Line_arc_2<CK> &A2, const bool b2)
//{
//  return CK().compare_x_2_object()(A1, b1, A2, b2);
//}

//template< class CK >
//inline
//CGAL::Comparison_result 
//compare_x(const Circular_arc_point_2<CK> &p, const Circular_arc_point_2<CK> &q)
//{
//  return CK().compare_x_2_object()(p, q);
//}
//
//template< class CK >
//inline
//CGAL::Comparison_result 
//compare_y(const Circular_arc_point_2<CK> &p, const Circular_arc_point_2<CK> &q)
//{
//  return CK().compare_y_2_object()(p, q);
//}
//
//template< class CK >
//inline
//CGAL::Comparison_result 
//compare_xy(const Circular_arc_point_2<CK> &p, const Circular_arc_point_2<CK> &q)
//{
//  return CK().compare_xy_2_object()(p, q);
//}

template< class CK >
inline
CGAL::Comparison_result 
compare_y_to_right(const Line_arc_2<CK> &A1,
		   const Line_arc_2<CK> &A2,
		   const bool b)
{
  return CK().compare_y_to_right_2_object()(A1, A2, b);
}

template < class CK >
inline
bool
point_in_range(const Line_arc_2<CK> &A, const Circular_arc_point_2<CK> &p) 
{
  return CK().in_range_2_object()(A, p);
}

template < class CK >
CGAL::Comparison_result       
compare_y_at_x(const Circular_arc_point_2<CK> &p, const Line_arc_2<CK> &a)
{
  return CK().compare_y_at_x_2_object()(p, a);
}

template < class CK, class OutputIterator >
OutputIterator
make_x_monotone(const Line_arc_2<CK> &A, OutputIterator it)
{
  return CK().make_x_monotone_2_object()(A, it);
}

} // namespace CGAL

#endif // CGAL_GLOBAL_FUNCTIONS_ON_LINE_ARCS_2_H
