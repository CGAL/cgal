// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
//
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado, 
//             Sebastien Loriot, Julien Hazebrouck, Damien Leroy

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_GLOBAL_FUNCTIONS_CIRCULAR_KERNEL_3_H
#define CGAL_SPHERICAL_KERNEL_GLOBAL_FUNCTIONS_CIRCULAR_KERNEL_3_H

// global functions
namespace CGAL {

template <class SK>
Circular_arc_point_3<SK>
x_extremal_point(const Circle_3<SK> & c, bool i)
{
  return SphericalFunctors::x_extremal_point<SK>(c,i);
}

template <class SK, class OutputIterator>
OutputIterator
x_extremal_points(const Circle_3<SK> & c, OutputIterator res)
{
  return SphericalFunctors::x_extremal_points<SK>(c,res);
}

template <class SK>
Circular_arc_point_3<SK>
y_extremal_point(const Circle_3<SK> & c, bool i)
{
  return SphericalFunctors::y_extremal_point<SK>(c,i);
}

template <class SK, class OutputIterator>
OutputIterator
y_extremal_points(const Circle_3<SK> & c, OutputIterator res)
{
  return SphericalFunctors::y_extremal_points<SK>(c,res);
}

template <class SK>
Circular_arc_point_3<SK>
z_extremal_point(const Circle_3<SK> & c, bool i)
{
  return SphericalFunctors::z_extremal_point<SK>(c,i);
}

template <class SK, class OutputIterator>
OutputIterator
z_extremal_points(const Circle_3<SK> & c, OutputIterator res)
{
  return SphericalFunctors::z_extremal_point<SK>(c,res);
}

template <class SK>
Circular_arc_point_3<SK>
x_extremal_point(const Sphere_3<SK> & c, bool i)
{
  return SphericalFunctors::x_extremal_point<SK>(c,i);
}

template <class SK, class OutputIterator>
OutputIterator
x_extremal_points(const Sphere_3<SK> & c, OutputIterator res)
{
  return SphericalFunctors::x_extremal_points<SK>(c,res);
}

template <class SK>
Circular_arc_point_3<SK>
y_extremal_point(const Sphere_3<SK> & c, bool i)
{
  return SphericalFunctors::y_extremal_point<SK>(c,i);
}

template <class SK, class OutputIterator>
OutputIterator
y_extremal_points(const Sphere_3<SK> & c, OutputIterator res)
{
  return SphericalFunctors::y_extremal_points<SK>(c,res);
}

template <class SK>
Circular_arc_point_3<SK>
z_extremal_point(const Sphere_3<SK> & c, bool i)
{
  return SphericalFunctors::z_extremal_point<SK>(c,i);
}

template <class SK, class OutputIterator>
OutputIterator
z_extremal_points(const Sphere_3<SK> & c, OutputIterator res)
{
  return SphericalFunctors::z_extremal_points<SK>(c,res);
}

template< class CK >
inline
CGAL::Comparison_result 
compare_x(const Circular_arc_point_3<CK> &p, const Circular_arc_point_3<CK> &q)
{
  return CK().compare_x_3_object()(p, q);
}

template< class CK >
inline
CGAL::Comparison_result 
compare_y(const Circular_arc_point_3<CK> &p, const Circular_arc_point_3<CK> &q)
{
  return CK().compare_y_3_object()(p, q);
}

template< class CK >
inline
CGAL::Comparison_result 
compare_z(const Circular_arc_point_3<CK> &p, const Circular_arc_point_3<CK> &q)
{
  return CK().compare_z_3_object()(p, q);
}

template< class CK >
inline
CGAL::Comparison_result 
compare_xy(const Circular_arc_point_3<CK> &p, const Circular_arc_point_3<CK> &q)
{
  return CK().compare_xy_3_object()(p, q);
}

template< class CK >
inline
CGAL::Comparison_result 
compare_xyz(const Circular_arc_point_3<CK> &p, const Circular_arc_point_3<CK> &q)
{
  return CK().compare_xyz_3_object()(p, q);
}

template <class SK>
CGAL::Circle_type
classify(const Circle_3<SK>& c,const Sphere_3<SK> & s)
{
  return SphericalFunctors::classify_circle_3<SK>(c,s);
}

template <class SK>
bool
is_theta_monotone(const Circular_arc_3<SK>& arc,const Sphere_3<SK> & s)
{
  return SphericalFunctors::is_theta_monotone_3<SK>(arc,s);
}

template <class SK>
CGAL::Comparison_result
compare_theta(const Circular_arc_point_3<SK>& pt1,const Circular_arc_point_3<SK>& pt2,const Sphere_3<SK>& sphere)
{
  return SphericalFunctors::compare_theta_of_pts<SK>(pt1,pt2,sphere);
}

template <class SK>
CGAL::Comparison_result
compare_theta(const Circular_arc_point_3<SK>& pt,const Vector_3<SK>& v,const Sphere_3<SK>& sphere)
{
  return SphericalFunctors::compare_theta_pt_vector<SK>(pt,v,sphere);
}

template <class SK>
CGAL::Comparison_result
compare_theta(const Vector_3<SK>& v,const Circular_arc_point_3<SK>& pt,const Sphere_3<SK>& sphere)
{
  return CGAL::opposite(SphericalFunctors::compare_theta_pt_vector<SK>(pt,v,sphere));
}


template <class SK>
CGAL::Comparison_result
compare_theta(const Vector_3<SK>&m1,const Vector_3<SK>&m2)
{ return SphericalFunctors::compare_theta_vectors<SK>(m1,m2); }

template <class SK>
CGAL::Comparison_result
compare_theta_z(const Circular_arc_point_3<SK>& pt1,const Circular_arc_point_3<SK>& pt2,const Sphere_3<SK>& sphere)
{
  return SphericalFunctors::compare_theta_z<SK>(pt1,pt2,sphere);
}


template <class SK>
typename SK::Circular_arc_point_3
theta_extremal_point(const Circle_3<SK>& circle,const Sphere_3<SK>& sphere,bool is_smallest)
{
  return SphericalFunctors::theta_extremal_point<SK>(circle,sphere,is_smallest);
}

template <class SK,class OutputIterator>
OutputIterator
theta_extremal_points(const Circle_3<SK>& circle,const Sphere_3<SK>& sphere,OutputIterator out_it)
{
  return SphericalFunctors::theta_extremal_points<SK>(circle,sphere,out_it);
}

} //namespace CGAL

#endif // CGAL_SPHERICAL_KERNEL_GLOBAL_FUNCTIONS_CIRCULAR_KERNEL_3_H
