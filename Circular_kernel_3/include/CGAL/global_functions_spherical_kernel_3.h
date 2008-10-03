// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://pmachado@scm.gforge.inria.fr/svn/cgal/trunk/Circular_kernel_2/include/CGAL/global_functions_on_line_arcs_2.h $
// $Id: global_functions_on_line_arcs_2.h 45945 2008-10-01 12:25:23Z pmachado $
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Julien Hazebrouck, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_GLOBAL_FUNCTIONS_CIRCULAR_KERNEL_3_H
#define CGAL_SPHERICAL_KERNEL_GLOBAL_FUNCTIONS_CIRCULAR_KERNEL_3_H

// global functions
CGAL_BEGIN_NAMESPACE

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
  return SphericalFunctors::x_extremal_point<SK>(c,res);
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
  return SphericalFunctors::y_extremal_point<SK>(c,res);
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
  return SphericalFunctors::x_extremal_point<SK>(c,res);
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
  return SphericalFunctors::y_extremal_point<SK>(c,res);
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
  return SphericalFunctors::z_extremal_point<SK>(c,res);
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

CGAL_END_NAMESPACE

#endif // CGAL_SPHERICAL_KERNEL_GLOBAL_FUNCTIONS_CIRCULAR_KERNEL_3_H
