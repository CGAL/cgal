// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_INTERNAL_KERNEL_CARTESIAN_KINETIC_KERNEL_BASE_H
#define CGAL_KINETIC_INTERNAL_KERNEL_CARTESIAN_KINETIC_KERNEL_BASE_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/internal/Kernel/Cartesian_moving_point_1.h>
#include <CGAL/Kinetic/internal/Kernel/Cartesian_moving_point_2.h>
#include <CGAL/Kinetic/internal/Kernel/Cartesian_moving_point_3.h>
#include <CGAL/Kinetic/internal/Kernel/Cartesian_moving_weighted_point_3.h>
#include <CGAL/Kinetic/internal/Kernel/cartesian_predicates_2.h>
#include <CGAL/Kinetic/internal/Kernel/cartesian_predicates_3.h>
#include <CGAL/Kinetic/internal/Kernel/Reverse_time.h>
#include <CGAL/Kinetic/internal/Kernel/Delaunay_lifting.h>
#include <CGAL/Kinetic/Certificate_generator.h>
#include <CGAL/Kinetic/internal/Kernel/Certificate.h>
#include <CGAL/Kinetic/internal/Kernel/Center.h>

CGAL_KINETIC_BEGIN_INTERNAL_NAMESPACE

//! A kinetic kernel using cartesian coordinates
/*!
  It takes a PolynomialKernel as a template parameter. The PolynomialKernel is used to define the Motion_function and the Certificate_function.
*/
template <class Function_kernel_k, class This>
class Cartesian_kinetic_kernel_base
{
public:
  typedef Cartesian_kinetic_kernel_base<Function_kernel_k, This> Base;
  Cartesian_kinetic_kernel_base(Function_kernel_k pk): k_(pk){};
  Cartesian_kinetic_kernel_base(){};

  //! The type of function used to represent coordinates.
  typedef typename Function_kernel_k::Function Motion_function;
  //! The type of function used to represent the results of certificate 
  typedef typename Function_kernel_k::Function Certificate_function;

  typedef CGAL::Kinetic::internal::Certificate<Function_kernel_k> Certificate;

  //! I am not sure if I want to expose this.
  typedef Function_kernel_k Function_kernel;

  //! A 1d Point
  typedef Cartesian_moving_point_1<Motion_function> Point_1;

  //! A 2d Point
  typedef Cartesian_moving_point_2<Motion_function> Point_2;

  //! A 3d Point
  typedef Cartesian_moving_point_3<Motion_function> Point_3;

  //! A 3d weighted Point
  typedef Cartesian_moving_weighted_point_3<Motion_function> Weighted_point_3;

  //! A 3d lifted Point
  // typedef CGALi::Cartesian_moving_lifted_point_3<Motion_function> Moving_lifted_point_3;


  struct Is_constant {
    template <class T>
    bool operator()(const T&t) const {
      return t.is_constant();
    }
  };

  Is_constant is_constant_object() const {
    return Is_constant();
  }

  //! 2D orientation
  /*!
    Takes 3 Point_2.
  */
  typedef Certificate_generator<This, Cartesian_orientation_2<This> > Positive_orientation_2;
  Positive_orientation_2 positive_orientation_2_object() const
  {
    return Positive_orientation_2(k_);
  }
  //! 3D orientation
  /*!
    Takes 4 Point_3.
  */
  typedef Certificate_generator<This, Cartesian_orientation_3<This> > Positive_orientation_3;
  Positive_orientation_3 positive_orientation_3_object() const
  {
    return Positive_orientation_3(k_);
  }

  typedef  Cartesian_orientation_3<This> Orientation_3;
  Orientation_3 orientation_3_object() const
  {
    return Orientation_3();
  }

  //! The in_circle test.
  typedef Certificate_generator<This, Cartesian_side_of_oriented_circle_2<This> > Positive_side_of_oriented_circle_2;
  Positive_side_of_oriented_circle_2 positive_side_of_oriented_circle_2_object() const
  {
    return Positive_side_of_oriented_circle_2(k_);
  }

  //! The 3D in_circle test.
  typedef Certificate_generator<This, Cartesian_side_of_oriented_sphere_3<This> > Positive_side_of_oriented_sphere_3;
  Positive_side_of_oriented_sphere_3 positive_side_of_oriented_sphere_3_object() const
  {
    return Positive_side_of_oriented_sphere_3(k_);
  }

  //! The power test for weighted points.
  typedef Certificate_generator<This, Cartesian_power_test_3<This> > Positive_power_test_3;
  Positive_power_test_3 positive_power_test_3_object() const
  {
    return Positive_power_test_3(k_);
  }

  //! An orientation test for weighted points.
  typedef Certificate_generator<This, Cartesian_weighted_orientation_3<This> > Weighted_positive_orientation_3;
  Weighted_positive_orientation_3 weighted_positive_orientation_3_object() const
  {
    return Weighted_positive_orientation_3(k_);
  }

  typedef  Cartesian_weighted_orientation_3<This> Weighted_orientation_3;
  Weighted_orientation_3 weighted_orientation_3_object() const
  {
    return Weighted_orientation_3();
  }


  //! Compare the x coordinates of two points
  typedef Certificate_generator<This, Cartesian_less_x_1<This> > Is_less_x_1;
  Is_less_x_1 is_less_x_1_object() const {return Is_less_x_1(k_);}

  //! Compare the x coordinates of two points
  typedef Certificate_generator<This, Cartesian_less_x_2<This> > Is_less_x_2;
  Is_less_x_2 is_less_x_2_object() const {return Is_less_x_2(k_);}

  //! Compare the y coordinate of two points
  typedef Certificate_generator<This, Cartesian_less_y_2<This> > Is_less_y_2;
  Is_less_y_2 is_less_y_2_object() const {return Is_less_y_2(k_);}

  //! Compare the x coordinate of two points
  typedef Certificate_generator<This, Cartesian_less_x_3<This> > Is_less_x_3;
  Is_less_x_3 is_less_x_3_object() const {return Is_less_x_3(k_);}

  //! Compare the y coordinate of two points
  typedef Certificate_generator<This, Cartesian_less_y_3<This> > Is_less_y_3;
  Is_less_y_3 is_less_y_3_object() const {return Is_less_y_3(k_);}

  //! Compare the z coordinate of two points
  typedef Certificate_generator<This, Cartesian_less_z_3<This> > Is_less_z_3;
  Is_less_z_3 is_less_z_3_object() const {return Is_less_z_3(k_);}

  //! computes the lifted coordinate under the lifting map
  typedef Delaunay_lifting<This> Delaunay_lifting_3;
  Delaunay_lifting_3 Delaunay_lifting_3_object() const
  {
    return Delaunay_lifting_3();
  }
  //! Finds the center of an object
  typedef Center<This> Center_3;
  Center_3 center_3_object() const
  {
    return Center_3();
  }
  //! Return the PolynomialKernel
  const Function_kernel &function_kernel_object() const
  {
    return k_;
  }

  typedef internal::Reverse_time<This> Reverse_time;
  Reverse_time reverse_time_object() const
  {
    return Reverse_time(k_.negate_variable_object());
  }

protected:
  Function_kernel k_;
};

CGAL_KINETIC_END_INTERNAL_NAMESPACE

//#include <CGAL/Kinetic_internals/kernel_undefs.h>
#endif
