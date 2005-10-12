#ifndef CGAL_KDS_INTERNAL_KERNEL_CARTESIAN_KINETIC_KERNEL_BASE_H
#define CGAL_KDS_INTERNAL_KERNEL_CARTESIAN_KINETIC_KERNEL_BASE_H
#include <CGAL/KDS/basic.h>
#include <CGAL/KDS/internal/Kernel/Cartesian_moving_point_1.h>
#include <CGAL/KDS/internal/Kernel/Cartesian_moving_point_2.h>
#include <CGAL/KDS/internal/Kernel/Cartesian_moving_point_3.h>
#include <CGAL/KDS/internal/Kernel/Cartesian_moving_weighted_point_3.h>
#include <CGAL/KDS/internal/Kernel/cartesian_predicates_2.h>
#include <CGAL/KDS/internal/Kernel/cartesian_predicates_3.h>
#include <CGAL/KDS/internal/Kernel/Reverse_time.h>
#include <CGAL/KDS/internal/Kernel/Delaunay_lifting.h>
#include <CGAL/KDS/internal/Kernel/Center.h>

CGAL_KDS_BEGIN_INTERNAL_NAMESPACE

//! A kinetic kernel using cartesian coordinates
/*!
  It takes a PolynomialKernel as a template parameter. The PolynomialKernel is used to define the Motion_function and the Certificate_function.
*/
template <class Polynomial_k, class This>
class Cartesian_kinetic_kernel_base {
public:
  typedef Cartesian_kinetic_kernel_base<Polynomial_k, This> Base;
  Cartesian_kinetic_kernel_base(Polynomial_k pk): k_(pk){};
  Cartesian_kinetic_kernel_base(){};


  //! The type of function used to represent coordinates.
  typedef typename Polynomial_k::Function Motion_function;
 //! The type of function used to represent the results of certificate computations.
  typedef typename Polynomial_k::Function Certificate_function;

	//! I am not sure if I want to expose this.
  typedef Polynomial_k Polynomial_kernel;
    
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

  //! 2D orientation
  /*!
    Takes 3 Point_2.
  */
  typedef Cartesian_orientation_2<This> Orientation_2;
  Orientation_2 orientation_2_object() const {
    return Orientation_2();
  }
  //! 3D orientation
  /*!
    Takes 4 Point_3.
  */
  typedef Cartesian_orientation_3<This> Orientation_3;
  Orientation_3 orientation_3_object() const {
    return Orientation_3();
  }

  //! The in_circle test.
  typedef Cartesian_side_of_oriented_circle_2<This> Side_of_oriented_circle_2;
  Side_of_oriented_circle_2 side_of_oriented_circle_2_object() const {
    return Side_of_oriented_circle_2();
  }

  //! The 3D in_circle test.
  typedef Cartesian_side_of_oriented_sphere_3<This> Side_of_oriented_sphere_3;
  Side_of_oriented_sphere_3 side_of_oriented_sphere_3_object() const {
    return Side_of_oriented_sphere_3();
  }

  //! The power test for weighted points.
  typedef Cartesian_power_test_3<This> Power_test_3;
  Power_test_3 power_test_3_object() const {
    return Power_test_3();
  }

  //! An orientation test for weighted points.
  typedef Cartesian_weighted_orientation_3<This> Weighted_orientation_3;
  Weighted_orientation_3 weighted_orientation_3_object() const {
    return Weighted_orientation_3();
  }
  
  //! Compare the x coordinates of two points
  typedef Cartesian_less_x_1<This> Less_x_1;
  Less_x_1 less_x_1_object() const {return Less_x_1();}

  //! Compare the x coordinates of two points
  typedef Cartesian_less_x_2<This> Less_x_2;
  Less_x_2 less_x_2_object() const {return Less_x_2();}

  //! Compare the y coordinate of two points
  typedef Cartesian_less_y_2<This> Less_y_2;
  Less_y_2 less_y_2_object() const {return Less_y_2();}

  //! Compare the x coordinate of two points
  typedef Cartesian_less_x_3<This> Less_x_3;
  Less_x_3 less_x_3_object() const {return Less_x_3();}

  //! Compare the y coordinate of two points
  typedef Cartesian_less_y_3<This> Less_y_3;
  Less_y_3 less_y_3_object() const {return Less_y_3();}

  //! Compare the z coordinate of two points
  typedef Cartesian_less_z_3<This> Less_z_3;
  Less_z_3 less_z_3_object() const {return Less_z_3();}

  //! computes the lifted coordinate under the lifting map
  typedef Delaunay_lifting<This> Delaunay_lifting_3;
  Delaunay_lifting_3 Delaunay_lifting_3_object() const {
    return Delaunay_lifting_3();
  }
  //! Finds the center of an object
  typedef Center<This> Center_3;
  Center_3 center_3_object() const {
    return Center_3();
  }
  //! Return the PolynomialKernel
  const Polynomial_kernel &polynomial_kernel_object() const{
    return k_;
  }
  
  typedef Reverse_time<This> Reverse_time;
  Reverse_time reverse_time_object() const {
    return Reverse_time(k_.negate_variable_object());
  }

protected:
  Polynomial_kernel k_;
};

CGAL_KDS_END_INTERNAL_NAMESPACE

//#include <CGAL/KDS_internals/kernel_undefs.h>
#endif
