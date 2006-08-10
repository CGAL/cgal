// Copyright (c) 2005  INRIA Sophia-Antipolis (France) 
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//           Julien Hazebrouck
//           Damien Leroy
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_TYPE_EQUALITY_WRAPPER_H
#define CGAL_SPHERICAL_KERNEL_TYPE_EQUALITY_WRAPPER_H

#include <CGAL/user_classes.h>
#include <CGAL/Kernel/Type_equality_wrapper.h>


namespace CGAL {

template < typename K_base, typename Kernel >
struct Spherical_kernel_type_equality_wrapper
  : public Type_equality_wrapper<K_base, Kernel>
{
    typedef K_base                                  Kernel_base;
    typedef CGAL::Circle_3<Kernel>                  Circle_3;
    typedef CGAL::Circular_arc_point_3<Kernel>      Circular_arc_point_3;
    typedef CGAL::Circular_arc_3<Kernel>            Circular_arc_3;
    typedef CGAL::Line_arc_3<Kernel>                Line_arc_3;
    typedef CGAL::Root_of_2<typename Kernel_base::FT>  Root_of_2;
};

}

#endif // CGAL_SPHERICAL_KERNEL_TYPE_EQUALITY_WRAPPER_H
