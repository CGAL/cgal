Author: Constantinos Tsirogiannis

#ifndef CGAL_CURVED_KERNEL_TYPE_EQUALITY_WRAPPER_H
#define CGAL_CURVED_KERNEL_TYPE_EQUALITY_WRAPPER_H

#include <CGAL/user_classes.h>
#include <CGAL/Kernel/Type_equality_wrapper.h>
#include <CGAL/Circular_arc_2.h>
#include <CGAL/Circular_arc_endpoint_2.h>
#include <CGAL/Line_arc_2.h>
#include <CGAL/Root_of_2.h>


CGAL_BEGIN_NAMESPACE

// This is a kernel wrapper which provides the type equality between
// Kernel::Point_2 and CGAL::Point_2<Kernel>, by deriving from
// K_base::Point_2 (and similar for the other types).

template < typename K_base, typename Kernel >
struct Curved_kernel_type_equality_wrapper
  : public K_base
{

    typedef K_base                                  Kernel_base;

 
    typedef CGAL::Circular_arc_2<Kernel>               Circular_arc_2;     
    typedef CGAL::Line_arc_2<Kernel>                   Line_arc_2;
    typedef CGAL::Circular_arc_point_2<Kernel>         Circular_arc_endpoint_2;
    typedef CGAL::Root_of_2<typename Kernel_base::FT>  Root_of_2;
    
    
    
	//Something has to be done with these 3, maybe a lazy Algebraic kernel?
	   
    //typedef Polynomial_for_circles_2_2<Kernel>   Polynomial_for_circles_2_2;
    //typedef Polynomial_1_2<Kernel>               Polynomial_1_2;
    //typedef Root_of_2<Kernel>                    Root_of_2;   

    



};

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_TYPE_EQUALITY_WRAPPER_H
