#ifndef CGAL_KERNEL_D_CARTESIAN_COMPLETE_H
#define CGAL_KERNEL_D_CARTESIAN_COMPLETE_H

#include <CGAL/Kernel_d/function_objects_cartesian.h>
#include <CGAL/Kernel_d/Cartesian_per_dimension.h>
#include <CGAL/Kernel_d/Types/Segment.h>
#include <CGAL/Kernel_d/Types/Sphere.h>
#include <CGAL/Kernel_d/Types/Hyperplane.h>
#include <CGAL/Kernel_d/Types/Aff_transformation.h>
#include <CGAL/Kernel_d/Types/Line.h>
#include <CGAL/Kernel_d/Types/Ray.h>
#include <CGAL/Kernel_d/Types/Iso_box.h>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/has_xxx.hpp>


namespace CGAL {
template<class R_,bool force_=false,class Derived_=Default> struct Cartesian_complete_types 
: public R_
{
  CGAL_CONSTEXPR Cartesian_complete_types(){}
  CGAL_CONSTEXPR Cartesian_complete_types(int d):R_(d){}
};

template<class R_,bool force_=false,class Derived_=Default> struct Cartesian_complete_constructors 
: public R_
{
  CGAL_CONSTEXPR Cartesian_complete_constructors(){}
  CGAL_CONSTEXPR Cartesian_complete_constructors(int d):R_(d){}
};

template<class R_,bool force_=false,class Derived_=Default> struct Cartesian_complete_predicates 
: public R_
{
  CGAL_CONSTEXPR Cartesian_complete_predicates(){}
  CGAL_CONSTEXPR Cartesian_complete_predicates(int d):R_(d){}
};

template<class R_,bool force_=false,class Derived_=Default> struct Cartesian_complete_computes 
: public R_
{
  CGAL_CONSTEXPR Cartesian_complete_computes(){}
  CGAL_CONSTEXPR Cartesian_complete_computes(int d):R_(d){}
};

}

#endif // CGAL_KERNEL_D_CARTESIAN_COMPLETE_H
