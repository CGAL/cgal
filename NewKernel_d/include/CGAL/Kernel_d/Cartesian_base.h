#ifndef CGAL_KERNEL_D_CARTESIAN_BASE_H
#define CGAL_KERNEL_D_CARTESIAN_BASE_H

#include <CGAL/basic.h>
#include <CGAL/Kernel_d/Cartesian_complete.h>
#include <CGAL/Kernel_d/Cartesian_LA_base.h>

namespace CGAL {

template < typename FT_, typename Dim_>
struct Cartesian_base_d : public
Cartesian_complete_predicates<
Cartesian_complete_constructors<
Cartesian_complete_computes<
Cartesian_complete_types<
	Cartesian_LA_base_d<FT_,Dim_>
>, false, Cartesian_base_d<FT_,Dim_>
>, false, Cartesian_base_d<FT_,Dim_>
>, false, Cartesian_base_d<FT_,Dim_>
>
{
    CGAL_CONSTEXPR Cartesian_base_d(){}
    CGAL_CONSTEXPR Cartesian_base_d(int d):Dimension_base<Dim_>(d){}
};

} //namespace CGAL

#endif // CGAL_KERNEL_D_CARTESIAN_BASE_H
