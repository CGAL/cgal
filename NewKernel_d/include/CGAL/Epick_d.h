#ifndef CGAL_EPICK_D_H
#define CGAL_EPICK_D_H
#include <CGAL/Kernel_d/Cartesian_base.h>
#include <CGAL/Kernel_d/Cartesian_static_filters.h>
#include <CGAL/Kernel_d/Cartesian_filter_K.h>
#include <CGAL/Kernel_d/Wrapper/Cartesian_wrap.h>
#include <CGAL/Kernel_d/Kernel_d_interface.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Interval_nt.h>


namespace CGAL {
#define CGAL_BASE \
  Cartesian_filter_K< \
    Cartesian_base_d<double, Dim>, \
    Cartesian_base_d<Interval_nt_advanced, Dim>, \
    Cartesian_base_d<Gmpq, Dim> \
  >
template<class Dim>
struct Epick_d_help1
: CGAL_BASE
{
  CGAL_CONSTEXPR Epick_d_help1(){}
  CGAL_CONSTEXPR Epick_d_help1(int d):CGAL_BASE(d){}
};
#undef CGAL_BASE
#define CGAL_BASE \
  Cartesian_static_filters<Dim,Epick_d_help1<Dim>,Epick_d_help2<Dim> >
template<class Dim>
struct Epick_d_help2
: CGAL_BASE
{
  CGAL_CONSTEXPR Epick_d_help2(){}
  CGAL_CONSTEXPR Epick_d_help2(int d):CGAL_BASE(d){}
};
#undef CGAL_BASE
#define CGAL_BASE \
  Kernel_d_interface< Cartesian_wrap< Epick_d_help2< \
    typename boost::conditional< d_==UNKNOWN_DIMENSION, \
                                 Dynamic_dimension_tag, \
				 Dimension_tag<d_> \
    >::type \
  > > >
template<int d_>
struct Epick_d
: CGAL_BASE
{
  CGAL_CONSTEXPR Epick_d(){}
  CGAL_CONSTEXPR Epick_d(int d):CGAL_BASE(d){}
};
#undef CGAL_BASE
}
#endif
