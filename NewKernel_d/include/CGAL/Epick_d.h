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
template<int d_>
struct Epick_d_help1
: Cartesian_filter_K<
    Cartesian_base_d<double, Dimension_tag<d_> >,
    Cartesian_base_d<Interval_nt_advanced, Dimension_tag<d_> >,
    Cartesian_base_d<Gmpq, Dimension_tag<d_> >
  >
{};
template<int d_>
struct Epick_d_help2
: Cartesian_static_filters<Dimension_tag<d_>,Epick_d_help1<d_>,Epick_d_help2<d_> >
{};
template<int d_>
struct Epick_d
: Kernel_d_interface<Cartesian_wrap<Epick_d_help2<d_> > >
{};
}
#endif
