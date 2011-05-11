#ifndef CGAL_KERNEL_D_LAZY_CARTESIAN_H
#define CGAL_KERNEL_D_LAZY_CARTESIAN_H

#include <CGAL/basic.h>
#include <CGAL/Lazy.h>
#include <CGAL/Filtered_predicate.h>

namespace CGAL {

template <class EK_, class AK_, class E2A_, class Kernel>
struct Lazy_cartesian
{
    //CGAL_CONSTEXPR Lazy_cartesian(){}
    //CGAL_CONSTEXPR Lazy_cartesian(int d):Base_(d){}

    typedef AK_   Approximate_kernel;
    typedef EK_   Exact_kernel;
    typedef E2A_  E2A;
    typedef Approx_converter<Kernel, Approximate_kernel>   C2A;
    typedef Exact_converter<Kernel, Exact_kernel>    C2E;
    typedef CGAL::Lazy_exact_nt<typename Exact_kernel::FT>  FT;
    typedef FT RT;

    typedef typename Exact_kernel::Rep_tag Rep_tag;
    typedef typename Exact_kernel::Kernel_tag Kernel_tag;

    typedef typename Same_uncertainty_nt<bool, FT>::type
	    Boolean;
    typedef typename Same_uncertainty_nt<CGAL::Sign, FT>::type
	    Sign;
    typedef typename Same_uncertainty_nt<CGAL::Comparison_result, FT>::type
	    Comparison_result;
    typedef typename Same_uncertainty_nt<CGAL::Orientation, FT>::type
	    Orientation;
    typedef typename Same_uncertainty_nt<CGAL::Oriented_side, FT>::type
	    Oriented_side;
    typedef typename Same_uncertainty_nt<CGAL::Bounded_side, FT>::type
	    Bounded_side;
    typedef typename Same_uncertainty_nt<CGAL::Angle, FT>::type
	    Angle;

#define CGAL_Kernel_obj(X) \
    typedef Lazy<typename Approximate_kernel::X, typename Exact_kernel::X, typename Exact_kernel::FT, E2A>  X;

#include <CGAL/Kernel_d/interface_macros.h>

    template<class T,int i=0> struct Predicate {
	    typedef typename Approximate_kernel::template Predicate<T>::type FA;
	    typedef typename Exact_kernel::template Predicate<T>::type FE;
	    typedef Filtered_predicate<FE,FA,C2E,C2A> type;
    };
    template<class T,int i=0> struct Compute {
	    typedef typename Approximate_kernel::template Compute<T>::type FA;
	    typedef typename Exact_kernel::template Compute<T>::type FE;
	    typedef Lazy_construction_nt<Kernel,FA,FE> type;
    };
    template<class T,int i=0> struct Construct {
	    typedef typename Approximate_kernel::template Construct<T>::type FA;
	    typedef typename Exact_kernel::template Construct<T>::type FE;
	    typedef Lazy_construction<Kernel,FA,FE> type;
    };

    //TODO:
    //typedef ????????? Cartesian_const_iterator;
    //typedef ????????? Construct_cartesian_const_iterator
//
//    template<int i> struct Compute<Compute_cartesian_coordinate_tag,i> {
//	    typedef Compute_cartesian_coordinate type;
//    };
//    template<int i> struct Construct<Construct_cartesian_const_iterator_tag,i> {
//	    typedef Construct_cartesian_const_iterator type;
//    };
};


} //namespace CGAL

#endif // CGAL_KERNEL_D_LAZY_CARTESIAN_H
