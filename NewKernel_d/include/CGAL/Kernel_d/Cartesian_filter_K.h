#ifndef CGAL_KERNEL_D_CARTESIAN_FILTER_K_H
#define CGAL_KERNEL_D_CARTESIAN_FILTER_K_H

#include <CGAL/basic.h>
#include <CGAL/Kernel_d/Cartesian_converter.h>
#include <CGAL/Filtered_predicate.h>

namespace CGAL {

template < typename Base_, typename AK_, typename EK_ >
struct Cartesian_filter_K : public Base_
{
    CGAL_CONSTEXPR Cartesian_filter_K(){}
    CGAL_CONSTEXPR Cartesian_filter_K(int d):Base_(d){}
    typedef Base_ Kernel_base;
    typedef AK_ AK;
    typedef EK_ EK;
    typedef CartesianD_converter<Kernel_base,AK> C2A;
    typedef CartesianD_converter<Kernel_base,EK> C2E;

    template<class T,class D=void,class=typename map_functor_type<T>::type> struct Functor :
	    Kernel_base::template Functor<T,D> {};
    template<class T,class D> struct Functor<T,D,Predicate_tag> {
	    typedef typename AK::template Functor<T>::type AP;
	    typedef typename EK::template Functor<T>::type EP;
	    typedef Filtered_predicate<EP,AP,C2E,C2A> type;
    };
// TODO:
//    template<class T> struct Functor<T,No_filter_tag,Predicate_tag> :
//	    Kernel_base::template Functor<T,No_filter_tag> {};
// TODO:
// detect when Less_cartesian_coordinate doesn't need filtering
};

} //namespace CGAL

#endif // CGAL_KERNEL_D_CARTESIAN_FILTER_K_H
