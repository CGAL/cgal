#ifndef CGAL_KERNEL_D_CARTESIAN_FILTER_K_H
#define CGAL_KERNEL_D_CARTESIAN_FILTER_K_H

#include <CGAL/basic.h>
#include <CGAL/Kernel_d/KernelD_converter.h>
#include <CGAL/Filtered_predicate.h>
#include <boost/mpl/if.hpp>

namespace CGAL {

template < typename Base_, typename AK_, typename EK_ >
struct Cartesian_filter_K : public Base_,
	private Store_kernel<AK_>, private Store_kernel2<EK_>
{
    CGAL_CONSTEXPR Cartesian_filter_K(){}
    CGAL_CONSTEXPR Cartesian_filter_K(int d):Base_(d){}
    //FIXME: or do we want an instance of AK and EK belonging to this kernel,
    //instead of a reference to external ones?
    CGAL_CONSTEXPR Cartesian_filter_K(AK_ const&a,EK_ const&b):Base_(),Store_kernel<AK_>(a),Store_kernel2<EK_>(b){}
    CGAL_CONSTEXPR Cartesian_filter_K(int d,AK_ const&a,EK_ const&b):Base_(d),Store_kernel<AK_>(a),Store_kernel2<EK_>(b){}
    typedef Base_ Kernel_base;
    typedef AK_ AK;
    typedef EK_ EK;
    typedef typename Store_kernel<AK_>::reference_type AK_rt;
    AK_rt approximate_kernel()const{return this->kernel();}
    typedef typename Store_kernel2<EK_>::reference2_type EK_rt;
    EK_rt exact_kernel()const{return this->kernel2();}

    //TODO: C2A/C2E could be able to convert *this into this->kernel() or this->kernel2().
    typedef KernelD_converter<Kernel_base,AK> C2A;
    typedef KernelD_converter<Kernel_base,EK> C2E;

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
