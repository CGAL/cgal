#ifndef CGAL_KERNEL_D_LAZY_CARTESIAN_H
#define CGAL_KERNEL_D_LAZY_CARTESIAN_H

#include <CGAL/basic.h>
#include <CGAL/Lazy.h>
#include <CGAL/Default.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/iterator_from_indices.h>

namespace CGAL {

template <class EK_, class AK_, class E2A_/*, class Kernel_=Default*/>
struct Lazy_cartesian : Dimension_base<typename EK_::Default_ambient_dimension>
{
    //CGAL_CONSTEXPR Lazy_cartesian(){}
    //CGAL_CONSTEXPR Lazy_cartesian(int d):Base_(d){}

    //TODO: Do we want to store an AK and an EK? Or just references?
    //FIXME: references would be better I guess.
    //TODO: In any case, make sure that we don't end up storing this kernel for
    //nothing (it is not empty but references empty kernels or something)
    AK_ ak; EK_ ek;
    AK_ const& approximate_kernel()const{return ak;}
    EK_ const& exact_kernel()const{return ek;}

    typedef Lazy_cartesian<EK_,AK_,E2A_/*,Kernel_*/> Self;
    //typedef typename Default::Get<Kernel_,Self>::type Kernel;
    typedef Self  Kernel;
    typedef AK_   Approximate_kernel;
    typedef EK_   Exact_kernel;
    typedef E2A_  E2A;
    typedef Approx_converter<Kernel, Approximate_kernel>   C2A;
    typedef Exact_converter<Kernel, Exact_kernel>    C2E;
    typedef CGAL::Lazy_exact_nt<typename Exact_kernel::FT>  FT;
    typedef FT RT;

    typedef typename Exact_kernel::Rep_tag Rep_tag;
    typedef typename Exact_kernel::Kernel_tag Kernel_tag;
    typedef typename Exact_kernel::Default_ambient_dimension Default_ambient_dimension;
    typedef typename Exact_kernel::Max_ambient_dimension Max_ambient_dimension;

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

    // Doesn't look like we need an explicit list.
    template <class T,class=void> struct Type {
	    typedef Lazy<
		    typename Approximate_kernel::template Type<T>::type,
		    typename Exact_kernel::template Type<T>::type,
		    typename Exact_kernel::FT, E2A> type;
    };
    template <class D> struct Type<FT_tag,D> {
      typedef FT type;
    };
    typedef typename typeset_intersection<
      typename Approximate_kernel::Object_list,
      typename Exact_kernel::Object_list
	>::type Object_list;

    template<class T,class D=void,class=typename map_functor_type<T>::type> struct Functor {
	    typedef Null_functor type;
    };
	    //FIXME: what do we do with D here?
    template<class T,class D> struct Functor<T,D,Predicate_tag> {
	    typedef typename Approximate_kernel::template Functor<T>::type FA;
	    typedef typename Exact_kernel::template Functor<T>::type FE;
	    typedef Filtered_predicate<FE,FA,C2E,C2A> type;
    };
    template<class T,class D> struct Functor<T,D,Compute_tag> {
	    typedef typename Approximate_kernel::template Functor<T>::type FA;
	    typedef typename Exact_kernel::template Functor<T>::type FE;
	    typedef Lazy_construction_nt<Kernel,FA,FE> type;
    };
    template<class T,class D> struct Functor<T,D,Construct_tag> {
	    typedef typename Approximate_kernel::template Functor<T>::type FA;
	    typedef typename Exact_kernel::template Functor<T>::type FE;
	    typedef Lazy_construction<Kernel,FA,FE> type;
    };

    typedef typename typeset_intersection<
      typename Approximate_kernel::Iterator_list,
      typename Exact_kernel::Iterator_list
	>::type Iterator_list;


    //TODO: handle the case without nth_element
#if 0
    template<class T>struct Default_nth_element : private Store_kernel<Kernel> {
      Default_nth_element(){}
      Default_nth_element(Kernel const&k):Store_kernel<Kernel>(k){}
      typedef /*???*/ result_type;
      template<class U> result_type operator()(CGAL_FORWARDABLE(U) u, int i) {
	typename /*???*/ ci(this->kernel());
	std::advance(ci, i);
	return *i;
      }
    };
#endif

    template <class T> struct Iterator {
      //WARNING: this fails because it is not lazy enough:
      //typedef typename Read_tag_type<Self,typename iterator_tag_traits<T>::value_tag>::type V;
      typedef typename Type<typename iterator_tag_traits<T>::value_tag>::type V;
      typedef Iterator_from_indices<
	const typename Type<typename iterator_tag_traits<T>::container>::type,
	const V, V,
	typename Functor<typename iterator_tag_traits<T>::nth_element>::type
      > type;
    };
    //typedef typename Iterator<Point_cartesian_const_iterator_tag>::type Point_cartesian_const_iterator;
    //typedef typename Iterator<Vector_cartesian_const_iterator_tag>::type Vector_cartesian_const_iterator;

    template<class U>
    struct Construct_iter : private Store_kernel<Kernel> {
	    Construct_iter(){}
	    Construct_iter(Kernel const&k):Store_kernel<Kernel>(k){}
	    //FIXME: pass the kernel to the functor in the iterator
	    typedef U result_type;
	    template<class T>
	    result_type operator()(T const& t,Begin_tag)const{
		    return result_type(t,0,this->kernel());
	    }
	    template<class T>
	    result_type operator()(T const& t,End_tag)const{
		    return result_type(t,Self().dimension(),this->kernel());
	    }
    };
    template<class T,class D> struct Functor<T,D,Construct_iterator_tag> {
	    typedef Construct_iter<typename Iterator<typename map_result_tag<T>::type>::type> type;
    };


    //TODO: what about other functors of the Misc category?
    // for Point_dimension, we should apply it to the approximate point
    // for printing, we should??? just not do printing this way?
};


} //namespace CGAL

#endif // CGAL_KERNEL_D_LAZY_CARTESIAN_H
