#ifndef CGAL_KERNEL_D_LAZY_CARTESIAN_H
#define CGAL_KERNEL_D_LAZY_CARTESIAN_H

#include <CGAL/basic.h>
#include <CGAL/algorithm.h>
#include <CGAL/Lazy.h>
#include <CGAL/Default.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/iterator_from_indices.h>
#include <CGAL/Kernel_d/Define_kernel_types.h>

namespace CGAL {

template<class K,class T>
struct Nth_iterator_element : private Store_kernel<K> {
  Nth_iterator_element(){}
  Nth_iterator_element(K const&k):Store_kernel<K>(k){}
  typedef typename Read_tag_type<K, typename iterator_tag_traits<T>::value_tag>::type result_type;
  template<class U> result_type operator()(CGAL_FORWARDABLE(U) u, int i) const {
    typename K::template Functor<Construct_ttag<T> >::type ci(this->kernel());
    return *cpp0x::next(ci(CGAL_FORWARD(U,u),Begin_tag()),i);
  }
};
      //typedef typename Functor<typename iterator_tag_traits<T>::nth_element>::type nth_elem;
template<class K, class T, bool = iterator_tag_traits<T>::has_nth_element>
struct Select_nth_element_functor {
  typedef Nth_iterator_element<K, T> type;
};
template<class K, class T>
struct Select_nth_element_functor <K, T, true> :
  K::template Functor<typename iterator_tag_traits<T>::nth_element> {};

namespace internal {
  template<class A,class B,class C,bool/*is_NT=false*/>
    struct Lazy_construction_maybe_nt {
      typedef Lazy_construction<A,B,C> type;
    };
  template<class A,class B,class C>
    struct Lazy_construction_maybe_nt<A,B,C,true> {
      typedef Lazy_construction_nt<A,B,C> type;
    };
}

template <class EK_, class AK_, class E2A_, class Kernel_>
struct Lazy_cartesian_types
{
    typedef typename typeset_intersection<
      typename AK_::Object_list,
      typename EK_::Object_list
	>::type Object_list;

    typedef typename typeset_intersection<
      typename AK_::Iterator_list,
      typename EK_::Iterator_list
	>::type Iterator_list;

    template <class T,class=void> struct Type {
	    typedef Lazy<
		    typename Read_tag_type<AK_,T>::type,
		    typename Read_tag_type<EK_,T>::type,
		    typename EK_::FT, E2A_> type;
    };
    template <class D> struct Type<FT_tag,D> {
      typedef CGAL::Lazy_exact_nt<typename EK_::FT>  type;
    };
    template <class D> struct Type<RT_tag,D> {
      typedef CGAL::Lazy_exact_nt<typename EK_::RT>  type;
    };

    template <class T> struct Iterator {
      typedef typename iterator_tag_traits<T>::value_tag Vt;
      typedef typename Type<Vt>::type V;
      typedef typename Select_nth_element_functor<AK_,T>::type AF;
      typedef typename Select_nth_element_functor<EK_,T>::type EF;

      typedef typename internal::Lazy_construction_maybe_nt<
	Kernel_, AF, EF, is_NT_tag<Vt>::value
	>::type nth_elem;

      typedef Iterator_from_indices<
	const typename Type<typename iterator_tag_traits<T>::container>::type,
	const V, V, nth_elem
      > type;
    };
};

template <class EK_, class AK_, class E2A_/*, class Kernel_=Default*/>
struct Lazy_cartesian : Dimension_base<typename EK_::Default_ambient_dimension>,
  Define_kernel_types<Lazy_cartesian_types<EK_,AK_,E2A_,Lazy_cartesian<EK_,AK_,E2A_> > >
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

    typedef Lazy_cartesian Self;
    typedef Define_kernel_types<Lazy_cartesian_types<EK_,AK_,E2A_,Self> > Base;
    //typedef typename Default::Get<Kernel_,Self>::type Kernel;
    typedef Self  Kernel;
    typedef AK_   Approximate_kernel;
    typedef EK_   Exact_kernel;
    typedef E2A_  E2A;
    typedef Approx_converter<Kernel, Approximate_kernel>   C2A;
    typedef Exact_converter<Kernel, Exact_kernel>    C2E;
    typedef CGAL::Lazy_exact_nt<typename Exact_kernel::FT>  FT;
    typedef CGAL::Lazy_exact_nt<typename Exact_kernel::RT>  RT;

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
	    typedef Construct_iter<typename Read_tag_type<Base,typename map_result_tag<T>::type>::type> type;
    };


    //TODO: what about other functors of the Misc category?
    // for Point_dimension, we should apply it to the approximate point
    // for printing, we should??? just not do printing this way?
};


} //namespace CGAL

#endif // CGAL_KERNEL_D_LAZY_CARTESIAN_H
