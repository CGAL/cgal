// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Marc Glisse

#ifndef CGAL_KERNEL_D_LAZY_CARTESIAN_H
#define CGAL_KERNEL_D_LAZY_CARTESIAN_H

#include <CGAL/basic.h>
#include <CGAL/algorithm.h>
#include <CGAL/Lazy.h>
#include <CGAL/Default.h>
#include <CGAL/NewKernel_d/Filtered_predicate2.h>
#include <CGAL/iterator_from_indices.h>
#include <CGAL/NewKernel_d/Define_kernel_types.h>
#include <boost/function_output_iterator.hpp>

namespace CGAL {

template<class K,class T>
struct Nth_iterator_element : private Store_kernel<K> {
  Nth_iterator_element(){}
  Nth_iterator_element(K const&k):Store_kernel<K>(k){}
  typedef typename Get_type<K, typename iterator_tag_traits<T>::value_tag>::type result_type;
  template<class U> result_type operator()(U&& u, int i) const {
    typename Get_functor<K, Construct_ttag<T> >::type ci(this->kernel());
    return *std::next(ci(std::forward<U>(u),Begin_tag()),i);
  }
};
      //typedef typename Functor<typename iterator_tag_traits<T>::nth_element>::type nth_elem;
template<class K, class T, bool = iterator_tag_traits<T>::has_nth_element>
struct Select_nth_element_functor {
  typedef Nth_iterator_element<K, T> type;
};
template<class K, class T>
struct Select_nth_element_functor <K, T, true> :
  Get_functor<K, typename iterator_tag_traits<T>::nth_element> {};

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

// Whenever a construction takes iterator pairs as input, whether they point to double of Lazy objects, copy the ranges inside the lazy result so they are available for update_exact(). We analyze the input to try and guess where iterator pairs are. I would prefer if each functor had a specific signature (no overload in this layer) so we wouldn't have to guess.
namespace Lazy_internal {
template<class...>struct typelist{};
template<std::size_t>struct arg_i{};
template<std::size_t>struct arg_i_begin{};
template<std::size_t>struct arg_i_end{};
template<std::size_t>struct arg_i_ip1_range{};
template<class,class,class,class=void>struct analyze_args;
template<class T,class U>struct analyze_args<T,U,typelist<>> {
  typedef T creator;
  typedef U reader;
};
template<class...T,class...U,class V,class...W>
struct analyze_args<typelist<T...>,typelist<U...>,typelist<V,W...>,std::enable_if_t<!is_iterator_type<V,std::input_iterator_tag>::value>> :
analyze_args<typelist<T...,arg_i<sizeof...(U)>>,typelist<U...,arg_i<sizeof...(T)>>,typelist<W...>> {};
template<class...T,class...U,class It,class...W>
struct analyze_args<typelist<T...>,typelist<U...>,typelist<It,It,W...>,std::enable_if_t<is_iterator_type<It,std::input_iterator_tag>::value>> :
analyze_args<typelist<T...,arg_i_ip1_range<sizeof...(U)>>,typelist<U...,arg_i_begin<sizeof...(T)>,arg_i_end<sizeof...(T)>>,typelist<W...>> {};
template<class...T> using analyze_args_for_lazy = analyze_args<typelist<>,typelist<>,typelist<T...>>;
template<class,class>struct extract1;
template<std::size_t i,class T>struct extract1<arg_i<i>,T>:std::tuple_element<i,T>{};
template<std::size_t i,class T>struct extract1<arg_i_ip1_range<i>,T>{
  typedef std::tuple_element_t<i,T> E;
  typedef std::remove_cv_t<std::remove_reference_t<E>> It;
  typedef typename std::iterator_traits<It>::value_type element_type;
  // TODO: find a way to use an array of the right size, at least for the most frequent constructions
  typedef std::vector<element_type> type;
};
template<std::size_t i,class...T>decltype(auto)
do_extract(arg_i<i>,std::tuple<T...>const&t)
{return std::get<i>(t);}
template<std::size_t i,class...T>decltype(auto)
do_extract(arg_i_begin<i>,std::tuple<T...>const&t)
{return std::begin(std::get<i>(t));}
template<std::size_t i,class...T>decltype(auto)
do_extract(arg_i_end<i>,std::tuple<T...>const&t)
{return std::end(std::get<i>(t));}
template<std::size_t i,class...T>decltype(auto)
do_extract(arg_i_ip1_range<i>,std::tuple<T...>const&t)
{
  typedef std::tuple<T...> L;
  typedef std::tuple_element_t<i,L> E;
  typedef std::remove_cv_t<std::remove_reference_t<E>> It;
  typedef typename std::iterator_traits<It>::value_type element_type;
  typedef std::vector<element_type> type;
  return type(std::get<i>(t),std::get<i+1>(t));
}
template<class,class>struct data_from_input;
template<class...T,class U>struct data_from_input<typelist<T...>,U> {
  typedef std::tuple<typename extract1<T,U>::type...> type;
};
}
template<typename AT, typename ET, typename AC, typename EC, typename E2A, typename...L>
class Lazy_rep_XXX :
  public Lazy_rep< AT, ET, E2A >, private EC
{
  // `default_construct<T>()` is the same as `T{}`. But, this is a
  // workaround to a MSVC-2015 bug (fixed in MSVC-2017): its parser
  // seemed confused by `T{}` somewhere below.
  template <typename T>
  static T default_construct() { return T(); }

  // Lazy_rep_0 does not inherit from EC or take a parameter AC. It has different constructors.
  static_assert(sizeof...(L)>0, "Use Lazy_rep_0 instead");
  template <class Ei, class Ai, class E2Ai, class Ki> friend class Lazy_kernel_base;
  typedef Lazy_internal::analyze_args_for_lazy<L...> Args;
  // How to go from l to Lazy_rep's data
  typedef typename Args::creator Creator;
  // How to go back
  typedef typename Args::reader Reader;
  // what Lazy_rep should store
  typedef typename Lazy_internal::data_from_input<Creator, std::tuple<L...>>::type LL;
  mutable LL l; // L...l; is not yet allowed.
  const EC& ec() const { return *this; }
  template<class...T>
  void update_exact_helper(Lazy_internal::typelist<T...>) const {
    this->et = new ET(ec()( CGAL::exact( Lazy_internal::do_extract(T{},l) ) ... ) );
    this->at = E2A()(*(this->et));
    l = LL(); // There should be a nicer way to clear. Destruction for instance. With this->et as a witness of whether l has already been destructed.
  }
  public:
  void update_exact() const {
    update_exact_helper(Reader{});
  }
  template<class...LL>
  Lazy_rep_XXX(const AC& ac, const EC& ec, LL const&...ll) :
    Lazy_rep_XXX(Creator{},ac,ec,std::forward_as_tuple(ll...),ll...){};
  private:
  // Currently we construct the vectors, then move them into the tuple. It would be nicer to construct them in their final destination, because eventually we will also have arrays instead of vectors.
  template<class...T,class LLL,class...LL>
  Lazy_rep_XXX(Lazy_internal::typelist<T...>, const AC& ac, const EC& ec, LLL const&lll, LL const&...ll) :
    Lazy_rep<AT, ET, E2A>(ac(CGAL::approx(ll)...)), EC(ec), l(Lazy_internal::do_extract(default_construct<T>(),lll)...)
  {
    //this->set_depth(std::max({ -1, (int)CGAL::depth(ll)...}) + 1);
    this->set_depth(1); // FIXME: now that we have ranges, we could actually compute the depth if we cared...
  }
  // TODO: print_dag needs a specific implementation for Lazy_rep_XXX
};
template<class Tag, class LK, class ET>struct Select_converter { typedef typename LK::E2A type; };
template<class LK, class ET>struct Select_converter<Compute_tag,LK,ET> { typedef To_interval<ET> type; };
template<typename T, typename LK>
struct Lazy_construction2 {
  static const bool Protection = true;

  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename Get_functor<AK, T>::type AC;
  typedef typename Get_functor<EK, T>::type EC;
  typedef typename map_result_tag<T>::type result_tag;
  typedef typename Get_type<AK, result_tag>::type AT;
  typedef typename Get_type<EK, result_tag>::type ET;
  typedef typename Get_type<LK, result_tag>::type result_type;
  // same as Handle = Lazy< AT, ET, E2A>
  typedef typename Select_converter<typename Get_functor_category<LK,T>::type, LK, ET>::type E2A;

  Lazy_construction2(){}
  Lazy_construction2(LK const&k):ac(k.approximate_kernel()),ec(k.exact_kernel()){}
  CGAL_NO_UNIQUE_ADDRESS AC ac;
  CGAL_NO_UNIQUE_ADDRESS EC ec;

  template<class...L>
  std::enable_if_t<(sizeof...(L)>0), result_type> operator()(L const&...l) const {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      return new Lazy_rep_XXX<AT, ET, AC, EC, E2A, L...>(ac, ec, l...);
    } catch (Uncertain_conversion_exception&) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      return new Lazy_rep_0<AT,ET,E2A>(ec(CGAL::exact(l)...));
    }
  }
  // FIXME: this forces us to have default constructors for all types, try to make its instantiation lazier
  // Actually, that may be the clearing in update_exact().
  result_type operator()() const
  {
    return new Lazy_rep_0<AT,ET,E2A>();
  }
};

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

    template <class T,class=typename Get_type_category<Kernel_,T>::type> struct Type {};
    template <class T> struct Type<T,Object_tag> {
            typedef Lazy<
                    typename Get_type<AK_,T>::type,
                    typename Get_type<EK_,T>::type,
                    E2A_> type;
    };
    template <class T> struct Type<T,Number_tag> {
      typedef CGAL::Lazy_exact_nt<typename Get_type<EK_,T>::type>  type;
    };

    template <class T> struct Iterator {
      typedef typename iterator_tag_traits<T>::value_tag Vt;
      typedef typename Type<Vt>::type V;
      typedef typename Select_nth_element_functor<AK_,T>::type AF;
      typedef typename Select_nth_element_functor<EK_,T>::type EF;

      // TODO: we should use Lazy_construction2, but this seems ok for now, we never construct iterators from iterators.
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
struct Lazy_cartesian :
  Lazy_cartesian_types<EK_,AK_,E2A_,Lazy_cartesian<EK_,AK_,E2A_> >
{
    constexpr Lazy_cartesian(){}
    constexpr Lazy_cartesian(int d):ak(d),ek(d){}

    //TODO: Do we want to store an AK and an EK? Or just references?
    //FIXME: references would be better I guess.
    //TODO: In any case, make sure that we don't end up storing this kernel for
    //nothing (it is not empty but references empty kernels or something)
    CGAL_NO_UNIQUE_ADDRESS AK_ ak;
    CGAL_NO_UNIQUE_ADDRESS EK_ ek;
    AK_ const& approximate_kernel()const{return ak;}
    EK_ const& exact_kernel()const{return ek;}

    int dimension()const{return ak.dimension();}
    void set_dimension(int dim){ak.set_dimension(dim);ek.set_dimension(dim);}

    // For compilers that do not handle [[no_unique_address]]
    typedef boost::mpl::and_<
      internal::Do_not_store_kernel<AK_>,
      internal::Do_not_store_kernel<EK_> > Do_not_store_kernel;

    typedef typename EK_::Dimension Dimension; // ?
    typedef Lazy_cartesian Self;
    typedef Lazy_cartesian_types<EK_,AK_,E2A_,Self> Base;
    //typedef typename Default::Get<Kernel_,Self>::type Kernel;
    typedef Self  Kernel;
    typedef AK_   Approximate_kernel;
    typedef EK_   Exact_kernel;
    typedef E2A_  E2A;
    //typedef Approx_converter<Kernel, Approximate_kernel>   C2A;
    //typedef Exact_converter<Kernel, Exact_kernel>    C2E;
    struct C2A {
      C2A(){}
      C2A(Kernel const&, Approximate_kernel const&){}
      template<class T>decltype(auto)operator()(T const&t)const{return CGAL::approx(t);}
    };
    struct C2E {
      C2E(){}
      C2E(Kernel const&, Exact_kernel const&){}
      template<class T>decltype(auto)operator()(T const&t)const{return  CGAL::exact(t);}
    };

    typedef typename Exact_kernel::Rep_tag Rep_tag;
    typedef typename Exact_kernel::Kernel_tag Kernel_tag;
    typedef typename Exact_kernel::Default_ambient_dimension Default_ambient_dimension;
    typedef typename Exact_kernel::Max_ambient_dimension Max_ambient_dimension;
    //typedef typename Exact_kernel::Flat_orientation Flat_orientation;
    // Check that Approximate_kernel agrees with all that...

    template<class T,class D=void,class=typename Get_functor_category<Lazy_cartesian,T,D>::type> struct Functor {
            typedef Null_functor type;
    };
            //FIXME: what do we do with D here?
    template<class T,class D> struct Functor<T,D,Predicate_tag> {
            typedef typename Get_functor<Approximate_kernel, T>::type FA;
            typedef typename Get_functor<Exact_kernel, T>::type FE;
            typedef Filtered_predicate2<Lazy_cartesian,FE,FA,C2E,C2A> type;
    };
    template<class T,class D> struct Functor<T,D,Compute_tag> {
            typedef Lazy_construction2<T,Kernel> type;
    };
    template<class T,class D> struct Functor<T,D,Construct_tag> {
            typedef Lazy_construction2<T,Kernel> type;
    };
    template<class D> struct Functor<Point_dimension_tag,D,Misc_tag> {
            typedef typename Get_functor<Approximate_kernel, Point_dimension_tag>::type FA;
            struct type {
              FA fa;
              type(){}
              type(Kernel const&k):fa(k.approximate_kernel()){}
              template<class P>
              int operator()(P const&p)const{return fa(CGAL::approx(p));}
            };
    };
    template<class D> struct Functor<Vector_dimension_tag,D,Misc_tag> {
            typedef typename Get_functor<Approximate_kernel, Vector_dimension_tag>::type FA;
            struct type {
              FA fa;
              type(){}
              type(Kernel const&k):fa(k.approximate_kernel()){}
              template<class V>
              int operator()(V const&v)const{return fa(CGAL::approx(v));}
            };
    };
    template<class D> struct Functor<Linear_base_tag,D,Misc_tag> {
      // Don't filter that one, as there is no guarantee that the interval
      // basis would be in any way related to the exact basis, the most obvious
      // difference being the order of the vectors.
      // Don't try to be generic until we have more than just this one.
      typedef typename Get_functor<Exact_kernel, Linear_base_tag>::type FE;
      typedef typename Get_type<Approximate_kernel, Vector_tag>::type AT;
      typedef typename Get_type<Exact_kernel, Vector_tag>::type ET;
      typedef typename Base::template Type<Vector_tag>::type V; // Lazy<AT, ET, E2A>
      struct type {
        FE fe;
        type(){}
        type(Kernel const&k):fe(k.exact_kernel()){}
        template<class Iter, class Oter>
        void operator()(Iter i, Iter e, Oter o)const{
          fe(CGAL::exact(i), CGAL::exact(e),
              boost::make_function_output_iterator(
                [&o](ET const&v){
                  *o++ = V(new Lazy_rep_0<AT,ET,E2A>(v));
                }
              )
          );
        }
      };
    };

    typedef typename Base::template Iterator<Point_cartesian_const_iterator_tag>::type Point_cartesian_const_iterator;
    typedef typename Base::template Iterator<Vector_cartesian_const_iterator_tag>::type Vector_cartesian_const_iterator;

    // This is really specific to point/vector coordinate iterators
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
                    typedef typename Get_functor<Approximate_kernel, Point_dimension_tag>::type PD;
                    return result_type(t,PD(this->kernel().approximate_kernel())(CGAL::approx(t)),this->kernel());
            }
    };
    template<class T,class D> struct Functor<T,D,Construct_iterator_tag> {
            typedef Construct_iter<typename Base::template Iterator<typename map_result_tag<T>::type>::type> type;
    };


    //TODO: what about other functors of the Misc category?
    // for Point_dimension, we should apply it to the approximate point
    // for printing, we should??? just not do printing this way?
};


} //namespace CGAL

#endif // CGAL_KERNEL_D_LAZY_CARTESIAN_H
