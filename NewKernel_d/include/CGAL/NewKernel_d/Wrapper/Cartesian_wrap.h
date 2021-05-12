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

#ifndef CGAL_KERNEL_D_CARTESIAN_WRAP_H
#define CGAL_KERNEL_D_CARTESIAN_WRAP_H

#include <CGAL/basic.h>
#include <CGAL/is_iterator.h>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4003) // not enough actual parameters for macro 'BOOST_PP_EXPAND_I'
                                // https://lists.boost.org/boost-users/2014/11/83291.php
#endif
#include <CGAL/NewKernel_d/Wrapper/Point_d.h>
#include <CGAL/NewKernel_d/Wrapper/Vector_d.h>
#include <CGAL/NewKernel_d/Wrapper/Segment_d.h>
#include <CGAL/NewKernel_d/Wrapper/Sphere_d.h>
#include <CGAL/NewKernel_d/Wrapper/Hyperplane_d.h>
#include <CGAL/NewKernel_d/Wrapper/Weighted_point_d.h>

#include <CGAL/NewKernel_d/Wrapper/Ref_count_obj.h>

#include <boost/mpl/or.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/mpl/vector.hpp>

//TODO: do we want to store the kernel ref in the Object wrappers? It would allow for additions and operator[] and things like that to work, but objects would still need to be created by functors.

namespace CGAL {
namespace internal {
BOOST_MPL_HAS_XXX_TRAIT_DEF(Is_wrapper)
template<class T,bool=has_Is_wrapper<T>::value> struct Is_wrapper {
        enum { value=false };
        typedef Tag_false type;
};
template<class T> struct Is_wrapper<T,true> {
        typedef typename T::Is_wrapper type;
        enum { value=type::value };
};

template<class T,bool=is_iterator_type<T,std::input_iterator_tag>::value> struct Is_wrapper_iterator {
        enum { value=false };
        typedef Tag_false type;
};
template<class T> struct Is_wrapper_iterator<T,true> :
        Is_wrapper<typename std::iterator_traits<typename CGAL::decay<T>::type>::value_type>
{ };

struct Forward_rep {
//TODO: make a good C++0X version with perfect forwarding
//#ifdef CGAL_CXX11
//template <class T,class=typename std::enable_if<!Is_wrapper<typename std::decay<T>::type>::value&&!Is_wrapper_iterator<typename std::decay<T>::type>::value>::type>
//T&& operator()(typename std::remove_reference<T>::type&& t) const {return static_cast<T&&>(t);};
//template <class T,class=typename std::enable_if<!Is_wrapper<typename std::decay<T>::type>::value&&!Is_wrapper_iterator<typename std::decay<T>::type>::value>::type>
//T&& operator()(typename std::remove_reference<T>::type& t) const {return static_cast<T&&>(t);};
//
//template <class T,class=typename std::enable_if<Is_wrapper<typename std::decay<T>::type>::value>::type>
//typename Type_copy_cvref<T,typename std::decay<T>::type::Rep>::type&&
//operator()(T&& t) const {
//        return static_cast<typename Type_copy_cvref<T,typename std::decay<T>::type::Rep>::type&&>(t.rep());
//};
//
//template <class T,class=typename std::enable_if<Is_wrapper_iterator<typename std::decay<T>::type>::value>::type>
//transforming_iterator<Forward_rep,typename std::decay<T>::type>
//operator()(T&& t) const {
//        return make_transforming_iterator(std::forward<T>(t),Forward_rep());
//};
//#else
template <class T,bool=Is_wrapper<T>::value,bool=Is_wrapper_iterator<T>::value> struct result_;
template <class T> struct result_<T,false,false>{typedef T const& type;};
template <class T> struct result_<T,true,false>{typedef typename decay<T>::type::Rep const& type;};
template <class T> struct result_<T,false,true>{typedef transforming_iterator<Forward_rep,typename decay<T>::type> type;};
template<class> struct result;
template<class T> struct result<Forward_rep(T)> : result_<T> {};

template <class T> typename boost::disable_if<boost::mpl::or_<Is_wrapper<T>,Is_wrapper_iterator<T> >,T>::type const& operator()(T const& t) const {return t;}
template <class T> typename boost::disable_if<boost::mpl::or_<Is_wrapper<T>,Is_wrapper_iterator<T> >,T>::type& operator()(T& t) const {return t;}

template <class T> typename T::Rep const& operator()(T const& t, typename boost::enable_if<Is_wrapper<T> >::type* = 0) const {return t.rep();}

template <class T> transforming_iterator<Forward_rep,typename boost::enable_if<Is_wrapper_iterator<T>,T>::type> operator()(T const& t) const {return make_transforming_iterator(t,Forward_rep());}
//#endif
};
}

template <class B, class K, class T, bool = Provides_type<B, T>::value>
struct Map_wrapping_type : Get_type<B, T> {};
#define CGAL_REGISTER_OBJECT_WRAPPER(X) \
  template <class B, class K> \
  struct Map_wrapping_type <B, K, X##_tag, true> { \
    typedef Wrap::X##_d<K> type; \
  }
CGAL_REGISTER_OBJECT_WRAPPER(Point);
CGAL_REGISTER_OBJECT_WRAPPER(Vector);
CGAL_REGISTER_OBJECT_WRAPPER(Segment);
CGAL_REGISTER_OBJECT_WRAPPER(Sphere);
CGAL_REGISTER_OBJECT_WRAPPER(Hyperplane);
CGAL_REGISTER_OBJECT_WRAPPER(Weighted_point);
#undef CGAL_REGISTER_OBJECT_WRAPPER

// Note: this tends to be an all or nothing thing currently, wrapping
// only some types breaks, probably because we don't check whether the
// return type is indeed wrapped.
template < typename Base_ , typename Derived_ = Default >
struct Cartesian_wrap : public Base_
{
    constexpr Cartesian_wrap(){}
    constexpr Cartesian_wrap(int d):Base_(d){}
    typedef Base_ Kernel_base;
    typedef Cartesian_wrap Self;
    // TODO: pass the 2 types Self and Derived to the wrappers, they can use Self for most purposes and Derived only for Kernel_traits' typedef R.
    typedef typename Default::Get<Derived_, Self>::type Derived;
    // FIXME: The list doesn't belong here.
    typedef boost::mpl::vector<Point_tag,Segment_tag,Sphere_tag,Vector_tag,Hyperplane_tag> Wrapped_list;

    template <class T>
    struct Type : Map_wrapping_type<Base_, Derived, T> {};

    //Translate the arguments
    template <class T, class D = void,
      class=typename Get_functor_category<Derived,T>::type,
      bool=Provides_functor<Kernel_base, T>::value,
      bool=boost::mpl::contains<Wrapped_list,typename map_result_tag<T>::type>::type::value>
    struct Functor {
            typedef typename Get_functor<Kernel_base, T>::type B;
            struct type {
                    B b;
                    type(){}
                    type(Self const&k):b(k){}
                    template<class...U> decltype(auto) operator()(U&&...u)const{
                            return b(internal::Forward_rep()(u)...);
                    }
            };
    };

    // Preserve the difference between Null_functor and nothing.
    template <class T, class D, class C, bool b>
    struct Functor <T, D, C, false, b>
      : Get_functor <Kernel_base, T> {};

    //Translate both the arguments and the result
    //TODO: Check Is_wrapper instead of relying on map_result_tag?
    template<class T,class D> struct Functor<T,D,Construct_tag,true,true> {
            typedef typename Get_functor<Kernel_base, T>::type B;
            struct type {
                    B b;
                    type(){}
                    type(Self const&k):b(k){}
                    typedef typename map_result_tag<T>::type result_tag;
                    // FIXME: Self or Derived?
                    typedef typename Get_type<Self,result_tag>::type result_type;
                    template<class...U> result_type operator()(U&&...u)const{
                            return result_type(Eval_functor(),b,internal::Forward_rep()(u)...);
                    }
            };
    };

};

template < typename Base_ >
struct Cartesian_refcount : public Base_
{
    constexpr Cartesian_refcount(){}
    constexpr Cartesian_refcount(int d):Base_(d){}
    typedef Base_ Kernel_base;
    typedef Cartesian_refcount Self;

    // FIXME: Use object_list, or a list passed as argument, or anything
    // automatic.
    template <class T, class=void> struct Type : Get_type<Base_, T> {};
#define CGAL_Kernel_obj(X,Y) \
    template <class D> struct Type<X##_tag, D> { typedef Ref_count_obj<Cartesian_refcount, X##_tag> type; };

    CGAL_Kernel_obj(Point,point)
    CGAL_Kernel_obj(Vector,vector)
#undef CGAL_Kernel_obj

    template<class T> struct Dispatch {
            //typedef typename map_functor_type<T>::type f_t;
            typedef typename map_result_tag<T>::type r_t;
            enum {
                    is_nul = boost::is_same<typename Get_functor<Kernel_base, T>::type,Null_functor>::value,
                    ret_rcobj = boost::is_same<r_t,Point_tag>::value || boost::is_same<r_t,Vector_tag>::value
            };
    };

    //Translate the arguments
    template<class T,class D=void,bool=Dispatch<T>::is_nul,bool=Dispatch<T>::ret_rcobj> struct Functor {
            typedef typename Get_functor<Kernel_base, T>::type B;
            struct type {
                    B b;
                    type(){}
                    type(Self const&k):b(k){}
                    typedef typename B::result_type result_type;
                    template<class...U> result_type operator()(U&&...u)const{
                            return b(internal::Forward_rep()(u)...);
                    }
            };
    };

    //Translate both the arguments and the result
    template<class T,class D,bool b> struct Functor<T,D,true,b> {
            typedef Null_functor type;
    };

    template<class T,class D> struct Functor<T,D,false,true> {
            typedef typename Get_functor<Kernel_base, T>::type B;
            struct type {
                    B b;
                    type(){}
                    type(Self const&k):b(k){}
                    typedef typename map_result_tag<T>::type result_tag;
                    typedef typename Get_type<Self,result_tag>::type result_type;
                    template<class...U> result_type operator()(U&&...u)const{
                            return result_type(Eval_functor(),b,internal::Forward_rep()(u)...);
                    }
            };
    };

};

} //namespace CGAL

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_KERNEL_D_CARTESIAN_WRAP_H
