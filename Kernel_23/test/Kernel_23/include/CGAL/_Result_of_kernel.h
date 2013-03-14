// Copyright (c) 2011 GeometryFactory (France). All rights reserved.
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Philipp Moeller

#ifndef CGAL_RESULT_OF_KERNEL_H
#define CGAL_RESULT_OF_KERNEL_H

#include <CGAL/config.h>

#if !defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES) && !defined(CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE) && !defined(CGAL_CFG_NO_CPP0X_STATIC_ASSERT)

#define CGAL_RESULT_OF_KERNEL 1

#include <type_traits>

#include <boost/mpl/has_xxx.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/type_traits/is_scalar.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/result_of.hpp>

#include <CGAL/Cartesian/Cartesian_base.h>
#include <CGAL/Homogeneous/Homogeneous_base.h>

#include <CGAL/Kernel/Type_equality_wrapper.h>

#include <boost/tr1/functional.hpp>


// Unfortunately this test is only an approximation. At this point
// cpp11::result_of behaves as C++11 result_of and we cannot force it
// back into an old mode thanks to the include guards. We instead use
// a TR1 implementation of result_of to compare the types. This is the
// best we can go for.

namespace CGAL {
  // avoid crashes with what we already have in namespace internal
  namespace result_of_kernel {
    // required for smooth wrapping of the functors
    BOOST_MPL_HAS_XXX_TRAIT_DEF(result_type)
    BOOST_MPL_HAS_XXX_TRAIT_DEF(Rep)

    template<typename A, typename B>
    struct Lazy_is_same {
      typedef boost::is_same<A, typename B::Rep> type;
    };

    // trickery to get rid of the inequality that appears with the
    // return_base_tag versions of construct calls
    template<typename A, typename B, bool t = boost::mpl::or_< 
                                       boost::is_same<A, B>, 
                                       typename
                                       boost::mpl::eval_if< has_Rep<B>,
                                                            Lazy_is_same<A, B>,
                                                            boost::false_type >::type
                                       >::value >
    struct Rep_equal;
    
    template<typename A, typename B>
    struct Rep_equal<A, B, true> : boost::true_type {};
    template<typename A, typename B>
    struct Rep_equal<A, B, false> : boost::false_type {};
  }

  // This functor can wrap any DefaultConstructible functor. Iff there
  // is a result_type typedef it needs to be forwarded. In all other
  // cases cpp11::result_of is necessary to determine the return type.
  template<typename F, bool result_type = result_of_kernel::has_result_type<F>::value >
  struct AnyFunctor;

  template<typename F>
  struct AnyFunctor<F, true> {
    typedef typename F::result_type result_type;

    template<typename... Args>
    auto operator()(Args&&... args) const -> typename std::result_of<F(Args...)>::type {
      F f;
      // check the equality of a c++03 std::tr1::result_of and a c++11 result_of
      typedef typename std::result_of<F(Args...)>::type c11_return_type;
      typedef typename F::result_type c03_return_type;

      static_assert((result_of_kernel::Rep_equal<c11_return_type, c03_return_type>::value), 
                    "Type difference between actual return type and cpp11::result_of<>::type");
      
      return f(std::forward<Args>(args)...);
    }
  };

  template<typename F>
  struct AnyFunctor<F, false> {
    template<typename>
    struct result;
    
    template<typename Func, typename... Args>
    struct result<Func(Args...)> {
      typedef typename cpp11::result_of<F(Args...)>::type type;
    };

    // same as above
    template<typename... Args>
    auto operator()(Args&&... args) const -> typename std::result_of<F(Args...)>::type {
      F f;
      typedef typename std::result_of<F(Args...)>::type c11_return_type;
      typedef typename std::tr1::result_of<F( 
        typename
        std::remove_cv< 
          typename
          std::remove_reference< 
            Args&&
            >::type
          >::type ...
        )>::type c03_return_type;

      static_assert((result_of_kernel::Rep_equal<c11_return_type, c03_return_type>::value), 
                    "Type difference between actual return type and cpp11::result_of<>::type");

      return f(std::forward<Args>(args)...);
    }
  };

  // usual copy pasta from simple_cartesian and cartesian plus an
  // template template parameter to inject the true base
  template < typename RT_, typename FT_, typename Kernel_, template<typename, typename, typename> class Inject>
  struct Result_of_base
    : public Inject< RT_, FT_, Kernel_>
  {
    typedef RT_                                           RT;
    typedef FT_                                           FT;

    // The mechanism that allows to specify reference-counting or not.
    template < typename T >
    struct Handle { typedef T   type; };

    template < typename Kernel2 >
    struct Base { typedef Result_of_base<RT_, FT_, Kernel2, Inject>  Type; };
  };

  // alias cartesian to something similar to Homogeneous. Requires gcc 4.7.
  // NB. Kernel is either last of first argument to add a little more confusion.

  // template<typename FT, typename, typename K>
  // using Cartesian_base_3 = Cartesian_base<K, FT>;

  // hack instead
  template<typename FT, typename, typename K>
  struct Cartesian_base_3 : public Cartesian_base<K, FT> { };
  
  template < typename FT_, typename Kernel_ >
  struct Result_of_cartesian_base : public Result_of_base<FT_, FT_, Kernel_, Cartesian_base_3>
  {
    typedef Kernel_ K;

    #define CGAL_Kernel_pred(Y,Z) typedef CartesianKernelFunctors::Y<K> Y; \
    Y Z() const { return Y(); }
    #define CGAL_Kernel_cons(Y,Z) CGAL_Kernel_pred(Y,Z)

    #include <CGAL/Kernel/interface_macros.h>
  };

  template < typename RT_, typename FT_, typename Kernel_ >
  struct Result_of_homogeneous_base : public Result_of_base<RT_, FT_, Kernel_, Homogeneous_base> 
  {
    typedef Kernel_ K;

    #define CGAL_Kernel_pred(Y,Z) typedef HomogeneousKernelFunctors::Y<K> Y; \
    Y Z() const { return Y(); }
    #define CGAL_Kernel_cons(Y,Z) CGAL_Kernel_pred(Y,Z)

    #include <CGAL/Kernel/interface_macros.h>
  };

  template < typename FT_ >
  struct Result_of_cartesian
    : public Type_equality_wrapper<
    Result_of_cartesian_base<FT_, Result_of_cartesian<FT_> >,
    Result_of_cartesian<FT_> >
  {
    // this has to be delayed until here as AnyFunctor will
    // instantiate its arguments and only here lookup for all the
    // typedefs inside the functor will be possible
    typedef Result_of_cartesian K;
    #define CGAL_Kernel_pred(Y,Z) typedef AnyFunctor< CartesianKernelFunctors::Y<K> > Y; \
    Y Z() const { return Y(); }
    #define CGAL_Kernel_cons(Y,Z) CGAL_Kernel_pred(Y,Z)
    
    #include <CGAL/Kernel/interface_macros.h>
  };

  // same as above
  template < typename RT_, typename FT_ >
  struct Result_of_homogeneous
    : public Type_equality_wrapper<
    Result_of_homogeneous_base<RT_, FT_, Result_of_homogeneous<RT_, FT_> >,
    Result_of_homogeneous<RT_, FT_ > >
  {
    typedef Result_of_homogeneous K;
    #define CGAL_Kernel_pred(Y,Z) typedef AnyFunctor< HomogeneousKernelFunctors::Y<K> > Y; \
    Y Z() const { return Y(); }
    #define CGAL_Kernel_cons(Y,Z) CGAL_Kernel_pred(Y,Z)

    #include <CGAL/Kernel/interface_macros.h>
  };
}

#endif /* C++11 GUARD */
#endif /* CGAL_RESULT_OF_KERNEL_H */
