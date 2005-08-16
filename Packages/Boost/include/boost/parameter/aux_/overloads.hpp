// Copyright David Abrahams, Daniel Wallin 2003. Use, modification and 
// distribution is subject to the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)

// This file generates overloads in this format:
//
//     template<class A0, class A1>
//     typename aux::make_arg_list<
//         PS0,A0
//       , aux::make_arg_list<
//             PS1,A1
//           , mpl::identity<aux::empty_arg_list>
//         >
//     >::type
//     operator()(A0 const& a0, A1 const& a1) const
//     {
//         typedef typename aux::make_arg_list<
//             PS0,A0
//           , aux::make_arg_list<
//                 PS1,A1
//               , mpl::identity<aux::empty_arg_list>
//             >
//         >::type arg_tuple;
//
//         return arg_tuple(
//             a0
//           , a1
//           , aux::void_()
//             ...
//         );
//     }
//

#if !defined(BOOST_PP_IS_ITERATING)
# error Boost.Parameters - do not include this file!
#endif

#define N BOOST_PP_ITERATION()

#define BOOST_PARAMETER_open_list(z, n, text) \
    aux::make_arg_list< \
        BOOST_PP_CAT(PS, n), BOOST_PP_CAT(A, n) \

#define BOOST_PARAMETER_close_list(z, n, text) > 

#define BOOST_PARAMETER_arg_list(n) \
    BOOST_PP_ENUM(N, BOOST_PARAMETER_open_list, _) \
  , mpl::identity<aux::empty_arg_list> \
    BOOST_PP_REPEAT(N, BOOST_PARAMETER_close_list, _) 

template<BOOST_PP_ENUM_PARAMS(N, class A)>
typename BOOST_PARAMETER_arg_list(N)::type
operator()(BOOST_PP_ENUM_BINARY_PARAMS(N, A, const& a)) const
{
    typedef typename BOOST_PARAMETER_arg_list(N)::type arg_tuple;

    return arg_tuple(
        BOOST_PP_ENUM_PARAMS(N, a)
        BOOST_PP_ENUM_TRAILING_PARAMS(
            BOOST_PP_SUB(BOOST_PARAMETER_MAX_ARITY, N)
          , aux::void_() BOOST_PP_INTERCEPT
        ));
}

#undef BOOST_PARAMETER_arg_list
#undef BOOST_PARAMETER_open_list
#undef BOOST_PARAMETER_close_list
#undef N

