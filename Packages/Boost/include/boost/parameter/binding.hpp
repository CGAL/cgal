// Copyright David Abrahams 2005. Distributed under the Boost
// Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_PARAMETER_BINDING_DWA200558_HPP
# define BOOST_PARAMETER_BINDING_DWA200558_HPP

# include <boost/mpl/apply.hpp>
# include <boost/parameter/aux_/result_of0.hpp>

#  include <boost/parameter/aux_/void.hpp>

# if BOOST_WORKAROUND(BOOST_MSVC, < 1300)
#  include <boost/parameter/aux_/void.hpp>
#  include <boost/type_traits/is_same.hpp>
# endif 

namespace boost { namespace parameter { 

// A metafunction that, given an argument pack, returns the type of
// the parameter identified by the given keyword.  If no such
// parameter has been specified, returns Default
# if !BOOST_WORKAROUND(BOOST_MSVC, < 1300)
template <class Parameters, class Keyword, class Default = void>
struct binding
  : mpl::apply_wrap2<
        typename Parameters::binding,Keyword,Default
    >
{};
# else
template <class Parameters, class Keyword, class Default = aux::void_>
struct binding
{
    typedef typename mpl::apply_wrap2<
        typename Parameters::binding,Keyword,
        typename mpl::if_<is_same<Default,aux::void_>,void,Default>::type
    >::type type;
};
# endif

// A metafunction that, given an argument pack, returns the type of
// the parameter identified by the given keyword.  If no such
// parameter has been specified, returns the type returned by invoking
// DefaultFn
template <class Parameters, class Keyword, class DefaultFn>
struct lazy_binding
{
  typedef typename mpl::apply_wrap2<
      typename Parameters::binding
    , Keyword
    , typename aux::result_of0<DefaultFn>::type
  >::type type;
};


}} // namespace boost::parameter

#endif // BOOST_PARAMETER_BINDING_DWA200558_HPP
