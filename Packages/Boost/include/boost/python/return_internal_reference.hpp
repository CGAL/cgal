// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef RETURN_INTERNAL_REFERENCE_DWA2002131_HPP
# define RETURN_INTERNAL_REFERENCE_DWA2002131_HPP

# include <boost/python/detail/prefix.hpp>

# include <boost/python/default_call_policies.hpp>
# include <boost/python/reference_existing_object.hpp>
# include <boost/python/with_custodian_and_ward.hpp>
# include <boost/mpl/if.hpp>

namespace boost { namespace python { 

namespace detail
{
  template <std::size_t>
  struct return_internal_reference_owner_arg_must_be_greater_than_zero
# if defined(__GNUC__) && __GNUC__ >= 3 || defined(__EDG__)
  {}
# endif
  ;
}

template <std::size_t owner_arg = 1, class BasePolicy_ = default_call_policies>
struct return_internal_reference
    : with_custodian_and_ward_postcall<0, owner_arg, BasePolicy_>
{
 private:
    BOOST_STATIC_CONSTANT(bool, legal = owner_arg > 0);
 public:
    typedef typename mpl::if_c<
        legal
        , reference_existing_object
        , detail::return_internal_reference_owner_arg_must_be_greater_than_zero<owner_arg>
    >::type result_converter;
};

}} // namespace boost::python

#endif // RETURN_INTERNAL_REFERENCE_DWA2002131_HPP
