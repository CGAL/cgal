// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef DEFAULT_CALL_POLICIES_DWA2002131_HPP
# define DEFAULT_CALL_POLICIES_DWA2002131_HPP

# include <boost/python/detail/prefix.hpp>
# include <boost/mpl/if.hpp>
# include <boost/python/to_python_value.hpp>
# include <boost/type_traits/transform_traits.hpp>

namespace boost { namespace python { 

template <class T> struct to_python_value;

namespace detail
{
// for "readable" error messages
  template <class T> struct specify_a_return_value_policy_to_wrap_functions_returning
# if defined(__GNUC__) && __GNUC__ >= 3 || defined(__EDG__)
  {}
# endif 
  ;
}

struct default_result_converter;

struct default_call_policies
{
    // Ownership of this argument tuple will ultimately be adopted by
    // the caller.
    template <class ArgumentPackage>
    static bool precall(ArgumentPackage const&)
    {
        return true;
    }

    // Pass the result through
    template <class ArgumentPackage>
    static PyObject* postcall(ArgumentPackage const&, PyObject* result)
    {
        return result;
    }

    typedef default_result_converter result_converter;
    typedef PyObject* argument_package;
};

struct default_result_converter
{
    template <class R>
    struct apply
    {
        BOOST_STATIC_CONSTANT(bool, is_illegal = is_reference<R>::value || is_pointer<R>::value);
        
        typedef typename mpl::if_c<
            is_illegal
            , detail::specify_a_return_value_policy_to_wrap_functions_returning<R>
            , boost::python::to_python_value<
                typename add_reference<typename add_const<R>::type>::type
                >
        >::type type;
    };
};

// Exceptions for c strings an PyObject*s
template <>
struct default_result_converter::apply<char const*>
{
    typedef boost::python::to_python_value<char const*const&> type;
};

template <>
struct default_result_converter::apply<PyObject*>
{
    typedef boost::python::to_python_value<PyObject*const&> type;
};

}} // namespace boost::python

#endif // DEFAULT_CALL_POLICIES_DWA2002131_HPP
