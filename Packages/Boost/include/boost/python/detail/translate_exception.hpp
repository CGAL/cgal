// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef TRANSLATE_EXCEPTION_DWA2002810_HPP
# define TRANSLATE_EXCEPTION_DWA2002810_HPP

# include <boost/python/detail/exception_handler.hpp>

# include <boost/call_traits.hpp>
# include <boost/type_traits/add_const.hpp>

# include <boost/function/function0.hpp>

namespace boost { namespace python { namespace detail { 

// A ternary function object used to translate C++ exceptions of type
// ExceptionType into Python exceptions by invoking an object of type
// Translate. Typically the translate function will be curried with
// boost::bind().
template <class ExceptionType, class Translate>
struct translate_exception
{
    typedef typename add_reference<
        typename add_const<ExceptionType>::type
    >::type exception_cref;
    
    inline bool operator()(
        exception_handler const& handler
      , function0<void> const& f
      , typename call_traits<Translate>::param_type translate) const
    {
        try
        {
            return handler(f);
        }
        catch(exception_cref e)
        {
            translate(e);
            return true;
        }
    }
};

}}} // namespace boost::python::detail

#endif // TRANSLATE_EXCEPTION_DWA2002810_HPP
