// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef EXCEPTION_TRANSLATOR_DWA2002810_HPP
# define EXCEPTION_TRANSLATOR_DWA2002810_HPP

# include <boost/python/detail/prefix.hpp>

# include <boost/bind.hpp>
# include <boost/type.hpp>
# include <boost/python/detail/translate_exception.hpp>
# include <boost/python/detail/exception_handler.hpp>

namespace boost { namespace python { 

template <class ExceptionType, class Translate>
void register_exception_translator(Translate const& translate, boost::type<ExceptionType>* = 0)
{
    detail::register_exception_handler(
        bind<bool>(detail::translate_exception<ExceptionType,Translate>(), _1, _2, translate)
        );
}

}} // namespace boost::python

#endif // EXCEPTION_TRANSLATOR_DWA2002810_HPP
