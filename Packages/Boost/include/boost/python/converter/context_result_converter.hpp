// Copyright David Abrahams 2003. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef CONTEXT_RESULT_CONVERTER_DWA2003917_HPP
# define CONTEXT_RESULT_CONVERTER_DWA2003917_HPP

namespace boost { namespace python { namespace converter { 

// A ResultConverter base class used to indicate that this result
// converter should be constructed with the original Python argument
// list.
struct context_result_converter {};

}}} // namespace boost::python::converter

#endif // CONTEXT_RESULT_CONVERTER_DWA2003917_HPP
