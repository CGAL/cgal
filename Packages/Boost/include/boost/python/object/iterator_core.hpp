// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef ITERATOR_CORE_DWA2002512_HPP
# define ITERATOR_CORE_DWA2002512_HPP

# include <boost/python/object_fwd.hpp>

namespace boost { namespace python { namespace objects {

BOOST_PYTHON_DECL object const& identity_function();
BOOST_PYTHON_DECL void stop_iteration_error();

}}} // namespace boost::python::object

#endif // ITERATOR_CORE_DWA2002512_HPP
