// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef RAW_PYOBJECT_DWA2002628_HPP
# define RAW_PYOBJECT_DWA2002628_HPP

namespace boost { namespace python { namespace detail { 

//
// Define some types which we can use to get around the vagaries of
// PyObject*. We will use these to initialize object instances, and
// keep them in namespace detail to make sure they stay out of the
// hands of users. That is much simpler than trying to grant
// friendship to all the appropriate parties.
//

// New references are normally checked for null
struct new_reference_t;
typedef new_reference_t* new_reference;

// Borrowed references are assumed to be non-null
struct borrowed_reference_t;
typedef borrowed_reference_t* borrowed_reference;

// New references which aren't checked for null
struct new_non_null_reference_t;
typedef new_non_null_reference_t* new_non_null_reference;

}}} // namespace boost::python::detail

#endif // RAW_PYOBJECT_DWA2002628_HPP
