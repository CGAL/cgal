// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef SLICE_NIL_DWA2002620_HPP
# define SLICE_NIL_DWA2002620_HPP

# include <boost/python/detail/prefix.hpp>

namespace boost { namespace python { namespace api {

class object;

enum slice_nil
{
# ifndef _ // Watch out for GNU gettext users, who #define _(x)
      _
# endif 
};

template <class T>
struct slice_bound
{
    typedef object type;
};

template <>
struct slice_bound<slice_nil>
{
    typedef slice_nil type;
};

}

using api::slice_nil;
# ifndef _ // Watch out for GNU gettext users, who #define _(x)
using api::_;
# endif 

}} // namespace boost::python

#endif // SLICE_NIL_DWA2002620_HPP
