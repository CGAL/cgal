// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef OVERLOADS_FWD_DWA2002101_HPP
# define OVERLOADS_FWD_DWA2002101_HPP

namespace boost { namespace python { namespace detail { 

// forward declarations
struct overloads_base;
  
template <class OverloadsT, class NameSpaceT, class SigT>
inline void define_with_defaults(char const* name, OverloadsT const&, NameSpaceT&, SigT const&);

}}} // namespace boost::python::detail

#endif // OVERLOADS_FWD_DWA2002101_HPP
