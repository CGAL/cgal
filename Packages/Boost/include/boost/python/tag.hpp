// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef TAG_DWA2002720_HPP
# define TAG_DWA2002720_HPP

# include <boost/python/detail/prefix.hpp>

namespace boost { namespace python { 

// used only to prevent argument-dependent lookup from finding the
// wrong function in some cases. Cheaper than qualification.
enum tag_t { tag };

}} // namespace boost::python

#endif // TAG_DWA2002720_HPP
