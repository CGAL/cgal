// Copyright David Abrahams 2003. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef PREFIX_DWA2003531_HPP
# define PREFIX_DWA2003531_HPP

// The rule is that <Python.h> must be included before any system
// headers (so it can get control over some awful macros).
// Unfortunately, Boost.Python needs to #include <limits.h> first, at
// least... but this gets us as close as possible.
# include <boost/python/detail/wrap_python.hpp>
# include <boost/python/detail/config.hpp>

#endif // PREFIX_DWA2003531_HPP
