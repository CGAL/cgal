// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef ARGS_FWD_DWA2002927_HPP
# define ARGS_FWD_DWA2002927_HPP

# include <boost/python/detail/prefix.hpp>

# include <boost/python/handle.hpp>
# include <boost/config.hpp>
# include <cstddef>
# include <utility>

namespace boost { namespace python { 

namespace detail
{
  struct keyword
  {
      keyword(char const* name_=0)
       : name(name_)
      {}
      
      char const* name;
      handle<> default_value;
  };
  
  template <std::size_t nkeywords = 0> struct keywords;
  
  typedef std::pair<keyword const*, keyword const*> keyword_range;
  
  template <>
  struct keywords<0>
  {
      BOOST_STATIC_CONSTANT(std::size_t, size = 0);
      static keyword_range range() { return keyword_range(); }
  };

  namespace error
  {
    template <int keywords, int function_args>
    struct more_keywords_than_function_arguments
    {
        typedef char too_many_keywords[keywords > function_args ? -1 : 1];
    };
  }
}

}} // namespace boost::python

#endif // ARGS_FWD_DWA2002927_HPP
