// Copyright David Abrahams 2004. Use, modification and distribution is
// subject to the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#ifndef INDIRECT_REFERENCE_DWA200415_HPP
# define INDIRECT_REFERENCE_DWA200415_HPP

// dereferenceable_traits provides access to the value_type and
// reference of a Dereferenceable type.

# include <boost/detail/is_incrementable.hpp>
# include <boost/iterator/iterator_traits.hpp>
# include <boost/type_traits/remove_cv.hpp>
# include <boost/mpl/eval_if.hpp>
# include <boost/pointee.hpp>

namespace boost { 

namespace detail
{
  template <class P>
  struct smart_ptr_reference
  {
      typedef typename boost::pointee<P>::type& type;
  };
}

template <class P>
struct indirect_reference
  : mpl::eval_if<
        detail::is_incrementable<P>
      , iterator_reference<P>
      , detail::smart_ptr_reference<P>
    >
{
};
  
} // namespace boost

#endif // INDIRECT_REFERENCE_DWA200415_HPP
