/* boost random/detail/iterator_mixin.hpp header file
 *
 * Copyright Jens Maurer 2000-2001
 * Permission to use, copy, modify, sell, and distribute this software
 * is hereby granted without fee provided that the above copyright notice
 * appears in all copies and that both that copyright notice and this
 * permission notice appear in supporting documentation,
 *
 * Jens Maurer makes no representations about the suitability of this
 * software for any purpose. It is provided "as is" without express or
 * implied warranty.
 *
 * See http://www.boost.org for most recent version including documentation.
 *
 * Revision history
 */

#ifndef BOOST_ITERATOR_MIXIN_HPP
#define BOOST_ITERATOR_MIXIN_HPP

#include <boost/operators.hpp>

namespace boost {

// must be in boost namespace, otherwise the inline friend trick fails
template<class Generator, class ResultType>
class generator_iterator_mixin_adapter
  : incrementable<Generator>, equality_comparable<Generator>
{
public:
  typedef std::input_iterator_tag iterator_category;
  typedef ResultType value_type;
  typedef std::ptrdiff_t difference_type;
  typedef const value_type * pointer;
  typedef const value_type & reference;
  Generator& operator++() { v = cast()(); return cast(); }
  const value_type& operator*() const { return v; }

protected:
  // instantiate from derived classes only
  generator_iterator_mixin_adapter() { }
  void iterator_init() { operator++(); }
private:
  Generator & cast() { return static_cast<Generator&>(*this); }
  value_type v;
};

} // namespace boost

#endif // BOOST_ITERATOR_MIXIN_HPP
