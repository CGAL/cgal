// (C) Copyright Toon Knapen    2001.
// (C) Copyright David Abrahams 2003.
// (C) Copyright Roland Richter 2003.
// Permission to copy, use, modify, sell and distribute this software
// is granted provided this copyright notice appears in all copies.
// This software is provided "as is" without express or implied
// warranty, and with no claim as to its suitability for any purpose.

#ifndef BOOST_PERMUTATION_ITERATOR_HPP
#define BOOST_PERMUTATION_ITERATOR_HPP

#include <iterator>

#include <boost/iterator/iterator_adaptor.hpp>


namespace boost
{

template< class ElementIterator
        , class IndexIterator>
class permutation_iterator
  : public iterator_adaptor< 
             permutation_iterator<ElementIterator, IndexIterator>
           , IndexIterator, typename detail::iterator_traits<ElementIterator>::value_type
           , use_default, typename detail::iterator_traits<ElementIterator>::reference>
{
  typedef iterator_adaptor< 
            permutation_iterator<ElementIterator, IndexIterator>
          , IndexIterator, typename detail::iterator_traits<ElementIterator>::value_type
          , use_default, typename detail::iterator_traits<ElementIterator>::reference> super_t;

  friend class iterator_core_access;

public:
  permutation_iterator() : m_elt_iter() {}

  explicit permutation_iterator(ElementIterator x, IndexIterator y) 
      : super_t(y), m_elt_iter(x) {}

  template<class OtherElementIterator, class OtherIndexIterator>
  permutation_iterator(
      permutation_iterator<OtherElementIterator, OtherIndexIterator> const& r
      , typename enable_if_convertible<OtherElementIterator, ElementIterator>::type* = 0
      , typename enable_if_convertible<OtherIndexIterator, IndexIterator>::type* = 0
      )
    : super_t(r.base()), m_elt_iter(r.m_elt_iter)
  {}

private:
    typename super_t::reference dereference() const
        { return *(m_elt_iter + *this->base()); }

    ElementIterator m_elt_iter;
};


template <class ElementIterator, class IndexIterator>
permutation_iterator<ElementIterator, IndexIterator> 
make_permutation_iterator( ElementIterator e, IndexIterator i )
{
    return permutation_iterator<ElementIterator, IndexIterator>( e, i );
}


} // namespace boost

#endif
