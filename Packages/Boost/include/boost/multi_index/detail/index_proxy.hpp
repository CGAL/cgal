/* Copyright 2003-2004 Joaquín M López Muñoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org/libs/multi_index for library home page.
 */

#ifndef BOOST_MULTI_INDEX_DETAIL_INDEX_PROXY_HPP
#define BOOST_MULTI_INDEX_DETAIL_INDEX_PROXY_HPP

#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
#include <boost/config.hpp> /* keep it first to prevent nasty warns in MSVC */
#include <boost/detail/workaround.hpp>

#if BOOST_WORKAROUND(BOOST_MSVC,<1300)
#include <algorithm>
#include <boost/multi_index/detail/index_iterator_fwd.hpp>
#include <boost/multi_index/detail/safe_mode.hpp>

namespace boost{

namespace multi_index{

namespace detail{

/* In safe mode, index iterators are derived from safe_iterator<Index>,
 * where Index is the type of the index where the iterator belongs. Due
 * to the long symbol names of indices, MSVC++ 6.0 often issues a
 * LNK1179 (duplicate comdat) error. To workaround this problem,
 * index_proxy is used instead. index_proxy<Node> acts as an index
 * over nodes of type Node in all aspects relevant to safe_iterator, and
 * its shorter symbol name makes life easier for MSVC++ 6.0.
 */

template<typename Node>
class index_proxy:public safe_container<index_proxy<Node> >
{
protected:
  index_proxy(Node* header_):header(header_){}

  void swap(index_proxy<Node>& x)
  {
    std::swap(header,x.header);
    safe_container<index_proxy<Node> >::swap(x);
  }

public:
  typedef index_iterator<Node> iterator;
  typedef index_iterator<Node> const_iterator;

  index_iterator<Node> begin()const
  {
    return index_iterator<Node>(
      Node::begin(header),const_cast<index_proxy*>(this));
  }

  index_iterator<Node> end()const
  {
    return index_iterator<Node>(
      Node::end(header),const_cast<index_proxy*>(this));
  }

private:
  Node* header;
};

} /* namespace multi_index::detail */

} /* namespace multi_index */

} /* namespace boost */

#endif /* workaround */

#endif /* BOOST_MULTI_INDEX_ENABLE_SAFE_MODE */

#endif
