/* Copyright 2003-2004 Joaquín M López Muñoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org/libs/multi_index for library home page.
 */

#ifndef BOOST_MULTI_INDEX_DETAIL_INDEX_ITERATOR_HPP
#define BOOST_MULTI_INDEX_DETAIL_INDEX_ITERATOR_HPP

#include <boost/config.hpp> /* keep it first to prevent nasty warns in MSVC */
#include <boost/detail/workaround.hpp>
#include <boost/multi_index/detail/index_iterator_fwd.hpp>
#include <boost/multi_index/detail/index_proxy.hpp>
#include <boost/multi_index/detail/safe_mode.hpp>
#include <boost/operators.hpp>

namespace boost{

namespace multi_index{

namespace detail{

/* An iterator template for nodes of multi_index::detail::index.
 * Built with the aid boost::bidirectional_iterator_helper from
 * boost/operators.hpp.
 */

#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
#if BOOST_WORKAROUND(BOOST_MSVC,<1300)
template<typename Node>
class index_iterator:
  public boost::bidirectional_iterator_helper<
    index_iterator<Node>,
    typename Node::value_type,
    std::ptrdiff_t,
    const typename Node::value_type*,
    const typename Node::value_type&>,
  public safe_iterator<index_proxy<Node> >
#else
template<typename Node,typename Container>
class index_iterator:
  public boost::bidirectional_iterator_helper<
    index_iterator<Node,Container>,
    typename Node::value_type,
    std::ptrdiff_t,
    const typename Node::value_type*,
    const typename Node::value_type&>,
  public safe_iterator<Container>
#endif
#else
template<typename Node>
class index_iterator:
  public boost::bidirectional_iterator_helper<
    index_iterator<Node>,
    typename Node::value_type,
    std::ptrdiff_t,
    const typename Node::value_type*,
    const typename Node::value_type&>
#endif

{
#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
public:

#if BOOST_WORKAROUND(BOOST_MSVC,<1300)
  typedef index_proxy<Node> container_type;
#else
  typedef Container         container_type;
#endif

private:
  typedef safe_iterator<container_type> safe_super;

public:
  index_iterator():node(0){}
  index_iterator(Node* node_,container_type* cont_):
    safe_super(cont_),node(node_){}

  index_iterator& operator=(const index_iterator& x)
  {
    BOOST_MULTI_INDEX_CHECK_VALID_ITERATOR(x);
    safe_super::operator=(x);
    node=x.node;
    return *this;
  }

#else
public:
  index_iterator(){}
  index_iterator(Node* node_):node(node_){}
#endif

  const typename Node::value_type& operator*()const
  {
    BOOST_MULTI_INDEX_CHECK_VALID_ITERATOR(*this);
    BOOST_MULTI_INDEX_CHECK_DEREFERENCEABLE_ITERATOR(*this);
    return node->value;
  }

  index_iterator& operator++()
  {
    BOOST_MULTI_INDEX_CHECK_VALID_ITERATOR(*this);
    BOOST_MULTI_INDEX_CHECK_INCREMENTABLE_ITERATOR(*this);
    Node::increment(node);
    return *this;
  }

  index_iterator& operator--()
  {
    BOOST_MULTI_INDEX_CHECK_VALID_ITERATOR(*this);
    BOOST_MULTI_INDEX_CHECK_DECREMENTABLE_ITERATOR(*this);
    Node::decrement(node);
    return *this;
  }

  friend bool operator==(const index_iterator& x,const index_iterator& y)
  {
    BOOST_MULTI_INDEX_CHECK_VALID_ITERATOR(x);
    BOOST_MULTI_INDEX_CHECK_VALID_ITERATOR(y);
    BOOST_MULTI_INDEX_CHECK_SAME_OWNER(x,y);
    return x.node==y.node;
  }

  /* get_node is not to be used by the user */

  Node* get_node()const{return node;}

private:
  Node* node;
};

} /* namespace multi_index::detail */

} /* namespace multi_index */

} /* namespace boost */

#endif
