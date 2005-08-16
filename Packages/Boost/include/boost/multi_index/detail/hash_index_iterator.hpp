/* Copyright 2003-2005 Joaquín M López Muñoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org/libs/multi_index for library home page.
 */

#ifndef BOOST_MULTI_INDEX_DETAIL_HASH_INDEX_ITERATOR_HPP
#define BOOST_MULTI_INDEX_DETAIL_HASH_INDEX_ITERATOR_HPP

#if defined(_MSC_VER)&&(_MSC_VER>=1200)
#pragma once
#endif

#include <boost/config.hpp> /* keep it first to prevent nasty warns in MSVC */
#include <boost/detail/workaround.hpp>
#include <boost/multi_index/detail/hash_index_proxy.hpp>
#include <boost/multi_index/detail/safe_mode.hpp>
#include <boost/operators.hpp>

#if !defined(BOOST_MULTI_INDEX_DISABLE_SERIALIZATION)
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/split_member.hpp>
#endif

namespace boost{

namespace multi_index{

namespace detail{

/* An iterator template for nodes of multi_index::detail::hashed_index.
 * Built with the aid boost::forward_iterator_helper from
 * boost/operators.hpp.
 */

#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
#if BOOST_WORKAROUND(BOOST_MSVC,<1300)
template<typename Node,typename BucketArray>
class hashed_index_iterator:
  public boost::forward_iterator_helper<
    hashed_index_iterator<Node,BucketArray>,
    typename Node::value_type,
    std::ptrdiff_t,
    const typename Node::value_type*,
    const typename Node::value_type&>,
  public safe_iterator<hashed_index_proxy<Node,BucketArray> >
#else
template<typename Node,typename BucketArray,typename Container>
class hashed_index_iterator:
  public boost::forward_iterator_helper<
    hashed_index_iterator<Node,BucketArray,Container>,
    typename Node::value_type,
    std::ptrdiff_t,
    const typename Node::value_type*,
    const typename Node::value_type&>,
  public safe_iterator<Container>
#endif
#else
template<typename Node,typename BucketArray>
class hashed_index_iterator:
  public boost::forward_iterator_helper<
    hashed_index_iterator<Node,BucketArray>,
    typename Node::value_type,
    std::ptrdiff_t,
    const typename Node::value_type*,
    const typename Node::value_type&>
#endif

{
#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
public:

#if BOOST_WORKAROUND(BOOST_MSVC,<1300)
  typedef hashed_index_proxy<Node,BucketArray> container_type;
#else
  typedef Container                            container_type;
#endif

private:
  typedef safe_iterator<container_type> safe_super;

public:
  hashed_index_iterator():node(0),buckets(0){}
  hashed_index_iterator(
    Node* node_,BucketArray& buckets_,container_type* cont_):
    safe_super(cont_),node(node_),buckets(&buckets_){}

  hashed_index_iterator& operator=(const hashed_index_iterator& x)
  {
    BOOST_MULTI_INDEX_CHECK_VALID_ITERATOR(x);
    safe_super::operator=(x);
    node=x.node;
    buckets=x.buckets;
    return *this;
  }

#else
public:
  hashed_index_iterator(){}
  hashed_index_iterator(Node* node_,BucketArray& buckets_):
    node(node_),buckets(&buckets_){}
#endif

  const typename Node::value_type& operator*()const
  {
    BOOST_MULTI_INDEX_CHECK_VALID_ITERATOR(*this);
    BOOST_MULTI_INDEX_CHECK_DEREFERENCEABLE_ITERATOR(*this);
    return node->value;
  }

  hashed_index_iterator& operator++()
  {
    BOOST_MULTI_INDEX_CHECK_VALID_ITERATOR(*this);
    BOOST_MULTI_INDEX_CHECK_INCREMENTABLE_ITERATOR(*this);
    Node::increment(node,buckets->begin(),buckets->end());
    return *this;
  }

  friend bool operator==(const hashed_index_iterator& x,const hashed_index_iterator& y)
  {
    BOOST_MULTI_INDEX_CHECK_VALID_ITERATOR(x);
    BOOST_MULTI_INDEX_CHECK_VALID_ITERATOR(y);
    BOOST_MULTI_INDEX_CHECK_SAME_OWNER(x,y);
    return x.node==y.node;
  }

  /* get_node is not to be used by the user */

  Node* get_node()const{return node;}

private:
  Node*        node;
  BucketArray* buckets;

#if !defined(BOOST_MULTI_INDEX_DISABLE_SERIALIZATION)
  /* serialization */

  friend class boost::serialization::access;

  BOOST_SERIALIZATION_SPLIT_MEMBER()

  typedef typename Node::base_type node_base_type;

  template<class Archive>
  void save(Archive& ar,const unsigned int version)const
  {
    BOOST_MULTI_INDEX_CHECK_VALID_ITERATOR(*this);
    node_base_type* bnode=node;
    ar<<serialization::make_nvp("pointer",bnode);
    ar<<serialization::make_nvp("pointer",buckets);
  }

  template<class Archive>
  void load(Archive& ar,const unsigned int version)
  {
    node_base_type* bnode;
    ar>>serialization::make_nvp("pointer",bnode);
    node=static_cast<Node*>(bnode);
    ar>>serialization::make_nvp("pointer",buckets);

#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
    safe_super::uncheck();
#endif

  }
#endif
};

} /* namespace multi_index::detail */

} /* namespace multi_index */

} /* namespace boost */

#endif
