/* Copyright 2003-2004 Joaquín M López Muñoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org/libs/multi_index for library home page.
 */

#ifndef BOOST_MULTI_INDEX_DETAIL_HASH_hashed_index_proxy_HPP
#define BOOST_MULTI_INDEX_DETAIL_HASH_hashed_index_proxy_HPP

#if defined(_MSC_VER)&&(_MSC_VER>=1200)
#pragma once
#endif

#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
#include <boost/config.hpp> /* keep it first to prevent nasty warns in MSVC */
#include <boost/detail/workaround.hpp>

#if BOOST_WORKAROUND(BOOST_MSVC,<1300)
#include <algorithm>
#include <boost/multi_index/detail/hash_index_iterator_fwd.hpp>
#include <boost/multi_index/detail/safe_mode.hpp>

namespace boost{

namespace multi_index{

namespace detail{

/* Analogous of index_proxy for interoperability with
 * hashed_index_iterator.
 */

template<typename Node,typename BucketArray>
class hashed_index_proxy:
  public safe_container<hashed_index_proxy<Node,BucketArray> >
{
protected:
  hashed_index_proxy(Node* header_,BucketArray& buckets_):
    header(header_),buckets(&buckets_){}

  void swap(hashed_index_proxy<Node,BucketArray>& x)
  {
    std::swap(header,x.header);
    std::swap(buckets,x.buckets);
    safe_container<hashed_index_proxy<Node,BucketArray> >::swap(x);
  }

public:
  typedef hashed_index_iterator<Node,BucketArray> iterator;
  typedef hashed_index_iterator<Node,BucketArray> const_iterator;

  hashed_index_iterator<Node,BucketArray> end()const
  {
    return hashed_index_iterator<Node,BucketArray>(
      Node::end(header),*buckets,const_cast<hashed_index_proxy*>(this));
  }

private:
  Node*        header;
  BucketArray* buckets;
};

} /* namespace multi_index::detail */

} /* namespace multi_index */

} /* namespace boost */

#endif /* workaround */

#endif /* BOOST_MULTI_INDEX_ENABLE_SAFE_MODE */

#endif
