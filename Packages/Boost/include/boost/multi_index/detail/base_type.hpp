/* Copyright 2003-2004 Joaquín M López Muñoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org/libs/multi_index for library home page.
 */

#ifndef BOOST_MULTI_INDEX_DETAIL_BASE_TYPE_HPP
#define BOOST_MULTI_INDEX_DETAIL_BASE_TYPE_HPP

#include <boost/config.hpp> /* keep it first to prevent nasty warns in MSVC */
#include <boost/detail/workaround.hpp>
#include <boost/mpl/bind.hpp>
#include <boost/mpl/reverse_iter_fold.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/multi_index_container_fwd.hpp>
#include <boost/multi_index/detail/header_holder.hpp>
#include <boost/multi_index/detail/index_base.hpp>
#include <boost/multi_index/detail/is_index_list.hpp>
#include <boost/multi_index/detail/msvc_index_specifier.hpp>
#include <boost/multi_index/ordered_index_fwd.hpp>
#include <boost/multi_index/detail/prevent_eti.hpp>
#include <boost/static_assert.hpp>

namespace boost{

namespace multi_index{

namespace detail{

/* MPL machinery to construct a linear hierarchy of indices out of
 * a index list.
 */

#if BOOST_WORKAROUND(BOOST_MSVC,<1310)
struct index_applier
{
  template<typename IndexSpecifierIterator,typename Super>
  struct apply:
    msvc_index_specifier< mpl::deref<IndexSpecifierIterator>::type>::
      template result_index_class<Super>
  {
  }; 
};
#else
struct index_applier
{
  template<typename IndexSpecifierIterator,typename Super>
  struct apply
  {
    typedef typename mpl::deref<IndexSpecifierIterator>::type index_specifier;
    typedef typename index_specifier::
      BOOST_NESTED_TEMPLATE index_class<Super>::type type;
  }; 
};
#endif

template<typename Value,typename IndexSpecifierList,typename Allocator>
struct multi_index_base_type
{
  BOOST_STATIC_ASSERT(detail::is_index_list<IndexSpecifierList>::value);

  typedef typename prevent_eti<
    multi_index_container<Value,IndexSpecifierList,Allocator>,
    typename mpl::reverse_iter_fold<
      IndexSpecifierList,
      index_base<Value,IndexSpecifierList,Allocator>,
      mpl::bind2<
        index_applier,
        mpl::_2,
        mpl::_1
      >
    >::type
  >::type type;
};

} /* namespace multi_index::detail */

} /* namespace multi_index */

} /* namespace boost */

#endif
