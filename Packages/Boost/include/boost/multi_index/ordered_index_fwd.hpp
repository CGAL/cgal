/* Copyright 2003-2004 Joaquín M López Muñoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org/libs/multi_index for library home page.
 */

#ifndef BOOST_MULTI_INDEX_ORDERED_INDEX_FWD_HPP
#define BOOST_MULTI_INDEX_ORDERED_INDEX_FWD_HPP

#include <boost/multi_index/detail/ord_index_args.hpp>

namespace boost{

namespace multi_index{

namespace detail{

template<
  typename KeyFromValue,typename Compare,
  typename Super,typename TagList,typename Category
>
class index;

template<
  typename KeyFromValue1,typename Compare1,
  typename Super1,typename TagList1,typename Category1,
  typename KeyFromValue2,typename Compare2,
  typename Super2,typename TagList2,typename Category2
>
bool operator==(
  const index<KeyFromValue1,Compare1,Super1,TagList1,Category1>& x,
  const index<KeyFromValue2,Compare2,Super2,TagList2,Category2>& y);

template<
  typename KeyFromValue1,typename Compare1,
  typename Super1,typename TagList1,typename Category1,
  typename KeyFromValue2,typename Compare2,
  typename Super2,typename TagList2,typename Category2
>
bool operator<(
  const index<KeyFromValue1,Compare1,Super1,TagList1,Category1>& x,
  const index<KeyFromValue2,Compare2,Super2,TagList2,Category2>& y);

template<
  typename KeyFromValue1,typename Compare1,
  typename Super1,typename TagList1,typename Category1,
  typename KeyFromValue2,typename Compare2,
  typename Super2,typename TagList2,typename Category2
>
bool operator!=(
  const index<KeyFromValue1,Compare1,Super1,TagList1,Category1>& x,
  const index<KeyFromValue2,Compare2,Super2,TagList2,Category2>& y);

template<
  typename KeyFromValue1,typename Compare1,
  typename Super1,typename TagList1,typename Category1,
  typename KeyFromValue2,typename Compare2,
  typename Super2,typename TagList2,typename Category2
>
bool operator>(
  const index<KeyFromValue1,Compare1,Super1,TagList1,Category1>& x,
  const index<KeyFromValue2,Compare2,Super2,TagList2,Category2>& y);

template<
  typename KeyFromValue1,typename Compare1,
  typename Super1,typename TagList1,typename Category1,
  typename KeyFromValue2,typename Compare2,
  typename Super2,typename TagList2,typename Category2
>
bool operator>=(
  const index<KeyFromValue1,Compare1,Super1,TagList1,Category1>& x,
  const index<KeyFromValue2,Compare2,Super2,TagList2,Category2>& y);

template<
  typename KeyFromValue1,typename Compare1,
  typename Super1,typename TagList1,typename Category1,
  typename KeyFromValue2,typename Compare2,
  typename Super2,typename TagList2,typename Category2
>
bool operator<=(
  const index<KeyFromValue1,Compare1,Super1,TagList1,Category1>& x,
  const index<KeyFromValue2,Compare2,Super2,TagList2,Category2>& y);

template<
  typename KeyFromValue,typename Compare,
  typename Super,typename TagList,typename Category
>
void swap(
  index<KeyFromValue,Compare,Super,TagList,Category>& x,
  index<KeyFromValue,Compare,Super,TagList,Category>& y);

} /* namespace multi_index::detail */

/* index specifiers */

template<
  typename Arg1,typename Arg2=detail::null_arg,typename Arg3=detail::null_arg
>
struct ordered_unique;

template<
  typename Arg1,typename Arg2=detail::null_arg,typename Arg3=detail::null_arg
>
struct ordered_non_unique;

} /* namespace multi_index */

} /* namespace boost */

#endif
