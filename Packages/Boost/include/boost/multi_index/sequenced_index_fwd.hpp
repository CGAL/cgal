/* Copyright 2003-2004 Joaquín M López Muñoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org/libs/multi_index for library home page.
 */

#ifndef BOOST_MULTI_INDEX_SEQUENCED_INDEX_FWD_HPP
#define BOOST_MULTI_INDEX_SEQUENCED_INDEX_FWD_HPP

#include <boost/multi_index/tag.hpp>

namespace boost{

namespace multi_index{

namespace detail{

template<typename Super,typename TagList>
class sequenced_index;

template<
  typename Super1,typename TagList1,
  typename Super2,typename TagList2
>
bool operator==(
  const sequenced_index<Super1,TagList1>& x,
  const sequenced_index<Super2,TagList2>& y);

template<
  typename Super1,typename TagList1,
  typename Super2,typename TagList2
>
bool operator<(
  const sequenced_index<Super1,TagList1>& x,
  const sequenced_index<Super2,TagList2>& y);

template<
  typename Super1,typename TagList1,
  typename Super2,typename TagList2
>
bool operator!=(
  const sequenced_index<Super1,TagList1>& x,
  const sequenced_index<Super2,TagList2>& y);

template<
  typename Super1,typename TagList1,
  typename Super2,typename TagList2
>
bool operator>(
  const sequenced_index<Super1,TagList1>& x,
  const sequenced_index<Super2,TagList2>& y);

template<
  typename Super1,typename TagList1,
  typename Super2,typename TagList2
>
bool operator>=(
  const sequenced_index<Super1,TagList1>& x,
  const sequenced_index<Super2,TagList2>& y);

template<
  typename Super1,typename TagList1,
  typename Super2,typename TagList2
>
bool operator<=(
  const sequenced_index<Super1,TagList1>& x,
  const sequenced_index<Super2,TagList2>& y);

template<typename Super,typename TagList>
void swap(
  sequenced_index<Super,TagList>& x,
  sequenced_index<Super,TagList>& y);

} /* namespace multi_index::detail */

/* index specifiers */

template <typename TagList=tag<> >
struct sequenced;

} /* namespace multi_index */

} /* namespace boost */

#endif
