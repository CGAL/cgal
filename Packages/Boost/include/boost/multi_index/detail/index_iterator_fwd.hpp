/* Copyright 2003-2004 Joaquín M López Muñoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org/libs/multi_index for library home page.
 */

#ifndef BOOST_MULTI_INDEX_DETAIL_INDEX_ITERATOR_FWD_HPP
#define BOOST_MULTI_INDEX_DETAIL_INDEX_ITERATOR_FWD_HPP

#include <boost/config.hpp>
#include <boost/detail/workaround.hpp>

namespace boost{

namespace multi_index{

namespace detail{

#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
#if BOOST_WORKAROUND(BOOST_MSVC,<1300)
template<typename Node>
class index_iterator;
#else
template<typename Node,typename Container>
class index_iterator;
#endif
#else
template<typename Node>
class index_iterator;
#endif

} /* namespace multi_index::detail */

} /* namespace multi_index */

} /* namespace boost */

#endif
