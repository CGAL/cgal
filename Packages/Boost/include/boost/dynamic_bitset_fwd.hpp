// (C) Copyright Chuck Allison and Jeremy Siek 2001, 2002.
//
// Permission to copy, use, modify, sell and distribute this software
// is granted provided this copyright notice appears in all
// copies. This software is provided "as is" without express or
// implied warranty, and with no claim as to its suitability for any
// purpose.

// See http://www.boost.org/libs/dynamic_bitset for documentation.

#ifndef BOOST_DYNAMIC_BITSET_FWD_HPP
#define BOOST_DYNAMIC_BITSET_FWD_HPP

#include <memory>

namespace boost {

template <typename Block = unsigned long,
          typename Allocator = std::allocator<Block> >
class dynamic_bitset;

} // namespace boost

#endif // BOOST_DYNAMIC_BITSET_FWD_HPP
