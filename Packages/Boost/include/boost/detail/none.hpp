// Copyright (C) 2003, Fernando Luis Cacciola Carballal.
//
// Use, modification, and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/lib/optional for documentation.
//
// You are welcome to contact the author at:
//  fernando_cacciola@hotmail.com
//
#ifndef BOOST_UTILITY_NONE_17SEP2003_HPP
#define BOOST_UTILITY_NONE_17SEP2003_HPP

#include "boost/detail/none_t.hpp"

// NOTE: Borland users have to include this header outside any precompiled headers
// (bcc<=5.64 cannot include instance data in a precompiled header)

namespace boost {

namespace {

detail::none_t const none = ((detail::none_t)0) ;

}

} // namespace boost

#endif

