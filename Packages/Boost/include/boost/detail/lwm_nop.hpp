#ifndef BOOST_DETAIL_LWM_NOP_HPP_INCLUDED
#define BOOST_DETAIL_LWM_NOP_HPP_INCLUDED

// MS compatible compilers support #pragma once

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

//
//  boost/detail/lwm_nop.hpp
//
//  Copyright (c) 2002 Peter Dimov and Multi Media Ltd.
//
//  Permission to copy, use, modify, sell and distribute this software
//  is granted provided this copyright notice appears in all copies.
//  This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.
//

namespace boost
{

namespace detail
{

class lightweight_mutex
{
public:

    typedef lightweight_mutex scoped_lock;
};

} // namespace detail

} // namespace boost

#endif // #ifndef BOOST_DETAIL_LWM_NOP_HPP_INCLUDED
