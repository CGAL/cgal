//  Boost noncopyable.hpp header file  --------------------------------------//

//  (C) Copyright Boost.org 1999-2003. Permission to copy, use, modify, sell
//  and distribute this software is granted provided this copyright
//  notice appears in all copies. This software is provided "as is" without
//  express or implied warranty, and with no claim as to its suitability for
//  any purpose.

//  See http://www.boost.org/libs/utility for documentation.

#ifndef BOOST_NONCOPYABLE_HPP_INCLUDED
#define BOOST_NONCOPYABLE_HPP_INCLUDED

namespace boost {

//  Private copy constructor and copy assignment ensure classes derived from
//  class noncopyable cannot be copied.

//  Contributed by Dave Abrahams

class noncopyable
{
 protected:
    noncopyable() {}
    ~noncopyable() {}
 private:  // emphasize the following members are private
    noncopyable( const noncopyable& );
    const noncopyable& operator=( const noncopyable& );
};

} // namespace boost

#endif  // BOOST_NONCOPYABLE_HPP_INCLUDED
