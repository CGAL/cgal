// Copyright (C) 2001
// Mac Murrett
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee,
// provided that the above copyright notice appear in all copies and
// that both that copyright notice and this permission notice appear
// in supporting documentation.  Mac Murrett makes no representations
// about the suitability of this software for any purpose.  It is
// provided "as is" without express or implied warranty.
//
// See http://www.boost.org for most recent version including documentation.

#ifndef BOOST_FORCE_CAST_MJM012402_HPP
#define BOOST_FORCE_CAST_MJM012402_HPP

namespace boost {

namespace detail {

namespace thread {


// force_cast will convert anything to anything.

// general case
template<class Return_Type, class Argument_Type>
inline Return_Type &force_cast(Argument_Type &rSrc)
{   return(*reinterpret_cast<Return_Type *>(&rSrc));    }

// specialization for const
template<class Return_Type, class Argument_Type>
inline const Return_Type &force_cast(const Argument_Type &rSrc)
{   return(*reinterpret_cast<const Return_Type *>(&rSrc));  }


} // namespace thread

} // namespace detail

} // namespace boost


#endif // BOOST_FORCE_CAST_MJM012402_HPP
