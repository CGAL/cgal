//-----------------------------------------------------------------------------
// boost variant/detail/cast_storage.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2003
// Eric Friedman
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_VARIANT_DETAIL_CAST_STORAGE_HPP
#define BOOST_VARIANT_DETAIL_CAST_STORAGE_HPP

#include "boost/config.hpp"

namespace boost {
namespace detail { namespace variant {

///////////////////////////////////////////////////////////////////////////////
// (detail) function template cast_storage
//
// Casts the given storage to the specified type, but with qualification.
//

template <typename T>
inline T& cast_storage(
      void* storage
      BOOST_APPEND_EXPLICIT_TEMPLATE_TYPE(T)
    )
{
    return *static_cast<T*>(storage);
}

template <typename T>
inline const T& cast_storage(
      const void* storage
      BOOST_APPEND_EXPLICIT_TEMPLATE_TYPE(T)
    )
{
    return *static_cast<const T*>(storage);
}

}} // namespace detail::variant
} // namespace boost

#endif // BOOST_VARIANT_DETAIL_CAST_STORAGE_HPP
