// Copyright (C) 2000, 2001 Stephen Cleary
//
// This file can be redistributed and/or modified under the terms found
//  in "copyright.html"
// This software and its documentation is provided "as is" without express or
//  implied warranty, and with no claim as to its suitability for any purpose.
//
// See http://www.boost.org for updates, documentation, and revision history.

#ifndef BOOST_POOLFWD_HPP
#define BOOST_POOLFWD_HPP

#include <boost/config.hpp> // for workarounds

// std::size_t
#include <cstddef>

// boost::details::pool::default_mutex
#include <boost/pool/detail/mutex.hpp>

namespace boost {

//
// Location: <boost/pool/simple_segregated_storage.hpp>
//
template <typename SizeType = std::size_t>
class simple_segregated_storage;

//
// Location: <boost/pool/pool.hpp>
//
struct default_user_allocator_new_delete;
struct default_user_allocator_malloc_free;

template <typename UserAllocator = default_user_allocator_new_delete>
class pool;

//
// Location: <boost/pool/object_pool.hpp>
//
template <typename T, typename UserAllocator = default_user_allocator_new_delete>
class object_pool;

//
// Location: <boost/pool/singleton_pool.hpp>
//
template <typename Tag, unsigned RequestedSize,
    typename UserAllocator = default_user_allocator_new_delete,
    typename Mutex = details::pool::default_mutex,
    unsigned NextSize = 32>
struct singleton_pool;

//
// Location: <boost/pool/pool_alloc.hpp>
//
struct pool_allocator_tag;

template <typename T,
    typename UserAllocator = default_user_allocator_new_delete,
    typename Mutex = details::pool::default_mutex,
    unsigned NextSize = 32>
class pool_allocator;

struct fast_pool_allocator_tag;

template <typename T,
    typename UserAllocator = default_user_allocator_new_delete,
    typename Mutex = details::pool::default_mutex,
    unsigned NextSize = 32>
class fast_pool_allocator;

} // namespace boost

#endif
