// Copyright (C) 2000 Stephen Cleary
//
// This file can be redistributed and/or modified under the terms found
//  in "copyright.html"
// This software and its documentation is provided "as is" without express or
//  implied warranty, and with no claim as to its suitability for any purpose.
//
// See http://www.boost.org for updates, documentation, and revision history.

#ifndef BOOST_POOL_GUARD_HPP
#define BOOST_POOL_GUARD_HPP

// Extremely Light-Weight guard glass

namespace boost {

namespace details {
namespace pool {

template <typename Mutex>
class guard
{
  private:
    Mutex & mtx;

    guard(const guard &);
    void operator=(const guard &);

  public:
    explicit guard(Mutex & nmtx)
    :mtx(nmtx) { mtx.lock(); }

    ~guard() { mtx.unlock(); }
};

} // namespace pool
} // namespace details

} // namespace boost

#endif
