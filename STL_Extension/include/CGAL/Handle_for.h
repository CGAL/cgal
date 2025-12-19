// Copyright (c) 1999,2001,2003
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stefan Schirra, Sylvain Pion

#ifndef CGAL_HANDLE_FOR_H
#define CGAL_HANDLE_FOR_H

#include <CGAL/disable_warnings.h>

#include <CGAL/config.h>
#include <CGAL/assertions.h> // for CGAL_assume

#include <boost/config.hpp>
#include <CGAL/memory.h>
#include <algorithm>
#include <cstddef>
#include <atomic>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4345) // Avoid warning https://learn.microsoft.com/en-us/previous-versions/wewb47ee(v=vs.120)
#endif
namespace CGAL {

template <class T, class Alloc = CGAL_ALLOCATOR(T) >
class Handle_for
{
    // Wrapper that adds the reference counter.
    struct RefCounted {
        T t;
        std::atomic_uint count;
        template <class... U>
        RefCounted(U&&...u ) : t{std::forward<U>(u)...}, count(1) {}
    };


    typedef std::allocator_traits<Alloc> Alloc_traits;
    typedef typename Alloc_traits::template rebind_alloc<RefCounted>           Allocator;
    typedef std::allocator_traits<Allocator> Allocator_traits;
    typedef typename Alloc_traits::template rebind_traits<RefCounted>::pointer pointer;

    static Allocator   allocator;
    pointer            ptr_;

public:

    typedef T element_type;

    typedef std::ptrdiff_t Id_type ;

    Handle_for()
    {
        pointer p = allocator.allocate(1);
        ptr_ = new (p) RefCounted();
    }

    Handle_for(const element_type& t)
    {
        pointer p = allocator.allocate(1);
        ptr_ = new (p) RefCounted(t);
    }

    Handle_for(element_type && t)
    {
        pointer p = allocator.allocate(1);
        ptr_ = new (p) RefCounted(std::move(t));
    }

/* I comment this one for now, since it's preventing the automatic conversions
   to take place.  We'll see if it's a problem later.
    template < typename T1 >
    Handle_for(const T1& t1)
    {
        pointer p = allocator.allocate(1);
        ptr_ = new (p) RefCounted(t1);
    }
*/

    template < typename T1, typename T2, typename... Args >
    Handle_for(T1 && t1, T2 && t2, Args && ... args)
    {
        pointer p = allocator.allocate(1);
        ptr_ = new (p) RefCounted(std::forward<T1>(t1), std::forward<T2>(t2), std::forward<Args>(args)...);
    }

    Handle_for(const Handle_for& h) noexcept(!CGAL_ASSERTIONS_ENABLED)
      : ptr_(h.ptr_)
    {
        // CGAL_assume (ptr_->count > 0);
        if (is_currently_single_threaded())
          ptr_->count.store(ptr_->count.load(std::memory_order_relaxed) + 1, std::memory_order_relaxed);
        else
          ptr_->count.fetch_add(1, std::memory_order_relaxed);
    }

    Handle_for&
    operator=(const Handle_for& h) noexcept(!CGAL_ASSERTIONS_ENABLED)
    {
        Handle_for tmp = h;
        swap(tmp);
        return *this;
    }

    Handle_for&
    operator=(const element_type &t)
    {
        if (is_shared())
            *this = Handle_for(t);
        else
            ptr_->t = t;

        return *this;
    }

    // Note : I don't see a way to make a useful move constructor, apart
    //        from e.g. using nullptr as a ptr value, but this is drastic.

    Handle_for&
    operator=(Handle_for && h) noexcept
    {
        swap(h);
        return *this;
    }

    Handle_for&
    operator=(element_type && t)
    {
        if (is_shared())
            *this = Handle_for(std::move(t));
        else
            ptr_->t = std::move(t);

        return *this;
    }

    ~Handle_for()
    {
      if (is_currently_single_threaded()) {
        auto c = ptr_->count.load(std::memory_order_relaxed);
        if (c == 1) {
          Allocator_traits::destroy(allocator, ptr_);
          allocator.deallocate(ptr_, 1);
        } else {
          ptr_->count.store(c - 1, std::memory_order_relaxed);
        }
      } else {
      // TSAN does not support fences :-(
#if !defined __SANITIZE_THREAD__ && !__has_feature(thread_sanitizer)
        if (ptr_->count.load(std::memory_order_relaxed) == 1
            || ptr_->count.fetch_sub(1, std::memory_order_release) == 1) {
          std::atomic_thread_fence(std::memory_order_acquire);
#else
        if (ptr_->count.fetch_sub(1, std::memory_order_acq_rel) == 1) {
#endif
          Allocator_traits::destroy(allocator, ptr_);
          allocator.deallocate(ptr_, 1);
        }
      }
    }

    void
    initialize_with(const element_type& t)
    {
        // kept for backward compatibility.  Use operator=(t) instead.
        *this = t;
    }

    Id_type id() const noexcept { return Ptr() - static_cast<T const*>(0); }

    bool identical(const Handle_for& h) const noexcept { return Ptr() == h.Ptr(); }


    // Ptr() is the "public" access to the pointer to the object.
    // The non-const version asserts that the instance is not shared.
    const element_type *
    Ptr() const noexcept
    {
       return &(ptr_->t);
    }

    /*
    // The assertion triggers in a couple of places, so I comment it for now.
    T *
    Ptr()
    {
      CGAL_assertion(!is_shared());
      return &(ptr_->t);
    }
    */

    bool
    is_shared() const noexcept
    {
        return ptr_->count.load(std::memory_order_relaxed) > 1;
    }

    bool
    unique() const noexcept
    {
        return !is_shared();
    }

    long
    use_count() const noexcept
    {
        return ptr_->count.load(std::memory_order_relaxed);
    }

    void
    swap(Handle_for& h) noexcept
    {
      std::swap(ptr_, h.ptr_);
    }

protected:

    void
    copy_on_write()
    {
      if ( is_shared() ) Handle_for(ptr_->t).swap(*this);
    }

    // ptr() is the protected access to the pointer.  Both const and non-const.
    // Redundant with Ptr().
    element_type *
    ptr() noexcept
    { return &(ptr_->t); }

    const element_type *
    ptr() const noexcept
    { return &(ptr_->t); }
};


template <class T, class Allocator>
typename Handle_for<T, Allocator>::Allocator
Handle_for<T, Allocator>::allocator;

template <class T, class Allocator>
inline
void
swap(Handle_for<T, Allocator> &h1, Handle_for<T, Allocator> &h2)
{
    h1.swap(h2);
}

template <class T, class Allocator>
inline
bool
identical(const Handle_for<T, Allocator> &h1,
          const Handle_for<T, Allocator> &h2)
{
    return h1.identical(h2);
}

template <class T> inline bool identical(const T &t1, const T &t2) { return &t1 == &t2; }

template <class T, class Allocator>
inline
const T&
get_pointee_or_identity(const Handle_for<T, Allocator> &h)
{
    return *(h.Ptr());
}

template <class T>
inline
const T&
get_pointee_or_identity(const T &t)
{
    return t;
}

} //namespace CGAL

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#include <CGAL/enable_warnings.h>

#endif // CGAL_HANDLE_FOR_H
