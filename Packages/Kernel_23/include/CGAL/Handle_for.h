// Copyright (c) 1999,2001,2003  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Stefan Schirra, Sylvain Pion
 
#ifndef CGAL_HANDLE_FOR_H
#define CGAL_HANDLE_FOR_H

#include <CGAL/memory.h>
#include <algorithm>

CGAL_BEGIN_NAMESPACE

template <class T, class Alloc = CGAL_ALLOCATOR(T) >
class Handle_for
{
    // Wrapper that adds the reference counter.
    struct RefCounted {
        T t;
        unsigned int count;
    };

    typedef typename Alloc::template rebind<RefCounted>::other  Allocator;
    typedef typename Allocator::pointer                         pointer;

    static Allocator   allocator;
    pointer            ptr_;

public:

    typedef T element_type;

    Handle_for()
    {
        // Use a unique static instance to speed up default construction.
        // It's a static variable of a function instead of the class to
        // avoid the requirement of a default constructor for T().
        static const Handle_for def = Handle_for(T());
        ptr_ = def.ptr_;
        ++(ptr_->count);
    }

    // TODO :
    // We should also think about providing template constructors in
    // order to forward the functionality of T to Handle_for<T> without
    // the need to an intermediate copy.
    // Currently it's not working, because some places use conversions.

    Handle_for(const T& t)
      : ptr_(allocator.allocate(1))
    {
        new (&(ptr_->t)) T(t);
        ptr_->count = 1;
    }

    Handle_for(const Handle_for& h)
      : ptr_(h.ptr_)
    {
        ++(ptr_->count);
    }

    ~Handle_for()
    {
      if (! is_shared() ) {
          allocator.destroy( ptr_);
          allocator.deallocate( ptr_, 1);
      }
      else
	  --(ptr_->count);
    }

    Handle_for&
    operator=(const Handle_for& h)
    {
        Handle_for tmp = h;
        swap(tmp);
        return *this;
    }

    Handle_for&
    operator=(const T &t)
    {
        if (is_shared())
            *this = Handle_for(t);
        else
            ptr_->t = t;

        return *this;
    }

    void
    initialize_with(const T& t)
    {
        // kept for backward compatibility.  Use operator=(t) instead.
        *this = t;
    }

    bool
    identical(const Handle_for& h) const
    { return ptr_ == h.ptr_; }

    long int
    id() const
    { return reinterpret_cast<long int>(&*ptr_); }

    // Ptr() is the "public" access to the pointer to the object.
    // The non-const version asserts that the instance is not shared.
    const T *
    Ptr() const
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
    is_shared() const
    {
	return ptr_->count > 1;
    }

    void
    swap(Handle_for& h)
    {
      std::swap(ptr_, h.ptr_);
    }

protected:

    void
    copy_on_write()
    {
      if ( is_shared() )
      {
        pointer tmp_ptr = allocator.allocate(1);
        new (&(tmp_ptr->t)) T(ptr_->t);
        tmp_ptr->count = 1;
        --(ptr_->count);
        ptr_ = tmp_ptr;
      }
    }

    // ptr() is the protected access to the pointer.  Both const and non-const.
    // Redundant with Ptr().
    T *
    ptr()
    { return &(ptr_->t); }

    const T *
    ptr() const
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

CGAL_END_NAMESPACE

#endif // CGAL_HANDLE_FOR_H
