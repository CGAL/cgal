// Copyright (c) 1999,2001,2003  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stefan Schirra, Sylvain Pion
 
#ifndef CGAL_HANDLE_FOR_H
#define CGAL_HANDLE_FOR_H

#include <CGAL/config.h>
#include <CGAL/assertions.h> // for CGAL_assume

#include <boost/config.hpp>
#include <CGAL/memory.h>
#include <algorithm>
#include <cstddef>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4345) // Avoid warning  http://msdn.microsoft.com/en-us/library/wewb47ee(VS.80).aspx
#endif
namespace CGAL {

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
    
    typedef std::ptrdiff_t Id_type ;

    Handle_for()
    {
        pointer p = allocator.allocate(1);
        new (&(p->t)) element_type(); // we get the warning here
        p->count = 1;
        ptr_ = p;
    }

    Handle_for(const element_type& t)
    {
        pointer p = allocator.allocate(1);
        new (&(p->t)) element_type(t);
        p->count = 1;
        ptr_ = p;
    }

#ifndef CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE
    Handle_for(element_type && t)
    {
        pointer p = allocator.allocate(1);
        new (&(p->t)) element_type(std::move(t));
        p->count = 1;
        ptr_ = p;
    }
#endif

/* I comment this one for now, since it's preventing the automatic conversions
   to take place.  We'll see if it's a problem later.
    template < typename T1 >
    Handle_for(const T1& t1)
    {
        pointer p = allocator.allocate(1);
        new (&(p->t)) T(t1);
        p->count = 1;
        ptr_ = p;
    }
*/

#if !defined CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES && !defined CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE
    template < typename T1, typename T2, typename... Args >
    Handle_for(T1 && t1, T2 && t2, Args && ... args)
    {
        pointer p = allocator.allocate(1);
        new (&(p->t)) element_type(std::forward<T1>(t1), std::forward<T2>(t2), std::forward<Args>(args)...);
        p->count = 1;
        ptr_ = p;
    }
#else
    template < typename T1, typename T2 >
    Handle_for(const T1& t1, const T2& t2)
    {
        pointer p = allocator.allocate(1);
        new (&(p->t)) element_type(t1, t2);
        p->count = 1;
        ptr_ = p;
    }

    template < typename T1, typename T2, typename T3 >
    Handle_for(const T1& t1, const T2& t2, const T3& t3)
    {
        pointer p = allocator.allocate(1);
        new (&(p->t)) element_type(t1, t2, t3);
        p->count = 1;
        ptr_ = p;
    }

    template < typename T1, typename T2, typename T3, typename T4 >
    Handle_for(const T1& t1, const T2& t2, const T3& t3, const T4& t4)
    {
        pointer p = allocator.allocate(1);
        new (&(p->t)) element_type(t1, t2, t3, t4);
        p->count = 1;
        ptr_ = p;
    }
#endif // CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

    Handle_for(const Handle_for& h)
      : ptr_(h.ptr_)
    {
	CGAL_assume (ptr_->count > 0);
        ++(ptr_->count);
    }

    Handle_for&
    operator=(const Handle_for& h)
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

#ifndef CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE
    // Note : I don't see a way to make a useful move constructor, apart
    //        from e.g. using NULL as a ptr value, but this is drastic.

    Handle_for&
    operator=(Handle_for && h)
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
#endif

    ~Handle_for()
    {
      if (--(ptr_->count) == 0) {
          allocator.destroy( ptr_);
          allocator.deallocate( ptr_, 1);
      }
    }

    void
    initialize_with(const element_type& t)
    {
        // kept for backward compatibility.  Use operator=(t) instead.
        *this = t;
    }

    Id_type id() const { return Ptr() - static_cast<T const*>(0); }
    
    bool identical(const Handle_for& h) const { return Ptr() == h.Ptr(); }


    // Ptr() is the "public" access to the pointer to the object.
    // The non-const version asserts that the instance is not shared.
    const element_type *
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

    bool
    unique() const
    {
	return !is_shared();
    }

    long
    use_count() const
    {
	return ptr_->count;
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
      if ( is_shared() ) Handle_for(ptr_->t).swap(*this);
    }

    // ptr() is the protected access to the pointer.  Both const and non-const.
    // Redundant with Ptr().
    element_type *
    ptr()
    { return &(ptr_->t); }

    const element_type *
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

#endif // CGAL_HANDLE_FOR_H
