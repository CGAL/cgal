// ======================================================================
//
// Copyright (c) 1999,2001,2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : Handle_for.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra, Sylvain Pion
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 
#ifndef CGAL_HANDLE_FOR_H
#define CGAL_HANDLE_FOR_H

#include <CGAL/memory.h>

CGAL_BEGIN_NAMESPACE

template <class T, class Alloc = CGAL_ALLOCATOR(T) >
class Handle_for
{
    // Wrapper that adds the reference counter.
    struct RefCounted
    {
        RefCounted(const T& t) : t_(t) {}
        RefCounted(const RefCounted& r) : t_(r.t_), count(1) {}

        T* base_ptr() { return &t_; }
        const T& t() const { return t_; }

	// Speeds things up with LEDA and New_delete_allocator<>.
        CGAL_MEMORY(RefCounted)

        void add_reference()    { ++count; }
        void remove_reference() { --count; }
        bool is_shared() const  { return count > 1; }

    private:

        T t_;
        unsigned int count;
    };

    typedef typename Alloc::template rebind<RefCounted>::other  Allocator;

  public:

    typedef T element_type;

    Handle_for()
    {
        // Use a unique static instance to speed up default construction.
        // It's a static variable of a function instead of the class to
        // avoid the requirement of a default constructor for T().
        static const Handle_for def = Handle_for(T());
        ptr_ = def.ptr_;
        ptr_->add_reference();
    }

    Handle_for(const T& t)
      : ptr_(allocator.allocate(1))
    {
        allocator.construct(ptr_, RefCounted(t));
    }

    Handle_for(const Handle_for& h)
      : ptr_(h.ptr_)
    {
        ptr_->add_reference();
    }

    ~Handle_for()
    {
	remove_reference();
    }

    Handle_for&
    operator=(const Handle_for& h)
    {
        Handle_for tmp(h);
        swap(tmp);
        return *this;
    }

    Handle_for&
    operator=(const T &t)
    {
        if (is_shared())
            *this = Handle_for(t);
        else
            *ptr() = t;

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

    // Ptr() is the "public" access to the pointer.  Both const and non-const.
    // non-const does copy-on-write.
    const T *
    Ptr() const
    { return ptr_->base_ptr(); }

    /*
    T *
    Ptr()
    {
	copy_on_write();
	return ptr_;
    }
    */

    bool
    is_shared() const
    {
	return ptr_->is_shared();
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
        RefCounted* tmp_ptr = allocator.allocate(1);
        allocator.construct( tmp_ptr, *ptr_);
        ptr_->remove_reference();
        ptr_ = tmp_ptr;
      }
    }

    // ptr() is the protected access to the pointer.  Both const and non-const.
    T *
    ptr()
    { return ptr_->base_ptr(); }

    const T *
    ptr() const
    { return ptr_->base_ptr(); }

private:

    void
    remove_reference()
    {
      if (! is_shared() ) {
          allocator.destroy( ptr_);
          allocator.deallocate( ptr_, 1);
      }
      else
	  ptr_->remove_reference();
    }

    static Allocator                 allocator;
    typename Allocator::pointer      ptr_;
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
