// ======================================================================
//
// Copyright (c) 1999,2001 The CGAL Consortium
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

// There are basically 2 ways of constructing an object deriving from
// Handle_for :
// - call the default constructor of Handle_for<T>, then eventually use
//   initialize_with(const T&), which uses the assignment.
// - call the constructor
//   Handle_for(Handle_for::TO_BE_USED_ONLY_WITH_CONSTRUCT_WITH), then call
//   construct_with(const T&).  You HAVE to call it.


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

    struct TO_BE_USED_ONLY_WITH_CONSTRUCT_WITH {};

    Handle_for()
    {
	CGAL_assertion_code(ptr_ = NULL;)
	construct_with(T());
    }

    Handle_for(TO_BE_USED_ONLY_WITH_CONSTRUCT_WITH)
    {
	CGAL_assertion_code(ptr_ = NULL;)
    }

    Handle_for(const T& t)
    {
	CGAL_assertion_code(ptr_ = NULL;)
	construct_with(t);
    }

    Handle_for(const Handle_for& h)
    {
        ptr_ = h.ptr_;
        ptr_->add_reference();
    }

    ~Handle_for()
    {
	remove_reference();
    }

    Handle_for&
    operator=(const Handle_for& h)
    {
        h.ptr_->add_reference();
	remove_reference();
        ptr_ = h.ptr_;
        return *this;
    }

    void
    initialize_with(const T& t)
    {
	*ptr() = t;
    }

    void
    construct_with(const T& t)
    {
	CGAL_assertion(ptr_ == NULL);
        ptr_ = allocator.allocate(1);
        allocator.construct(ptr_, RefCounted(t));
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

CGAL_END_NAMESPACE

#endif // CGAL_HANDLE_FOR_H
