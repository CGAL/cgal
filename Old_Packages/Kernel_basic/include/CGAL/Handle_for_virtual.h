// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : Handle_for_virtual.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 

#ifndef CGAL_HANDLE_FOR_VIRTUAL_H
#define CGAL_HANDLE_FOR_VIRTUAL_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

class Ref_counted_virtual
{
  public:
    Ref_counted_virtual() : count(1) {}
    Ref_counted_virtual(const Ref_counted_virtual&) : count(1) {}

    void  add_reference() { ++count; }
    void  remove_reference() { --count; }
    bool  is_referenced() { return (count != 0); }
    bool  is_shared() { return (count > 1); }

    virtual ~Ref_counted_virtual() {}

  protected:
    unsigned int count;
};


template <class RefCounted>
// RefCounted must provide
// add_reference()
// remove_reference()
// bool is_referenced()
// bool is_shared()
// and initialize count to 1 in default and copy constructor
class Handle_for_virtual
{
  public:

    Handle_for_virtual(const RefCounted& rc)
    {
      ptr = new RefCounted(rc);
    }

    Handle_for_virtual()
    {
      ptr = NULL;
    }

    Handle_for_virtual( const Handle_for_virtual& h)
    {
      ptr = h.ptr;
      ptr->add_reference();
    }

    ~Handle_for_virtual()
    {
      ptr->remove_reference();
      if ( !ptr->is_referenced() )
	  delete ptr;
    }

    Handle_for_virtual&
    operator=( const Handle_for_virtual& h)
    {
      h.ptr->add_reference();
      ptr->remove_reference();
      if ( !ptr->is_referenced() )
	  delete ptr;
      ptr = h.ptr;
      return *this;
    }

// protected:
    typedef RefCounted element_type;

    template <class T>
    void
    initialize_with( const T& rc)
    {
	ptr = new T(rc);
    }

    bool
    identical( const Handle_for_virtual& h) const
    { return ptr == h.ptr; }

    long int
    id() const
    { return reinterpret_cast<long int>(&*ptr); }

    const RefCounted *
    Ptr() const
    { return ptr; }

    /*
    T *
    Ptr()
    {
	copy_on_write();
	return ptr;
    }
    */

    /*
// private:
    void
    copy_on_write()
    {
      if ( ptr->is_shared() )
      {
        RefCounted* tmp_ptr = allocator.allocate(1);
        allocator.construct( tmp_ptr, *ptr);
        ptr->remove_reference();
        ptr = tmp_ptr;
      }
    }
    */
protected:

    /*
    RefCounted * ptr() const
    {
	return ptr;
    }
    */

private:

    RefCounted * ptr;
};

CGAL_END_NAMESPACE

#endif // CGAL_HANDLE_FOR_VIRTUAL_H
