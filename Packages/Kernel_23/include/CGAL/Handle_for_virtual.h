// Copyright (c) 1999  Utrecht University (The Netherlands),
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
// Author(s)     : Stefan Schirra
 

#ifndef CGAL_HANDLE_FOR_VIRTUAL_H
#define CGAL_HANDLE_FOR_VIRTUAL_H

#include <CGAL/basic.h>
#include <typeinfo>

CGAL_BEGIN_NAMESPACE

class Ref_counted_virtual
{
  public:
    Ref_counted_virtual() : count(1) {}
    Ref_counted_virtual(const Ref_counted_virtual&) : count(1) {}

    void  add_reference() { ++count; }
    void  remove_reference() { --count; }
    bool  is_referenced() const { return (count != 0); }
    bool  is_shared() const { return (count > 1); }

    virtual const std::type_info & type() const
    { return typeid(void); }

    virtual const void * object_ptr() const
    { return NULL; }

    virtual ~Ref_counted_virtual() {}

  protected:
    unsigned int count;
};


template <class RefCounted>
// RefCounted must provide
// add_reference()
// remove_reference()
// bool is_referenced() const
// bool is_shared() const
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

    const void * object_ptr() const
    {
      return ptr->object_ptr();
    }

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
