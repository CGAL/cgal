// Copyright (c) 1999  
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
// Author(s)     : Stefan Schirra
 

#ifndef CGAL_HANDLE_FOR_VIRTUAL_H
#define CGAL_HANDLE_FOR_VIRTUAL_H

#include <CGAL/config.h>
#include <typeinfo>
#include <cstddef>

namespace CGAL {

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

    typedef std::ptrdiff_t Id_type ;
    
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

#ifndef CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE
    Handle_for_virtual&
    operator=( Handle_for_virtual && h)
    {
      swap(h);
      return *this;
    }
#endif

// protected:
    typedef RefCounted element_type;

    template <class T>
    void
    initialize_with( const T& rc)
    {
	ptr = new T(rc);
    }

    Id_type id() const { return Ptr() - static_cast<RefCounted const*>(0); }
    
    bool identical( const Handle_for_virtual& h) const { return Ptr() == h.Ptr(); }


    void
    swap(Handle_for_virtual & h)
    { std::swap(h.ptr, ptr); }

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

    RefCounted * ptr;
};

template <class RefCounted>
inline bool identical(const Handle_for_virtual<RefCounted> &h1, const Handle_for_virtual<RefCounted> &h2) { return h1.identical(h2); }

} //namespace CGAL

#endif // CGAL_HANDLE_FOR_VIRTUAL_H
