// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  :         :
//
// file          : include/CGAL/Pointer.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Olivier Devillers, Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================



#ifndef CGAL_POINTER_H
#define CGAL_POINTER_H

#include <cassert>
#include <iterator>
#include <CGAL/circulator.h>

CGAL_BEGIN_NAMESPACE 

template<class T>
class Pointer
{

    private:
    T* _pointer;
    
    public:
    typedef T value_type;
  //typedef Pointer<T> Pointer;

    inline Pointer<T>()
        : _pointer(NULL)
    {}
    
    inline Pointer<T>(const T* p)
        : _pointer((T*)p)
    {}
    
    inline Pointer<T>& operator=(const T*& p)
    {
        ptr() = p ;
        return *this;
    }
    
    inline Pointer<T>& operator=(const Pointer<T>& p)
    {
        ptr() = p.ptr();
        return *this;
    }
    inline T& operator*() const
    {
        return *ptr();
    }
    
    inline T* operator->() const
    {
        return ptr();
    }
    inline void Delete()
    {
        delete ptr();
        clear();
    }
    
    inline void clear()
    {
        ptr() = NULL;
    }
    
    
    inline bool operator==(const Pointer<T>& p) const
    {
        return ( ptr() == p.ptr() );
    }
    
    inline bool operator!=(const Pointer<T>& p) const
    {
        return !(*this == p);
    }
    
  inline bool is_null()
  {
    return (ptr() == NULL );
  }


  inline bool operator==(CGAL_NULL_TYPE n) const
    {
        assert( n == 0);
        return ( ptr() == NULL );
    }
    
    inline bool operator!=(CGAL_NULL_TYPE n) const
    {
        return !(*this == n);
    }

 
    public:
    inline T*& ptr()       {return _pointer;}
    inline T*  ptr() const {return _pointer;}
    

};

CGAL_END_NAMESPACE

#endif // CGAL_POINTER_H
