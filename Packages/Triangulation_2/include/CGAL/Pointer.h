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
// release       : $CGAL_Revision: CGAL $
// date          :
//
// file          : include/CGAL/Pointer.h
// source        : ~yvinec/Cgal/Triangulation/include/CGAL/RCS/Pointer.h,v
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Olivier Devillers, Mariette Yvinec
//
// coordinator   : Mariette Yvinec  < Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================



#ifndef CGAL_POINTER_H
#define CGAL_POINTER_H

#include <assert.h>
#include <iterator.h>
#include <CGAL/circulator.h>

template<class T>
class CGAL_Pointer
{

    private:
    T* _pointer;
    
    public:
    typedef T value_type;
    typedef CGAL_Pointer<T> Pointer;
    inline CGAL_Pointer()
        : _pointer(NULL)
    {}
    
    inline CGAL_Pointer(const T* p)
        : _pointer((T*)p)
    {}
    
    inline Pointer& operator=(const T*& p)
    {
        ptr() = p ;
        return *this;
    }
    
    inline Pointer& operator=(const Pointer& p)
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
    
    
    inline bool operator==(const Pointer& p) const
    {
        return ( ptr() == p.ptr() );
    }
    
    inline bool operator!=(const Pointer& p) const
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

#endif // CGAL_POINTER_H
