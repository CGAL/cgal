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
// file          : d_tuple.h
// package       : _d
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sven Schoenherr
//                 Bernd Gaertner
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL__D_TUPLE_H
#define CGAL__D_TUPLE_H

CGAL_BEGIN_NAMESPACE

template < class T >
class _d_tuple : public Rep
{
    public:
    const int d;
    T* e;

    _d_tuple(int dim = 0, bool cartesian = true) : d(dim)
    {
    if (cartesian)
        e = new T[d];
    else
        e = new T[d+1];
    }

    ~_d_tuple ()
    {
        delete[] e;
    }
};
CGAL_END_NAMESPACE


#endif // CGAL__D_TUPLE_H
