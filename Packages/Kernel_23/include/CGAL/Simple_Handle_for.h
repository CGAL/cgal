// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
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
// file          : include/CGAL/Simple_Handle_for.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
// coordinator   : MPI, Saarbruecken  (<Susan.Hert@mpi-sb.mpg.de>)
//
// ======================================================================

#ifndef CGAL_SIMPLE_HANDLE_FOR_H
#define CGAL_SIMPLE_HANDLE_FOR_H

CGAL_BEGIN_NAMESPACE

template < class Stored >
class Simple_Handle_for
{
public:

    Simple_Handle_for()
    {}

    Simple_Handle_for(const Stored& rc)
	: _s(rc) {}

    typedef Stored element_type;

    void
    initialize_with(const Stored& rc)
    {
      _s = rc;
    }

    long int
    id() const
    { return reinterpret_cast<long int>(&_s); }

    bool
    identical(const Simple_Handle_for& h) const
    { return id() == h.id(); } // Or should it always return false ?

    const Stored * Ptr() const
    { return &_s; }

    Stored * Ptr()
    { return &_s; }

    const Stored * ptr() const
    { return &_s; }

    Stored * ptr()
    { return &_s; }

    bool
    is_shared() const
    {
	return false;
    }

private:
    void
    copy_on_write()
    {}

    Stored _s;
};

CGAL_END_NAMESPACE

#endif // CGAL_SIMPLE_HANDLE_FOR_H
