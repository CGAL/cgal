// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Handle_hash_function.h
// package       : Hash_map
// chapter       : STL_Extensions
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
//                 Lutz Kettner <kettner@inf.ethz.ch>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: Hash Map by chained hashing, default hash function
// ============================================================================

#ifndef CGAL_HANDLE_HASH_FUNCTION_H
#define CGAL_HANDLE_HASH_FUNCTION_H

#include <CGAL/basic.h>
#include <cstddef>

CGAL_BEGIN_NAMESPACE

struct Handle_hash_function {
    typedef std::size_t result_type;
    template <class H> 
    std::size_t operator() (const H& h) const { 
        return std::size_t(&*h) / 
            sizeof( typename std::iterator_traits<H>::value_type);
    } 
};

CGAL_END_NAMESPACE

#endif // CGAL_HANDLE_HASH_FUNCTION_H
// EOF

