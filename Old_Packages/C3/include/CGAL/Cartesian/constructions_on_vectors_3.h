// ==========================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// --------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/constructions_on_vectors_3.h
// source        : include/CGAL/Cartesian/constructions_on_vectors_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ==========================================================================


#ifndef CGAL_CARTESIAN_CONSTRUCTIONS_ON_VECTORS_3_H
#define CGAL_CARTESIAN_CONSTRUCTIONS_ON_VECTORS_3_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#include <CGAL/Cartesian/redefine_names_3.h>
#endif

#ifndef CGAL_CARTESIAN_VECTOR_3_H
#include <CGAL/Cartesian/Vector_3.h>
#endif // CGAL_CARTESIAN_VECTOR_3_H
#ifndef CGAL_CONSTRUCTIONS_KERNEL_FTC3_H
#include <CGAL/constructions/kernel_ftC3.h>
#endif // CGAL_CONSTRUCTIONS_KERNEL_FTC3_H

CGAL_BEGIN_NAMESPACE

template < class R >
VectorC3<R CGAL_CTAG>
cross_product(const VectorC3<R CGAL_CTAG>& v,
              const VectorC3<R CGAL_CTAG>& w)
{
    return VectorC3<R CGAL_CTAG>( v.y() * w.z() - v.z() * w.y() ,
                         v.z() * w.x() - v.x() * w.z() ,
                         v.x() * w.y() - v.y() * w.x() );
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_CONSTRUCTIONS_ON_VECTORS_3_H
