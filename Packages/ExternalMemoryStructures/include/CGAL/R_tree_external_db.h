// ======================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 2000, August 11
//
// file          : include/CGAL/ExternalMemoryStructures/R_tree_external_db.h
// package       : ExternalMemoryStructures (0.631)
// maintainer    : Philipp Kramer <kramer@inf.ethz.ch>
// chapter       : $CGAL_Chapter: Basic / External Data Structures $
// source        : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Gabriele Neyer<neyer@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Peter Widmayer <widmayer@inf.ethz.ch>)
//
// Implementation of the R_tree_external_db.
// ======================================================================
#ifndef __R_tree_external_db_H__
#define __R_tree_external_db_H__
#include<CGAL/IO_tree_traits.h>
#include<CGAL/cache.h>

CGAL_BEGIN_NAMESPACE

template<const unsigned int buffersize>
class R_tree_external_db:public CGAL::IO_tree_traits
                        <CGAL::cache<buffersize> >{};

CGAL_END_NAMESPACE
#endif
