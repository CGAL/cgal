// ============================================================================
//
// Copyright (c) 1999,2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : test_kernel_archetype_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Matthias Baesken
//
//
// coordinator   : MPI, Saarbruecken
// ============================================================================
 

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

#define CGAL_NO_DEPRECATED_CODE
#define CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS

#include <CGAL/Kernel_archetype.h>

// needed in kernel testsuite ...
CGAL::Test_vector_2 operator-(const CGAL::Test_vector_2& v)
{ return v; }  

#include "CGAL/_test_new_2.h"

typedef CGAL::Kernel_archetype   Kernel;

int main()
{
  test_new_2( Kernel() );
  return 0;
}

