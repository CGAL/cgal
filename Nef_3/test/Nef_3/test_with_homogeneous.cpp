// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : test/Nef_3/nef_3.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Id$
// revision_date : $Date$
//
// author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
// maintainer    : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// ============================================================================

#define CGAL_NEF3_SORT_OUTPUT 1

#include <CGAL/basic.h>
#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer NT;
#else
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz NT;
#endif
#include <CGAL/Homogeneous.h>
#include <CGAL/Timer.h>
#include <CGAL/test_Nef_3.h>

int main() {
  typedef CGAL::Homogeneous<NT>              H_kernel;
  
#ifdef CGAL_CFG_ISTREAM_INT_BUG
  std::locale::global(std::locale("C")); 
#endif

  CGAL::test_Nef_3<H_kernel>  test_H;

  CGAL::Timer t;
  t.start();
  test_H.run_test();
  t.stop();
  std::cout << "Homogeneous kernel successful in: " << t.time() 
            << " seconds " << std::endl;
}

