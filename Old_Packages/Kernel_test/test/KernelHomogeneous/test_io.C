// ============================================================================
//
// Copyright (c) 2001 The CGAL Consortium
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
// file          : test_io.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
//
// coordinator   : MPI, Saarbruecken
// ============================================================================
 

#include <CGAL/basic.h>
#include <cassert>

#include <CGAL/Homogeneous.h>
#include <CGAL/Quotient.h>
#include "../Kernel/include/CGAL/_test_io.h"

int
main()
{
  typedef   CGAL::Homogeneous<double>     Cls;
  std::cout << "Testing IO with Homogeneous<double> :" << std::endl;
  _test_io( Cls() );
  return 0;
}
