
// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : Ieee_float.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#include <iomanip>

union IEEE_fpst_float {
  float    a;
  struct { unsigned sign : 1;
           unsigned exp  : 8;
           unsigned high :23;
         } b;
  struct { unsigned L    :32;
         } c;
};         

void
show( IEEE_fpst_float* p)
{ 
  std::cout << std::endl;
  std::cout << std::hex << std::setw(8) << std::setfill('0') << p->c.L;
  std::cout << std::endl;
}
