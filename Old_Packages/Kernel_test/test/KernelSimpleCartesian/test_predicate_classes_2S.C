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
// source        :
// file          : test_predicate_classes_2S.C
// revision      : 3.8
// revision_date : 08 Oct 2000 
// author(s)     : Stefan Schirra
//
// maintainer    : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de> 
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <CGAL/basic.h>
#include <cassert>
#ifdef CGAL_USE_GMP
# include <CGAL/Gmpz.h>
typedef CGAL::Gmpz    Precise_integer;
#else
# ifdef CGAL_USE_LEDA
#  include <CGAL/leda_integer.h>
typedef leda_integer  Precise_integer;
# endif // CGAL_USE_LEDA
#endif // CGAL_USE_GMP


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/Point_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/_test_fctobj_point_line_2.C>

int
main()
{
  typedef   CGAL::Simple_cartesian<CGAL::Quotient<Precise_integer> >     C_Cls;
  typedef   CGAL::Point_2<C_Cls>  C_Point;
  typedef   CGAL::Line_2<C_Cls>   C_Line;
  std::cout << "Testing 2d with Simple_cartesian<Q<Precise_integer> > :";
  std::cout << std::endl;
  _test_fctobj_point_line_2( C_Point(), C_Line());

  return 0;
}
