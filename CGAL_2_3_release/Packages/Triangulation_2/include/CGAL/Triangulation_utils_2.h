// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
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
// file          : include/CGAL/Triangulation_utils_2.h
// revision      : $Revision$
// author(s)     : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//                 Sylvain Pion    <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================
#ifndef CGAL_TRIANGULATION_UTILS_2_H
#define CGAL_TRIANGULATION_UTILS_2_H

CGAL_BEGIN_NAMESPACE 

class Triangulation_cw_ccw_2
{
public:
  int ccw(const int i) const
    {
      CGAL_triangulation_precondition( ((unsigned int) i) < 3 );
      return (i==2) ? 0 : i+1;
    }

  int cw(const int i) const
    {
      CGAL_triangulation_precondition( ((unsigned int) i) < 3 );
      return (i==0) ? 2 : i-1;
    }
};

CGAL_END_NAMESPACE

#endif //CGAL_TRIANGULATION_UTILS_2_H
