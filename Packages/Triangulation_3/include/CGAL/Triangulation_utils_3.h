// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// file          : include/CGAL/Triangulation_cw_2.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================
#ifndef CGAL_TRIANGULATION_CW_2_H
#define CGAL_TRIANGULATION_CW_2_H

inline unsigned int ccw(const unsigned int i) const
{
  CGAL_precondition( i== 0 || i == 1 || i==2 );
  static const unsigned int tab_ccw[] = {1,2,0};
  return tab_ccw[i];
}

inline unsigned int cw(const unsigned int i) const
{
  CGAL_precondition( i== 0 || i == 1 || i==2 );
  static const unsigned int tab_cw[] = {2,0,1};
  return tab_cw[i];
}

#endif CGAL_TRIANGULATION_CW_2_H

