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

inline unsigned int ccw(const unsigned int i) 
{
  CGAL_precondition( i < 3 );
  static const unsigned int tab_ccw[] = {1,2,0};
  return tab_ccw[i];
}

inline unsigned int cw(const unsigned int i) 
{
  CGAL_precondition( i < 3 );
  static const unsigned int tab_cw[] = {2,0,1};
  return tab_cw[i];
}

inline unsigned int nextposaroundij(const unsigned int i, 
				    const unsigned int j)
{
  // index of the next cell when turning around the
  // oriented edge vertex(i) vertex(j) in 3d
  CGAL_precondition( ( i < 4 ) && ( j < 4 ) && ( i != j ) );
  static const unsigned int tab_nextposaroundij[4][4] = {
    5, 2, 3, 1,
    3, 5, 0, 2,
    1, 3, 5, 0,
    2, 0, 1, 5 }; // the diagonal has no meaning
  return tab_nextposaroundij[i][j];
}

#endif CGAL_TRIANGULATION_CW_2_H

