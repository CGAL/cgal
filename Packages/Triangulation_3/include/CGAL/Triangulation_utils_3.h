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
// file          : include/CGAL/Triangulation_utils_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis 
//                 (Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_UTILS_3_H
#define CGAL_TRIANGULATION_UTILS_3_H

CGAL_BEGIN_NAMESPACE

class Triangulation_utils_3
{
public:
//   inline unsigned int ccw(const unsigned int i) const
//     {
//       CGAL_precondition( i < 3 );
//       static const unsigned int tab_ccw[] = {1,2,0};
//       return tab_ccw[i];
//     }
  
//   inline unsigned int cw(const unsigned int i) const
//     {
//       CGAL_precondition( i < 3 );
//       static const unsigned int tab_cw[] = {2,0,1};
//       return tab_cw[i];
//     }

  inline unsigned int ccw(const unsigned int i) const
    {
      CGAL_triangulation_precondition( ((unsigned int) i) < 3 );
      return (i==2) ? 0 : i+1;
    }
  
  inline unsigned int cw(const unsigned int i) const
    {
      CGAL_triangulation_precondition( ((unsigned int) i) < 3 );
      return (i==0) ? 2 : i-1;
    }

  inline unsigned int nextposaround(const unsigned int i, 
				    const unsigned int j) const
{
  // index of the next cell when turning around the
  // oriented edge vertex(i) vertex(j) in 3d
  CGAL_precondition( ( i < 4 ) && ( j < 4 ) && ( i != j ) );
  static const unsigned int tab_nextposaroundij[4][4] = {
    {5, 2, 3, 1},
    {3, 5, 0, 2},
    {1, 3, 5, 0},
    {2, 0, 1, 5} }; // the diagonal has no meaning
  return tab_nextposaroundij[i][j];
}

};

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_UTILS_3_H
