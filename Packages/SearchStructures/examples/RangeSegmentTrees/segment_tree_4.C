// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of the
// Computational Geometry Algorithms Library (CGAL).
//
// Every use of CGAL requires a license. Licenses come in three kinds:
//
// - For academic research and teaching purposes, permission to use and
//   copy the software and its documentation is hereby granted free of  
//   charge, provided that
//   (1) it is not a component of a commercial product, and
//   (2) this notice appears in all copies of the software and
//       related documentation.
// - Development licenses grant access to the source code of the library 
//   to develop programs. These programs may be sold to other parties as 
//   executable code. To obtain a development license, please contact
//   the CGAL Consortium (at cgal@cs.uu.nl).
// - Commercialization licenses grant access to the source code and the
//   right to sell development licenses. To obtain a commercialization 
//   license, please contact the CGAL Consortium (at cgal@cs.uu.nl).
//
// This software and documentation is provided "as-is" and without
// warranty of any kind. In no event shall the CGAL Consortium be
// liable for any damage of any kind.
//
// The CGAL Consortium consists of Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Free University of Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany) Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
//
// release       : CGAL-1.2
// release_date  : 1999, January 18
//
// file          : src/examples/RangeSegmentTrees/segment_tree_4.C
// source        : src/examples/RangeSegmentTrees/segment_tree_4.C
// revision      : $Revision$
// revision_date : $Date$/
// author(s)     : Gabriele Neyer
//
// coordinator   : Peter Widmayer, ETH Zurich
//
// Implementation: Testprogram for 4-dimensional Segment Trees
// A four dimensional Segment Tree is defined in this class.
// Ti is the type of each dimension of the tree.
//
// email         : cgal@cs.uu.nl
//
// ======================================================================
#include <CGAL/basic.h>
#include <iostream>
#include <CGAL/Segment_tree_k.h>
#include "include/Tree_Traits.h"
#include <vector>
#include <iterator>
#include <list>

typedef CGAL::Segment_tree_4<CGAL::Tree_traits_4> Segment_tree_4_type;

int main()
{

  typedef CGAL::Tree_traits_4::Interval Interval;
  typedef CGAL::Tree_traits_4::Key Key;
  // definition of the two-dimensional segment tree
  std::list<Interval> InputList, OutputList, N;

  // insertion of the tree elements into the sequence container
  InputList.push_back(Interval(Key(1, 4, 5, 2.5), Key(2,5,7, 4.7)));
  InputList.push_back(Interval(Key(2,6.7,5, 3.5), Key(4,6.9, 8, 7.9)));
  InputList.push_back(Interval(Key(2,4.55, 8, 5.5), Key(5, 7.88, 10, 9.9)));
  InputList.push_back(Interval(Key(2, 4.66, 5, 4.6), Key(6, 8.99, 8, 8.4)));
 
  // creation of the segment tree
  std::list<Interval>::iterator first = InputList.begin();
  std::list<Interval>::iterator last = InputList.end();

  Segment_tree_4_type Segment_tree_4(first,last);

  // perform a window query
  Interval a(Key(3,4.88, 6, 5.4), Key(6, 8.999, 9, 8.6));
  Segment_tree_4.window_query(a,std::back_inserter(OutputList));

  // output of the querey elements on stdout
  std::list<Interval>::iterator j = OutputList.begin();
  std::cerr << "\n window_query (3,4.88, 6, 5.4),(6, 8.999, 9, 8.6)\n";
  while(j!=OutputList.end())
  {
    std::cerr << (*j).first.key_1 << "," << (*j).first.key_2 << ", " 
	 << (*j).first.key_3 << ", " << (*j).first.key_4 << "-" 
	 << (*j).second.key_1 << "," << (*j).second.key_2 << ", " 
	 << (*j).second.key_3 << " , " <<  (*j).second.key_4 << std::endl;
    j++;
  }
  Interval b(Key(2,6.8,9, 5.9),Key(3,7,10, 6.7));
  Segment_tree_4.enclosing_query(b,std::back_inserter(N));
  j = N.begin();
  std::cerr << "\n enclosing_query (2,6.8,9, 5.9), (3,7,10, 6.7)\n";
  while(j!=N.end())
  {
    std::cerr << (*j).first.key_1 << "," << (*j).first.key_2 << ", " 
	 << (*j).first.key_3 << ", " << (*j).first.key_4 << "-" 
	 << (*j).second.key_1 << "," << (*j).second.key_2 << ", " 
	 << (*j).second.key_3 << " , " <<  (*j).second.key_4 << std::endl;
    j++;
  }
  if(Segment_tree_4.CSegment_tree_4->is_valid())
    std::cerr << "Tree  is valid\n";
  else
    std::cerr << "Tree is not valid\n";
  return 0; 
}


