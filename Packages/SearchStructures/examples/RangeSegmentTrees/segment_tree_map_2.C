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
// file          : src/examples/RangeSegmentTrees/segment_tree_map_2.C
// source        : src/examples/RangeSegmentTrees/segment_tree_map_2.C
// revision      : $Revision$
// revision_date : $Date$/
// author(s)     : Gabriele Neyer
//
// coordinator   : Peter Widmayer, ETH Zurich
//
// Implementation: Testprogram for 2-dimensional Segment Trees
// A two dimensional Segment Tree is defined in this class.
// Ti is the type of each dimension of the tree.
//
// email         : cgal@cs.uu.nl
//
// ======================================================================
#include <CGAL/basic.h>
#include <iostream>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <utility>
#include <vector>
#include <iterator>
//#include <tempbuf.h>
#include <list>
#include <CGAL/Segment_tree_k.h>
#include <CGAL/Range_segment_tree_traits.h>



typedef CGAL::Cartesian<double> Representation;
typedef CGAL::Segment_tree_map_traits_2<Representation, char> Traits;
typedef CGAL::Segment_tree_2<Traits > Segment_tree_2_type;

int main()
{
  typedef Traits::Interval Interval;
  typedef Traits::Pure_interval Pure_interval;
  typedef Traits::Key Key;
  // definition of the two-dimensional segment tree
  std::list<Interval> InputList, OutputList, N;



  // insertion of the tree elements into the sequence container
  InputList.push_back(Interval(Pure_interval(Key(1,5), Key(2,7)),'a'));
  InputList.push_back(Interval(Pure_interval(Key(2,7), Key(3,8)),'b'));
  InputList.push_back(Interval(Pure_interval(Key(6,9), Key(9,13)),'c'));
  InputList.push_back(Interval(Pure_interval(Key(1,3), Key(3,9)),'d'));
 
  // creation of the segment tree
  std::list<Interval>::iterator first = InputList.begin();
  std::list<Interval>::iterator last = InputList.end();

  Segment_tree_2_type Segment_tree_2(first,last);

  // perform a window query
  Interval a=Interval(Pure_interval(Key(3,6), Key(7,12)),'e');
  Segment_tree_2.window_query(a,std::back_inserter(OutputList));

  // output of the querey elements on stdout
  std::list<Interval>::iterator j = OutputList.begin();
  std::cerr << "\n window_query (3,6),(7,12)\n";
  while(j!=OutputList.end())
  {
    std::cerr << (*j).first.first.x() << "-" << (*j).first.second.x() << " " 
	 << (*j).first.first.y() << "-" << (*j).first.second.y() << std::endl; 
    j++;
  }


  
  Interval b=Interval(Pure_interval(Key(6,10),Key(7,11)), 'f');
  Segment_tree_2.enclosing_query(b,std::back_inserter(N));
  j = N.begin();
  std::cerr << "\n enclosing_query (6,10),(7,11)\n";
  while(j!=N.end())
  {
    std::cerr << (*j).first.first.x() << "-" << (*j).first.second.x() << " " 
	 << (*j).first.first.y() << "-" << (*j).first.second.y() << std::endl; 
    j++;
  }
  if(Segment_tree_2.CSegment_tree_2->is_valid())
    std::cerr << "Tree  is valid\n";
  else
    std::cerr << "Tree is not valid\n";
  return 0; 
}


