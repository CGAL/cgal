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
// file          : src/examples/RangeSegmentTrees/range_tree_1.C
// source        : src/examples/RangeSegmentTrees/range_tree_1.C
// revision      : $Revision$
// revision_date : $Date$/
// author(s)     : Gabriele Neyer
//
// coordinator   : Peter Widmayer, ETH Zurich
//
// Implementation: Testprogram for 1-dimensional Range Trees
// A two dimensional Range Tree is defined in this class.
// Ti is the type of each dimension of the tree.
//
// email         : cgal@cs.uu.nl
//
// ======================================================================

#include <CGAL/basic.h> 

#include <iostream>
#include <CGAL/Range_tree_k.h>
#include "include/Tree_Traits.h"
#include <vector>
#include <iterator>
#include <utility>



typedef  CGAL::Range_tree_1< CGAL::Tree_traits_1> Range_tree_1_type;

int main()
{
  typedef  CGAL::Tree_traits_1::Key Key;
  typedef  CGAL::Tree_traits_1::Interval Interval;
  std::vector<Key> InputList, OutputList;
  std::vector<Key>::iterator first, last, current;

  InputList.push_back(8.0);
  InputList.push_back(1.0);
  InputList.push_back(3.9);
  InputList.push_back(2.0);
  InputList.push_back(5.0);
  InputList.push_back(4.8);
  InputList.push_back(7.7);
  InputList.push_back(6.8);

  first = InputList.begin();
  last = InputList.end();

  Range_tree_1_type Range_tree_1(first,last);

  Interval win(Interval(4.6, 6.8));

  std::cerr << "\n Window Query (4.6, 6.8) \n";
  Range_tree_1.window_query(win, std::back_inserter(OutputList));
  current=OutputList.begin();

  while(current!=OutputList.end())
  {
    std::cerr << (*current) << std::endl;
    current++;
  }

  if(Range_tree_1.CRange_tree_1->is_valid())
    std::cerr << "Tree is valid\n";
  else
    std::cerr << "Tree is not valid\n";
  return 0; 
}


