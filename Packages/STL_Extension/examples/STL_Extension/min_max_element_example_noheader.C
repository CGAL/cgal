// Copyright (c) 2003  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of an example program for CGAL. This example
// program may be used, distributed and modified without limitation.

#include <CGAL/algorithm.h>
#include <vector>
#include <iostream>

using std::vector;
using std::pair;
using std::cout;
using std::endl;
using CGAL::min_max_element;

int main()
{
  vector< int > v;
  v.push_back(3);
  v.push_back(6);
  v.push_back(5);
  typedef std::vector< int >::iterator iterator;
  pair< iterator, iterator > p = min_max_element(v.begin(), v.end());
  cout << "min = " << *p.first << ", max = " << *p.second << endl;
  return 0;
}
