// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of an example program for CGAL.  This example
// program may be used, distributed and modified without limitation.
//


// small example for compilation
// check of Delaunay_triangulation_2 using the kernel concept
// archetype

#include <CGAL/basic.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Kernel_archetype.h>

typedef CGAL::Kernel_archetype  K;
typedef K::Point_2              Point;

typedef CGAL::Delaunay_triangulation_2<K>  Delaunay_triang_2;

Delaunay_triang_2 dt;

int main()
{
  std::list<Point> input;  
  Point act;

  dt.insert(act);
  return 0;  
}
