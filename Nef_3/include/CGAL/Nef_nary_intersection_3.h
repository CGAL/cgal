// Copyright (c) 1997-2002,2005,2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Peter Hachenberger

#ifndef CGAL_NEF_NARY_INTERSECTION_3_H
#define CGAL_NEF_NARY_INTERSECTION_3_H

#include <CGAL/license/Nef_3.h>

#include <CGAL/disable_warnings.h>

#include <list>

namespace CGAL {

template<class Polyhedron>
class Nef_nary_intersection_3 {

  int inserted;
  std::list<Polyhedron> queue;
  typedef typename std::list<Polyhedron>::iterator pit;
  Polyhedron empty;

 public:
  Nef_nary_intersection_3() : inserted(0) {}

  void intersect() {
    pit i1(queue.begin()), i2(i1);
    ++i2;

    Polyhedron tmp(*i1 * *i2);

    queue.pop_front();
    queue.pop_front();
    queue.push_front(tmp);
  }

  void add_polyhedron(const Polyhedron& P) {
    queue.push_front(P);
    ++inserted;
    for(int i=2;(inserted%i) == 0; i*=2) {
      intersect();
    }
  }

  Polyhedron get_intersection() {
    if (queue.empty())
      return empty;
    while(queue.size() > 1)
      intersect();
    inserted = 0;
    return queue.front();
  }
};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_NEF_NARY_INTERSECTION_3_H
