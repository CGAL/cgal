// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France) and
// Notre Dame University (U.S.A.).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>


#ifndef CGAL_POINT_CONTAINER_H
#define CGAL_POINT_CONTAINER_H

#include <CGAL/basic.h>
#include <list>
#include <vector>

CGAL_BEGIN_NAMESPACE

template<class P, class C = std::list<P> >
class Point_container
{
public:
  typedef C  Container;
  typedef P  Point_2;
  typedef typename Container::iterator   Point_handle;
  typedef typename Container::size_type  size_type;

private:
  typedef Point_container<Point_2,Container> Self;

public:
  Point_container() {}

  Point_handle insert(const Point_2& p)
  {
    c.push_back(p);
    return --c.end();
  }

  void remove(Point_handle handle)
  {
    c.erase(handle);
  }

  void swap(const Self& other)
  {
    c.swap(other.c);
  }

  void clear() {
    c.clear();
  }

  size_type size() const { return c.size(); }

private:
  Container c;
};

#if 1
template<class P, unsigned long S>
class Array_point_container
{
public:
  //  typedef C  Container;
  typedef P  Point_2;
  enum { Size = S };
  typedef Point_2*      Point_handle;
  typedef unsigned long size_type;

  Array_point_container() {
    last = 0;
    c = new Point_2[Size];
  }

  Point_handle insert(const Point_2& p)
  {
    CGAL_precondition( last < Size );

    c[last] = p;
    Point_handle h = &c[last];
    last++;
    return h;
  }

  std::size_t size() const { return c.size(); }

  void clear() {
    last = 0;
  }

private:
  unsigned long last;
  Point_2* c;
};
#endif

CGAL_END_NAMESPACE

#endif // CGAL_POINT_CONTAINER_H
