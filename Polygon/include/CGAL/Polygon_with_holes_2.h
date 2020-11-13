// Copyright (c) 2005
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_POLYGON_WITH_HOLES_2_H
#define CGAL_POLYGON_WITH_HOLES_2_H

#include <CGAL/Polygon_2.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <vector>
#include <algorithm>

namespace CGAL {

/*!
\ingroup PkgPolygon2Ref

The class `Polygon_with_holes_2` models the concept `GeneralPolygonWithHoles_2`.
It represents a linear polygon with holes. It is parameterized with two
types (`Kernel` and `Container`) that are used to instantiate
the type `Polygon_2<Kernel,Container>`. The latter is used to
represents the outer boundary and the boundary of the holes (if any exist).

\cgalModels `GeneralPolygonWithHoles_2`

*/
template <class Kernel,
          class Containter = std::vector<typename Kernel::Point_2> >
class Polygon_with_holes_2 :
  public General_polygon_with_holes_2<CGAL::Polygon_2<Kernel, Containter> >
{
public:

  typedef CGAL::Polygon_2<Kernel, Containter>        Polygon_2;
  typedef General_polygon_with_holes_2<Polygon_2>    Base;
  typedef typename Base::Hole_const_iterator         Hole_const_iterator;
  typedef typename Base::Size                        Size;

  /*! %Default constructor. */
  Polygon_with_holes_2 () :
    Base()
  {}

  /*! Constructor from the base class. */
  Polygon_with_holes_2 (const Base& base) :
    Base (base)
  {}

  /*! Constructor from a polygon. */
  explicit Polygon_with_holes_2 (const Polygon_2& pgn_boundary) :
    Base (pgn_boundary)
  {}

  /*! Constructor from a polygon (outer boundary) and hole polygons. */
  template <class HolesInputIterator>
  Polygon_with_holes_2 (const Polygon_2& pgn_boundary,
                        HolesInputIterator h_begin,
                        HolesInputIterator h_end) :
    Base (pgn_boundary, h_begin, h_end)
  {}

  /*! Obtain the bounding box of the polygon with holes */
  Bbox_2 bbox() const { return this->outer_boundary().bbox(); }
};

//-----------------------------------------------------------------------//
//                          operator<<
//-----------------------------------------------------------------------//

/*!
This operator exports a polygon with holes to the output stream `out`.

An ASCII and a binary format exist. The format can be selected with
the \cgal modifiers for streams, `set_ascii_mode()` and `set_binary_mode()`
respectively. The modifier `set_pretty_mode()` can be used to allow for (a
few) structuring comments in the output. Otherwise, the output would
be free of comments. The default for writing is ASCII without
comments.

The number of points of the outer boundary is exported followed by the
points themselves in counterclockwise order. Then, the number of holes
is exported, and for each hole, the number of points on its outer
boundary is exported followed by the points themselves in clockwise
order.

\relates Polygon_with_holes_2
*/
template <class Kernel_, class Container_>
std::ostream& operator<<(std::ostream &os,
                         const Polygon_with_holes_2<Kernel_, Container_>& p)
{
  typename Polygon_with_holes_2<Kernel_,Container_>::Hole_const_iterator i;

  switch(get_mode(os)) {
    case IO::ASCII :
      os << p.outer_boundary() << ' ' << p.number_of_holes()<<' ';
      for (i = p.holes_begin(); i != p.holes_end(); ++i) {
        os << *i << ' ';
      }
      return os;

    case IO::BINARY :
       os << p.outer_boundary() << p.number_of_holes();
      for (i = p.holes_begin(); i != p.holes_end(); ++i) {
        os << *i ;
      }
      return os;

    default:
      os << "Polygon_with_holes_2(" << std::endl;
      if(p.is_unbounded())
        os << "No outer bounary" << std::endl;
      else
      {
        os << "Boundary(" << std::endl;
        os << p.outer_boundary() << std::endl;
      }

      os << "Holes" << std::endl;
      os << p.number_of_holes()<<std::endl;
      for (i = p.holes_begin(); i != p.holes_end(); ++i) {
        os <<" "<< *i << std::endl;
      }

      os << ")" << std::endl;
      return os;
  }
}

//-----------------------------------------------------------------------//
//                          operator>>
//-----------------------------------------------------------------------//

/*!
This operator imports a polygon with holes from the input stream `in`.

An ASCII and a binary format exist. The stream detects the format
automatically and can read both.

The format consists of the number of points of the outer boundary followed
by the points themselves in counterclockwise order, followed by the number of holes,
and for each hole, the number of points of the outer boundary is followed
by the points themselves in clockwise order.

\relates Polygon_with_holes_2
*/
template <class Kernel_, class Container_>
std::istream &operator>>(std::istream &is,
                         Polygon_with_holes_2<Kernel_, Container_>& p)
{
  typedef CGAL::Polygon_2<Kernel_, Container_> Polygon_2;
  p.clear();
  is >> p.outer_boundary();

  unsigned int n; // number of holes;
  is >> n;
  if(is)
  {
     for (unsigned int i=0; i<n; i++)
     {
       Polygon_2 hole;
       is >> hole;
       p.add_hole(hole);
     }
  }

  return is;
}


//-----------------------------------------------------------------------//
//                          operator==
//-----------------------------------------------------------------------//
template <class Kernel_, class Container_>
bool operator==(const Polygon_with_holes_2<Kernel_, Container_>& p1,
                const Polygon_with_holes_2<Kernel_, Container_>& p2)
{
  typedef typename
    Polygon_with_holes_2<Kernel_, Container_>::Hole_const_iterator HCI;
  typedef CGAL::Polygon_2<Kernel_, Container_> Polygon_2;
  if(&p1 == &p2)
    return (true);

  if(p1.number_of_holes() != p2.number_of_holes())
    return (false);

  if(p1.outer_boundary() != p2.outer_boundary())
    return (false);

  std::list<Polygon_2> tmp_list(p2.holes_begin(), p2.holes_end());

  HCI i = p1.holes_begin();
  for(; i!= p1.holes_end(); ++i)
  {
    typename std::list<Polygon_2>::iterator j =
      (std::find(tmp_list.begin(), tmp_list.end(), *i));

    if(j == tmp_list.end())
      return (false);

    tmp_list.erase(j);
  }


  CGAL_assertion(tmp_list.empty());
  return (true);
}

//-----------------------------------------------------------------------//
//                          operator!=
//-----------------------------------------------------------------------//
template <class Kernel_, class Container_>
inline bool operator!=(const Polygon_with_holes_2<Kernel_, Container_>& p1,
                       const Polygon_with_holes_2<Kernel_, Container_>& p2)
{
  return (!(p1==p2));
}

//-----------------------------------------------------------------------//
//                          operator==
//-----------------------------------------------------------------------//
template <class Kernel_, class Container_>
bool operator==(const Polygon_with_holes_2<Kernel_, Container_>& p1,
                const Polygon_2<Kernel_, Container_>& p2)
{
  return (p1.outer_boundary() == p2  &&  !p1.number_of_holes());
}

template <class Kernel_, class Container_>
inline bool operator==(const Polygon_2<Kernel_, Container_>& p1,
                       const Polygon_with_holes_2<Kernel_, Container_>& p2)
{
  return (p2 == p1);
}

template <class Kernel_, class Container_>
inline bool operator!=(const Polygon_with_holes_2<Kernel_, Container_>& p1,
                       const Polygon_2<Kernel_, Container_>& p2)
{
  return (!(p1==p2));
}

template <class Kernel_, class Container_>
inline bool operator!=(const Polygon_2<Kernel_, Container_>& p1,
                       const Polygon_with_holes_2<Kernel_, Container_>& p2)
{
  return (!(p1==p2));
}



} //namespace CGAL

#endif
