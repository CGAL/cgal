// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_GENERAL_POLYGON_SET_2_H
#define CGAL_GENERAL_POLYGON_SET_2_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/basic.h>

#include <CGAL/Boolean_set_operations_2/Gps_default_dcel.h>
#include <CGAL/General_polygon_set_on_surface_2.h>
#include <CGAL/Arrangement_2/Arr_default_planar_topology.h>
#include <CGAL/Arrangement_2.h>

namespace CGAL {

// General_polygon_set_2
template <class Traits_, class Dcel_ = Gps_default_dcel<Traits_> >
class General_polygon_set_2 : public General_polygon_set_on_surface_2
  <Traits_, typename Default_planar_topology<Traits_, Dcel_>::Traits>
{
protected:
  typedef General_polygon_set_2<Traits_, Dcel_>           Self;

public:
  typedef Traits_                                         Traits_2;
  typedef Dcel_                                           Dcel;

  typedef General_polygon_set_on_surface_2 <Traits_2,
    typename Default_planar_topology<Traits_2, Dcel >::Traits>
                                                          Base;

  typedef CGAL::Arrangement_2<Traits_2, Dcel>             Arrangement_2;

  typedef typename Base::Polygon_2                        Polygon_2;
  typedef typename Base::Polygon_with_holes_2             Polygon_with_holes_2;

  // default costructor
  General_polygon_set_2() : Base() {}

  // constructor from a traits object
  General_polygon_set_2(const Traits_2& traits) : Base(traits) {}

  // constructor from a polygon
  explicit General_polygon_set_2(const Polygon_2& pgn) : Base(pgn) {}

  // constructor from a polygon with holes
  explicit General_polygon_set_2(const Polygon_with_holes_2& pwh) : Base(pwh) {}

  // constructor from a polygon and a traits object
  explicit General_polygon_set_2(const Polygon_2& pgn, const Traits_2& traits) :
    Base(pgn, traits)
  {}

  // constructor from a polygon with holes and a traits object
  explicit General_polygon_set_2(const Polygon_with_holes_2& pwh,
                                 const Traits_2& traits) :
    Base(pwh, traits)
  {}

  // For some reason the below functions (the ones that we call "using" for)
  // are hidden by the function in this class and are not found in the parent's
  // class (General_polygon_set_on_surface_2) when they are called on an
  // object of type General_polygon_set_2.
  // Check in the Vandervoorde / Stroustrup books what is the exact reason.
  // (There may be a better and more correct solution.)
  using Base::intersection;
  using Base::join;
  using Base::symmetric_difference;

  inline void intersection(const Self& ps1, const Self& ps2)
  {
    Base::intersection(static_cast<const Base&>(ps1),
                       static_cast<const Base&>(ps2));
  }

  inline void join(const Self& ps1, const Self& ps2)
  {
    Base::join(static_cast<const Base&>(ps1), static_cast<const Base&>(ps2));
  }

  inline void symmetric_difference(const Self& ps1, const Self& ps2)
  {
    Base::symmetric_difference(static_cast<const Base&>(ps1),
                               static_cast<const Base&>(ps2));
  }

  //@{

  /*! Obtain a const reference to the underlying arrangement
   * \return the underlying arrangement.
   */
  const Arrangement_2& arrangement() const
  { return *(static_cast<const Arrangement_2*>(this->m_arr)); }

  /*! Obtain a reference to the underlying arrangement
   * \return the underlying arrangement.
   */
  Arrangement_2& arrangement()
  { return *(static_cast<Arrangement_2*>(this->m_arr)); }

  //@}
};



//-----------------------------------------------------------------------//
//                          operator<<
//-----------------------------------------------------------------------//
/*!
This operator exports a General_polygon_with_holes_2 to the output stream `out`.

An \ascii and a binary format exist. The format can be selected with
the \cgal modifiers for streams, `set_ascii_mode(0` and `set_binary_mode()`
respectively. The modifier `set_pretty_mode()` can be used to allow for (a
few) structuring comments in the output. Otherwise, the output would
be free of comments. The default for writing is \ascii without comments.

The number of curves of the outer boundary is exported followed by the
curves themselves. Then, the number of holes
is exported, and for each hole, the number of curves on its outer
boundary is exported followed by the curves themselves.

\relates General_polygon_with_holes_2
*/
template <class Traits_, class Dcel_ = Gps_default_dcel<Traits_> >
std::ostream
&operator<<(std::ostream &os, const General_polygon_set_2<Traits_, Dcel_>& g)
{
  typedef typename General_polygon_set_2<Traits_, Dcel_>::Polygon_with_holes_2 Polygon_with_holes_2;
  std::vector<Polygon_with_holes_2> v;
  const unsigned int n = g.number_of_polygons_with_holes();

  v.reserve(n);
  g.polygons_with_holes(std::back_inserter(v));

  switch(IO::get_mode(os)) {
    case IO::ASCII :
      os << n;
      for (const Polygon_with_holes_2 &H: v) {
        os << ' ' << H;
      }
      return os;

    case IO::BINARY :
      write(os, n);
      for (const Polygon_with_holes_2 &H: v) {
        os << H;
      }
      return os;


    default:
      os << "General_polygon_set_2( " << std::endl;
      os << n << " ";
      for (const Polygon_with_holes_2 &H: v) {
        os << H;
      }
      return os;
  }
}

//-----------------------------------------------------------------------//
//                          operator>>
//-----------------------------------------------------------------------//

/*!
This operator imports a General_polygon_with_holes_2 from the input stream `in`.

Both ASCII and binary formats are supported, and the format is automatically detected.

The format consists of the number of curves of the outer boundary
followed by the curves themselves, followed
by the number of holes, and for each hole, the number of curves on its
outer boundary is followed by the curves themselves.

\relates General_polygon_with_holes_2
*/
template <class Traits_, class Dcel_ = Gps_default_dcel<Traits_> >
std::istream
&operator>>(std::istream &is, General_polygon_set_2<Traits_, Dcel_>& g)
{
  typedef typename General_polygon_set_2<Traits_, Dcel_>::Polygon_with_holes_2 Polygon_with_holes_2;

  g.clear();

  unsigned int n;

  switch(IO::get_mode(is)) {
    case IO::ASCII :
      is >> n;
      break;
    case IO::BINARY :
      read(is, n);
      break;
  }

  for (unsigned int i = 0; i < n; i++) {
    Polygon_with_holes_2 pgn_hole;
    is >> pgn_hole;
    g.insert(pgn_hole);
  }

  return is;
}


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
