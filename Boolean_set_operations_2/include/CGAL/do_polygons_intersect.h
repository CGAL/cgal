// Copyright (c) 2021  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Baruch Zukerman <baruchzu@post.tau.ac.il>
//            Ron Wein        <wein@post.tau.ac.il>
//            Efi Fogel       <efif@post.tau.ac.il>
//            Simon Giraudot  <simon.giraudot@geometryfactory.com>

#ifndef CGAL_DO_POLYGONS_INTERSECT_H
#define CGAL_DO_POLYGONS_INTERSECT_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <boost/utility/enable_if.hpp>

#include <CGAL/disable_warnings.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2/Bso_internal_functions.h>

namespace CGAL {

/// \name unregularized do_intersect() functions.
//@{

template <class Kernel, class Container>
inline bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                         const Polygon_2<Kernel, Container>& pgn2)
{
  return (_do_intersect(pgn1, pgn2));
}

template <class Kernel, class Container, class Traits>
inline bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                         const Polygon_2<Kernel, Container>& pgn2,
                         Traits& tr)
{
  return (_do_intersect(pgn1, pgn2, tr));
}

template <class Kernel, class Container>
inline bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                         const Polygon_with_holes_2<Kernel, Container>& pgn2)
{
  return (_do_intersect(pgn1, pgn2));
}

template <class Kernel, class Container, class Traits>
inline bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                         const Polygon_with_holes_2<Kernel, Container>& pgn2,
                         Traits& tr)
{
  return (_do_intersect(pgn1, pgn2, tr));
}

template <class Kernel, class Container>
inline bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                         const Polygon_2<Kernel, Container>& pgn2)
{
  return (_do_intersect(pgn1, pgn2));
}

template <class Kernel, class Container, class Traits>
inline bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                         const Polygon_2<Kernel, Container>& pgn2,
                         Traits& tr)
{
  return (_do_intersect(pgn1, pgn2, tr));
}

template <class Kernel, class Container>
inline bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                         const Polygon_with_holes_2<Kernel, Container>& pgn2)
{
  return (_do_intersect(pgn1, pgn2));
}

template <class Kernel, class Container, class Traits>
inline bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                         const Polygon_with_holes_2<Kernel, Container>& pgn2,
                         Traits& tr)
{
  return (_do_intersect(pgn1, pgn2, tr));
}

//@}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif //CGAL_DO_POLYGONS_INTERSECT_H
