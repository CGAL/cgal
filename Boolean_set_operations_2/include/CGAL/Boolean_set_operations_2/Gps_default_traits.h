// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_GPS_DEFAULT_TRAITS_H
#define CGAL_GPS_DEFAULT_TRAITS_H

#include <CGAL/license/Boolean_set_operations_2.h>


#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <CGAL/General_polygon_2.h>
#include <CGAL/Gps_segment_traits_2.h>
#include <CGAL/Gps_traits_2.h>

namespace CGAL {

template <class Polygon>
struct Gps_default_traits
{};


template <class Kernel, class Container>
struct Gps_default_traits<CGAL::Polygon_2<Kernel, Container> >
{
  typedef Gps_segment_traits_2<Kernel,
                               Container,
                               Arr_segment_traits_2<Kernel> >    Traits;
};

template <class Kernel, class Container>
struct Gps_default_traits<CGAL::Polygon_with_holes_2<Kernel, Container> >
{
  typedef Gps_segment_traits_2<Kernel,
                               Container,
                               Arr_segment_traits_2<Kernel> >    Traits;
};

template <class Polygon>
struct Gps_default_traits<CGAL::General_polygon_with_holes_2<Polygon> >
{
  typedef typename Gps_default_traits<Polygon>::Traits Traits;
};

template <class Arr_traits>
struct Gps_default_traits<CGAL::General_polygon_2<Arr_traits> >
{
  typedef Gps_traits_2<Arr_traits>    Traits;
};

} //namespace CGAL

#endif
