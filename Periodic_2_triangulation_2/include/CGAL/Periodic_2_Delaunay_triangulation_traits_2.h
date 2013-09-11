// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Nico Kruithof <Nico@nghk.nl>

#ifndef CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_TRAITS_2_H
#define CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_TRAITS_2_H

#include <CGAL/Periodic_2_triangulation_traits_2.h>

namespace CGAL
{

/// The Periodic_2_Delaunay_triangulation_traits_2 is equal to
/// Periodic_2_triangulation_traits_2
template < typename K, typename Off = CGAL::Periodic_2_offset_2 >
class Periodic_2_Delaunay_triangulation_traits_2 :
  public Periodic_2_triangulation_traits_2<K, Off>
{
};



} //namespace CGAL

#endif // CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_TRAITS_2_H
