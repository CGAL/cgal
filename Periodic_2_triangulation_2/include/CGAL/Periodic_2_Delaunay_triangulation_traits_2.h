// Copyright (c) 2011   INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// $URL: svn+ssh://cgal/svn/cgal/branches/experimental-packages/Periodic_2_triangulation_2/include/CGAL/Periodic_2_triangulation_traits_2.h $
// $Id: Periodic_2_triangulation_traits_2.h 60448 2010-12-21 15:40:31Z nicokruithof $
//
//
// Author(s)     : Nico Kruithof <Nico@nghk.nl>

#ifndef CGAL_PERIODIC_2_TRIANGULATION_TRAITS_2_H
#define CGAL_PERIODIC_2_TRIANGULATION_TRAITS_2_H

#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>

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

#endif // CGAL_PERIODIC_2_TRIANGULATION_TRAITS_2_H
