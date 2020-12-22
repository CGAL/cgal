// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_ADAPTATION_TRAITS_FUNCTORS_H
#define CGAL_VORONOI_DIAGRAM_2_ADAPTATION_TRAITS_FUNCTORS_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

//=========================================================================
//=========================================================================

struct Null_functor
{
  Null_functor() {}
  template<typename T> Null_functor(T /*t*/) {}
};

//=========================================================================

template<class Functor>
struct Functor_exists
{
  typedef Tag_true  Value;
};

template<>
struct Functor_exists<Null_functor>
{
  typedef Tag_false Value;
};

template<class AT, class SI> class Default_caching_site_inserter;

template<class AT>
struct Functor_exists< Default_caching_site_inserter<AT,Null_functor> >
{
  typedef Tag_false Value;
};

//=========================================================================
//=========================================================================


} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL


#endif // CGAL_VORONOI_DIAGRAM_2_ADAPTATION_TRAITS_FUNCTORS_H
