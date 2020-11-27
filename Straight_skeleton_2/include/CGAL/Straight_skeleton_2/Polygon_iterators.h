// Copyright (c) 2007 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>

#ifndef CGAL_STRAIGHT_SKELETON_POLYGON_ITERATOR_H
#define CGAL_STRAIGHT_SKELETON_POLYGON_ITERATOR_H 1

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/Polygon_with_holes_2.h>

#include <boost/shared_ptr.hpp>

namespace CGAL {
namespace CGAL_SS_i {

// to distinguish between SequenceContainers of points, and GeneralPolygonWithHoles_2
BOOST_MPL_HAS_XXX_TRAIT_DEF(Hole_const_iterator)

// Generic container
template<class Poly>
inline typename Poly::const_iterator vertices_begin ( Poly const& aPoly )
{ return aPoly.begin() ; }

template<class Poly>
inline typename Poly::const_iterator vertices_end ( Poly const& aPoly )
{ return aPoly.end() ; }

template<class Poly>
inline typename Poly::const_iterator vertices_begin ( boost::shared_ptr<Poly> const& aPoly )
{ return aPoly->begin() ; }

template<class Poly>
inline typename Poly::const_iterator vertices_end ( boost::shared_ptr<Poly> const& aPoly )
{ return aPoly->end() ; }

// Polygon_2
template<class K, class C>
inline typename Polygon_2<K,C>::Vertex_const_iterator
vertices_begin ( Polygon_2<K,C> const& aPoly )
{ return aPoly.vertices_begin() ; }

template<class K, class C>
inline typename Polygon_2<K,C>::Vertex_const_iterator
vertices_end( Polygon_2<K,C> const& aPoly )
{ return aPoly.vertices_end() ; }

template<class K, class C>
inline typename Polygon_2<K,C>::Vertex_const_iterator vertices_begin ( boost::shared_ptr<Polygon_2<K,C> > const& aPoly )
{ return aPoly->vertices_begin() ; }

template<class K, class C>
inline typename Polygon_2<K,C>::Vertex_const_iterator vertices_end( boost::shared_ptr<Polygon_2<K,C> > const& aPoly )
{ return aPoly->vertices_end() ; }

// Polygon_with_holes_2
template<class K, class C>
inline typename Polygon_with_holes_2<K,C>::Polygon_2::Vertex_const_iterator
vertices_begin ( Polygon_with_holes_2<K,C> const& aPoly )
{ return aPoly.outer_boundary().vertices_begin() ; }

template<class K, class C>
inline typename Polygon_with_holes_2<K,C>::Polygon_2::Vertex_const_iterator
vertices_end( Polygon_with_holes_2<K,C> const& aPoly )
{ return aPoly.outer_boundary().vertices_end() ; }

template<class K, class C>
inline typename Polygon_with_holes_2<K,C>::Polygon_2::Vertex_const_iterator
vertices_begin ( boost::shared_ptr<Polygon_with_holes_2<K,C> > const& aPoly )
{ return aPoly->outer_boundary().vertices_begin() ; }

template<class K, class C>
inline typename Polygon_with_holes_2<K,C>::Polygon_2::Vertex_const_iterator
vertices_end( boost::shared_ptr<Polygon_with_holes_2<K,C> > const& aPoly )
{ return aPoly->outer_boundary().vertices_end() ; }

} // namespace CGAL_SS_i
} // namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_POLYGON_ITERATOR_H //
