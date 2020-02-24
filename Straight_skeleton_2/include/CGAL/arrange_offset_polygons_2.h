// Copyright (c) 2006-2008 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//

// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar> 
//
#ifndef CGAL_ARRANGE_OFFSET_POLYGONS_2_H
#define CGAL_ARRANGE_OFFSET_POLYGONS_2_H

#include <CGAL/license/Straight_skeleton_2.h>


#include <boost/shared_ptr.hpp>
#include <CGAL/Polygon_with_holes_2.h>

namespace CGAL {

//
// This should only be used to arrange the polygons coming from Polygon_offset_builder
// as it uses their known properties:
//
//  Polygons are simple
//  Outer polygons are CCW while holes are CW
//  Outer polygons do not contain other outer polygons, only holes
//  Every hole is contained in one and only one outer polygon
//
template<class K, class InputPolygonPtrIterator, class OutputPolygonWithHolesPtrIterator>
bool arrange_offset_polygons_2 ( InputPolygonPtrIterator           aBegin
                               , InputPolygonPtrIterator           aEnd
                               , OutputPolygonWithHolesPtrIterator rOut
                               , K const&
                               )
{
  typedef typename std::iterator_traits<InputPolygonPtrIterator>::difference_type difference_type ;
  
  typedef Polygon_2<K>            Polygon ;
  typedef Polygon_with_holes_2<K> PolygonWithHoles ;
  
  typedef boost::shared_ptr<Polygon>          PolygonPtr ;
  typedef boost::shared_ptr<PolygonWithHoles> PolygonWithHolesPtr ;
  
  typedef typename Polygon::Vertex_const_iterator Vertex_const_iterator ;
  
  difference_type lSize = std::distance(aBegin,aEnd);
  
  std::vector<PolygonWithHolesPtr> lTable(lSize);
  
  for ( InputPolygonPtrIterator it = aBegin ; it != aEnd ; ++ it )
  {
    difference_type lIdx = std::distance(aBegin,it);
    
    PolygonPtr lPoly = *it ;
    Orientation lOrient = lPoly->orientation();
    
    // It's an outer boundary
    if ( lOrient == COUNTERCLOCKWISE )
    {
      PolygonWithHolesPtr lOuter( new PolygonWithHoles(*lPoly) );
      *rOut ++ = lOuter ; 
      lTable[lIdx] = lOuter ;
    }  
  }  
  
  for ( InputPolygonPtrIterator it = aBegin ; it != aEnd ; ++ it )
  {
    PolygonPtr lPoly = *it ;
    
    difference_type lIdx = std::distance(aBegin,it);
    
    // Is a hole?
    if ( !lTable[lIdx] )
    {
      PolygonWithHolesPtr lParent ;
      
      for ( difference_type j = 0 ; j < lSize && !lParent ; ++ j )
      {
        PolygonWithHolesPtr lOuter = lTable[j];
        if ( lOuter )
          for ( Vertex_const_iterator v = lPoly->vertices_begin(), ve = lPoly->vertices_end(); v != ve && !lParent ; ++ v )
            if ( lOuter->outer_boundary().bounded_side(*v) == ON_BOUNDED_SIDE )
              lParent = lOuter ;
      }
      
      if (lParent == nullptr)
        return false;
      
      lParent->add_hole(*lPoly);
    }
  }  
  return true;
}

template<class K, class C>
std::vector< boost::shared_ptr< Polygon_with_holes_2<K,C> > >
inline
arrange_offset_polygons_2 ( std::vector< boost::shared_ptr< Polygon_2<K,C> > > const& aPolygons, bool& no_error)
{
  std::vector< boost::shared_ptr< Polygon_with_holes_2<K,C> > > rResult ;
  
  no_error = arrange_offset_polygons_2(aPolygons.begin(), aPolygons.end(), std::back_inserter(rResult), K() ) ;
  
  return rResult ;
}

template<class K, class C>
std::vector< boost::shared_ptr< Polygon_with_holes_2<K,C> > >
inline
arrange_offset_polygons_2 ( std::vector< boost::shared_ptr< Polygon_2<K,C> > > const& aPolygons)
{
  bool no_error;
  return arrange_offset_polygons_2(aPolygons, no_error);
}

} // end namespace CGAL


#endif 
// EOF //
