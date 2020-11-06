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

#include <CGAL/Straight_skeleton_2/Polygon_iterators.h>

#include <CGAL/assertions.h>
#include <CGAL/Polygon_with_holes_2.h>

#include <boost/range/value_type.hpp>
#include <boost/shared_ptr.hpp>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>

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
template<class K, class InputPolygonPtrIterator, class OutputPolygonWithHolesPtrIterator,
         class PolygonWithHoles = Polygon_with_holes_2<K> >
bool arrange_offset_polygons_2 ( InputPolygonPtrIterator           aBegin
                               , InputPolygonPtrIterator           aEnd
                               , OutputPolygonWithHolesPtrIterator rOut
                               , const K&
                               )
{
  bool bk_poly_assert_mode = get_use_polygon_assertions();
  set_use_polygon_assertions(false); // disable assertions in Polygon_2 function as we may manipulate strictly simple polygons

  typedef typename std::iterator_traits<InputPolygonPtrIterator>::difference_type difference_type ;
  typedef typename std::iterator_traits<InputPolygonPtrIterator>::value_type PolygonPtr ;

  typedef boost::shared_ptr<PolygonWithHoles> PolygonWithHolesPtr ;

  difference_type lSize = std::distance(aBegin,aEnd);

  std::vector<PolygonWithHolesPtr> lTable(lSize);

  for ( InputPolygonPtrIterator it = aBegin ; it != aEnd ; ++ it )
  {
    difference_type lIdx = std::distance(aBegin,it);

    const PolygonPtr lPoly = *it ;

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
    const PolygonPtr lPoly = *it ;

    difference_type lIdx = std::distance(aBegin,it);

    // Is a hole?
    if ( !lTable[lIdx] )
    {
      PolygonWithHolesPtr lParent ;

      for ( difference_type j = 0 ; j < lSize && !lParent ; ++ j )
      {
        PolygonWithHolesPtr lOuter = lTable[j];
        if ( lOuter )
          for ( auto v = CGAL_SS_i::vertices_begin(lPoly), ve = CGAL_SS_i::vertices_end(lPoly); v != ve && !lParent ; ++ v )
          {
            if ( lOuter->outer_boundary().bounded_side(*v) == ON_BOUNDED_SIDE )
              lParent = lOuter ;

          }
      }

      if (lParent == nullptr)
      {
        set_use_polygon_assertions(bk_poly_assert_mode);
        return false;
      }

      lParent->add_hole(*lPoly);
    }
  }
  set_use_polygon_assertions(bk_poly_assert_mode);
  return true;
}

template<class PolygonWithHoles, class Polygon>
std::vector< boost::shared_ptr<PolygonWithHoles> >
inline
arrange_offset_polygons_2 ( std::vector<boost::shared_ptr<Polygon> > const& aPolygons,
                            bool& no_error)
{
  typedef std::vector< boost::shared_ptr<PolygonWithHoles> > result_type ;
  typedef std::back_insert_iterator<result_type> Back_inserter;

  typedef typename PolygonWithHoles::General_polygon_2 Polygon_2 ;
  typedef typename Kernel_traits<typename boost::range_value<Polygon_2>::type>::Kernel K ;

  typedef typename std::vector<boost::shared_ptr<Polygon> >::const_iterator PolygonIterator ;

  result_type rResult ;
  no_error = arrange_offset_polygons_2<K, PolygonIterator, Back_inserter, PolygonWithHoles>(
               aPolygons.begin(), aPolygons.end(), std::back_inserter(rResult), K() ) ;

  return rResult ;
}

template<class PolygonWithHoles, class Polygon>
std::vector< boost::shared_ptr<PolygonWithHoles> >
inline
arrange_offset_polygons_2 ( std::vector<boost::shared_ptr<Polygon> > const& aPolygons)
{
  bool no_error;
  return arrange_offset_polygons_2<PolygonWithHoles>(aPolygons, no_error);
}

} // end namespace CGAL

#endif // CGAL_ARRANGE_OFFSET_POLYGONS_2_H
