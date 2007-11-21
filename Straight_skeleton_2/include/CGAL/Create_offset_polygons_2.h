// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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

// $URL: svn+ssh://CGAL/svn/cgal/trunk/Straight_skeleton_2/include/CGAL/Straight_skeleton_builder_2.h $
// $Id: Straight_skeleton_builder_2.h 40951 2007-11-19 16:33:25Z fcacciola $
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_CREATE_OFFSET_POLYGONS_2_H
#define CGAL_CREATE_OFFSET_POLYGONS_2_H

#include <CGAL/Create_straight_skeleton_2.h>
#include <CGAL/Polygon_offset_builder_2.h>
#include <CGAL/compute_outer_frame_margin.h>
#include <CGAL/Polygon_2.h>

CGAL_BEGIN_NAMESPACE

template<class Kernel, class OutPolygon>
std::vector< boost::shared_ptr<OutPolygon> > 
create_offset_polygons_2 ( typename Kernel::FT aOffset , boost::shared_ptr< Straight_skeleton_2<Kernel> > aSs )
{
  typedef boost::shared_ptr<OutPolygon> OutPolygonPtr ; 
  typedef std::vector<OutPolygonPtr>    OutPolygonPtrVector ;
   
  typedef Straight_skeleton_2<Kernel>                                 Ss ; 
  typedef Polygon_offset_builder_traits_2<Kernel>                     OffsetBuilderTraits;
  typedef Polygon_offset_builder_2<Ss,OffsetBuilderTraits,OutPolygon> OffsetBuilder;
  
  OutPolygonPtrVector rR ;
  
  CGAL_precondition(aSs);
  
  if ( aSs )
  {
    OffsetBuilder ob(*aSs);
    ob.construct_offset_contours(aOffset, std::back_inserter(rR) ) ;
  }
    
  return rR ;
}

template<class Kernel>
std::vector< boost::shared_ptr< Polygon_2<Kernel> > > 
create_offset_polygons_2 ( typename Kernel::FT aOffset , boost::shared_ptr< Straight_skeleton_2<Kernel> > aSs )
{
  return create_offset_polygons_2<Kernel, Polygon_2<Kernel> >(aOffset,aSs);
}

template<class Kernel, class Polygon, class HoleIterator>
std::vector< boost::shared_ptr<Polygon> >
create_interior_skeleton_and_offset_polygons_2 ( typename Kernel::FT aOffset
                                               , Polygon const&      aOuterBoundary
                                               , HoleIterator        aHolesBegin
                                               , HoleIterator        aHolesEnd
                                               )
{
  return create_offset_polygons_2<Kernel,Polygon>(aOffset
                                                 ,create_interior_straight_skeleton_2<Kernel>(aOuterBoundary, aHolesBegin, aHolesEnd )
                                                 );
    
}

template<class Kernel, class Polygon, class HoleIterator>
std::vector< boost::shared_ptr<Polygon> >
create_exterior_skeleton_and_offset_polygons_2 ( typename Kernel::FT aOffset
                                               , Polygon const&      aOuterBoundary
                                               , HoleIterator        aHolesBegin
                                               , HoleIterator        aHolesEnd
                                               , bool                aDontReverseOrientation = false
                                               )
{
  return create_offset_polygons_2<Kernel,Polygon>(aOffset
                                                 ,create_exterior_straight_skeleton_2<Kernel>(aOffset
                                                                                             ,aOuterBoundary
                                                                                             ,aHolesBegin
                                                                                             ,aHolesEnd
                                                                                             ,aDontReverseOrientation 
                                                                                             )
                                                 );
}

template<class Kernel>
std::vector< boost::shared_ptr< Polygon_2<Kernel> > >
create_interior_skeleton_and_offset_polygons_2 ( typename Kernel::FT      aOffset
                                               , Polygon_2<Kernel> const& aOuterBoundary
                                               )
{
  std::vector< Polygon_2<Kernel> > no_holes ;
  return create_interior_skeleton_and_offset_polygons_2<Kernel>(aOffset
                                                               ,aOuterBoundary
                                                               ,no_holes.begin()
                                                               ,no_holes.end()
                                                               );
}

template<class Kernel>
std::vector< boost::shared_ptr< Polygon_2<Kernel> > >
create_exterior_skeleton_and_offset_polygons_2 ( typename Kernel::FT      aOffset
                                               , Polygon_2<Kernel> const& aOuterBoundary
                                               )
{
  std::vector< Polygon_2<Kernel> > no_holes ;
  return create_exterior_skeleton_and_offset_polygons_2<Kernel>(aOffset
                                                               ,aOuterBoundary
                                                               ,no_holes.begin()
                                                               ,no_holes.end()
                                                               );
}


CGAL_END_NAMESPACE


#endif
// EOF //
