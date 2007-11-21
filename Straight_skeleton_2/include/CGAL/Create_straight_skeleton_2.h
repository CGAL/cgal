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
#ifndef CGAL_CREATE_STRAIGHT_SKELETON_2_H
#define CGAL_CREATE_STRAIGHT_SKELETON_2_H

#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/compute_outer_frame_margin.h>
#include <CGAL/Polygon_2.h>

CGAL_BEGIN_NAMESPACE


template<class Poly>
inline typename Poly::const_iterator vertices_begin ( Poly const& aPoly ) { return aPoly.begin() ; }

template<class Poly>
inline typename Poly::const_iterator vertices_end ( Poly const& aPoly ) { return aPoly.end() ; }


template<class Kernel>
inline typename Polygon_2<Kernel>::Vertex_const_iterator vertices_begin ( Polygon_2<Kernel> const& aPoly ) 
{ return aPoly.vertices_begin() ; }

template<class Kernel>
inline typename Polygon_2<Kernel>::Vertex_const_iterator vertices_end( Polygon_2<Kernel> const& aPoly ) 
{ return aPoly.vertices_end() ; }

template<class Poly>
inline Poly reverse_polygon ( Poly const& aPoly ) { return Poly( aPoly.rbegin(), aPoly.rend() ) ; }

template<class Kernel>
inline Polygon_2<Kernel> reverse_polygon( Polygon_2<Kernel> const& aPoly ) 
{
  Polygon_2<Kernel> r ( aPoly ) ;
  r.reverse_orientation();
  return r ;  
}

template<class Kernel, class Contour, class HoleIterator>
boost::shared_ptr< Straight_skeleton_2<Kernel> >
create_interior_straight_skeleton_2 ( Contour const& aOuterContour
                                    , HoleIterator   aHolesBegin
                                    , HoleIterator   aHolesEnd
                                    )
{
  typedef Straight_skeleton_2<Kernel> Ss ;
  typedef boost::shared_ptr<Ss> SsPtr ;

  typedef Straight_skeleton_builder_traits_2<Kernel>      SsBuilderTraits;
  typedef Straight_skeleton_builder_2<SsBuilderTraits,Ss> SsBuilder;
  
  SsBuilder ssb ;
  
  ssb.enter_contour( vertices_begin(aOuterContour), vertices_end(aOuterContour) ) ;
  
  for ( HoleIterator hi = aHolesBegin ; hi != aHolesEnd ; ++ hi )
    ssb.enter_contour( vertices_begin(*hi), vertices_end(*hi) ) ;
  
  return ssb.construct_skeleton();
}

template<class Kernel>
boost::shared_ptr< Straight_skeleton_2<Kernel> >
create_interior_straight_skeleton_2 ( Polygon_2<Kernel> const&  aOuterContour )
{
  std::vector< Polygon_2<Kernel> > no_holes ;
  return create_interior_straight_skeleton_2< Kernel>(aOuterContour, no_holes.begin(), no_holes.end() );
}

template<class Kernel, class Contour, class HoleIterator>
boost::shared_ptr< Straight_skeleton_2<Kernel> >
create_exterior_straight_skeleton_2 ( typename Kernel::FT aMaxOffset
                                    , Contour const&      aOuterContour
                                    , HoleIterator        aHolesBegin
                                    , HoleIterator        aHolesEnd
                                    , bool                aDontReverseOrientation = false
                                    )
{
  typedef typename Kernel::FT      FT ;
  typedef typename Kernel::Point_2 Point_2 ;
  
  typedef Straight_skeleton_2<Kernel> Ss ;
  typedef boost::shared_ptr<Ss> SsPtr ;
  
  SsPtr rSkeleton ;
  
  boost::optional<FT> margin = compute_outer_frame_margin( vertices_begin(aOuterContour)
                                                         , vertices_end  (aOuterContour)
                                                         , aMaxOffset
                                                         );

  if ( margin )
  {
    
    Bbox_2 bbox = bbox_2(vertices_begin(aOuterContour), vertices_end(aOuterContour));

    FT fxmin = bbox.xmin() - *margin ;
    FT fxmax = bbox.xmax() + *margin ;
    FT fymin = bbox.ymin() - *margin ;
    FT fymax = bbox.ymax() + *margin ;

    Contour frame ;
    frame.push_back( Point_2(fxmin,fymin) ) ;
    frame.push_back( Point_2(fxmax,fymin) ) ;
    frame.push_back( Point_2(fxmax,fymax) ) ;
    frame.push_back( Point_2(fxmin,fymax) ) ;

    std::vector<Contour> holes ;
    
    holes.push_back( aDontReverseOrientation ? aOuterContour : reverse_polygon(aOuterContour) ) ;
        
    for ( HoleIterator hi = aHolesBegin ; hi != aHolesEnd ; ++ hi )
      holes.push_back( aDontReverseOrientation ? *hi : reverse_polygon(*hi) ) ;
      
    rSkeleton = create_interior_straight_skeleton_2<Kernel>(frame, holes.begin(), holes.end() ) ;  
  }
  
  return rSkeleton ;
}


template<class Kernel>
boost::shared_ptr< Straight_skeleton_2<Kernel> >
create_exterior_straight_skeleton_2 ( typename Kernel::FT      aMaxOffset
                                    , Polygon_2<Kernel> const& aOuterContour
                                    , bool                     aDontReverseOrientation = false
                                    )
{
  std::vector< Polygon_2<Kernel> > no_holes ;
  return create_exterior_straight_skeleton_2<Kernel>(aMaxOffset
                                                    ,aOuterContour
                                                    ,no_holes.begin()
                                                    ,no_holes.end() 
                                                    ,aDontReverseOrientation
                                                    );
}

CGAL_END_NAMESPACE


#endif // CGAL_STRAIGHT_SKELETON_BUILDER_2_H //
// EOF //
