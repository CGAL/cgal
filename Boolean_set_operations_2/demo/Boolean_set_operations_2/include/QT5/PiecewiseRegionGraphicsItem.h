// Copyright (c) 2009  GeometryFactory Sarl (France).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s) : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//             Ronnie Gandhi <ronniegandhi19999@gmail.com>

#ifndef CGAL_QT_PIECEWISE_REGION_GRAPHICS_ITEM_BEZIER_H
#define CGAL_QT_PIECEWISE_REGION_GRAPHICS_ITEM_BEZIER_H

#include "QT5/PiecewiseBoundaryGraphicsItem.h"

namespace CGAL {

namespace Qt {

template <class Piecewise_region_, class Draw_piece_, class Piece_bbox_>
class Piecewise_region_graphics_item_bezier : public Piecewise_boundary_graphics_item_bezier< typename Piecewise_region_::General_polygon_2, Draw_piece_, Piece_bbox_ > 
{
  typedef Piecewise_boundary_graphics_item_bezier< typename Piecewise_region_::General_polygon_2, Draw_piece_, Piece_bbox_> Base ;
  
  typedef Piecewise_region_ Piecewise_region ;
  typedef Draw_piece_       Draw_piece ;
  typedef Piece_bbox_       Piece_bbox ;
  
  typedef typename Piecewise_region::Hole_const_iterator Hole_const_itertator ;
  
public:

  Piecewise_region_graphics_item_bezier( Piecewise_region* aRegion, Draw_piece const& aPieceDrawer = Draw_piece(), Piece_bbox const& aPieceBBox = Piece_bbox() )
    :
     Base(aPieceDrawer, aPieceBBox)
    ,mRegion(aRegion)
  {}  

public:

  virtual bool isModelEmpty() const { return !mRegion || mRegion->outer_boundary().size() ; }
  
protected:
  
  Piecewise_region_graphics_item_bezier( Draw_piece const& aPieceDrawer = Draw_piece(), Piece_bbox const& aPieceBBox = Piece_bbox() )
    :
     Base(aPieceDrawer, aPieceBBox)
  {}  
  
  virtual void update_bbox( Piecewise_graphics_item_base_bezier::Bbox_builder& aBboxBuilder)
  {
    if ( mRegion ) 
      update_region_bbox(*mRegion, aBboxBuilder ) ;
  }    

  virtual void draw_model ( QPainterPath& aPath ) 
  {
    if ( mRegion )
      draw_region(*mRegion,aPath);  
  }

  void update_region_bbox( Piecewise_region const& aRegion, Piecewise_graphics_item_base_bezier::Bbox_builder& aBboxBuilder ) ;
  void draw_region       ( Piecewise_region const& aRegion, QPainterPath& aPath ) ;
  
protected:

  Piecewise_region* mRegion;
};

template <class R, class D, class P>
void Piecewise_region_graphics_item_bezier<R,D,P>::update_region_bbox( Piecewise_region const& aRegion, Piecewise_graphics_item_base_bezier::Bbox_builder& aBboxBuilder )
{
  this->update_boundary_bbox( aRegion.outer_boundary(), aBboxBuilder ) ;//"This" added for qt5 version !
  
  for( Hole_const_itertator hit = aRegion.holes_begin(); hit != aRegion.holes_end(); ++ hit )
    this->update_boundary_bbox(*hit,aBboxBuilder);//"This" added for qt5 version !
}

template <class R, class D, class P>
void Piecewise_region_graphics_item_bezier<R,D,P>::draw_region( Piecewise_region const& aRegion, QPainterPath& aPath )
{
  this->draw_boundary( aRegion.outer_boundary(), aPath ) ;//This added for qt5 version !
  
  for( Hole_const_itertator hit = aRegion.holes_begin(); hit != aRegion.holes_end(); ++ hit )
    this->draw_boundary(*hit,aPath);//"This" added for qt5 version !
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_PIECEWISE_REGION_GRAPHICS_ITEM_BEZIER_H
