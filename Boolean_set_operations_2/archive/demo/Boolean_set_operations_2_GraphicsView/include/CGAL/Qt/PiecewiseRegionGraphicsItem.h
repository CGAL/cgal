// Copyright (c) 2009  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Fernando Cacciola <fernando.cacciola@geometryfactory.com>

#ifndef CGAL_QT_PIECEWISE_REGION_GRAPHICS_ITEM_H
#define CGAL_QT_PIECEWISE_REGION_GRAPHICS_ITEM_H

#include <CGAL/Qt/PiecewiseBoundaryGraphicsItem.h>

namespace CGAL {

namespace Qt {

template <class Piecewise_region_, class Draw_piece_, class Piece_bbox_>
class Piecewise_region_graphics_item : public Piecewise_boundary_graphics_item< typename Piecewise_region_::General_polygon_2, Draw_piece_, Piece_bbox_ >
{
  typedef Piecewise_boundary_graphics_item< typename Piecewise_region_::General_polygon_2, Draw_piece_, Piece_bbox_> Base ;

  typedef Piecewise_region_ Piecewise_region ;
  typedef Draw_piece_       Draw_piece ;
  typedef Piece_bbox_       Piece_bbox ;

  typedef typename Piecewise_region::Hole_const_iterator Hole_const_itertator ;

public:

  Piecewise_region_graphics_item( Piecewise_region* aRegion, Draw_piece const& aPieceDrawer = Draw_piece(), Piece_bbox const& aPieceBBox = Piece_bbox() )
    :
     Base(aPieceDrawer, aPieceBBox)
    ,mRegion(aRegion)
  {}

public:

  virtual bool isModelEmpty() const { return !mRegion || mRegion->outer_boundary().size() ; }

protected:

  Piecewise_region_graphics_item( Draw_piece const& aPieceDrawer = Draw_piece(), Piece_bbox const& aPieceBBox = Piece_bbox() )
    :
     Base(aPieceDrawer, aPieceBBox)
  {}

  virtual void update_bbox( Piecewise_graphics_item_base::Bbox_builder& aBboxBuilder)
  {
    if ( mRegion )
      update_region_bbox(*mRegion, aBboxBuilder ) ;
  }

  virtual void draw_model ( QPainterPath& aPath )
  {
    if ( mRegion )
      draw_region(*mRegion,aPath);
  }

  void update_region_bbox( Piecewise_region const& aRegion, Piecewise_graphics_item_base::Bbox_builder& aBboxBuilder ) ;
  void draw_region       ( Piecewise_region const& aRegion, QPainterPath& aPath ) ;

protected:

  Piecewise_region* mRegion;
};

template <class R, class D, class P>
void Piecewise_region_graphics_item<R,D,P>::update_region_bbox( Piecewise_region const& aRegion, Piecewise_graphics_item_base::Bbox_builder& aBboxBuilder )
{
  this->update_boundary_bbox( aRegion.outer_boundary(), aBboxBuilder ) ;//"This" added for qt5 version !

  for( Hole_const_itertator hit = aRegion.holes_begin(); hit != aRegion.holes_end(); ++ hit )
    this->update_boundary_bbox(*hit,aBboxBuilder);//"This" added for qt5 version !
}

template <class R, class D, class P>
void Piecewise_region_graphics_item<R,D,P>::draw_region( Piecewise_region const& aRegion, QPainterPath& aPath )
{
  this->draw_boundary( aRegion.outer_boundary(), aPath ) ;//This added for qt5 version !

  for( Hole_const_itertator hit = aRegion.holes_begin(); hit != aRegion.holes_end(); ++ hit )
    this->draw_boundary(*hit,aPath);//"This" added for qt5 version !
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_PIECEWISE_REGION_GRAPHICS_ITEM_H
