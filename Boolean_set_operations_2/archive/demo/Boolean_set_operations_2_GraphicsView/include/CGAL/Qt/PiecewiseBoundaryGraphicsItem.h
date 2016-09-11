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
// 
//
// Author(s) : Fernando Cacciola <fernando.cacciola@geometryfactory.com>

#ifndef CGAL_QT_PIECEWISE_BOUNDARY_GRAPHICS_ITEM_H
#define CGAL_QT_PIECEWISE_BOUNDARY_GRAPHICS_ITEM_H

#include <CGAL/Qt/PiecewiseGraphicsItemBase.h>

namespace CGAL {

namespace Qt {

template <class Piecewise_boundary_, class Draw_piece_, class Piece_bbox_>
class Piecewise_boundary_graphics_item : public Piecewise_graphics_item_base
{
  typedef Piecewise_boundary_ Piecewise_boundary ;
  typedef Draw_piece_         Draw_piece ;
  typedef Piece_bbox_         Piece_bbox ;
  
  typedef typename Piecewise_boundary::Curve_const_iterator Curve_piece_const_iterator ;

public:

  Piecewise_boundary_graphics_item( Piecewise_boundary* aBoundary
                                  , Draw_piece   const& aPieceDrawer = Draw_piece()
                                  , Piece_bbox   const& aPieceBBox   = Piece_bbox()
                                  )
    :
    mBoundary   (aBoundary)
   ,mPieceDrawer(aPieceDrawer)
   ,mPieceBBox  (aPieceBBox)
  {}  

public:

  virtual bool isModelEmpty() const { return !mBoundary || mBoundary->is_empty() ; }
  
protected:
  
  Piecewise_boundary_graphics_item( Draw_piece const& aPieceDrawer = Draw_piece() 
                                  , Piece_bbox const& aPieceBBox   = Piece_bbox()
                                  )
    :
    mBoundary   (0)
   ,mPieceDrawer(aPieceDrawer)
   ,mPieceBBox  (aPieceBBox)
  {}  
  
  virtual void update_bbox( Bbox_builder& aBboxBuilder)
  {
    if ( mBoundary ) 
      update_boundary_bbox(*mBoundary, aBboxBuilder ) ;
  }    

  virtual void draw_model ( QPainterPath& aPath ) 
  {
    if ( mBoundary )
      draw_boundary(*mBoundary,aPath);  
  }

  void update_boundary_bbox( Piecewise_boundary const& aBoundary, Bbox_builder& aBboxBuilder )
  {
    for( Curve_piece_const_iterator pit = aBoundary.curves_begin(); pit != aBoundary.curves_end(); ++ pit )
      aBboxBuilder.add(mPieceBBox(*pit));
  }
  
  void draw_boundary( Piecewise_boundary const& aBoundary, QPainterPath& aPath ) ;
  
protected:

  Piecewise_boundary* mBoundary;
  Draw_piece          mPieceDrawer ;
  Piece_bbox          mPieceBBox ;    
};

template <class B, class D, class P>
void Piecewise_boundary_graphics_item<B,D,P>::draw_boundary( Piecewise_boundary const& aBoundary, QPainterPath& aPath )
{
  int c = 0 ;
  for( Curve_piece_const_iterator pit = aBoundary.curves_begin(); pit != aBoundary.curves_end(); ++ pit, ++c )
    mPieceDrawer(*pit,aPath,c);
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_PIECEWISE_BOUNDARY_GRAPHICS_ITEM_H
