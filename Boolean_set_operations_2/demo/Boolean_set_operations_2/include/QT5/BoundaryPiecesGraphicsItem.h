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
//             Ronnie Gandhi<ronniegandhi19999@gmail.com>

#ifndef CGAL_QT_BOUNDARY_PIECES_GRAPHICS_ITEM_BEZIER_H
#define CGAL_QT_BOUNDARY_PIECES_GRAPHICS_ITEM_BEZIER_H

#include "QT5/PiecewiseGraphicsItemBase.h"

namespace CGAL {

namespace Qt {




template <class Boundary_pieces_, class Draw_piece_, class Piece_bbox_>
class Boundary_pieces_graphics_item_bezier : public Piecewise_graphics_item_base_bezier
{
  typedef Boundary_pieces_ Boundary_pieces ;
  typedef Draw_piece_      Draw_piece ;
  typedef Piece_bbox_      Piece_bbox ;
  
  typedef typename Boundary_pieces::const_iterator Piece_const_iterator ;

public:

  Boundary_pieces_graphics_item_bezier( Boundary_pieces*    aBoundary
                               , Draw_piece   const& aPieceDrawer = Draw_piece()
                               , Piece_bbox   const& aPieceBBox   = Piece_bbox()
                               )
    :
     mBoundary   (aBoundary)
    ,mPieceDrawer(aPieceDrawer)
    ,mPieceBBox  (aPieceBBox)
  {}  

public:

  virtual bool isModelEmpty() const { return !mBoundary || mBoundary->size() == 0 ; }
  
protected:
  
  virtual void update_bbox( Bbox_builder& aBboxBuilder)
  {
    if ( mBoundary ) 
    {
      for( Piece_const_iterator pit = mBoundary->begin(); pit != mBoundary->end(); ++ pit )
        aBboxBuilder.add(mPieceBBox(*pit));
    }  
  }    

  virtual void draw_model ( QPainterPath& aPath ) 
  {
    if ( mBoundary )
    {
      int c = 0 ;
      for( Piece_const_iterator pit = mBoundary->begin(); pit != mBoundary->end(); ++ pit, ++c )
        mPieceDrawer(*pit,aPath,c);
    }  
  }

protected:

  Boundary_pieces* mBoundary;
  Draw_piece       mPieceDrawer ;
  Piece_bbox       mPieceBBox ;    
};


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_BOUNDARY_PIECES_GRAPHICS_ITEM_BEZIER_H
