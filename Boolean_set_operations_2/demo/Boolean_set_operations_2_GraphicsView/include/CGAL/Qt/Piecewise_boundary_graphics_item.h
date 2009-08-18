// Copyright (c) 2009  GeometryFactory Sarl (France).
// All rights reserved.
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
// $URL$
// $Id$
// 
//
// Author(s) : Fernando Cacciola <fernando.cacciola@geometryfactory.com>

#ifndef CGAL_QT_PIECEWISE_BOUNDARY_GRAPHICS_ITEM_H
#define CGAL_QT_PIECEWISE_BOUNDARY_GRAPHICS_ITEM_H

#include <CGAL/Qt/Piecewise_graphics_item_base.h>

namespace CGAL {

namespace Qt {

template <class Piecewise_boundary_, class Draw_piece_>
class Piecewise_boundary_graphics_item : public Piecewise_graphics_item_base
{
  typedef Piecewise_boundary_ Piecewise_boundary ;
  typedef Draw_piece_         Draw_piece ;
  
  typedef typename Piecewise_boundary::Curve_const_iterator Curve_piece_const_iterator ;

public:

  Piecewise_boundary_graphics_item( Piecewise_boundary* aBoundary, Draw_piece const& aPieceDrawer = Draw_piece() )
    :
    mBoundary      (aBoundary)
    mPieceDrawer(aPieceDrawer)
  {}  

public:

  virtual bool isModelEmpty() const { return !mBoundary || mBoundary->is_empty() ; }
  
protected:
  
  Piecewise_boundary_graphics_item( Draw_piece const& aPieceDrawer = Draw_piece() )
    :
    mBoundary   (0)
    mPieceDrawer(aPieceDrawer)
  {}  
  
  virtual void update_bbox( Bbox_builder& aBBoxBuilder)
  {
    if ( mBoundary ) 
      update_boundary_bbox(*mBoundary, aBBoxBuilder ) ;
  }    

  virtual void draw_model ( QPainterPath& aPath ) 
  {
    if ( mBoundary )
      draw_boundary(*mBoundary,aPath);  
  }

  void update_boundary_bbox( Piecewise_boundary const& aBoundary, Bbox_builder& aBBoxBuilder )
  {
    aBBoxBuilder.add( aBoundary.bbox() ) ;
  }
  
  void draw_boundary( Piecewise_boundary const& aBoundary, QPainterPath& aPath ) ;
  
protected:

  Piecewise_boundary*  mBoundary;
  Draw_piece           mDrawer ;
};

template <class B, class D>
void Piecewise_boundary_graphics_item<B,D>::draw_boundary( Piecewise_boundary const& aBoundary, QPainterPath& aPath )
{
  int c = 0 ;
  for( Curve_piece_const_iterator pit = aBoundary.curves_begin(); pit != aBoundary.curves_end(); ++ pit )
    mPieceDrawer(*pit,aPath,ToQtConverter(),c);
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_PIECEWISE_BOUNDARY_GRAPHICS_ITEM_H
