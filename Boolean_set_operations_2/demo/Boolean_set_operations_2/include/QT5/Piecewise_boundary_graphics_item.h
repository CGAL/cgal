// Copyright (c) 2012,2018  Tel-Aviv University (Israel).
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
// Author(s) : Apurva Bhatt <response2apurva@gmail.com>
//             Ronnie Gandhi <ronniegandhi19999@gmail.com>
//             Efi Fogel <efifogel@gmain.com>

#ifndef CGAL_QT_PIECEWISE_BOUNDARY_GRAPHICS_ITEM_H
#define CGAL_QT_PIECEWISE_BOUNDARY_GRAPHICS_ITEM_H

#include <CGAL/iterator.h>

#include "QT5/Piecewise_graphics_item_base.h"
#include "Typedefs.h"

//This class contains implementaion of drawing the polygons
namespace CGAL {
namespace Qt {

template <typename Gps_traits, class Draw_piece_, typename Piece_bbox_>
class Piecewise_boundary_graphics_item : public Piecewise_graphics_item_base {
  //for iterating the polygons
  typedef typename Gps_traits::Polygon_2 Piecewise_boundary;

  typedef Draw_piece_         Draw_piece;
  typedef Piece_bbox_         Piece_bbox;

  typedef typename Gps_traits::Curve_const_iterator  Curve_piece_const_iterator;

public:
  //to check if there are any polygon to draw or not
  virtual bool isModelEmpty() const
  { return !m_boundary || m_boundary->is_empty(); }

protected:
  //preotected constructor and default
  Piecewise_boundary_graphics_item(Draw_piece const& aPieceDrawer = Draw_piece(),
                                   Piece_bbox const& aPieceBBox = Piece_bbox()) :
    m_boundary(0),
    m_piece_drawer(aPieceDrawer),
    m_piece_BBox(aPieceBBox)
  {}

  //for updating bbox
  virtual void update_bbox(Bbox_builder& aBboxBuilder)
  {
    if (m_boundary) update_boundary_bbox(*m_boundary, aBboxBuilder);
  }

  //for drawing a boundary
  virtual void draw_model(QPainterPath& aPath)
  {
    if (m_boundary) draw_boundary(*m_boundary,aPath);
  }

  //updating the boundary of bbox
  void update_boundary_bbox(Piecewise_boundary const& aBoundary,
                            Bbox_builder& aBboxBuilder)
  {
    for (auto pit = aBoundary.curves_begin(); pit != aBoundary.curves_end();
         ++pit)
      aBboxBuilder.add(m_piece_BBox(*pit));
  }

  void draw_boundary(Piecewise_boundary const& aBoundary, QPainterPath& aPath);

protected:
  //boundary of the polygon
  Piecewise_boundary* m_boundary;
  Gps_traits* m_traits;

  Draw_piece m_piece_drawer;
  Piece_bbox m_piece_BBox;
};

//for drawing boundary of each polygon
template <class B, class D, class P>
void Piecewise_boundary_graphics_item<B,D,P>::
draw_boundary(Piecewise_boundary const& aBoundary, QPainterPath& aPath)
{
  int c = 0;
  for (auto pit = aBoundary.curves_begin(); pit != aBoundary.curves_end();
       ++pit, ++c)
    m_piece_drawer(*pit,aPath,c);
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_PIECEWISE_BOUNDARY_GRAPHICS_ITEM_H
