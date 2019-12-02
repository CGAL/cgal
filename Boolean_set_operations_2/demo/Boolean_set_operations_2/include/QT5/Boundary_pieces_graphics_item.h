// Copyright (c) 2012  Tel-Aviv University (Israel).
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
//             Ronnie Gandhi<ronniegandhi19999@gmail.com>
//             Efi Fogel <efifogel@gmain.com>

#ifndef CGAL_QT_BOUNDARY_PIECES_GRAPHICS_ITEM_H
#define CGAL_QT_BOUNDARY_PIECES_GRAPHICS_ITEM_H

#include <iostream>

#include "QT5/Piecewise_graphics_item_base.h"

using namespace std;
//This class is used by Graphics_view_circular_polygon and Graphics_view_linear_polygon. It helps them in drawing the outer boundary of a set of polygon_with_holes. It helps in GUI

namespace CGAL {
namespace Qt {

template <typename Boundary_pieces_, typename Draw_piece_, typename Piece_bbox_>
class Boundary_pieces_graphics_item : public Piecewise_graphics_item_base {
  //vector of Gps_traits::Curves_2
  typedef Boundary_pieces_ Boundary_pieces;

  typedef Draw_piece_      Draw_piece;
  typedef Piece_bbox_      Piece_bbox;

  typedef typename Boundary_pieces::const_iterator Piece_const_iterator;

public:
  //constructor
  Boundary_pieces_graphics_item(Boundary_pieces* aBoundary,
                                Draw_piece const& aPieceDrawer = Draw_piece(),
                                Piece_bbox const& aPieceBBox = Piece_bbox()) :
    m_boundary(aBoundary),
    m_piece_drawer(aPieceDrawer),
    m_piece_BBox(aPieceBBox)
  {}

public:
  virtual bool isModelEmpty() const
  { return !m_boundary || (m_boundary->size() == 0); }

protected:

  //for updating the bbox
  virtual void update_bbox( Bbox_builder& aBboxBuilder)
  {
    if (m_boundary) {
      //cout<<"update bbox of bgi"<<endl;
      for (auto pit = m_boundary->begin(); pit != m_boundary->end(); ++ pit)
        aBboxBuilder.add(m_piece_BBox(*pit));
    }
  }

  //for drawing the polygon
  virtual void draw_model ( QPainterPath& aPath )
  {
    if (m_boundary) {
      int c = 0;
      //cout<<"draw model of bgi"<<endl;
      for (auto pit = m_boundary->begin(); pit != m_boundary->end(); ++ pit, ++c)
        m_piece_drawer(*pit, aPath, c);
    }
  }

protected:
  //vector of Gps_traits::Curves_2
  Boundary_pieces* m_boundary;

  Draw_piece m_piece_drawer;
  Piece_bbox m_piece_BBox;
};


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_BOUNDARY_PIECES_GRAPHICS_ITEM_H
