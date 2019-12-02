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
//             Ronnie Gandhi <ronniegandhi19999@gmail.com>
//             Efi Fogel <efifogel@gmain.com>


#ifndef CGAL_QT_PIECEWISE_REGION_GRAPHICS_ITEM_H
#define CGAL_QT_PIECEWISE_REGION_GRAPHICS_ITEM_H

#include "Typedefs.h"
#include "QT5/Piecewise_boundary_graphics_item.h"

//This class contains the implementaion of polygon with holes.

namespace CGAL {
namespace Qt {

template <typename Gps_traits, typename Draw_piece_, typename Piece_bbox_>
class Piecewise_region_graphics_item :
    public Piecewise_boundary_graphics_item<Gps_traits, Draw_piece_,
                                            Piece_bbox_>
{
  typedef typename Gps_traits::Polygon_with_holes_2 Piecewise_region;

  // base is just used for initializing Piecewise_boundary_graphics_item no use of it in the current class
  typedef Piecewise_boundary_graphics_item<Gps_traits, Draw_piece_, Piece_bbox_>
                                Base;

  typedef Draw_piece_           Draw_piece;
  typedef Piece_bbox_           Piece_bbox;

  //for iterating the holes
  typedef typename Gps_traits::Hole_const_iterator Hole_const_itertator;

public:
  //cosntructor
  /*
  Piecewise_region_graphics_item(Gps_traits* aRegion, Draw_piece const& aPieceDrawer = Draw_piece(), Piece_bbox const& aPieceBBox = Piece_bbox())
    :
     Base(aPieceDrawer, aPieceBBox)
    ,m_traits(aRegion)
  {
      //m_traits=aRegion;
       m_region=m_traits::Polygon_with_holes_2;
  }
  */
public:
  virtual bool isModelEmpty() const
  { return !m_region || m_region->outer_boundary().size(); }

protected:

  Piecewise_region_graphics_item(Draw_piece const& aPieceDrawer = Draw_piece(),
                                 Piece_bbox const& aPieceBBox = Piece_bbox()) :
    Base(aPieceDrawer, aPieceBBox)
  {}

  //updating the boundary of bbox
  virtual void update_bbox(Piecewise_graphics_item_base::Bbox_builder& aBboxBuilder)
  {
    if (m_region) update_region_bbox(*m_region, aBboxBuilder);
  }

  //for drawing boundary of polygon with holes
  virtual void draw_model (QPainterPath& aPath)
  { if (m_region) draw_region(*m_region,aPath); }

  void update_region_bbox(Piecewise_region const& aRegion,
                          Piecewise_graphics_item_base::Bbox_builder& aBboxBuilder);
  void draw_region(Piecewise_region const& aRegion, QPainterPath& aPath);

protected:
  Piecewise_region* m_region;
  Gps_traits* m_traits;
};

//updating the boundary of bbox
template <typename R, typename D,typename  P>
void Piecewise_region_graphics_item<R, D, P>::
update_region_bbox(Piecewise_region const& aRegion,
                   Piecewise_graphics_item_base::Bbox_builder& aBboxBuilder)
{
  this->update_boundary_bbox(aRegion.outer_boundary(), aBboxBuilder);

  for (auto hit = aRegion.holes_begin(); hit != aRegion.holes_end(); ++ hit)
    this->update_boundary_bbox(*hit,aBboxBuilder);
}

//for drawing boundary of polygon with holes
template <typename R, typename D, typename P>
void Piecewise_region_graphics_item<R, D, P>::
draw_region(Piecewise_region const& aRegion, QPainterPath& aPath)
{
  this->draw_boundary(aRegion.outer_boundary(), aPath);

  for(auto hit = aRegion.holes_begin(); hit != aRegion.holes_end(); ++ hit)
    this->draw_boundary(*hit,aPath);
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_PIECEWISE_REGION_GRAPHICS_ITEM_H
