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

#ifndef CGAL_QT_PIECEWISE_SET_GRAPHICS_ITEM_H
#define CGAL_QT_PIECEWISE_SET_GRAPHICS_ITEM_H

#include "QT5/Piecewise_region_graphics_item.h"

//This class contains the implementaion of setting the classes which draws the polygons
namespace CGAL {
namespace Qt {

template <typename Piecewise_set_, typename Gps_traits, typename Draw_piece_,
          typename Piece_bbox_>
class Piecewise_set_graphics_item :
    public Piecewise_region_graphics_item< Gps_traits, Draw_piece_, Piece_bbox_>
{
  //Gps_traits class contains traits on which the operations are going to be
  // conducted

  //general polygon set
  typedef Piecewise_set_        Piecewise_set;

  typedef Draw_piece_           Draw_piece;
  typedef Piece_bbox_           Piece_bbox;

  //polygon with holes
  typedef typename Gps_traits::Polygon_with_holes_2 Region;

  // base is just used for initializing Piecewise_region_graphics_item no use
  // of it in the current class
  typedef Piecewise_region_graphics_item<Gps_traits, Draw_piece, Piece_bbox>
                                Base;

  //a container for polygon with holes
  typedef std::vector<Region> Region_vector;

  //an itertor for polygon with holes container
  typedef typename Region_vector::const_iterator Region_const_iterator;

public:
  //constructor
  Piecewise_set_graphics_item(Piecewise_set* aSet, Gps_traits Traits,
                              Draw_piece const& aPieceDrawer = Draw_piece(),
                              Piece_bbox const& aPieceBBox = Piece_bbox())
    :
    Base(aPieceDrawer,aPieceBBox),
    m_set(aSet),
    m_traits(Traits)
  {}

public:
  virtual bool isModelEmpty() const { return !m_set || m_set->is_empty(); }

protected:
  //updating the boundary of bbox
  virtual void update_bbox(Piecewise_graphics_item_base::Bbox_builder& aBboxBuilder)
  {
    if (m_set)
      update_set_bbox(*m_set, aBboxBuilder);
  }

  //for drawing polygon with holes set
  virtual void draw_model (QPainterPath& aPath)
  {
    if (m_set) draw_set(*m_set,aPath);
  }

  void update_set_bbox(Piecewise_set const& aSet,
                       Piecewise_graphics_item_base::Bbox_builder& aBboxBuilder);
  void draw_set(Piecewise_set const& aSet, QPainterPath& aPath);

protected:
  //general polygon set
  Piecewise_set* m_set;
  //Gps_traits
  Gps_traits m_traits;
};

//updating the boundary of bbox
template <typename S, typename P, typename D, typename F>
void Piecewise_set_graphics_item<S, P, D, F>::
update_set_bbox(Piecewise_set const& aSet,
                Piecewise_graphics_item_base::Bbox_builder& aBboxBuilder)
{
  Region_vector vec;

  aSet.polygons_with_holes(std::back_inserter(vec));

  //cout<<"inside update_set_bbox"<<endl;
  for (auto rit = vec.begin(); rit != vec.end(); ++rit)
    this->update_region_bbox(*rit,aBboxBuilder);
}

//for drawing polygon with holes set
template <typename S, typename P, typename D, typename F>
void Piecewise_set_graphics_item<S, P, D, F>::
draw_set(Piecewise_set const& aSet, QPainterPath& aPath)
{
  Region_vector vec;

  aSet.polygons_with_holes(std::back_inserter(vec));

  //cout<<"inside draw_set of psr"<<endl;
  for (auto rit = vec.begin(); rit != vec.end(); ++rit)
    this->draw_region(*rit,aPath);
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_PIECEWISE_SET_GRAPHICS_ITEM_H
