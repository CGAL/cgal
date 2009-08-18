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

#ifndef CGAL_QT_PIECEWISE_SET_GRAPHICS_ITEM_H
#define CGAL_QT_PIECEWISE_SET_GRAPHICS_ITEM_H

#include <CGAL/Qt/Piecewise_region_graphics_item.h>

namespace CGAL {

namespace Qt {

  namespace CGALi
  {
    template<class Piecewise_set> 
    struct Piecewise_set_traits 
    {
      typedef typename Piecewise_set::Base base ;
      
      typedef typename base::Polygon_with_holes_2 Region ;
    } ;
    
    template<class K, class C, class D>
    struct Piecewise_set_traits< Polygon_set_2<K,C,D> >
    {
      typedef Polygon_set_2<K,C,D> PS ;
      
      typedef typename PS::Polygon_with_holes_2 Region ;
    } ;
  }

template <class Piecewise_set_, class Draw_piece_>
class Piecewise_set_graphics_item : public Piecewise_region_graphics_item< typename CGALi::Piecewise_set_traits<Piecewise_set_>::Region, Draw_piece_ > 
{
  typedef Piecewise_set_ Piecewise_set ;
  typedef Draw_piece_    Draw_piece ;
  
  typedef typename CGALi::Piecewise_set_traits<Piecewise_set_>::Region Region ;
  
  typedef Piecewise_region_graphics_item<Region, Draw_piece> Base ;
 
  typedef std::vector<Region> Region_vector ;

  typedef typename Region_vector::const_iterator Region_const_iterator ;
  
public:

  Piecewise_set_graphics_item( Piecewise_set* aSet, Draw_piece const& aPieceDrawer = Draw_piece() )
    :
     Base(aPieceDrawer)
    ,mSet(aSet)
  {}  

public:

  virtual bool isModelEmpty() const { return !mSet || mSet->is_empty() ; }
  
protected:
  
  virtual void update_bbox( Bbox_builder& aBBoxBuilder)
  {
    if ( mSet ) 
      update_set_bbox(*mSet, aBBoxBuilder ) ;
  }    

  virtual void draw_model ( QPainterPath& aPath ) 
  {
    if ( mSet )
      draw_set(*mSet,aPath);  
  }

  void update_set_bbox( Piecewise_set const& aSet, Bbox_builder& aBBoxBuilder ) ;
  void draw_set       ( Piecewise_set const& aSet, QPainterPath& aPath ) ;
  
protected:

  Piecewise_set* mSet;
};

template <class S, class D>
void Piecewise_set_graphics_item<S,D>::update_set_bbox( Piecewise_set const& aSet, Bbox_builder& aBBoxBuilder )
{
  Region_vector vec ;
  
  aSet.polygons_with_holes( std::back_inserter(vec) ) ;
  
  for( Region_const_iterator rit = vec.begin(); rit != vec.end() ; ++ rit )
    update_region_bbox(*rit,aBBoxBuilder);
}

template <class S, class D>
void Piecewise_set_graphics_item<S,D>::draw_set( Piecewise_set const& aSet, QPainterPath& aPath )
{
  Region_vector vec ;
  
  aSet.polygons_with_holes( std::back_inserter(vec) ) ;
  
  for( Region_const_iterator rit = vec.begin(); rit != vec.end() ; ++ rit )
    draw_region(*rit,aPath);
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_PIECEWISE_SET_GRAPHICS_ITEM_H
