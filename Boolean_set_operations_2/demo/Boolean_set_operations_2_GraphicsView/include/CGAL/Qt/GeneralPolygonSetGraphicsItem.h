// Copyright (c) 2008  GeometryFactory Sarl (France).
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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/branches/experimental-packages/Polyline_simplification_2/demo/Polyline_simplification_2/include/CGAL/Qt/Polyline_simplification_2_graphics_item.h $
// $Id: Polyline_simplification_2_graphics_item.h 48710 2009-04-07 21:41:12Z fcacciola $
// 
//
// Author(s) : Fernando Cacciola <Fernando.Cacciola @geometryfactory.com>

#ifndef CGAL_QT_GENERAL_POLYGON_SET_GRAPHICS_ITEM_H
#define CGAL_QT_GENERAL_POLYGON_SET_GRAPHICS_ITEM_H

#include <vector>
#include <list>

#include <boost/optional.hpp>

#include <CGAL/function_objects.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>


#include <QPainter>
#include <QBrush>
#include <QPen>

namespace CGAL {

namespace Qt {

  namespace CGALi
  {
    template<class General_polygon_set> 
    class Polygon_set_traits 
    {
      typedef typename General_polygon_set::Base base ;
      
    public:
      
      typedef typename base::Polygon_with_holes_2    Polygon_with_holes_2 ;
      typedef typename base::Polygon_2               Polygon_2 ;
      typedef typename General_polygon_set::Traits_2 Traits_2;
    } ;
    
    template<class K, class C, class D>
    class Polygon_set_traits< Polygon_set_2<K,C,D> >
    {
      typedef Polygon_set_2<K,C,D> PS ;
      
    public:
      
      typedef typename PS::Polygon_with_holes_2 Polygon_with_holes_2 ;
      typedef typename PS::Polygon_2            Polygon_2 ;
      typedef typename PS::Traits_2             Traits_2;
    } ;
  }

template <class General_polygon_set_
         , class Linearizer_ = CGAL::Identity< typename General_polygon_set_::Polygon_with_holes_2 >
         >                                    
class GeneralPolygonSetGraphicsItem : public GraphicsItem
{
  typedef Linearizer_ Linearizer ;
  
  typedef CGALi::Polygon_set_traits<General_polygon_set_> GPS_traits ;
    
  typedef General_polygon_set_                      General_polygon_set ;
  typedef typename GPS_traits::Polygon_with_holes_2 General_polygon_with_holes ;
  typedef typename GPS_traits::Polygon_2            General_polygon ;
  typedef typename GPS_traits::Traits_2             General_traits;
  typedef typename General_traits::Point_2          General_point;
  typedef std::vector<General_polygon_with_holes>   General_polygon_with_holes_vector ;
   
  typedef typename Linearizer::result_type              Linear_polygon_with_holes ;
  typedef typename Linear_polygon_with_holes::Polygon_2 Linear_polygon ;
  typedef typename Linear_polygon::Point_2              Linear_point ;
  typedef typename Kernel_traits<Linear_point>::Kernel  Linear_kernel ;
  typedef std::vector<Linear_polygon_with_holes>        Linear_polygon_with_holes_vector ;
  
  typedef typename General_polygon_with_holes_vector::const_iterator General_pwh_const_iterator ; 
  typedef typename Linear_polygon_with_holes_vector::const_iterator  Linear_pwh_const_iterator ;
  typedef typename Linear_polygon_with_holes::Hole_const_iterator    Linear_hole_const_itertator ;
  typedef typename Linear_polygon::Vertex_const_iterator             Linear_vertex_const_iterator ;
  
  typedef Converter<Linear_kernel> ToQtConverter;
  
public:

  GeneralPolygonSetGraphicsItem( General_polygon_set* aSet,  Linearizer const& aL = Linearizer() );

  void modelChanged();

public:

  bool isModelEmpty() const { return !mBSet || mBSet->is_empty() ; }
  
  QRectF boundingRect() const { return mBounding_rect ; }
  
  void paint(QPainter* aPainter, const QStyleOptionGraphicsItem* aOption, QWidget* aWidget);
  
  const QBrush& brush() const { return mBrush; }
  
  void setBrush(const QBrush& aBrush ) { mBrush = aBrush; }

  const QPen& pen() const{ return mPen; }

  void setPen(const QPen& aPen) { mPen = aPen; }

protected:

  void updateSample();
  
  void updateBoundingBox();

  void dump_linear_polygon( Linear_polygon const& aPoly, QPainterPath& rPath ) ;
  
protected:

  General_polygon_set*             mBSet;
  Linearizer                       mLinearizer ;
  Linear_polygon_with_holes_vector mPList ;
  QRectF                           mBounding_rect;
  QBrush                           mBrush;
  QPen                             mPen;
  
  ToQtConverter to_Qt;
};


template <class General_polygon_with_holes, class Linearizer>
GeneralPolygonSetGraphicsItem<General_polygon_with_holes,Linearizer>::GeneralPolygonSetGraphicsItem(General_polygon_set* aSet, Linearizer const& aL )
  : mBSet(aSet)
  , mLinearizer(aL)
{
  updateBoundingBox();
}


template <class General_polygon_with_holes, class Linearizer>
void GeneralPolygonSetGraphicsItem<General_polygon_with_holes,Linearizer>::dump_linear_polygon( Linear_polygon const& aPoly, QPainterPath& rPath )
{
  QPointF lFirstP ;
  for ( Linear_vertex_const_iterator vit = aPoly.vertices_begin(); vit != aPoly.vertices_end() ; ++ vit )
  {
    QPointF lP = to_Qt(*vit);
    if ( vit == aPoly.vertices_begin() )
    {
      lFirstP = lP ;
      rPath.moveTo(lP);
    }
    else 
    {
      rPath.lineTo(lP);  
    }  
  }
  
  rPath.lineTo(lFirstP) ;
}


template <class General_polygon_with_holes, class Linearizer>
void GeneralPolygonSetGraphicsItem<General_polygon_with_holes,Linearizer>::paint( QPainter* aPainter, const QStyleOptionGraphicsItem* aOption, QWidget* aWidget )
{
  if ( ! isModelEmpty() )
  {
    QPainterPath lPath;
    
    for ( Linear_pwh_const_iterator rit = mPList.begin() ; rit != mPList.end() ; ++ rit )
    {
      dump_linear_polygon(rit->outer_boundary(), lPath);
      
      for ( Linear_hole_const_itertator hit = rit->holes_begin() ; hit != rit->holes_end() ; ++ hit )
      {
        dump_linear_polygon(*hit, lPath);
      }
    }
    
    aPainter->setPen  (mPen );
    aPainter->setBrush(mBrush);
    aPainter->drawPath(lPath);
  }
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template <class General_polygon_with_holes, class Linearizer>
void GeneralPolygonSetGraphicsItem<General_polygon_with_holes,Linearizer>::updateBoundingBox()
{
  if ( ! isModelEmpty() )
  {
    prepareGeometryChange();
    
    boost::optional<Bbox_2> lBBox ;
     
    for ( Linear_pwh_const_iterator rit = mPList.begin() ; rit != mPList.end() ; ++ rit )
    {
      if ( lBBox )
           lBBox = *lBBox + rit->outer_boundary().bbox();
      else lBBox =          rit->outer_boundary().bbox();
      
      for ( Linear_hole_const_itertator hit = rit->holes_begin() ; hit != rit->holes_end() ; ++ hit )
      {
        if ( lBBox )
             lBBox = *lBBox + hit->bbox();
        else lBBox =          hit->bbox();
      }  
    }
    
    if ( lBBox ) 
      mBounding_rect = to_Qt(*lBBox);
  }
}


template <class General_polygon_with_holes, class Linearizer>
void GeneralPolygonSetGraphicsItem<General_polygon_with_holes,Linearizer>::modelChanged()
{
  updateSample();
  updateBoundingBox();
  update();
}

template <class General_polygon_with_holes, class Linearizer>
void GeneralPolygonSetGraphicsItem<General_polygon_with_holes,Linearizer>::updateSample()
{
  if ( !isModelEmpty() )
  {
    mPList.clear();
    General_polygon_with_holes_vector vec ;
    mBSet->polygons_with_holes( std::back_inserter(vec) ) ;
    for( General_pwh_const_iterator lit = vec.begin(); lit != vec.end(); ++ lit )
      mPList.push_back(  mLinearizer(*lit) );
  }
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GENERAL_POLYGON_SET_GRAPHICS_ITEM_H
