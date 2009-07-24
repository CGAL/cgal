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

#ifndef CGAL_QT_BEZIER_POLYGON_WITH_HOLES_GRAPHICS_ITEM_H
#define CGAL_QT_BEZIER_POLYGON_WITH_HOLES_GRAPHICS_ITEM_H

#include <vector>
#include <list>

#include <boost/optional.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>

#include <CGAL/Bezier_polygonal_sampler_2.h>

#include <QPainter>
#include <QBrush>
#include <QPen>

namespace CGAL {

namespace Qt {

template <typename Bezier_polygon_with_holes_>
class BezierPolygonWithHolesGraphicsItem : public GraphicsItem
{
  typedef Bezier_polygon_with_holes_                            Bezier_polygon_with_holes ;
  typedef typename Bezier_polygon_with_holes::General_polygon_2 Bezier_polygon ;
  typedef std::list<Bezier_polygon_with_holes>                  Bezier_polygon_with_holes_list ;
  typedef typename Bezier_polygon::Traits_2                     Bezier_traits;
  typedef typename Bezier_traits::Point_2                       Bezier_point;
  
  typedef Simple_cartesian<double>                              Linear_kernel ;
  typedef Polygon_2<Linear_kernel>                              Linear_polygon ;
  typedef Polygon_with_holes_2<Linear_kernel>                   Linear_polygon_with_holes ;
  typedef std::list<Linear_polygon_with_holes>                  Linear_polygon_with_holes_list ;
  
  typedef typename Bezier_polygon_with_holes_list::const_iterator Bezier_pwh_list_const_iterator ;
  typedef typename Linear_polygon_with_holes_list::const_iterator Linear_pwh_list_const_iterator ;
  typedef typename Linear_polygon_with_holes::Hole_const_iterator Linear_hole_const_itertator ;
  typedef typename Linear_polygon::Vertex_const_iterator          Linear_vertex_const_iterator ;
  
  typedef Converter<Linear_kernel> ToQtConverter;
  
public:

  BezierPolygonWithHolesGraphicsItem( Bezier_polygon_with_holes_list* aList );

  void modelChanged();

public:

  bool isModelEmpty() const { return !mBList || mBList->size() == 0 || mBList->front().outer_boundary().size() == 0 ; }
  
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

  Bezier_polygon_with_holes_list* mBList;
  Linear_polygon_with_holes_list  mPList ;
  QRectF                          mBounding_rect;
  QBrush                          mBrush;
  QPen                            mPen;
  
  ToQtConverter to_Qt;
};


template <typename Bezier_polygon_with_holes>
BezierPolygonWithHolesGraphicsItem<Bezier_polygon_with_holes>::BezierPolygonWithHolesGraphicsItem(Bezier_polygon_with_holes_list* aList)
  : mBList(aList)
{
  if( ! isModelEmpty() )
       updateBoundingBox();
  else this->hide();
}


template <typename Bezier_polygon_with_holes>
void BezierPolygonWithHolesGraphicsItem<Bezier_polygon_with_holes>::dump_linear_polygon( Linear_polygon const& aPoly, QPainterPath& rPath )
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


template <typename Bezier_polygon_with_holes>
void BezierPolygonWithHolesGraphicsItem<Bezier_polygon_with_holes>::paint( QPainter*                       aPainter
                                                                         , const QStyleOptionGraphicsItem* aOption
                                                                         , QWidget*                        aWidget
                                                                         )
{
  if ( ! isModelEmpty() )
  {
    QPainterPath lPath;
    
    for ( Linear_pwh_list_const_iterator rit = mPList.begin() ; rit != mPList.end() ; ++ rit )
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
template <typename Bezier_polygon_with_holes>
void BezierPolygonWithHolesGraphicsItem<Bezier_polygon_with_holes>::updateBoundingBox()
{
  if ( ! isModelEmpty() )
  {
    prepareGeometryChange();
    
    boost::optional<Bbox_2> lBBox ;
     
    for ( Linear_pwh_list_const_iterator rit = mPList.begin() ; rit != mPList.end() ; ++ rit )
    {
      if ( !lBBox )
           lBBox =          rit->outer_boundary().bbox();
      else lBBox = *lBBox + rit->outer_boundary().bbox();
      
      for ( Linear_hole_const_itertator hit = rit->holes_begin() ; hit != rit->holes_end() ; ++ hit )
        lBBox = *lBBox + hit->bbox();
    }
    
    mBounding_rect = to_Qt(*lBBox);
  }
}


template <typename Bezier_polygon_with_holes>
void BezierPolygonWithHolesGraphicsItem<Bezier_polygon_with_holes>::modelChanged()
{
  if ( !isModelEmpty() /*&& this->isVisible()*/ )
       this->show();
  else this->hide();
  
  updateSample();
  updateBoundingBox();
  update();
}

template <typename Bezier_polygon_with_holes>
void BezierPolygonWithHolesGraphicsItem<Bezier_polygon_with_holes>::updateSample()
{
  if ( !isModelEmpty() )
  {
    mPList.clear();
    for( Bezier_pwh_list_const_iterator lit = mBList->begin(); lit != mBList->end(); ++ lit )
    {
      Bezier_polygon_with_holes const& lBezier_pwh = *lit ;
      Linear_polygon_with_holes lLinear_pwh ;
      sample_bezier_polygon_with_holes_2(lBezier_pwh, lLinear_pwh);
      mPList.push_back(lLinear_pwh);
    }  
  }
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_BEZIER_POLYGON_WITH_HOLES_GRAPHICS_ITEM_H
