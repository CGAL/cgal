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

#include <list>

#include <boost/shared_ptr.hpp>

#include <QBrush>
#include <QPen>

namespace CGAL {

namespace Qt {

template <typename Bezier_polygon_wiht_holes_>
class BezierPolygonWithHolesGraphicsItem : public GraphicsItem
{
  typedef Bezier_polygon_wiht_holes_ Bezier_polygon_wiht_holes ;
  
  typedef boost::shared_ptr<Bezier_polygon_wiht_holes> Bezier_polygon_wiht_holes_ptr ;
  
  typedef std::list<Bezier_polygon_wiht_holes_ptr> Bezier_polygon_wiht_holes_list ;
  
  typedef typename Bezier_polygon_wiht_holes_list::iterator list_iterator ;
  
  typedef typename Bezier_polygon_wiht_holes::General_polygon_2 Bezier_polygon ;
  
  typedef typename Bezier_polygon::Traits Traits;
  
  typedef typename Traits::Point_2 Point;
  
public:

  BezierPolygonWithHolesGraphicsItem( Bezier_polygon_wiht_holes_list* aList );

  void modelChanged();

public:

  QRectF boundingRect() const { return mBounding_rect ; }
  
  void paint(QPainter* aPainter, const QStyleOptionGraphicsItem* aOption, QWidget* aWidget);
  
  const QBrush& brush() const { return mBrush; }
  
  void setBrush(const QBrush& aBrush ) { mBrush = aBrush; }

  const QPen& pen() const{ return mPen; }

  void setPen(const QPen& aPen) { mPen = aPen; }

protected:

  void updateBoundingBox();

protected:

  Bezier_polygon_wiht_holes_list* mList;
  QPainter*                       mPainter;
  PainterOstream<Traits>          mPainter_ostream;
  Point                           mP;
  QRectF                          mBounding_rect;
  QBrush                          mBrush;
  QPen                            mPen;
};


template <typename Bezier_polygon_wiht_holes>
BezierPolygonWithHolesGraphicsItem<Bezier_polygon_wiht_holes>::BezierPolygonWithHolesGraphicsItem(Bezier_polygon_wiht_holes_list* aList)
  : mList(aList)
  , mPainter_ostream(0)
{
  if(mList->size() == 0)
    this->hide();
  
  updateBoundingBox();
}

template <class Bezier_polygon, class SampleOutputIterator>
SampleOutputIterator Sample_bezier_polygon( Bezier_polygon const& aBP, SampleOutputIterator aSampleOut )
{
}


template <class Bezier_polygon, class SampleOutputIterator>
SampleOutputIterator Sample_bezier_polygon( Bezier_polygon const& aBP, SampleOutputIterator aSampleOut )
{
}

template <class Bezier_polygon_wiht_holes, class SampleOutputIterator>
SampleOutputIterator Sample_bezier_polygon_with_holes( Bezier_polygon_with_holes const& aBPWH, SampleOutputIterator aSampleOut )
{
}

template <typename Bezier_polygon_wiht_holes>
void BezierPolygonWithHolesGraphicsItem<Bezier_polygon_wiht_holes>::paint( QPainter*                       aPainter
                                                                         , const QStyleOptionGraphicsItem* aOption
                                                                         , QWidget*                        aWidget
                                                                         )
{
  Converter<Traits> convert;
  
  QPainterPath lPath;
  
  std::vector<Point> lSample ;
   
  for ( list_iterator it = mList->begin() ; it != mList->end() ; ++ it )
    Sample_bezier_polygon_with_holes( *it, std::back_inserter(lSample) ) ;
      
  for ( typename std::vector<Point>::const_iterator it = lSample.begin(); it != lSample.end() ; ++ it )
  {
    Bezier_polygon_wiht_holes_ptr region = *sit ;
    QPointF firstPoint = convert(*it);
    border.moveTo(firstPoint);
    mPainter_ostream = PainterOstream<Traits>(painter);
  
    for(++it;
        it != poly->outer_boundary().vertices_end();
        ++it){
      border.lineTo(convert(*it)); 
    }
    border.lineTo(firstPoint);
  
  
    for(typename Bezier_polygon_wiht_holes::Hole_const_iterator hit = poly->holes_begin();
        hit != poly->holes_end();
        ++hit){
      typename Bezier_polygon_wiht_holes::General_polygon_2::Vertex_iterator it = hit->vertices_begin();
      QPointF firstPoint = convert(*it);
      border.moveTo(firstPoint);
      for(++it;
          it != hit->vertices_end();
          ++it){
        border.lineTo(convert(*it)); 
      }
      border.lineTo(firstPoint);
    }
    
  }
  
  aPainter->setPen  (this->Pen());
  aPainter->setBrush(this->mBrush());
  aPainter->drawPath(lPath);

  for ( typename std::vector<QPointF>::const_iterator vit = lVertices.begin() ; vit != lVertices.end(); ++ vit )
    aPainter->drawPoint(*vit);
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template <typename Bezier_polygon_wiht_holes>
void 
BezierPolygonWithHolesGraphicsItem<Bezier_polygon_wiht_holes>::updateBoundingBox()
{
  Converter<Traits> convert;
  prepareGeometryChange();
  if(poly->outer_boundary().size() == 0){
    return;
  }
  mBounding_rect = convert(poly->outer_boundary().bbox());
}


template <typename Bezier_polygon_wiht_holes>
void 
BezierPolygonWithHolesGraphicsItem<Bezier_polygon_wiht_holes>::modelChanged()
{
  if((poly->outer_boundary().size() == 0) ){
    this->hide();
  } else if((poly->outer_boundary().size() > 0) && (! this->isVisible())){
    this->show();
  }
  updateBoundingBox();
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_BEZIER_POLYGON_WITH_HOLES_GRAPHICS_ITEM_H
