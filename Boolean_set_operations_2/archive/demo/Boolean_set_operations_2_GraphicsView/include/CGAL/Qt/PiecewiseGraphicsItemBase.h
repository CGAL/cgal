// Copyright (c) 2009  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// 
//
// Author(s) : Fernando Cacciola <fernando.cacciola@geometryfactory.com>

#ifndef CGAL_QT_PIECEWISE_GRAPHICS_ITEM_BASE_H
#define CGAL_QT_PIECEWISE_GRAPHICS_ITEM_BASE_H

#include <boost/optional.hpp>
#include <boost/utility.hpp>

#include <CGAL/function_objects.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>

#include <QPainter>
#include <QBrush>
#include <QPen>

namespace CGAL {

namespace Qt {

class Piecewise_graphics_item_base : public GraphicsItem
{
protected:

  Piecewise_graphics_item_base() {}
  
public:

  void updateBoundingBox();
  
  void modelChanged()
  {
    updateBoundingBox();
    update();
  }
  
  QRectF boundingRect() const { return mBounding_rect ; }
  
  void paint(QPainter* aPainter, const QStyleOptionGraphicsItem* aOption, QWidget* aWidget);
  
  const QBrush& brush() const { return mBrush; }
  
  void setBrush(const QBrush& aBrush ) { mBrush = aBrush; }

  const QPen& pen() const{ return mPen; }

  void setPen(const QPen& aPen) { mPen = aPen; }

protected:

  typedef Converter< Simple_cartesian<double> > ToQtConverter;
  
  struct Bbox_builder
  {
    void add ( Bbox_2 const& aBbox ) 
    {
      if ( bbox )
           bbox = *bbox + aBbox;
      else bbox =         aBbox;
    }
    
    boost::optional<Bbox_2> bbox ;
  } ;

  virtual bool isModelEmpty() const = 0 ;
  
  virtual void draw_model ( QPainterPath& aPath ) = 0 ;
  
  virtual void update_bbox( Bbox_builder& aBBoxBuilder ) = 0 ;

protected:

  QRectF mBounding_rect;
  QBrush mBrush;
  QPen   mPen;
};


void Piecewise_graphics_item_base::paint( QPainter* aPainter, const QStyleOptionGraphicsItem* aOption, QWidget* aWidget )
{
  if ( ! isModelEmpty() )
  {
    QPainterPath lPath ;
    
    draw_model(lPath);
    
    aPainter->setPen  (mPen );
    aPainter->setBrush(mBrush);
    aPainter->drawPath(lPath);
  }
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
void Piecewise_graphics_item_base::updateBoundingBox()
{
  if ( ! isModelEmpty() )
  {
    prepareGeometryChange();
    
    Bbox_builder lBBoxBuilder ;
    
    update_bbox(lBBoxBuilder);
    
    if ( lBBoxBuilder.bbox ) 
    {
      ToQtConverter to_Qt ;
      mBounding_rect = to_Qt(*lBBoxBuilder.bbox);
    }  
  }
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_PIECEWISE_GRAPHICS_ITEM_BASE_H
