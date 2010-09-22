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
// $URL$
// $Id$
// 
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_POINTS_IN_KD_TREE_GRAPHICS_ITEM_H
#define CGAL_QT_POINTS_IN_KD_TREE_GRAPHICS_ITEM_H

#include <CGAL/Bbox_2.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Fuzzy_iso_box.h>

#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>

namespace CGAL {
namespace Qt {

template <typename KdTree>
class PointsInKdTreeGraphicsItem : public GraphicsItem
{
  typedef typename std::iterator_traits<typename KdTree::iterator>::value_type Point_2;
  typedef typename CGAL::Kernel_traits<Point_2>::Kernel Traits;
  typedef typename Traits::Iso_rectangle_2 Iso_rectangle_2;

  typedef CGAL::Fuzzy_iso_box<typename KdTree::Traits> Fuzzy_iso_box;

public:
  PointsInKdTreeGraphicsItem(KdTree* p_);

  void modelChanged();

public:
  QRectF boundingRect() const;
  
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  

  const QPen& verticesPen() const
  {
    return vertices_pen;
  }

  void setVerticesPen(const QPen& pen)
  {
    vertices_pen = pen;
  }

  bool drawVertices() const
  {
    return draw_vertices;
  }

  void setDrawVertices(const bool b)
  {
    draw_vertices = b;
    update();
  }

protected:
  void updateBoundingBox();

  KdTree * kdtree;
  QPainter* m_painter;
  PainterOstream<Traits> painterostream;


  QRectF bounding_rect;

  QPen vertices_pen;
  bool draw_vertices;
};


template <typename KdTree>
PointsInKdTreeGraphicsItem<KdTree>::PointsInKdTreeGraphicsItem(KdTree * p_)
  :  kdtree(p_), painterostream(0),  draw_vertices(true)   
{
  setVerticesPen(QPen(::Qt::red, 3.));
  if(kdtree->size() == 0){
    this->hide();
  }
  updateBoundingBox();
  setZValue(3);
  setFlag(QGraphicsItem::ItemUsesExtendedStyleOption, true);
}

template <typename KdTree>
QRectF 
PointsInKdTreeGraphicsItem<KdTree>::boundingRect() const
{
  return bounding_rect;
}




template <typename KdTree>
void 
PointsInKdTreeGraphicsItem<KdTree>::paint(QPainter *painter, 
                                    const QStyleOptionGraphicsItem *option,
                                    QWidget * widget)
{
  if(drawVertices()) {
    Converter<Traits> convert;

    std::vector<Point_2> visible_points;
    Iso_rectangle_2 isor = convert(option->exposedRect);
    std::cerr << "isor = " << isor << std::endl;
    std::cerr << "x = "  << option->exposedRect.x() << "  y = " << option->exposedRect.y()  << " w =  " << option->exposedRect.width() << " h =  " << option->exposedRect.height() << std::endl; 
    std::cerr << "bl = " << isor.vertex(0) << "  tr = " << isor.vertex(2)  << std::endl;
    Fuzzy_iso_box range(isor.vertex(0), isor.vertex(2));
    kdtree->search(std::back_inserter(visible_points), range);
    std::cerr << "draw " << visible_points.size() << " of " << kdtree->size() << " points" << std::endl;
    painter->setPen(verticesPen());
    QMatrix matrix = painter->matrix();
    painter->resetMatrix();
    for(typename std::vector<Point_2>::iterator it = visible_points.begin();
        it != visible_points.end();
        it++){
      QPointF point = matrix.map(convert(*it));
      painter->drawPoint(point);
    }
  }
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template <typename KdTree>
void 
PointsInKdTreeGraphicsItem<KdTree>::updateBoundingBox()
{
  Converter<Traits> convert;
  prepareGeometryChange();
  if(kdtree->size() == 0){
    return;
  }
  bounding_rect = convert(CGAL::bounding_box(kdtree->begin(), kdtree->end()));
}


template <typename KdTree>
void 
PointsInKdTreeGraphicsItem<KdTree>::modelChanged()
{
  if((kdtree->size() == 0) ){
    this->hide();
  } else if((kdtree->size() > 0) && (! this->isVisible())){
    this->show();
  }
  updateBoundingBox();
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_POINTS_IN_KD_TREE_GRAPHICS_ITEM_H
