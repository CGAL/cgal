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
//                 Nico Kruithof <Nico@nghk.nl>

#ifndef CGAL_QT_PERIODIC_TRIANGULATION_GRAPHICS_ITEM_H
#define CGAL_QT_PERIODIC_TRIANGULATION_GRAPHICS_ITEM_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>

#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>

namespace CGAL {
  namespace Qt {
    
    template <typename T>
    class PeriodicTriangulationGraphicsItem : public GraphicsItem
    {
      typedef typename T::Geom_traits Geom_traits;
    public:
      PeriodicTriangulationGraphicsItem(T* t_);
      
      void modelChanged();

      enum Iterator_type {
        STORED = 0,
        UNIQUE, // 1
        STORED_COVER_DOMAIN, // 2
        UNIQUE_COVER_DOMAIN, // 3
        NONE
      };
      void setEmphasizedSimplices(Iterator_type type) { this->type = type; }
      Iterator_type getEmphasizedSimplices() { return this->type; }
    public:
      
      QRectF boundingRect() const;
      
      void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
      
      virtual void operator()(typename T::Face_handle fh);
      
      const QPen& verticesPen() const
      {
        return vertices_pen;
      }
      
      const QPen& edgesPen() const
      {
        return edges_pen;
      }

      const QPen& facesPen() const
      {
        return faces_pen;
      }
      
      const QPen& domainPen() const
      {
        return domain_pen;
      }
      
      void setVerticesPen(const QPen& pen)
      {
        vertices_pen = pen;
      }
      
      void setEdgesPen(const QPen& pen)
      {
        edges_pen = pen;
      }

      void setFacesPen(const QPen& pen)
      {
        edges_pen = pen;
      }
      
      void setDomainPen(const QPen& pen)
      {
        domain_pen = pen;
      }
      
      bool visibleVertices() const
      {
        return visible_vertices;
      }
      
      void setVisibleVertices(const bool b)
      {
        visible_vertices = b;
        update();
      }
      
      bool visibleEdges() const
      {
        return visible_edges;
      }
      
      void setVisibleEdges(const bool b)
      {
        visible_edges = b;
        update();
      }
      
    protected:
      virtual void drawAll(QPainter *painter);
      void paintVertices(QPainter *painter);
      void paintOneVertex(const typename T::Point& point);
      virtual void paintVertex(typename T::Vertex_handle vh);
      void updateBoundingBox();
      
      T * t;
      QPainter* m_painter;
      PainterOstream<Geom_traits> painterostream;
      
      typename T::Vertex_handle vh;
      typename T::Point p;
      QRectF bounding_rect;
      
      QPen vertices_pen;
      QPen edges_pen;
      QPen faces_pen;
      QPen domain_pen;
      bool visible_edges;
      bool visible_vertices;

      Iterator_type type;
    };
    
    
    template <typename T>
    PeriodicTriangulationGraphicsItem<T>::PeriodicTriangulationGraphicsItem(T * t_)
    :  t(t_), painterostream(0),
       visible_edges(true), visible_vertices(true),
       type(NONE)
    {
      setVerticesPen(QPen(::Qt::red, 1.));
      setFacesPen(QPen(QColor(100,100,100)));
      setDomainPen(QPen(::Qt::blue, .01));
      if(t->number_of_vertices() == 0){
        this->hide();
      }
      updateBoundingBox();
      setZValue(3);
    }
    
    template <typename T>
    QRectF 
    PeriodicTriangulationGraphicsItem<T>::boundingRect() const
    {
      return bounding_rect;
    }
    
    
    template <typename T>
    void 
    PeriodicTriangulationGraphicsItem<T>::operator()(typename T::Face_handle fh)
    {
      if(visible_edges) {
        for (int i=0; i<3; i++) {
          if (fh < fh->neighbor(i)){
            m_painter->setPen(this->edgesPen());
            painterostream << t->segment(fh,i);
          }
        }
      }
      if(visible_vertices) {
        for (int i=0; i<3; i++) {
          paintVertex(fh->vertex(i));
        }
      }
    }
    
    template <typename T>
    void 
    PeriodicTriangulationGraphicsItem<T>::drawAll(QPainter *painter)
    {
      painterostream = PainterOstream<Geom_traits>(painter);

      if (type != NONE)
      {
        typename T::Iterator_type itype;
        switch (type) {
        case STORED:
          itype = T::STORED;
          break;
        case UNIQUE:
          itype = T::UNIQUE;
          break;
        case STORED_COVER_DOMAIN:
          itype = T::STORED_COVER_DOMAIN;
          break;
        case UNIQUE_COVER_DOMAIN:
          itype = T::UNIQUE_COVER_DOMAIN;
          break;
        case NONE:
        default:
          assert(false);
          itype = T::STORED;
          break;
        }
        
        Converter<Geom_traits> convert;

        QMatrix matrix = painter->matrix();
        painter->resetMatrix();

        { // Faces
          painter->setPen(QPen());
          painter->setBrush(QBrush(::Qt::green));
          for (typename T::Periodic_triangle_iterator tit = t->periodic_triangles_begin(itype);
               tit != t->periodic_triangles_end(itype); ++tit) {
            painter->drawConvexPolygon(matrix.map(convert(t->triangle(*tit))));
          }
          painter->setBrush(QBrush());
        }
        { // Edges
          QPen pen = edgesPen();
          pen.setWidth(pen.width() + 2);
          painter->setPen(pen);
          for (typename T::Periodic_segment_iterator sit = t->periodic_segments_begin(itype);
               sit != t->periodic_segments_end(itype); ++sit) {
            painter->drawLine(matrix.map(convert(t->segment(*sit))));
          }
        }
        { // Vertices
          QPen pen = verticesPen();
          pen.setWidth(pen.width() + 2);
          painter->setPen(pen);
          for (typename T::Periodic_point_iterator pit = t->periodic_points_begin(itype);
               pit != t->periodic_points_end(itype); ++pit) {
            painter->drawPoint(matrix.map(convert(t->point(*pit))));
          }
        }

        painter->setMatrix(matrix);
      }      

      if(visibleEdges()) {
        painter->setPen(this->edgesPen());
        t->draw_triangulation(painterostream);
      }
      
      paintVertices(painter);
    }
    
    template <typename T>
    void 
    PeriodicTriangulationGraphicsItem<T>::paintVertices(QPainter *painter)
    {
      if(visibleVertices()) {
        Converter<Geom_traits> convert;
        
        QMatrix matrix = painter->matrix();
        painter->resetMatrix();

        QPen pen = verticesPen();
        if (t->number_of_vertices() < 8) {
          int v_index=1;
          for (typename T::Unique_vertex_iterator vit = t->unique_vertices_begin();
               vit != t->unique_vertices_end(); ++vit) {
            pen.setColor(QColor(255*(v_index&1), 255*((v_index>>1)&1), 255*((v_index>>2)&1)));
            painter->setPen(pen);

            painter->drawPoint(matrix.map(convert(t->point(vit))));
            std::vector<typename T::Vertex_handle> copies = t->periodic_copies(vit);
            for (size_t i=0; i<copies.size(); ++i)
              painter->drawPoint(matrix.map(convert(t->point(copies[i]))));

            ++v_index;
          }
          
        } else {
          painter->setPen(verticesPen());
          for (typename T::Periodic_point_iterator ppit = t->periodic_points_begin();
               ppit != t->periodic_points_end(); ++ppit)
            {
              QPointF point = matrix.map(convert(t->point(*ppit)));
              painter->drawPoint(point);
            }
        }

        painter->setMatrix(matrix);
      }
    }
    
    template <typename T>
    void 
    PeriodicTriangulationGraphicsItem<T>::paintOneVertex(const typename T::Point& point)
    {
      Converter<Geom_traits> convert;
      
      m_painter->setPen(this->verticesPen());
      QMatrix matrix = m_painter->matrix();
      m_painter->resetMatrix();
      m_painter->drawPoint(matrix.map(convert(point)));
      m_painter->setMatrix(matrix);
    }
    
    template <typename T>
    void 
    PeriodicTriangulationGraphicsItem<T>::paintVertex(typename T::Vertex_handle vh)
    {
      Converter<Geom_traits> convert;
      
      m_painter->setPen(this->verticesPen());
      QMatrix matrix = m_painter->matrix();
      m_painter->resetMatrix();
      m_painter->drawPoint(matrix.map(convert(vh->point())));
      m_painter->setMatrix(matrix);
    }
    
    template <typename T>
    void 
    PeriodicTriangulationGraphicsItem<T>::paint(QPainter *painter, 
                                                const QStyleOptionGraphicsItem *,
                                                QWidget *)
    {
      drawAll(painter);

      painter->setPen(this->domainPen());
      const typename Geom_traits::Iso_rectangle_2 &domain = t->domain();
      double dx = domain.xmax()-domain.xmin();
      double dy = domain.ymax()-domain.ymin();
      typename T::Covering_sheets sheets = t->number_of_sheets();
      for (int x=0; x<sheets[0]; ++x) {
        for (int y=0; y<sheets[1]; ++y) {
          painter->drawRect((int)(domain.xmin() + x*dx),
                            (int)(domain.ymin() + y*dy),
                            (int)dx, (int)dy);
        }
      }
      m_painter = painter;
    }
    
    // We let the bounding box only grow, so that when vertices get removed
    // the maximal bbox gets refreshed in the GraphicsView
    template <typename T>
    void 
    PeriodicTriangulationGraphicsItem<T>::updateBoundingBox()
    {
      prepareGeometryChange();
      
      CGAL::Bbox_2 bb = t->domain().bbox();
      for (typename T::Periodic_triangle_iterator tit = t->periodic_triangles_begin(T::STORED_COVER_DOMAIN);
           tit != t->periodic_triangles_end(T::STORED_COVER_DOMAIN); ++tit) {
        bb = bb + t->triangle(*tit).bbox();
      }
      
      double xmin = bb.xmin();
      double ymin = bb.ymin();
      double dx = bb.xmax() - xmin;
      double dy = bb.ymax() - ymin;

      double delta = 0.05;
      xmin -= delta * dx;
      ymin -= delta * dy;
      dx += 2 * delta * dx;
      dy += 2 * delta * dy;

      bounding_rect = QRectF(xmin, ymin, dx, dy);
    }
    
    
    template <typename T>
    void 
    PeriodicTriangulationGraphicsItem<T>::modelChanged()
    {
      if((t->number_of_vertices() == 0) ){
        this->hide();
      } else if((t->number_of_vertices() > 0) && (! this->isVisible())){
        this->show();
      }
      updateBoundingBox();
      update();
    }
    
    
  } // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_PERIODIC_TRIANGULATION_GRAPHICS_ITEM_H
