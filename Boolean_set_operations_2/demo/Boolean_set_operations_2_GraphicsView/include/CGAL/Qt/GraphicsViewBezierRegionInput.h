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
// Author(s)     : Fernando Cacciola <Fernando.Cacciola@geometryfactory.com>

#ifndef CGAL_QT_GRAPHICS_VIEW_BEZIER_REGION_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_BEZIER_REGION_INPUT_H

#include <list>

#include <QGraphicsView>
#include <QGraphicsItem>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>

#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/GraphicsViewBezierBoundaryInput.h>
#include <CGAL/array.h>

namespace CGAL {
  namespace Qt {


    template <typename Traits_>
    class GraphicsViewBezierRegionInput : public GraphicsViewInput
    {
    public:

      typedef Traits_ Traits ;
      
      typedef CGAL::Gps_traits_2<Traits> Bezier_gps_traits;
      
      typedef typename Traits::Curve_2                      Curve;
      typedef typename Traits::X_monotone_curve_2           X_monotone_curve;
      typedef typename Traits::Point_2                      Point;
      typedef typename Bezier_gps_traits::General_polygon_2 Boundary;
          typedef Bezier_region_ Bezier_region ;

      typedef typename Bezier_region::General_polygon_2 Bezier_boundary;

      typedef Point_2< Simple_cartesian<double> > Linear_point ;

      GraphicsViewBezierRegionInput(QObject *parent, QGraphicsScene* s); 
      ~GraphicsViewBezierRegionInput();

      public slots:

        void processInput(CGAL::Object o);


    protected:

      virtual void keyPressEvent(QKeyEvent *event);

      bool eventFilter(QObject *obj, QEvent *event);

    private:

      Bezier_boundary bounday;
      std::vector<Bezier_boundary> boundaries; 
      Bezier_region br;  // this one collects the input polygons

      Bezier_region_graphics_item<Bezier_region> * brItem;
      GraphicsViewBezierBoundaryInput<Bezier_boundary> * bi;

      bool polygon_input;
      QGraphicsScene *scene_;  
    };


    template <class R>
    GraphicsViewBezierRegionInput<R>::GraphicsViewBezierRegionInput(QObject *parent, QGraphicsScene* s)
      : GraphicsViewInput(parent), scene_(s), polygon_input(false)
    {
      brItem = new Bezier_region_graphics_item<Bezier_region>(&br);
      brItem->setBrush(::Qt::yellow);
      scene_->addItem(brItem);
      brItem->hide();

      bi = new GraphicsViewBezierBoundaryInput<Bezier_boundary>(parent,s);
      QObject::connect(bi, SIGNAL(generate(CGAL::Object)),
        this, SLOT(processInput(CGAL::Object)));

      QObject::connect(this, SIGNAL(modelChanged()),
        brItem, SLOT(modelChanged()));

    }

    template <class R>
    GraphicsViewBezierRegionInput<R>::~GraphicsViewBezierRegionInput()
    {
      //delete brItem;
      //delete bi;
    }


    template <class R>
    void GraphicsViewBezierRegionInput<R>::processInput(CGAL::Object o)
    {
      std::vector<Linear_point> points;
      if(CGAL::assign(points, o))
      {
        if((points.size() == 1)&& polygon.size()>0)
        {

        } else 
        {
          polygon.clear();
          if(points.front() == points.back())
          {
            points.pop_back();
          }
          polygon.insert(polygon.vertices_begin(), points.begin(), points.end());
          if(boundaries.empty())
          {
            if(polygon.orientation() == CGAL::CLOCKWISE)
            {
              polygon.reverse_orientation();
            }
          } else 
          {
            if(polygon.orientation() == CGAL::COUNTERCLOCKWISE)
            {
              polygon.reverse_orientation();
            }
          }
          boundaries.push_back(polygon);
          typename std::list<Bezier_boundary>::iterator it = boundaries.begin();
          it++;
          pwh = Bezier_region(boundaries.front(), it, boundaries.end());
        }
        emit(modelChanged());
        polygon_input = false;
      } 
    }


    template <class R>
    void GraphicsViewBezierRegionInput<R>::keyPressEvent ( QKeyEvent * event ) 
    {
    }



    template <class R>
    bool 
      GraphicsViewBezierRegionInput<R>::eventFilter(QObject *obj, QEvent *event)
    {
      if(polygon_input){
        return bi->eventFilter(obj, event);
      } else {
        if (event->type() == QEvent::GraphicsSceneMousePress) {
          QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);

          if(mouseEvent->modifiers()  & ::Qt::ShiftModifier){
            return QObject::eventFilter(obj, event);;
          }
          if(mouseEvent->button() == ::Qt::LeftButton) {
            polygon_input = true;
            return bi->eventFilter(obj, event);
          } else if(mouseEvent->button() == ::Qt::RightButton) {
            emit(generate(CGAL::make_object(pwh)));
            pwh.clear();
            boundaries.clear();
            polygon_input = false;
            emit(modelChanged());
          }
          return true;
        } else if (event->type() == QEvent::KeyPress) {
          QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);
          keyPressEvent(keyEvent);
          return true;
        } else{
          // standard event processing
          return QObject::eventFilter(obj, event);
        }
      }
    } 

  } // namespace Qt

} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_VIEW_BEZIER_REGION_INPUT_H
