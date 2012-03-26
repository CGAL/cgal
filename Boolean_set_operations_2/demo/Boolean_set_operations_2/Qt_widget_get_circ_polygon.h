// Copyright (c) 2005  Tel-Aviv University (Israel).
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
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_QT_WIDGET_GET_CIRC_POLYGON_H
#define CGAL_QT_WIDGET_GET_CIRC_POLYGON_H

#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/General_polygon_2.h>
#include "Qt_widget_circle_segment_2.h"
#include "Qt_widget_X_monotone_circle_segment_2.h"
#include <qcursor.h>

namespace CGAL
{
  template <class Kernel>
  class  Qt_widget_get_circ_polygon : public Qt_widget_layer
  {
  protected:
    typedef Arr_circle_segment_traits_2<Kernel>          Traits_2;
    typedef General_polygon_2<Traits_2>                  Polygon;
    typedef typename Traits_2::Point_2                   Arc_point_2;

    typedef typename Traits_2::Curve_2                   Curve_2;
    typedef typename Traits_2::X_monotone_curve_2        X_monotone_curve_2;
    typedef typename Polygon::Curve_iterator             Curve_iterator;
    typedef typename Kernel::FT                          FT;
    typedef typename Kernel::Point_2                     Point_2;
    typedef typename Kernel::Segment_2                   Segment_2;
    typedef typename Kernel::Circle_2                    Circle_2;
    typedef Cartesian<double>                            Double_kernel;
    typedef Double_kernel::Point_2                       Double_point_2;
    typedef Double_kernel::Segment_2                     Double_segment_2;


    //Data members
    Polygon     m_pgn;
    bool        m_active;      //true if the first point was inserted
    bool        m_first_time; //true if it is the first time when
                              //draw the rubber band

    bool        m_is_circ_mode; // true if we are at circular arc mode
    Point_2     m_rubber;   //the new point of the rubber band
    Point_2     m_last_of_poly; //the last point of the polygon
    Point_2     m_rubber_old; //the old point of the rubber band

    Point_2     m_arc_source;  // the arc source (incase of circ arc)
    Point_2     m_arc_target; // the arc target (incase of circ arc)
    Curve_2     m_old_arc;

    QWidget::FocusPolicy  m_oldpolicy;
    QCursor               m_oldcursor;
    QCursor               m_cursor;

    bool m_ignore_move_event;


  public:

    Qt_widget_get_circ_polygon(const QCursor c=QCursor(Qt::crossCursor),
                               QObject* parent = 0,
                               const char* name = 0)
      : Qt_widget_layer(parent, name),
        m_active(false),
        m_first_time(true),
        m_is_circ_mode(false),
        m_cursor(c) ,
        m_ignore_move_event(false)
    {}

    void draw()
    {
      if(m_pgn.size() > 1)
      {
        widget->lock();
        RasterOp old_rasterop = widget->rasterOp();
        widget->get_painter().setRasterOp(XorROP);
        *widget << CGAL::GREEN;
        Curve_iterator before_end = m_pgn.curves_end();
        --before_end;

        for(Curve_iterator it = m_pgn.curves_begin(); it != before_end; ++it)
          *widget << *it;
        widget->setRasterOp(old_rasterop);
        widget->unlock();
      }
      return;
    };


  protected:

  bool is_pure(Qt::ButtonState s)
  {
    if((s & Qt::ControlButton) ||
       (s & Qt::ShiftButton) ||
       (s & Qt::AltButton))
      return false;

    return true;
  }


  void mousePressEvent(QMouseEvent *e)
  {
    if(!is_pure(e->state()) && (e->state() & Qt::ControlButton))
    {
      if(m_is_circ_mode)
        return;
      if(!m_active)
        return;

      FT x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);

       if (m_last_of_poly == Point_2(x,y))
          return;

       m_is_circ_mode = true;
      CircModeEvent(e);
      return;
    }
    if(e->button() == Qt::LeftButton && is_pure(e->state()))
    {
      FT x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);

      if(!m_active)
      {
        m_active=true;
        m_last_of_poly = Point_2(x, y);
      }
      else
      {
        if (m_last_of_poly == Point_2(x,y))
          return;
        m_rubber_old = Point_2(x, y);
        if(is_simple())
        {
          if(!m_is_circ_mode)
          {
            m_pgn.push_back(X_monotone_curve_2(m_last_of_poly, Point_2(x,y)));
            //show the last rubber as edge of the polygon
            widget->lock();
            RasterOp old_rasterop=widget->rasterOp();
            widget->get_painter().setRasterOp(XorROP);
            *widget << CGAL::WHITE;
            *widget << Segment_2(m_rubber, m_last_of_poly);
            *widget << CGAL::GREEN;
            *widget << Segment_2(m_rubber, m_last_of_poly);
            widget->setRasterOp(old_rasterop);
            widget->unlock();
            m_last_of_poly = Point_2(x, y);
          }
          else
          {
            //circ mode
             Traits_2 tr;
             typename Traits_2::Make_x_monotone_2 make_x = tr.make_x_monotone_2_object();
             std::vector<Object> xcurves_vec;
             xcurves_vec.reserve(2);
             make_x(m_old_arc, std::back_inserter(xcurves_vec));

             widget->lock();
             RasterOp old_rasterop=widget->rasterOp();
             widget->get_painter().setRasterOp(XorROP);
             *widget << CGAL::WHITE;
             *widget << m_old_arc;
             *widget << CGAL::GREEN;
              for(unsigned int i=0; i<xcurves_vec.size(); ++i)
             {
               X_monotone_curve_2 cv_arc;
               CGAL::assign(cv_arc, xcurves_vec[i]);
               m_pgn.push_back(cv_arc);
               *widget << (cv_arc);
             }
             widget->setRasterOp(old_rasterop);
             widget->unlock();
             m_last_of_poly = Point_2(x, y);
             m_is_circ_mode = false;
             m_last_of_poly = m_arc_target;
             m_rubber_old = m_arc_target;
             int xpixel = widget->x_pixel(CGAL::to_double(m_arc_target.x()));
             int ypixel = widget->y_pixel(CGAL::to_double(m_arc_target.y()));

             QPoint qp(xpixel, ypixel);
             QPoint qq = widget->mapToGlobal(qp);
             QCursor::setPos(qq);
          }
        }
      }
      return;
    }
    if(e->button() == Qt::RightButton && is_pure(e->state()))
    {
      if (m_active)
      {
        if(m_is_circ_mode)
          return; // if we are at circ mode, ignore
        if(is_simple(true))
        {
          if(m_pgn.is_empty())
            return;
          const Arc_point_2& first_point = m_pgn.curves_begin()->source();
          CGAL_assertion(!first_point.x().is_extended() && !first_point.y().is_extended());
          FT xs = first_point.x().alpha();
          FT ys = first_point.y().alpha();
          m_pgn.push_back(X_monotone_curve_2(m_last_of_poly, Point_2(xs, ys)));
          widget->new_object(make_object(m_pgn));
          m_active = false;
          m_first_time = true;
          m_pgn.clear();
        }
      }
    }
  };//end mousePressEvent


  void keyPressEvent(QKeyEvent *e)
  {
    switch ( e->key() )
    {
      case Key_Escape:

          if(m_is_circ_mode)
          {
            // special treatment if we are at circ mode
            widget->lock();
            RasterOp old_rasterop=widget->rasterOp();
            widget->get_painter().setRasterOp(XorROP);
            *widget << CGAL::WHITE;
            *widget << m_old_arc;
            *widget << Segment_2(m_rubber, m_arc_source);
            m_is_circ_mode = false;

            //move the cursor to m_rubber position
            int rubber_x_pixel = widget->x_pixel(CGAL::to_double(m_rubber.x()));
            int rubber_y_pixel = widget->y_pixel(CGAL::to_double(m_rubber.y()));

            QPoint qp(rubber_x_pixel, rubber_y_pixel);
            QPoint qq = widget->mapToGlobal(qp);
            QCursor::setPos(qq);
            m_ignore_move_event = true;

            widget->setRasterOp(old_rasterop);
            widget->unlock();

            return;
          }

         if(!m_pgn.is_empty())
         {
           // segment mode
           widget->lock();
           RasterOp old_rasterop=widget->rasterOp();
           widget->get_painter().setRasterOp(XorROP);
           *widget << CGAL::GREEN;

           Curve_iterator last = m_pgn.curves_end();
           --last;

           if(last->is_linear())
           {
            *widget << *last;
            *widget << CGAL::WHITE;
            *widget << Segment_2(m_rubber, m_last_of_poly);
            const Arc_point_2& last_point = last->source();

            CGAL_assertion(!last_point.x().is_extended() && !last_point.y().is_extended());
            FT xs = last_point.x().alpha();
            FT ys = last_point.y().alpha();

            *widget << Segment_2(m_rubber, Point_2(xs, ys));
            widget->setRasterOp(old_rasterop);
            widget->unlock();
            m_last_of_poly = Point_2(xs, ys);
            m_pgn.erase(last);
           }
           else
           {
             //circular arc, remove all original xcurves

             Curve_iterator curr = last;
             Curve_iterator prev = curr;
             while(curr != m_pgn.curves_begin())
             {
               --curr;
               if(curr->has_same_supporting_curve(*prev))
               {
                 prev = curr;
               }
               else
                 break;
             }
             Curve_iterator first;
             if(curr == m_pgn.curves_begin())
             {
               Curve_iterator next = curr;
               ++next;
               if(next == m_pgn.curves_end())
               {
                 first = curr;
               }
               else
                 if(curr->has_same_supporting_curve(*next))
                   first = curr;
                 else
                   first = ++curr;
             }
             else
             {
               first = ++curr;
             }
             Curve_iterator itr;
             for(itr = first; itr != m_pgn.curves_end(); ++itr)
             {
               *widget << *itr;
             }
             *widget << CGAL::WHITE;
             *widget << Segment_2(m_rubber, m_last_of_poly);
             const Arc_point_2& last_point = first->source();

             CGAL_assertion(!last_point.x().is_extended() && !last_point.y().is_extended());
             FT xs = last_point.x().alpha();
             FT ys = last_point.y().alpha();

             *widget << Segment_2(m_rubber, Point_2(xs, ys));
             widget->setRasterOp(old_rasterop);
             widget->unlock();
             m_last_of_poly = Point_2(xs, ys);
             itr = first;
             while(itr != m_pgn.curves_end())
             {
               Curve_iterator temp = itr;
               ++itr;
               m_pgn.erase(temp);
             }

           }
         }
         else
         {
           // the polygon is empty. if m_first_time is false,
           // erase the rubber segment and m_active becomes false
           if(!m_first_time)
           {
            widget->lock();
            RasterOp old_rasterop=widget->rasterOp();
            widget->get_painter().setRasterOp(XorROP);
            *widget << CGAL::WHITE;
            *widget << Segment_2(m_rubber, m_last_of_poly);
            widget->setRasterOp(old_rasterop);
            widget->unlock();
            m_first_time = true;
            m_active = false;
           }

         }
        break;
    }//endswitch
  }//end keyPressEvent


  void mouseMoveEvent(QMouseEvent *e)
  {
    if(!m_active)
      return;
    if(m_ignore_move_event)
    {
      m_ignore_move_event = false;
      return;
    }

    FT x, y;
    widget->x_real(e->x(), x);
    widget->y_real(e->y(), y);
    m_rubber = Point_2(x, y);

    if(m_is_circ_mode)
    {
      Curve_2 circ_arc(m_arc_source, m_rubber, m_arc_target);

      widget->lock();
      RasterOp old_rasterop=widget->rasterOp();
      widget->get_painter().setRasterOp(XorROP);
      *widget << CGAL::WHITE;
      if(!m_first_time)
      {
        *widget << m_old_arc;
      }
      *widget << circ_arc;
      m_old_arc = circ_arc;
      m_rubber_old = m_rubber;
      widget->setRasterOp(old_rasterop);
      m_first_time = false;
      widget->unlock();

      return;
    }



    widget->lock();
    RasterOp old_rasterop=widget->rasterOp();
    widget->get_painter().setRasterOp(XorROP);
    *widget << CGAL::WHITE;
    if(!m_first_time)
      *widget << Segment_2(m_rubber_old, m_last_of_poly);
    *widget << Segment_2(m_rubber, m_last_of_poly);
    m_first_time = false;
    m_rubber_old = m_rubber;
    widget->setRasterOp(old_rasterop);
    widget->unlock();
  };//end mouseMoveEvent


  void CircModeEvent(QMouseEvent *e)
  {
    if(!m_active)
      return;
     FT x, y;
     widget->x_real(e->x(), x);
     widget->y_real(e->y(), y);

     m_arc_target = Point_2(x, y);

     FT xs = m_last_of_poly.x();
     FT ys = m_last_of_poly.y();

     m_arc_source = Point_2(xs, ys);
     double last_pgn_pt_x = CGAL::to_double(xs);
     double last_pgn_pt_y = CGAL::to_double(ys);


     m_old_arc = Curve_2(m_arc_source, m_arc_target);
     widget->lock();
     *widget << CGAL::WHITE;
     RasterOp old_rasterop=widget->rasterOp();
     widget->get_painter().setRasterOp(XorROP);
     *widget << Segment_2(m_rubber, m_last_of_poly);
     *widget << m_old_arc;
     widget->setRasterOp(old_rasterop);
     widget->unlock();


     Double_point_2 ptt(last_pgn_pt_x, last_pgn_pt_y);

     double mid_x = (last_pgn_pt_x + CGAL::to_double(x)) /2;
     double mid_y = (last_pgn_pt_y + CGAL::to_double(y)) /2;

     int mid_point_x_pixel = widget->x_pixel(mid_x);
     int mid_point_y_pixel = widget->y_pixel(mid_y);

     QPoint qp(mid_point_x_pixel, mid_point_y_pixel);
     QPoint qq = widget->mapToGlobal(qp);
     QCursor::setPos(qq);
     m_ignore_move_event = true;

     m_is_circ_mode = true;

  };//end mouseDoubleClickEvent


  void activating()
  {
    m_oldcursor = widget->cursor();
    widget->setCursor(m_cursor);
    m_oldpolicy = widget->focusPolicy();
    widget->setFocusPolicy(QWidget::StrongFocus);
  };

  void deactivating()
  {
    m_pgn.clear();
    m_active = false;
    m_first_time = true;
    widget->setCursor(m_oldcursor);
    widget->setFocusPolicy(m_oldpolicy);
    widget->redraw();
  };

  private:

  bool is_simple(bool is_last_curve = false)
  {
    if(this->m_pgn.size() > 0)
    {
      X_monotone_curve_2  rubber_curve;
      if(!m_is_circ_mode)
      {
        if(is_last_curve)
        {
          const Arc_point_2& first_point = m_pgn.curves_begin()->source();
          CGAL_assertion(!first_point.x().is_extended() && !first_point.y().is_extended());
          FT xs = first_point.x().alpha();
          FT ys = first_point.y().alpha();
          rubber_curve = X_monotone_curve_2(m_last_of_poly, Point_2(xs, ys));
          return(does_curve_disjoint_interior(rubber_curve, is_last_curve));
        }
        else
        {
          rubber_curve = X_monotone_curve_2(m_last_of_poly, m_rubber_old);
          return(does_curve_disjoint_interior(rubber_curve, is_last_curve));
        }
      }
      else
      {
        // circ mode
        Traits_2 tr;
        typename Traits_2::Make_x_monotone_2 make_x = tr.make_x_monotone_2_object();
        std::vector<Object> xcurves_vec;
        xcurves_vec.reserve(3);
        make_x(m_old_arc, std::back_inserter(xcurves_vec));
        for(unsigned int i=0; i<xcurves_vec.size(); ++i)
        {
          X_monotone_curve_2 cv_arc;
          CGAL::assign(cv_arc, xcurves_vec[i]);
          if(!does_curve_disjoint_interior(cv_arc, is_last_curve))
            return false;
        }
        return true;
      }
    }
    return true;
  }

  bool does_curve_disjoint_interior(const X_monotone_curve_2& rubber_curve,
                                    bool is_last_curve)
  {

    Curve_iterator before_last_cv = this->m_pgn.curves_end();
      --before_last_cv;

      Traits_2 tr;
      typename Traits_2::Intersect_2 intersect_func = tr.intersect_2_object();
      Curve_iterator it = this->m_pgn.curves_begin();

      std::list<CGAL::Object> obj_list;

      intersect_func(rubber_curve, *it, std::back_inserter(obj_list));
      if(is_last_curve)
      {
        //the last curve will intersect the first one at their common
        // end point.

        if(m_pgn.size() == 1)
        {
          if(m_pgn.curves_begin()->is_linear())
            return false;
          //the polygon can have one circular arc (and now we close it)
          CGAL_assertion(obj_list.size() == 2);
          return true;
        }
        if(obj_list.size() > 1)
          return false;
        std::pair<Arc_point_2, unsigned int> inter_point;
        CGAL_assertion(obj_list.size() == 1);
        if(! CGAL::assign(inter_point, obj_list.front()))
          return false;
        obj_list.clear();
      }
      else
      {
        if(m_pgn.size() == 1)
        {
          // its the second curve,
          //can intersect the first one at the common end point.

          if(obj_list.empty())
            return true; // no intersections at all
          if(obj_list.size() == 1)
          {
            std::pair<Arc_point_2, unsigned int> inter_pt;
            bool succ = CGAL::assign(inter_pt, obj_list.front());
            if(!succ)
              return false; // overlap curves!!

            Traits_2 tr;
            // make sure that the intersection point is equal to the curve target
            return (tr.equal_2_object()(inter_pt.first, it->target()));
          }
          return false;
        }
        // its not the last curve (or the secind), cannot intersect the first curve.
        if(!obj_list.empty())
        {
         return false;
        }
      }
      ++it;
      for(; it != before_last_cv; ++it)
      {
        intersect_func(rubber_curve, *it, std::back_inserter(obj_list));
        if(!obj_list.empty())
          return false;
      }
      //if I'm out of this means that all the edges,
      //didn't intersect the last one
      intersect_func(rubber_curve, *it, std::back_inserter(obj_list));
      if(obj_list.empty())
        return true;
      if(obj_list.size() > 1)
        return false;

      std::pair<Arc_point_2, unsigned int> inter_point;
      if(CGAL::assign(inter_point, obj_list.front()))
      {
        Traits_2 tr;
        return (tr.equal_2_object()(rubber_curve.source(), it->target()));
      }
      return false;
  }

  };

} // namespace CGAL
#endif
