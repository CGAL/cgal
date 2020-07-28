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
//             Ronnie Gandhi <ronniegandhi19999@gmil.com>
//             Efi Fogel <efifogel@gmain.com>

#ifndef CGAL_QT_GRAPHICS_VIEW_CIRCULAR_POLYGON_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_CIRCULAR_POLYGON_INPUT_H

//#include <iostream>

#include <CGAL/auto_link/Qt.h>
#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Arr_circle_segment_traits_2.h>

#include <QPolygonF>
#include <QPointF>
#include <QGraphicsLineItem>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>

#include "QT5/Circular_polygons.h"
//#include "Typedefs.h"

namespace CGAL {
namespace Qt {

template <typename Kernel_>
class Graphics_view_circular_polygon_input : public GraphicsViewInput {
public:
  typedef Kernel_                                   Kernel;

  typedef CGAL::Gps_circle_segment_traits_2<Kernel> Gps_traits;
  typedef typename Gps_traits::Curve_2              Circular_curve;
  typedef typename Gps_traits::X_monotone_curve_2   Circular_X_monotone_curve;
  typedef typename Gps_traits::Polygon_2            Circular_polygon;
  typedef typename Circular_polygon::Point_2        Arc_point;
  typedef typename Kernel::FT                       FT;
  typedef typename Kernel::Vector_2                 Vector;
  typedef typename Kernel::Point_2                  Point;

  typedef std::vector<Circular_curve>               Circular_curve_vector;

  typedef typename Circular_curve_vector::const_iterator
    const_circular_curve_iterator;

  typedef Circular_boundary_pieces_graphics_item<Circular_curve_vector> GI;

  //constructor
  Graphics_view_circular_polygon_input(QObject* aParent,
                                       QGraphicsScene* aScene) :
    GraphicsViewInput(aParent),
    mScene(aScene),
    mOngoingPieceGI(new GI(&mOngoingPieceCtr)),
    mHandleGI(new QGraphicsLineItem()),
    mCircularPolygonPen(QColor(0, 255, 0)),
    mOngoingCurvePen(QColor(255, 215, 0)),
    mHandlePen(QColor(255, 165, 0)),
    mState(Start),
    m_bound_rect(true),
    m_last_circular(false),
    m_last(false)
  {
    mOngoingPieceGI->setPen(mOngoingCurvePen);
    mHandleGI->setPen(mHandlePen);

    mHandleGI->setLine(0,0,1,1);
    mHandleGI->hide();

    mCircularGI = new GI(&mCircularPolygonPieces);

    mCircularGI->setPen(mCircularPolygonPen);

    mScene->addItem(mOngoingPieceGI);
    mScene->addItem(mHandleGI);
    mScene->addItem(mCircularGI);
  }

  //destructor
  ~Graphics_view_circular_polygon_input() {}

  //manages function call to all mouse activities
  bool eventFilter(QObject* obj, QEvent* aEvent)
  {
    bool rHandled = false;

    if (aEvent->type() == QEvent::GraphicsSceneMousePress) {
      rHandled = mousePressEvent(static_cast<QGraphicsSceneMouseEvent*>(aEvent));
    }
    else if (aEvent->type() == QEvent::GraphicsSceneMouseRelease) {
      rHandled =
        mouseReleaseEvent(static_cast<QGraphicsSceneMouseEvent*>(aEvent));
    }
    else if (aEvent->type() == QEvent::GraphicsSceneMouseMove) {
      rHandled = mouseMoveEvent(static_cast<QGraphicsSceneMouseEvent*>(aEvent));
    }
    else if (aEvent->type() == QEvent::KeyPress) {
      rHandled = keyPressEvent(static_cast<QKeyEvent*>(aEvent));
    }

    if (!rHandled) rHandled = QObject::eventFilter(obj, aEvent);

    return rHandled;
  }

public:
  //a set of all states of drawing a circular polygon
  enum State {
    Start, PieceStarted, PieceOngoing, HandleOngoing, PieceEnded, CurveEnded
  };

  Point cvt(QPointF const& aP) const { return Point(aP.x(), aP.y()); }

  //All functions related to mouse activity
  bool mousePressEvent(QGraphicsSceneMouseEvent* aEvent)
  {
    bool rHandled = false;
    m_bound_rect = false;

    //Point lP = cvt(aEvent->QGraphicsSceneMouseEvent::scenePos());
    Point lP = cvt(aEvent->scenePos());
    if (aEvent->button() == ::Qt::LeftButton) {
      switch (mState) {
       case Start:
        mP0 = lP;
        mState = PieceStarted;
        rHandled = true;
        break;

       case PieceOngoing:
        mP1    = lP;
        mState = HandleOngoing;
        rHandled = true;
        break;

       default: break; //! \todo handle default case
      }
    }

    else  if (aEvent->button() == ::Qt::RightButton) {
      switch (mState) {
        case PieceOngoing:
          // allowing user to curve last piece as well
          m_last = true;
          mState = HandleOngoing;
          rHandled = true;
          break;

         default: break; //! \todo handle default case
       }
    }
    return rHandled;
  }

  bool mouseMoveEvent(QGraphicsSceneMouseEvent* aEvent)
  {
    bool rHandled = false;

    Point lP = cvt(aEvent->scenePos());

    switch (mState) {
     case PieceOngoing:
      mP1 = lP;
      UpdateOngoingPiece();
      rHandled = true;
      break;

     case HandleOngoing:
      if(m_last)
      {
        mP1 = cvt(mCircularPolygonPieces.front().source());
        m_last_circular = true;
      }
      UpdateHandle(lP);
      UpdateOngoingPiece();
      rHandled = true;
      break;

     case PieceEnded:
      mState = PieceOngoing;
      rHandled = true;
      break;

     default: break; //! \todo handle default case
    }

    return rHandled;
  }

  bool mouseReleaseEvent(QGraphicsSceneMouseEvent* aEvent)
  {
    bool rHandled = false;
    Point lP = cvt(aEvent->scenePos());
    if (aEvent->button() == ::Qt::LeftButton) {
      switch (mState) {
       case PieceStarted:
        mState = PieceOngoing;
        rHandled = true;
        break;

       case HandleOngoing:
        HideHandle();
        CommitOngoingPiece(lP);
        mState = PieceEnded;
        rHandled = true;
        break;

       default: break; //! \todo handle default case
      }
    }
    else if (aEvent->button() == ::Qt::RightButton) {
      switch (mState) {
       case HandleOngoing:
        //cout<<"hello in Graphics_view_circular_polygon"<<endl;
        if(m_last_circular)
        { 
          HideHandle();
          CommitOngoingPiece(lP);
        }
        m_bound_rect = false;
        CommitCurrCircularPolygon();
        ReStart();
        rHandled = true;
        //cout<<"right click over"<<endl;
        break;

       default: break; //! \todo handle default case
      }
    }
    return rHandled;
  }

  bool keyPressEvent(QKeyEvent* aEvent)
  {
    bool rHandled = false;

    if ((aEvent->key() == ::Qt::Key_Delete) ||
        (aEvent->key() == ::Qt::Key_Backspace))
    {
      RemoveLastPiece();
      mState = (mCircularPolygonPieces.size() > 0) ? PieceEnded : Start;
      rHandled = true;
    }
    else if (aEvent->key() == ::Qt::Key_Escape) {
      Reset();
      mState = Start;
      rHandled = true;
    }
    return rHandled;
  }

public:
  
  
  Circular_curve const* ongoing_piece() const
  { return (mOngoingPieceCtr.size() == 1) ? &mOngoingPieceCtr[0] : NULL; }

  void ReStart()
  {
    mH = boost::optional<Point>();
    mState = Start;
  }

  void Reset()
  {
    mCircularPolygonPieces.clear();
    mOngoingPieceCtr.clear();
    mCircularGI->modelChanged();
    mOngoingPieceGI->modelChanged();
    ReStart();
  }

  void HideHandle()
  {
    mHandleGI->hide();
  }

  Circular_curve CreatePiece()
  {
    if (mH) {
      Vector lD = *mH - mP1;
      Vector lU = lD * 1.5;
      Point  lH = mP1 - lU;
      return Circular_curve(mP0,lH,mP1);
    }
    else return Circular_curve(mP0,mP1);
  }

  void RemoveLastPiece()
  {
    mCircularPolygonPieces.pop_back();
    mOngoingPieceCtr.clear();
    mCircularGI->modelChanged();
    mOngoingPieceGI->modelChanged();
    if (mCircularPolygonPieces.size() > 0) {
      mP0 = cvt(mCircularPolygonPieces.back().target());
      UpdateOngoingPiece();
    }
    mH = boost::optional<Point>();
  }

  void UpdateOngoingPiece()
  {
    if (mOngoingPieceCtr.size() > 0) mOngoingPieceCtr.clear();
    mOngoingPieceCtr.push_back(CreatePiece());
    //cout<<"hi"<<endl;
    mOngoingPieceGI->modelChanged();
  }

  void CommitOngoingPiece(Point const& aP)
  {
    if (ongoing_piece()) {
      mCircularPolygonPieces.push_back(*ongoing_piece());
      mCircularGI->modelChanged();
      mOngoingPieceCtr.clear();
      mOngoingPieceGI->modelChanged();
      mP0 = mP1;
      mP1 = aP;
      mH = boost::optional<Point>();
    }
  }

  void UpdateHandle(Point const& aP)
  {
    if (squared_distance(mP1,aP) >= 4)
    {
      mH = aP;
      mHandleGI->setLine(to_double(mP1.x()), to_double(mP1.y()),
                         to_double(mH->x()), to_double(mH->y()));
      mHandleGI->show();
    }
    else
    {
      HideHandle();
      mH = boost::optional<Point>();
    }
  }

  Point cvt(typename Circular_curve::Point_2 const& aP)
  { return Point(to_double(aP.x()), to_double(aP.y())); }

  void CommitCurrCircularPolygon()
  {
    GenerateCircularPolygon();
    mOngoingPieceCtr.clear();
    //cout<<"mOngoingPieceCtr"<<endl;
    mOngoingPieceGI->modelChanged();

    mCircularPolygonPieces.clear();
    //cout<<"mCircularPolygonPieces"<<endl;
    mCircularGI->modelChanged();

    mH = boost::optional<Point>();

    HideHandle();
    //cout<<"polygon is comitted"<<endl;
  }

  void GenerateCircularPolygon()
  {
    if (mCircularPolygonPieces.size() > 0) {
      Gps_traits traits;
      auto make_x_monotone = traits.make_x_monotone_2_object();

      std::vector<Circular_X_monotone_curve> xcvs;
      for (auto it = mCircularPolygonPieces.begin();
           it != mCircularPolygonPieces.end(); ++ it)
      {
        std::vector<CGAL::Object> x_objs;
        std::vector<CGAL::Object>::const_iterator xoit;
        //cout<<"point 1"<<endl;
        make_x_monotone(*it, std::back_inserter(x_objs));
        //cout<<"add curves"<<endl;
        //cout<<"point 2"<<endl;
        //exception handling: if user draws a line and ends polygon
        Circular_X_monotone_curve xcv;
        xoit = x_objs.begin();
        CGAL::assign(xcv,*xoit);
        if (xcv.is_linear() && mCircularPolygonPieces.size() == 1) return;
        for (xoit = x_objs.begin(); xoit != x_objs.end(); ++xoit) {
          if (CGAL::assign(xcv, *xoit)) xcvs.push_back(xcv);
        }
        //cout<<"point 3"<<endl;
      }

      if (xcvs.size() > 0) {
        //cout<<"point 4"<<endl;

        if(!m_last_circular)
        {
          Arc_point const& first_point = xcvs.front().source();
          Arc_point const& last_point =  xcvs.back().target();

          CGAL_assertion(!first_point.x().is_extended() &&
                         !first_point.y().is_extended());
          CGAL_assertion(!last_point. x().is_extended() &&
                         !last_point .y().is_extended());
          FT fxs = first_point.x().alpha();
          FT fys = first_point.y().alpha();
          FT lxs = last_point .x().alpha();
          FT lys = last_point .y().alpha();
          xcvs.push_back(Circular_X_monotone_curve(Point(lxs,lys),
                                                   Point(fxs,fys)));
        }

        m_last = false;
        m_last_circular = false;

        //cout<<"add curves if circular"<<endl;
        Circular_polygon cp(xcvs.begin(), xcvs.end());
        //cout<<"point 5"<<endl;
        emit(generate(CGAL::make_object(cp)));
        //cout<<"point 6"<<endl;
      }
    }
    //cout<<"generated circular polygon"<<endl;
  }

  void get_BoundingRect()
  {
      // mOngoingPieceCtr.push_back(Linear_curve(Point(-10000000,-10000000),Point(-10000000,10000000)));
      // mOngoingPieceCtr.push_back(Linear_curve(Point(-10000000,10000000),Point(10000000,10000000)));
      // mOngoingPieceCtr.push_back(Linear_curve(Point(10000000,10000000),Point(10000000,-10000000)));

      // mLinearPolygonPieces.push_back(mOngoingPieceCtr[0]);
      // mLinearPolygonPieces.push_back(mOngoingPieceCtr[1]);
      // mLinearPolygonPieces.push_back(mOngoingPieceCtr[2]);

      // m_bound_rect = true;
      // CommitCurrLinearPolygon();
      // ReStart();

      m_bound_rect = true;

      mP0 = Point(-15500000,-10000000);
      mState = PieceStarted;

      mState = PieceOngoing;
      mP1 = Point(-15500000,10000000);
      UpdateOngoingPiece();

      mP1 = Point(-15500000,10000000);
      mState = HandleOngoing;
      HideHandle();
      CommitOngoingPiece(Point(-15500000,10000000));
      mState   = PieceEnded;

      mState   = PieceOngoing;
      mP1 = Point(15500000,10000000);
      UpdateOngoingPiece();

      mP1 = Point(15500000,10000000);
      mState = HandleOngoing;
      HideHandle();
      CommitOngoingPiece(Point(15500000,10000000));
      mState   = PieceEnded;

      mState   = PieceOngoing;
      mP1 = Point(15500000,-10000000);
      UpdateOngoingPiece();

      mP1 = Point(15500000,-10000000);
      mState = HandleOngoing;
      HideHandle();
      CommitOngoingPiece(Point(15500000,-10000000));
      mState   = PieceEnded;

      mState   = PieceOngoing;
      mP1 = Point(-9000000,-9000000);
      UpdateOngoingPiece();

      CommitCurrCircularPolygon();
      ReStart();
  }

  bool isboundingRect()
  {
    return m_bound_rect;
  }

public:
  QGraphicsScene* mScene;
  GI* mCircularGI;
  GI* mOngoingPieceGI;
  QGraphicsLineItem* mHandleGI;

  QPen mCircularPolygonPen;
  QPen mOngoingCurvePen;
  QPen mHandlePen;

  bool m_bound_rect;
  bool m_last_circular;
  bool m_last;

  Circular_curve_vector mCircularPolygonPieces;
  Circular_curve_vector mOngoingPieceCtr;

  int mState;

  Point mP0;
  Point mP1;

  boost::optional<Point> mH;
};

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_VIEW_CIRCULAR_POLYGON_INPUT_H
