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
//             Ronnie Gandhi <ronniegandhi19999@gmail.com>
//             Efi Fogel <efifogel@gmain.com>

#ifndef CGAL_QT_GRAPHICS_VIEW_LINEAR_POLYGON_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_LINEAR_POLYGON_INPUT_H

#include <limits>

#include <CGAL/auto_link/Qt.h>
#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>

#include <QtCore>
#include <QPolygonF>
#include <QPointF>
#include <QGraphicsLineItem>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <list>
#include "Typedefs.h"
#include "QT5/Linear_polygons.h"

namespace CGAL {
namespace Qt {

template <typename Kernel_>
class Graphics_view_linear_polygon_input : public GraphicsViewInput {
public:

  typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel1;
  typedef Kernel1::Point_2                            Point_2;
  typedef CGAL::Polygon_2<Kernel1>                    Polygon_2;
  typedef CGAL::Polygon_with_holes_2<Kernel1>         Polygon_with_holes_2;
  typedef std::list<Polygon_with_holes_2>             Pgn_with_holes_2_container;


  typedef Kernel_                                 Kernel;

  typedef Linear_traits                           Gps_traits;
  typedef typename Gps_traits::Curve_2            Linear_curve;
  typedef typename Gps_traits::X_monotone_curve_2 Linear_X_monotone_curve;
  typedef typename Gps_traits::Polygon_2          Linear_polygon;
  typedef typename Gps_traits::Point_2            Linear_point;
  typedef typename Kernel::FT                     FT;
  typedef typename Kernel::Vector_2               Vector;
  typedef typename Kernel::Point_2                Point;
  typedef std::vector<Linear_curve>               Linear_curve_vector;
  typedef typename Linear_curve_vector::const_iterator
    const_linear_curve_iterator;

  //where did it come from
  typedef Linear_boundary_pieces_graphics_item<Linear_curve_vector> GI;

  //constructor
  Graphics_view_linear_polygon_input(QObject* aParent, QGraphicsScene* aScene) :
    GraphicsViewInput(aParent),
    mScene(aScene),
    mLinearGI(new GI(&mLinearPolygonPieces)),
    mOngoingPieceGI(new GI(&mOngoingPieceCtr)),
    mHandleGI(new QGraphicsLineItem()),
    mLinearPolygonPen(QColor(0, 255, 0)),
    mOngoingCurvePen(QColor(255, 215, 0)),
    mHandlePen(QColor(255, 165, 0)),
    mState(Start),
    m_bound_rect(true),
    m_hole(false),
    m_is_mink(false)
  {
    mOngoingPieceGI->setPen(mOngoingCurvePen);
    mHandleGI->setPen(mHandlePen);

    mHandleGI->setLine(0,0,1,1);
    mHandleGI->hide();

    mLinearGI->setPen(mLinearPolygonPen);

    mScene->addItem(mOngoingPieceGI);
    mScene->addItem(mHandleGI);
    mScene->addItem(mLinearGI);
  }

  //destructor
  ~Graphics_view_linear_polygon_input() {}

  //decision making:wether to accept the event or reject it
  bool eventFilter(QObject* obj, QEvent* aEvent)
  {
    bool rHandled = false;

    if (aEvent->type() == QEvent::GraphicsSceneMousePress) {
      rHandled = mousePressEvent(static_cast<QGraphicsSceneMouseEvent*>(aEvent));
    }
    else if (aEvent->type() == QEvent::GraphicsSceneMouseRelease) {
      rHandled =
        mouseReleaseEvent( static_cast<QGraphicsSceneMouseEvent*>(aEvent));
    }
    else if (aEvent->type() == QEvent::GraphicsSceneMouseMove) {
      rHandled = mouseMoveEvent( static_cast<QGraphicsSceneMouseEvent*>(aEvent));
    }
    else if (aEvent->type() == QEvent::KeyPress) {
      rHandled = keyPressEvent( static_cast<QKeyEvent*>(aEvent));
    }

    if (! rHandled) rHandled = QObject::eventFilter(obj, aEvent);

    return rHandled;
  }

public:
  enum State {
    Start, PieceStarted, PieceOngoing, HandleOngoing, PieceEnded, CurveEnded
  };

  Point cvt(QPointF const& aP) const { return Point(aP.x(), aP.y()); }

  bool mousePressEvent(QGraphicsSceneMouseEvent *aEvent)
  {
    m_bound_rect = false;
    bool rHandled = false;
    Point lP = cvt(aEvent->scenePos());
    if (aEvent->button() == ::Qt::LeftButton) {
      switch (mState) {
       case Start:
        mP0 = lP;
        mState = PieceStarted;
        rHandled = true;
        break;

       case PieceOngoing:
        mP1 = lP;
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
      UpdateHandle(lP);
      UpdateOngoingPiece();
      rHandled = true;
      break;

     case PieceEnded:
      mState   = PieceOngoing;
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
    QPointF const& aP = aEvent ->scenePos();
    if (aEvent->button() == ::Qt::LeftButton) {
      switch (mState) {
       case PieceStarted:
        mState = PieceOngoing;
        rHandled = true;
        break;

       case HandleOngoing:
        HideHandle();
        CommitOngoingPiece(lP);
        mState   = PieceEnded;
        rHandled = true;
        break;

       default: break; //! \todo handle default case
      }
    }
    else if (aEvent->button() == ::Qt::RightButton) {
      switch (mState) {
       case PieceOngoing:
        //cout<<"hello in Graphics_view_linear_polygon"<<endl;
        m_bound_rect = false;
        CommitCurrLinearPolygon();
        ReStart();
        rHandled = true;
        //cout<<"right click over"<<endl;
        break;

       default: break; //! \todo handle default case
      }
    }
    return rHandled;
  }

  bool keyPressEvent(QKeyEvent *aEvent)
  {
    bool rHandled = false;

    if ((aEvent->key() == ::Qt::Key_Delete) ||
        (aEvent->key() == ::Qt::Key_Backspace))
    {
      RemoveLastPiece();
      mState = (mLinearPolygonPieces.size() > 0) ? PieceEnded : Start;
      rHandled = true;
    }
    else if (aEvent->key() == ::Qt::Key_Escape) {
      Reset();
      mState = Start;
      rHandled = true;
    }
    return rHandled;
  }

  Linear_curve const* ongoing_piece() const
  { return (mOngoingPieceCtr.size() == 1) ? &mOngoingPieceCtr[0] : NULL; }

  void ReStart()
  {
    mH = boost::optional<Point>();
    mState = Start;
  }

  void Reset()
  {
    mLinearPolygonPieces.clear();
    mOngoingPieceCtr.clear();
    mLinearGI->modelChanged();
    mOngoingPieceGI->modelChanged();
    ReStart();
  }

  void HideHandle() { mHandleGI->hide(); }

  Linear_curve CreatePiece() { return Linear_curve(mP0, mP1); }

  void RemoveLastPiece()
  {
    mLinearPolygonPieces.pop_back();
    mOngoingPieceCtr.clear();
    mLinearGI->modelChanged();
    mOngoingPieceGI->modelChanged();
    if (mLinearPolygonPieces.size() > 0) {
      mP0 = cvt(mLinearPolygonPieces.back().target());
      UpdateOngoingPiece();
    }
    mH = boost::optional<Point>();
  }

  void UpdateOngoingPiece()
  {
    if (mOngoingPieceCtr.size() > 0) mOngoingPieceCtr.clear();
    mOngoingPieceCtr.push_back(CreatePiece());
    mOngoingPieceGI->modelChanged();
  }

  void CommitOngoingPiece(Point const& aP)
  {
    if (ongoing_piece()) {
      mLinearPolygonPieces.push_back(*ongoing_piece());
      mLinearGI->modelChanged();
      mOngoingPieceCtr.clear();
      mOngoingPieceGI->modelChanged();
      mP0 = mP1;
      mP1 = aP;
      mH = boost::optional<Point>();
    }
  }

  void UpdateHandle(Point const& aP)
  {
    if (squared_distance(mP1,aP) >= 4) {
      mH = aP;
      mHandleGI->setLine(to_double(mP1.x()), to_double(mP1.y()),
                         to_double(mH->x()), to_double(mH->y()));
      mHandleGI->show();
    }
    else {
      HideHandle();
      mH = boost::optional<Point>();
    }
  }

  //Point cvt ( typename Linear_curve::Point_2 const& aP)
  //{ return Point( to_double(aP.x()), to_double(aP.y())); }
  Point cvt (Point const& aP)
  { return Point(to_double(aP.x()), to_double(aP.y())); }

  void CommitCurrLinearPolygon()
  {
    GenerateLinearPolygon();

    mOngoingPieceCtr.clear();
    mOngoingPieceGI->modelChanged();

    mLinearPolygonPieces.clear();
    mLinearGI->modelChanged();

    mH = boost::optional<Point>();

    HideHandle();
  }

  void GenerateLinearPolygon()
  {
    if (mLinearPolygonPieces.size() >  1) {
      Gps_traits traits;
      typename Gps_traits::Make_x_monotone_2 make_x_monotone =
        traits.make_x_monotone_2_object();

      std::vector<Linear_X_monotone_curve> xcvs;
      for (auto it = mLinearPolygonPieces.begin();
           it != mLinearPolygonPieces.end(); ++it)
      {
        std::vector<CGAL::Object> x_objs;
        std::vector<CGAL::Object>::const_iterator xoit;

        make_x_monotone(*it, std::back_inserter(x_objs));
        //cout<<"add curves"<<endl;

        for (xoit = x_objs.begin(); xoit != x_objs.end(); ++xoit) {
          Linear_X_monotone_curve xcv;
          if (CGAL::assign(xcv, *xoit)) xcvs.push_back (xcv);
        }
      }

      if (xcvs.size() > 0) {
        Linear_point const& first_point = xcvs.front().source();
        Linear_point const& last_point =  xcvs.back().target();
        FT fxs = first_point.x();
        FT fys = first_point.y();
        FT lxs = last_point .x();
        FT lys = last_point .y();
        xcvs.push_back(Linear_X_monotone_curve( Point(lxs,lys), Point(fxs,fys)));
        Linear_polygon lp(xcvs.begin(), xcvs.end());
        emit(generate(CGAL::make_object(lp)));
      }
    }
  }


  void get_BoundingRect()
  {
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

      CommitCurrLinearPolygon();
      ReStart();

  }

  void get_Minkowski_holes(Polygon_with_holes_2 polygon, Polygon_with_holes_2 p0)
  {
    typename Polygon_2::Vertex_const_iterator vit;
    m_bound_rect = false;
    Point_2 pt;
    Polygon_2 hole;

    int i;
    
    for(Polygon_with_holes_2::Hole_const_iterator hi = polygon.holes_begin(); hi != polygon.holes_end(); ++ hi )
    {
      hole =*hi;

      i=0;

      for (vit = hole.begin(); vit != hole.end(); ++vit)
      {
        if(i==0)
        {
          pt =  *vit;

          mP0 = Point(CGAL::to_double(pt.x()),CGAL::to_double(pt.y()));
          mState = PieceStarted;
        }
        else
        {
          pt =  *vit;

          mState = PieceOngoing;

          mP1 = Point(CGAL::to_double(pt.x()),CGAL::to_double(pt.y()));
          UpdateOngoingPiece();

          mP1 = Point(CGAL::to_double(pt.x()),CGAL::to_double(pt.y()));
          mState = HandleOngoing;
          HideHandle();
          CommitOngoingPiece(Point(CGAL::to_double(pt.x()),CGAL::to_double(pt.y())));
          mState   = PieceEnded;
        }
        i++;
      }

      mState   = PieceOngoing;

      mP1 = Point(-90,-90);
      UpdateOngoingPiece();

      CommitCurrLinearPolygon();
      ReStart();
    }
  }

  void get_Minkowski_result(Polygon_with_holes_2 polygon, Polygon_with_holes_2 p0)
  {
    typename Polygon_2::Vertex_const_iterator vit;
    m_bound_rect = false;
    Point_2 pt;

    int i=0;
    
    for (vit = polygon.outer_boundary().vertices_begin(); vit != polygon.outer_boundary().vertices_end(); ++vit)
    {
      if(i==0)
      {
        pt =  *vit;

        mP0 = Point(CGAL::to_double(pt.x()),CGAL::to_double(pt.y()));
        mState = PieceStarted;
      }
      else
      {
        pt =  *vit;

        mState = PieceOngoing;

        mP1 = Point(CGAL::to_double(pt.x()),CGAL::to_double(pt.y()));
        UpdateOngoingPiece();

        mP1 = Point(CGAL::to_double(pt.x()),CGAL::to_double(pt.y()));
        mState = HandleOngoing;
        HideHandle();
        CommitOngoingPiece(Point(CGAL::to_double(pt.x()),CGAL::to_double(pt.y())));
        mState   = PieceEnded;

      }

      i++;
    }

    mState   = PieceOngoing;

    mP1 = Point(-90,-90);
    UpdateOngoingPiece();

    CommitCurrLinearPolygon();
    ReStart();
  }

  bool isboundingRect()
  {
    return m_bound_rect;
  }

  bool ishole()
  {
    return m_hole;
  }

  bool is_mink()
  {
    return m_is_mink;
  }


public:
  QGraphicsScene* mScene;
  GI* mLinearGI;
  GI* mOngoingPieceGI;
  QGraphicsLineItem* mHandleGI;

  QPen mLinearPolygonPen;
  QPen mOngoingCurvePen;
  QPen mHandlePen;

  Linear_curve_vector mLinearPolygonPieces;
  Linear_curve_vector mLinearPolygonPieces_boundR;
  Linear_curve_vector mOngoingPieceCtr;

  int mState;

  bool m_bound_rect;

  bool m_is_mink;

  bool m_hole;

  Point mP0;
  Point mP1;
  //boost::optional returns optional return type
  //link to documentation https://theboostcpplibraries.com/boost.optional
  boost::optional<Point> mH;
};

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_VIEW_LINEAR_POLYGON_INPUT_H
