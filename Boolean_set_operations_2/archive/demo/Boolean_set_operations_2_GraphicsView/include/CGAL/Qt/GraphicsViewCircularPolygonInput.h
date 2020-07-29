// Copyright (c) 2010  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <Fernando.Cacciola@geometryfactory.com>
//

#ifndef CGAL_QT_GRAPHICS_VIEW_CIRCULAR_POLYGON_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_CIRCULAR_POLYGON_INPUT_H

#include <CGAL/auto_link/Qt.h>

#include <QPolygonF>
#include <QPointF>
#include <QGraphicsLineItem>
#include <QGraphicsScene>


#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/CircularPolygons.h>

namespace CGAL {

namespace Qt {

  template <class K>
  class GraphicsViewCircularPolygonInput : public GraphicsViewInput
  {
  public:

    typedef K Kernel ;

    typedef CGAL::Gps_circle_segment_traits_2<K> Gps_traits;

    typedef typename Gps_traits::Curve_2            Circular_curve;
    typedef typename Gps_traits::X_monotone_curve_2 Circular_X_monotone_curve;
    typedef typename Gps_traits::Polygon_2          Circular_polygon;
    typedef typename Circular_polygon::Point_2      Arc_point ;
    typedef typename Kernel::FT                     FT ;
    typedef typename Kernel::Vector_2               Vector ;
    typedef typename Kernel::Point_2                Point ;

    typedef std::vector<Circular_curve> Circular_curve_vector ;

    typedef typename Circular_curve_vector::const_iterator const_circular_curve_iterator ;

    typedef Circular_boundary_pieces_graphics_item<Circular_curve_vector> GI ;

    GraphicsViewCircularPolygonInput(QObject* aParent, QGraphicsScene* aScene)
      :
        GraphicsViewInput  ( aParent         )
      , mScene             ( aScene          )
      , mState             ( Start           )
      , mCircularPolygonPen( QColor(0,255,0) )
      , mOngoingCurvePen   ( QColor(255,0,0) )
      , mHandlePen         ( QColor(0,0,255) )
      , mCircularGI        ( 0               )
    {
      mOngoingPieceGI = new GI(&mOngoingPieceCtr) ;
      mHandleGI       = new QGraphicsLineItem();

      mOngoingPieceGI->setPen(mOngoingCurvePen);
      mHandleGI      ->setPen(mHandlePen);

      mHandleGI->setLine(0,0,1,1);
      mHandleGI->hide();

      mCircularGI = new GI(&mCircularPolygonPieces) ;

      mCircularGI->setPen(mCircularPolygonPen);

      mScene->addItem(mOngoingPieceGI);
      mScene->addItem(mHandleGI);
      mScene->addItem(mCircularGI);
    }

    ~GraphicsViewCircularPolygonInput()
    {
    }

    bool eventFilter(QObject *obj, QEvent *aEvent)
    {
      bool rHandled = false ;

      if (aEvent->type() == QEvent::GraphicsSceneMousePress)
      {
        rHandled = mousePressEvent( static_cast<QGraphicsSceneMouseEvent *>(aEvent) ) ;
      }
      else if (aEvent->type() == QEvent::GraphicsSceneMouseRelease)
      {
        rHandled = mouseReleaseEvent( static_cast<QGraphicsSceneMouseEvent *>(aEvent) ) ;
      }
      else if (aEvent->type() == QEvent::GraphicsSceneMouseMove)
      {
        rHandled = mouseMoveEvent( static_cast<QGraphicsSceneMouseEvent *>(aEvent) ) ;
      }
      else if (aEvent->type() == QEvent::KeyPress)
      {
        rHandled = keyPressEvent( static_cast<QKeyEvent *>(aEvent) ) ;
      }

      if ( !rHandled )
        rHandled = QObject::eventFilter(obj, aEvent);

      return rHandled ;
    }

  protected:

    enum State { Start, PieceStarted, PieceOngoing, HandleOngoing, PieceEnded, CurveEnded } ;

    Point cvt ( QPointF const& aP ) const { return Point(aP.x(),aP.y()) ; }

    bool mousePressEvent(QGraphicsSceneMouseEvent *aEvent)
    {
      bool rHandled = false ;

      Point lP = cvt(aEvent->scenePos());

      if ( aEvent->button() == ::Qt::LeftButton )
      {
        switch (mState)
        {
          case Start:
            mP0      = lP;
            mState   = PieceStarted;
            rHandled = true;
            break;

          case PieceOngoing:
            mP1      = lP;
            mState   = HandleOngoing;
            rHandled = true;
            break;
        }
      }

      return rHandled ;
    }


    bool mouseMoveEvent(QGraphicsSceneMouseEvent *aEvent)
    {
      bool rHandled = false ;

      Point lP = cvt(aEvent->scenePos());

      switch (mState)
      {
        case PieceOngoing:
          mP1 = lP;
          UpdateOngoingPiece();
          rHandled = true ;
          break;


        case HandleOngoing:

          UpdateHandle(lP);
          UpdateOngoingPiece();

          rHandled = true ;
          break;

        case PieceEnded:
          mState   = PieceOngoing;
          rHandled = true;
          break;
      }

      return rHandled ;
    }

    bool mouseReleaseEvent(QGraphicsSceneMouseEvent *aEvent)
    {
      bool rHandled = false ;

      Point lP = cvt(aEvent->scenePos());

      if ( aEvent->button() == ::Qt::LeftButton )
      {
        switch (mState)
        {
          case PieceStarted:
            mState   = PieceOngoing;
            rHandled = true;
            break;

          case HandleOngoing:
            HideHandle();
            CommitOngoingPiece(lP);
            mState   = PieceEnded;
            rHandled = true;
            break;
        }
      }
      else if ( aEvent->button() == ::Qt::RightButton )
      {
        switch (mState)
        {
          case PieceOngoing:
            CommitCurrCircularPolygon();
            ReStart();
            rHandled = true;
            break;
        }
      }

      return rHandled ;
    }

    bool keyPressEvent(QKeyEvent *aEvent)
    {
      bool rHandled = false ;

      if( aEvent->key() == ::Qt::Key_Delete || aEvent->key() == ::Qt::Key_Backspace )
      {
        RemoveLastPiece();
        mState   = mCircularPolygonPieces.size() > 0 ? PieceEnded : Start ;
        rHandled = true;
      }
      else if( aEvent->key() == ::Qt::Key_Escape)
      {
        Reset();
        mState   = Start;
        rHandled = true;
      }

      return rHandled ;
    }



  private:

    Circular_curve const* ongoing_piece() const { return mOngoingPieceCtr.size() == 1 ? &mOngoingPieceCtr[0] : NULL ; }

    void ReStart()
    {
      mH = boost::optional<Point>();
      mState = Start ;
    }

    void Reset()
    {
      mCircularPolygonPieces.clear();
      mOngoingPieceCtr      .clear();
      mCircularGI    ->modelChanged();
      mOngoingPieceGI->modelChanged();
      ReStart();
    }

    void HideHandle()
    {
      mHandleGI->hide();
    }

    Circular_curve CreatePiece()
    {
      if ( mH )
      {
        Vector lD = *mH - mP1 ;
        Vector lU = lD * 1.5 ;
        Point  lH = mP1 - lU ;
        return Circular_curve(mP0,lH,mP1);
      }
      else
      {
        return Circular_curve(mP0,mP1);
      }
    }


    void RemoveLastPiece()
    {
      mCircularPolygonPieces.pop_back();
      mOngoingPieceCtr      .clear();
      mCircularGI    ->modelChanged();
      mOngoingPieceGI->modelChanged();
      if ( mCircularPolygonPieces.size() > 0 )
      {
        mP0 = cvt(mCircularPolygonPieces.back().target());
        UpdateOngoingPiece();
      }
      mH = boost::optional<Point>();
    }

    void UpdateOngoingPiece()
    {
      if ( mOngoingPieceCtr.size() > 0 )
        mOngoingPieceCtr.clear();
      mOngoingPieceCtr.push_back(CreatePiece());
      mOngoingPieceGI->modelChanged();
    }

    void CommitOngoingPiece( Point const& aP )
    {
      if ( ongoing_piece() )
      {
        mCircularPolygonPieces.push_back( *ongoing_piece() ) ;
        mCircularGI->modelChanged();
        mOngoingPieceCtr.clear();
        mOngoingPieceGI->modelChanged();
        mP0 = mP1 ;
        mP1 = aP ;
        mH = boost::optional<Point>();
      }
    }

    void UpdateHandle(Point const& aP)
    {
      if ( squared_distance(mP1,aP) >= 4 )
      {
        mH = aP ;

        mHandleGI->setLine( to_double(mP1.x()), to_double(mP1.y()), to_double(mH->x()), to_double(mH->y()));
        mHandleGI->show();
      }
      else
      {
        HideHandle();
        mH = boost::optional<Point>();
      }
    }

    Point cvt ( typename Circular_curve::Point_2 const& aP ) { return Point( to_double(aP.x()), to_double(aP.y()) ) ; }

    void CommitCurrCircularPolygon()
    {
      GenerateCircularPolygon();

      mOngoingPieceCtr.clear();
      mOngoingPieceGI->modelChanged();

      mCircularPolygonPieces.clear();
      mCircularGI->modelChanged() ;

      mH = boost::optional<Point>();

      HideHandle();
    }

    void GenerateCircularPolygon()
    {
      if ( mCircularPolygonPieces.size() >  0 )
      {
        Gps_traits traits ;
        typename Gps_traits::Make_x_monotone_2 make_x_monotone = traits.make_x_monotone_2_object();

        std::vector<Circular_X_monotone_curve> xcvs;

        for ( const_circular_curve_iterator it = mCircularPolygonPieces.begin() ; it != mCircularPolygonPieces.end() ; ++ it )
        {
          std::vector<CGAL::Object>                 x_objs;
          std::vector<CGAL::Object>::const_iterator xoit;

          make_x_monotone ( *it, std::back_inserter (x_objs));

          for (xoit = x_objs.begin(); xoit != x_objs.end(); ++xoit)
          {
            Circular_X_monotone_curve xcv;
            if (CGAL::assign (xcv, *xoit))
              xcvs.push_back (xcv);
          }
        }

        if ( xcvs.size() > 0 )
        {
          Arc_point const& first_point = xcvs.front().source();
          Arc_point const& last_point =  xcvs.back ().target();
          CGAL_assertion(!first_point.x().is_extended() && !first_point.y().is_extended());
          CGAL_assertion(!last_point. x().is_extended() && !last_point .y().is_extended());
          FT fxs = first_point.x().alpha();
          FT fys = first_point.y().alpha();
          FT lxs = last_point .x().alpha();
          FT lys = last_point .y().alpha();
          xcvs.push_back(Circular_X_monotone_curve( Point(lxs,lys), Point(fxs,fys)));

          Circular_polygon cp(xcvs.begin(), xcvs.end());
          emit(generate(CGAL::make_object(cp)));
        }
      }
    }

  private:

    QGraphicsScene*    mScene ;
    GI*                mCircularGI ;
    GI*                mOngoingPieceGI ;
    QGraphicsLineItem* mHandleGI ;

    QPen mCircularPolygonPen ;
    QPen mOngoingCurvePen ;
    QPen mHandlePen ;

    Circular_curve_vector mCircularPolygonPieces ;
    Circular_curve_vector mOngoingPieceCtr ;

    int mState;

    Point mP0;
    Point mP1;

    boost::optional<Point> mH;

  }; // end class GraphicsViewCircularPolygonInput

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_VIEW_BEZIER_REGION_INPUT_H
