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
//

#ifndef CGAL_QT_GRAPHICS_VIEW_BEZIER_BOUNDARY_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_BEZIER_BOUNDARY_INPUT_H

#include <CGAL/auto_link/Qt4.h>

#include <QPolygonF>
#include <QPointF>
#include <QGraphicsLineItem> 

#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/BezierCurves.h>

namespace CGAL {

namespace Qt {

  template <class Traits_>
  class GraphicsViewBezierBoundaryInput : public GraphicsViewInput
  {
  public:

    typedef Traits_ Traits ;
    
    typedef CGAL::Gps_traits_2<Traits> Bezier_gps_traits;
    
    typedef typename Traits::Curve_2                      Curve;
    typedef typename Traits::X_monotone_curve_2           X_monotone_curve;
    typedef typename Bezier_gps_traits::General_polygon_2 Bezier_Boundary;
    typedef typename Traits::Rat_kernel::Vector_2         Vector ;
    typedef typename Traits::Rat_kernel::Point_2          Point ;
    
    typedef std::vector<Curve> Curve_vector ;
    
    struct Boundary
    {
      Boundary() : mIsClosed(false) {}
      
      Curve_vector mPieces ;
      bool         mIsClosed ; 
    } ;

    typedef boost::shared_ptr<Boundary> Boundary_ptr ;
    
    typedef std::vector<Boundary_ptr> Boundary_vector ;
    
    typedef typename Curve_vector::const_iterator const_curve_terator ;
    
    typedef Bezier_boundary_pieces_graphics_item<Curve_vector> GI ;

    GraphicsViewBezierBoundaryInput(QObject* aParent, QGraphicsScene* aScene)
      :
        GraphicsViewInput( aParent         )
      , mScene           ( aScene          )
      , mState           ( Start           )
      , mBoundaryPen     ( QColor(0,255,0) )
      , mOngoingCurvePen ( QColor(255,0,0) )
      , mHandlePen       ( QColor(0,0,255) )
      , mCurrBoundary    ( new Boundary()  ) 
      , mCurrBoundaryGI  ( 0               )
    {
      mOngoingPieceGI = new GI(&mOngoingPieceCtr) ;
      mHandle0GI      = new QGraphicsLineItem();
      mHandle1GI      = new QGraphicsLineItem();
      
      mOngoingPieceGI->setPen(mOngoingCurvePen);
      mHandle0GI     ->setPen(mHandlePen);
      mHandle1GI     ->setPen(mHandlePen);
      
      mHandle0GI->setLine(0,0,1,1);
      mHandle1GI->setLine(0,0,1,1);
      mHandle0GI->hide();
      mHandle1GI->hide();
      
      mScene->addItem(mOngoingPieceGI);
      mScene->addItem(mHandle0GI);
      mScene->addItem(mHandle1GI);
      
      SetupCurrBoundary();
    }
    
    ~GraphicsViewBezierBoundaryInput()
    {
      //mScene->removeItem(mGI);
      //delete mGI ;
      //RemoveHandleItems();      
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

    enum State { Start, PieceOrHandleStarted, PieceOngoing, HandleOngoing, PieceEnded, CurveEnded } ;
    
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
            mState   = PieceOrHandleStarted;
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
    
    bool mouseReleaseEvent(QGraphicsSceneMouseEvent *aEvent)
    {
      bool rHandled = false ;
      
      Point lP = cvt(aEvent->scenePos());
      
      if ( aEvent->button() == ::Qt::LeftButton )
      {
        switch (mState)
        {
          case PieceOrHandleStarted:
            mState   = PieceOngoing;
            rHandled = true;
            break;
            
          case HandleOngoing: 
            UpdateHandles(lP);
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
            CommitCurrBoundary();
            mState = Start ;     
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
        
          UpdateHandles(lP);
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
    
    bool keyPressEvent(QKeyEvent *aEvent)
    {
      bool rHandled = false ;
      
      if (aEvent->text() == "c" )
      {
        CloseCurrBundary();
        CommitCurrBoundary();
        mState = Start ;     
        rHandled = true;
      }
      
      return rHandled ;
    }

    
    virtual void generate_boundary() 
    {
/*
      Traits traits ;
      Traits::Make_x_monotone_2 make_x_monotone = traits.make_x_monotone_2_object();
      
      std::vector<X_monotone_curve> xcvs;

      for ( const_curve_terator it = mCurrBdryPieces.begin() ; it != mCurrBdryPieces.end() ; ++ it )
      {       
        std::vector<CGAL::Object>                 x_objs;
        std::vector<CGAL::Object>::const_iterator xoit;
        
        make_x_monotone ( *it, std::back_inserter (x_objs));
        
        for (xoit = x_objs.begin(); xoit != x_objs.end(); ++xoit) 
        {
          X_monotone_curve xcv;
          if (CGAL::assign (xcv, *xoit))
            xcvs.push_back (xcv);
        }    
      }
      emit(generate(CGAL::make_object( Boundary(xcvs.begin(), xcvs.begin()) )));
*/
    }
    
  private:

    Curve const* ongoing_piece() const { return mOngoingPieceCtr.size() == 1 ? &mOngoingPieceCtr[0] : NULL ; }
    
    void HideHandles()
    {
      mHandle0GI->hide();
      mHandle1GI->hide();
    }  

    Curve CreatePiece()
    {
      if ( mPrevH0 && mH1 && *mPrevH0 != *mH1 && *mPrevH0 != mP0 && *mH1 != mP1 )
      {
        Point lControlPoints[4] = { mP0
                                  , *mPrevH0
                                  , *mH1
                                  , mP1 
                                  } ;
        return Curve( lControlPoints, lControlPoints + 4 ) ;
      }
      else if ( mPrevH0 && !mH1 && *mPrevH0 != mP0 && *mPrevH0 != mP1 )
      {
        Point lControlPoints[3] = { mP0
                                  , *mPrevH0
                                  , mP1 
                                  } ;
        return Curve( lControlPoints, lControlPoints + 3 ) ;
      }
      else if ( !mPrevH0 && mH1 && *mH1 != mP0 && *mH1 != mP1 )
      {
        Point lControlPoints[3] = { mP0
                                  , *mH1
                                  , mP1 
                                  } ;
        return Curve( lControlPoints, lControlPoints + 3 ) ;
      }
      else
      {
        Point lControlPoints[2] = { mP0
                                  , mP1 
                                  } ;
        return Curve( lControlPoints, lControlPoints + 2 ) ;
      }
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
      mCurrBoundary->mPieces.push_back( *ongoing_piece() ) ;
      mCurrBoundaryGI->modelChanged();
      mOngoingPieceCtr.clear();
      mOngoingPieceGI->modelChanged();
      mP0 = mP1 ;
      mP1 = aP ;
      mPrevH0 = mH0 ;
      mH0 = mH1 = boost::optional<Point>();
    }      

    void UpdateHandles(Point const& aP)
    {
      if ( squared_distance(mP1,aP) >= 9 )
      {
        if ( ongoing_piece() )
        {
          mH0 = aP ;
          mH1 = mP1 - (aP - mP1);
          
          mHandle0GI->setLine( to_double(mP1.x()), to_double(mP1.y()), to_double(mH0->x()), to_double(mH0->y()));
          mHandle1GI->setLine( to_double(mP1.x()), to_double(mP1.y()), to_double(mH1->x()), to_double(mH1->y()));
          mHandle0GI->show();
          mHandle1GI->show();
        }
      }
      else
      {
        HideHandles();
        mH0 = mH1 = boost::optional<Point>();
      }
      
    }    

    void CloseCurrBundary()
    {
      if ( mCurrBoundary->mPieces.size() > 0 && ongoing_piece()!= NULL )
      {
        std::vector<Point> lControlPoints(ongoing_piece()->control_points_begin(),ongoing_piece()->control_points_end());
        
        lControlPoints.back() = mCurrBoundary->mPieces.front().control_point(0);
              
        mCurrBoundary->mPieces.push_back( Curve( lControlPoints.begin(), lControlPoints.end() ) ) ;
        
        mCurrBoundary->mIsClosed = true ;
      }
    }
        
    void CommitCurrBoundary()
    {
      mOngoingPieceCtr.clear();
      mOngoingPieceGI->modelChanged();
      
      if ( mCurrBoundary->mPieces.size() > 0 )
      {
        mComittedBoundaries.push_back(mCurrBoundary) ;
        
        mCommittedBoundariesGI.push_back(mCurrBoundaryGI);
      }
      
      SetupCurrBoundary();
      
      HideHandles();
    }
    
    void SetupCurrBoundary()
    {
      mCurrBoundary = Boundary_ptr( new Boundary() ) ;
      
      mCurrBoundaryGI = new GI(&mCurrBoundary->mPieces) ;
      
      mCurrBoundaryGI->setPen(mBoundaryPen);
      
      mScene->addItem(mCurrBoundaryGI);
    }
    
  private:
  
    QGraphicsScene*    mScene ;
    std::vector<GI*>   mCommittedBoundariesGI ;
    GI*                mCurrBoundaryGI ; 
    GI*                mOngoingPieceGI ; 
    QGraphicsLineItem* mHandle0GI ;          
    QGraphicsLineItem* mHandle1GI ;          

    QPen mBoundaryPen ;
    QPen mOngoingCurvePen ;
    QPen mHandlePen ;    
    
    Boundary_vector mComittedBoundaries ;  
    Boundary_ptr    mCurrBoundary ;  
    Curve_vector    mOngoingPieceCtr ;  

    int mState;
    
    Point mP0;
    Point mP1;
    
    boost::optional<Point> mPrevH0;
    boost::optional<Point> mH0;
    boost::optional<Point> mH1;
  
  }; // end class GraphicsViewBezierBoundaryInput

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_VIEW_BEZIER_BOUNDARY_INPUT_H
