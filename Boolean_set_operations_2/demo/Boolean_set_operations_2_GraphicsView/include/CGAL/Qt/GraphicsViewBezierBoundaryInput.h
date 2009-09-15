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
    typedef typename Bezier_gps_traits::General_polygon_2 Boundary;
    typedef typename Traits::Rat_kernel::Vector_2         Vector ;
    typedef typename Traits::Rat_kernel::Point_2          Point ;
    
    typedef std::vector<Curve> Curve_vector ;
    
    typedef typename Curve_vector::const_iterator const_curve_terator ;
    
    typedef Bezier_boundary_pieces_graphics_item<Curve_vector> GI ;

    GraphicsViewBezierBoundaryInput(QObject* parent, QGraphicsScene* s)
      : GraphicsViewInput(parent)
      , mScene(s)
      , mIsClosed(false)
      , mState(Start)
      , mOngoingPieceAdded(false)
      , mHandle0Item(NULL)
      , mHandle1Item(NULL)
      , mBezierPen( QColor(255,0,0) )
      , mHandlePen( QColor(0,0,255) )
    {
      mGI = new GI(&mPieces) ;
      mGI->setPen(mBezierPen);
      mScene->addItem(mGI);
    }
    
    ~GraphicsViewBezierBoundaryInput()
    {
      //mScene->removeItem(mGI);
      //delete mGI ;
      //RemoveHandleItems();      
    }
    
    bool eventFilter(QObject *obj, QEvent *event)
    {
      bool rHandled = false ;
      
      if (event->type() == QEvent::GraphicsSceneMousePress) 
      {
        rHandled = mousePressEvent( static_cast<QGraphicsSceneMouseEvent *>(event) ) ;
      } 
      else if (event->type() == QEvent::GraphicsSceneMouseRelease) 
      {
        rHandled = mouseReleaseEvent( static_cast<QGraphicsSceneMouseEvent *>(event) ) ;
      } 
      else if (event->type() == QEvent::GraphicsSceneMouseMove) 
      {
        rHandled = mouseMoveEvent( static_cast<QGraphicsSceneMouseEvent *>(event) ) ;
      } 
      else if (event->type() == QEvent::KeyPress) 
      {
        rHandled = keyPressEvent( static_cast<QKeyEvent *>(event) ) ;
      }
      
      if ( !rHandled )
        rHandled = QObject::eventFilter(obj, event);
              
      return rHandled ;  
    }
    
  protected:

    enum State
    {
      Start,
      CurveStarted,
      PieceStarted, PieceProceeding, PieceEnded,
      HandleStarted, HandleProceeding, HandleEnded
    } ;
    
    Point cvt ( QPointF const& aP ) const { return Point(aP.x(),aP.y()) ; }

    bool mousePressEvent(QGraphicsSceneMouseEvent *event)
    {
      bool rHandled = false ;
      
      Point lP = cvt(event->scenePos());
      
      if ( event->button() == ::Qt::LeftButton )
      {
        switch (mState)
        {
          case Start: 
            mState  = HandleStarted;
            mP1     = lP;
            rHandled = true;
            break;

          case PieceProceeding: 

            mState = PieceEnded;
            mP1 = lP;
            rHandled = true;
            break;
        }
      }
      
      return rHandled ;
    }
    
    bool mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
    {
      bool rHandled = false ;
      
      Point lP = cvt(event->scenePos());
      
      if ( event->button() == ::Qt::LeftButton )
      {
        switch (mState)
        {
          case HandleStarted:
            mState = PieceStarted;
            rHandled = true;
            break;
            
          case CurveStarted: 
            mState = PieceStarted;
            rHandled = true;
            break;

          case PieceEnded:
            CommitOngoingPiece();
            mState = PieceStarted;
            rHandled = true;
            break;

          case HandleProceeding: 
            CommitOngoingPiece();
            mState = HandleEnded;
            rHandled = true;
            break;
        }
      }
      else if ( event->button() == ::Qt::RightButton )
      {
        
      }
      
      return rHandled ;
    }

    bool mouseMoveEvent(QGraphicsSceneMouseEvent *event)
    {
      bool rHandled = false ;
      
      Point lP = cvt(event->scenePos());
      
      switch (mState)
      {
        case HandleStarted:
          mState = HandleProceeding ;
          mP0 = mP1;
          rHandled = true ;
          break;
          
        case PieceStarted: 
          mState = PieceProceeding;
          mP0 = mP1;
          rHandled = true ;
          break;

        case PieceProceeding: 
          mP1 = lP;
          UpdateOngoingPiece();
          rHandled = true ;
          break;

        case HandleProceeding:
        
          mState = HandleProceeding;
          
          if ( mPieces.size() > 0 )
          {
            mH0 = lP ;
            mH1 = mP1 - (lP - mP1);
            
            if ( !mHandle0Item )
            {
              mHandle0Item = new QGraphicsLineItem();
              mHandle0Item->setPen(mHandlePen);
              mScene->addItem(mHandle0Item);
            }
            
            if ( !mHandle1Item )
            {
              mHandle1Item = new QGraphicsLineItem();
              mHandle1Item->setPen(mHandlePen);
              mScene->addItem(mHandle1Item);
            }
            
            mHandle0Item->setLine( to_double(mP1.x()), to_double(mP1.y()), to_double(mH0->x()), to_double(mH0->y()));
            mHandle1Item->setLine( to_double(mP1.x()), to_double(mP1.y()), to_double(mH1->x()), to_double(mH1->y()));
          }
          else
          {
            mH1 = lP ;
            
            if ( !mHandle0Item )
            {
              mHandle0Item = new QGraphicsLineItem();
              mHandle0Item->setPen(mHandlePen);
              mScene->addItem(mHandle0Item);
            }
            
            mHandle0Item->setLine( to_double(mP0.x()), to_double(mP0.y()), to_double(mH1->x()), to_double(mH1->y()));
          }
          UpdateOngoingPiece();
          rHandled = true ;
          break;

        case HandleEnded: 
          mPrevH0 = mH0 ;
          RemoveHandleItems();
          mState = PieceStarted;
          rHandled = true ;
          break;
      }
      
      return rHandled ;
    }
    
    bool keyPressEvent(QKeyEvent *event)
    {
      return false ;
    }

    
    virtual void generate_boundary() 
    {
      Traits traits ;
      Traits::Make_x_monotone_2 make_x_monotone = traits.make_x_monotone_2_object();
      
      std::vector<X_monotone_curve> xcvs;

      for ( const_curve_terator it = mPieces.begin() ; it != mPieces.end() ; ++ it )
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
    }
    
  private:

    void RemoveHandleItems()
    {
      if ( mHandle0Item )
      {
        mScene->removeItem(mHandle0Item);
        delete mHandle0Item ;
        mHandle0Item = NULL ;
      }
      
      if ( mHandle1Item )
      {
        mScene->removeItem(mHandle1Item);
        delete mHandle1Item ;
        mHandle1Item = NULL ;
      }
    }  

    Curve CreatePiece()
    {
      if ( mPrevH0 && mH1 )
      {
        Point lControlPoints[4] = { mP0
                                  , midpoint(mP0,*mPrevH0)
                                  , midpoint(*mH1,mP1)
                                  , mP1 
                                  } ;
        return Curve( lControlPoints, lControlPoints + 4 ) ;
      }
      else if ( mPrevH0 && !mH1 )
      {
        Point lControlPoints[3] = { mP0
                                  , *mPrevH0
                                  , mP1 
                                  } ;
        return Curve( lControlPoints, lControlPoints + 3 ) ;
      }
      else if ( !mPrevH0 && mH1 )
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
    
    void AddPiece()
    {
      if ( mOngoingPieceAdded && mPieces.size() > 0 )
        mPieces.pop_back();
        
      mPieces.push_back( CreatePiece() ) ;
      
      mGI->modelChanged();
      
    }
    
    void UpdateOngoingPiece()
    {
      AddPiece();
      mOngoingPieceAdded = true ;
    }      
    
    void CommitOngoingPiece()
    {
      AddPiece();
      mOngoingPieceAdded = false ;
    }      
    
  private:
  
    QGraphicsScene*    mScene ;
    GI*                mGI ; 
    QGraphicsLineItem* mHandle0Item ;          
    QGraphicsLineItem* mHandle1Item ;          

    QPen mBezierPen ;
    QPen mHandlePen ;    
    
    Curve_vector mPieces ;  
    bool         mIsClosed;
    int          mState;
    bool         mOngoingPieceAdded;
    Point        mP0;
    Point        mP1;
    
    boost::optional<Point>  mPrevH0;
    boost::optional<Point>  mH0;
    boost::optional<Point>  mH1;
  
  }; // end class GraphicsViewBezierBoundaryInput

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_VIEW_BEZIER_BOUNDARY_INPUT_H
