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

#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/BezierCurves.h>

class QGraphicsScene;
class QGraphicsSceneMouseEvent;
class QGraphicsItem;
class QGraphicsPathItem;
class QKeyEvent;
class QEvent;
class QObject;

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
      , mBezierPhase(None)
      , mOngoingPieceAdded(false)
    {
      mGI = boost::shared_ptr<GI>( new GI(&mPieces)  ) ;
      mScene->addItem(mGI.get());
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
      HandleProceeding, HandleEnded
    } ;
    
    enum BezierPhase { None, Manual, Automatic } ;
    
    enum PieceType   { Unknown, Cubic } ;
    
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
            mState = CurveStarted;
            mP1 = lP;
            mPrevP1 = mP1;
            rHandled = true;
            break;

          case PieceProceeding: 

            mState = PieceEnded;
            mP1 = lP;
            bool lAutoClose = false;
            if (!mPrevP1 || (!!mPrevP1 && (*mPrevP1 != mP1)) )
            {
              mPrevP1 = mP1;

              AddPiece(Automatic);
            }

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
          case CurveStarted: 
            mState = PieceStarted;
            rHandled = true;
            break;

          case PieceEnded:
          
            mState = PieceStarted;
            mBezierPhase = None;
            rHandled = true;
            break;

          case HandleProceeding: 
            mState = HandleEnded;
            AddPiece(Cubic);
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
        case PieceStarted: 
          mState = PieceProceeding;
          mP0 = mP1;
          rHandled = true ;
          break;

        case PieceProceeding: 
          mState = PieceProceeding;
          mP1 = lP;
          rHandled = true ;
          break;

        case PieceEnded:
          if ( squared_distance(mP1, lP) > 4 )
          {
            mState = HandleProceeding;
            mBezierPhase = Manual;
            mPieces.pop_back();
            mH0 = lP;
            mH1 = mP1 - (mH0 - mP1);
            rHandled = true ;
          }
          break;

        case HandleProceeding: 
          mState = HandleProceeding;
          mH0 = lP;
          mH1 = mP1 - (mH0 - mP1);
          rHandled = true ;
          break;

        case HandleEnded: 
          mState = PieceStarted;
          mBezierPhase = Automatic;
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
  
    void TakebackOngoingPiece()
    {
      if (mOngoingPieceAdded)
      {
        mPieces.pop_back();
        mOngoingPieceAdded = false;
        mGI->modelChanged();
      }
    }

    // The end-points of the user-visible handles are not directly the knots of a bezier.
    void CalculateK(Point const& aP0, Point const& aP1, Point const& aP2, Point const& aP3,Point& aK0, Point& aK1 )
    {
      aK0 = midpoint(aP0, aP1);
      aK1 = midpoint(aP2, aP3);
    }

    // Adds a new piece.
    void AddPiece(int aType)
    {
      Point lControlPoints[4] = { mP0, mP0, mP1, mP1 } ;

      // This is null if the previous piece is not a bezier
      Curve const* lPreviousPiece = mPieces.size() > 0 ? & mPieces.back() : NULL ;

      // If "automatic" the piece is a bezier of degree 2 whenever 
      //   (a) the previous piece is also a bezier
      //   (b) and it is a manual bezier (so the current Bezier Phase is automatic)
      if (aType == Automatic)
      {
        bool lCurrIsCubic = !!lPreviousPiece && mBezierPhase == Automatic ;
        if (lCurrIsCubic)
        {
          // For an automatic degree 2 bezier the kont is taken as the reflexion
          // of the previous second knot around the current piece starting-point (mP0)
          Vector lV = (lPreviousPiece->control_point(2) - mP0);
          Point lH0 = mP0 - lV - lV;
          CalculateK(mP0, lH0, lH0, mP1, lControlPoints[1], lControlPoints[2]);
        }
        else
        {
          // The piece is a straight-line
          lControlPoints[1] = ORIGIN + ( ( 2.0 * ( mP0 - ORIGIN ) ) + ( mP1 - ORIGIN )           ) / 3.0;
          lControlPoints[1] = ORIGIN + ( ( mP0 - ORIGIN )           + ( 2.0 * ( mP1 - ORIGIN ) ) ) / 3.0;
        }  
      }
      else // aType == Cubic
      {
        // If the previous piece is a bezier this one is of degree 3.
        if ( !!lPreviousPiece)
        {
          // For a manual bezier of degree 3 the first kont is taken as the reflexion
          // of the previous second knot around the current piece starting-point (mP0)
          Vector lU = (lPreviousPiece->control_point(2) - mP0);
          Point lH0 = mP0 - lU - lU;
          CalculateK(mP0, lH0, mH1, mP1, lControlPoints[1], lControlPoints[2]);
        }
        else
        {
          CalculateK(mP0, mH1, mH1, mP1, lControlPoints[1], lControlPoints[2]);
        }
      }

      TakebackOngoingPiece();
      
      mPieces.push_back( Curve( lControlPoints, lControlPoints + 4 ) ) ;
      
      mGI->modelChanged();
    }

  private:
  
    QGraphicsScene*       mScene ;
    boost::shared_ptr<GI> mGI ;           
    
    Curve_vector mPieces ;  
    bool         mIsClosed;
    int          mState;
    int          mBezierPhase;
    bool         mOngoingPieceAdded;
    Point        mP0;
    Point        mP1;
    Point        mH0;
    Point        mH1;
    
    boost::optional<Point>  mPrevP1;
    
  
  }; // end class GraphicsViewBezierBoundaryInput

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_VIEW_BEZIER_BOUNDARY_INPUT_H
