// Copyright (c) 2010  GeometryFactory Sarl (France).
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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Boolean_set_operations_2/demo/Boolean_set_operations_2_GraphicsView/include/CGAL/Qt/GraphicsViewCircularPolygonInput.h $
// $Id: GraphicsViewCircularPolygonInput.h 53502 2009-12-18 16:47:24Z fcacciola $
// 
//
// Author(s)     : Fernando Cacciola <Fernando.Cacciola@geometryfactory.com>
//

#ifndef CGAL_QT_GRAPHICS_VIEW_CIRCULAR_POLYGON_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_CIRCULAR_POLYGON_INPUT_H

#include <CGAL/auto_link/Qt4.h>

#include <QPolygonF>
#include <QPointF>
#include <QGraphicsLineItem> 

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
            UpdateHandle(lP);
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
            CloseCurrBundary();
            CommitCurrCircularPolygon();
            mH = boost::optional<Point>();
            mState = Start ;     
            rHandled = true;
            break;
        }    
      }
      
      return rHandled ;
    }
    
    bool keyPressEvent(QKeyEvent *aEvent)
    {
      bool rHandled = false ;
      
      
      return rHandled ;
    }

    
    
  private:

    Circular_curve const* ongoing_piece() const { return mOngoingPieceCtr.size() == 1 ? &mOngoingPieceCtr[0] : NULL ; }
    
    void HideHandle()
    {
      mHandleGI->hide();
    }  

    Circular_curve CreatePiece()
    {
      return mH ? Circular_curve(mP0,*mH,mP1) : Circular_curve(mP0,mP1) ;
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
      if ( squared_distance(mP1,aP) >= 9 )
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
          

    void CloseCurrBundary()
    {
      if ( mCircularPolygonPieces.size() > 0 && ongoing_piece()!= NULL )
      {
       // mCircularPolygonPieces.push_back( Circular_curve( ongoing_piece()->target(), mCircularPolygonPieces.front().source() ) ) ;
        
        mCircularGI->modelChanged() ;
      }
    }
        
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
        Circular_polygon cp(xcvs.begin(), xcvs.end());
        emit(generate(CGAL::make_object(cp)));
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
