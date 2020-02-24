// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#ifndef ARRANGEMENT_CURVE_INPUT_CALLBACK_H
#define ARRANGEMENT_CURVE_INPUT_CALLBACK_H

#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>
#include <QEvent>
#include <QGraphicsLineItem>
#include <QGraphicsSceneMouseEvent>
#include <iostream>
#include "GraphicsViewCurveInput.h"
#include "Utils.h"

template <typename Arr_, typename ArrTraits = typename Arr_::Geometry_traits_2>
class ArrangementCurveInputCallback:
  public CGAL::Qt::GraphicsViewCurveInput< typename Arr_::Geometry_traits_2 >
{
public:
  typedef Arr_ Arrangement;
  typedef ArrTraits Traits;
  typedef CGAL::Qt::GraphicsViewCurveInput< Traits >    Superclass;
  typedef typename Arrangement::Vertex_iterator         Vertex_iterator;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename Kernel::Point_2                      Kernel_point_2;
  typedef typename ArrTraitsAdaptor< Traits >::Point_2  Point_2;
  typedef typename Kernel::Segment_2                    Segment_2;
  typedef typename Kernel::FT                           FT;

  ArrangementCurveInputCallback( Arrangement* arrangement_, QObject* parent ):
    Superclass( parent ),
    arrangement( arrangement_ )
  {
    this->snapToVertexStrategy.setArrangement( arrangement_ );

    QObject::connect( this, SIGNAL( generate( CGAL::Object ) ),
                      this, SLOT( processInput( CGAL::Object ) ) );
  }

  void processInput( CGAL::Object o )
  {
    Curve_2 curve;
    X_monotone_curve_2 xcurve;
    if ( CGAL::assign( curve, o ) )
    {
      CGAL::insert( *( this->arrangement ), curve );
    }
#if 0
    else if ( CGAL::assign( xcurve, o ) )
    {
      std::vector< X_monotone_curve_2 > box;
      box.push_back( xcurve );
      CGAL::insert( *( this->arrangement ), box.begin( ), box.end( ) );
    }
#endif
    
    Q_EMIT CGAL::Qt::GraphicsViewInput::modelChanged( );
  }

  void setScene( QGraphicsScene* scene )
  {
    this->Superclass::setScene( scene );
    this->snapToVertexStrategy.setScene( scene );
    this->snapToGridStrategy.setScene( scene );
  }

  void setArrangement( Arrangement* newArr )
  {
    this->arrangement = newArr;
  }

protected:
  Point_2 snapPoint( QGraphicsSceneMouseEvent* event )
  {
    if ( this->snapToGridEnabled )
    {
      return this->snapToGridStrategy.snapPoint( event );
    }
    else if ( this->snappingEnabled )
    {
      return this->snapToVertexStrategy.snapPoint( event );
    }
    else
    {
      Kernel_point_2 p = this->convert( event->scenePos( ) );
      Point_2 res = this->toArrPoint( p );
      return res;
    }
  }

  Arrangement* arrangement;
  SnapToArrangementVertexStrategy< Arrangement > snapToVertexStrategy;
  SnapToGridStrategy< typename Arrangement::Geometry_traits_2 >
    snapToGridStrategy;
  Arr_construct_point_2< Traits > toArrPoint;
}; // class ArrangementCurveInputCallback

#endif // ARRANGEMENT_SEGMENT_INPUT_CALLBACK_H
