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
#include "GraphicsViewCurveInput.h"

template <typename Arr_>
class ArrangementCurveInputCallback:
  public CGAL::Qt::GraphicsViewCurveInput< typename Arr_::Geometry_traits_2 >
{
public:
  typedef Arr_                                          Arrangement;
  typedef typename Arrangement::Geometry_traits_2       Traits;
  typedef CGAL::Qt::GraphicsViewCurveInput< Traits >    Superclass;
  typedef typename Traits::Curve_2                      Curve_2;

  ArrangementCurveInputCallback(Arrangement* arrangement_, QObject* parent, QGraphicsScene* scene):
    Superclass( parent, scene ),
    arrangement( arrangement_)
  {
    this->setScene(scene);

    QObject::connect( this, SIGNAL( generate( CGAL::Object ) ),
                      this, SLOT( processInput( CGAL::Object ) ) );
  }

  void processInput( CGAL::Object o )
  {
    Curve_2 curve;
    if ( CGAL::assign( curve, o ) )
    {
      CGAL::insert( *( this->arrangement ), curve );
    }

    Q_EMIT CGAL::Qt::GraphicsViewInput::modelChanged( );
  }

  void setScene( QGraphicsScene* scene )
  {
    this->Superclass::setScene( scene );
  }

  void setArrangement( Arrangement* newArr )
  {
    this->arrangement = newArr;
  }

protected:
  Arrangement* arrangement;
}; // class ArrangementCurveInputCallback

#endif // ARRANGEMENT_SEGMENT_INPUT_CALLBACK_H
