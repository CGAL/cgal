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

#ifndef ENVELOPE_CALLBACK_H
#define ENVELOPE_CALLBACK_H

#include "Callback.h"
#include "CurveGraphicsItem.h"

class EnvelopeCallbackBase : public CGAL::Qt::Callback
{
public:
  virtual void setEnvelopeEdgeColor( const QColor& color ) = 0;
  virtual const QColor& getEnvelopeEdgeColor( ) const = 0;
  virtual void setEnvelopeEdgeWidth( int width ) = 0;
  virtual int getEnvelopeEdgeWidth( ) const = 0;
  virtual void setEnvelopeVertexColor( const QColor& color ) = 0;
  virtual const QColor& getEnvelopeVertexColor( ) const = 0;
  virtual void setEnvelopeVertexRadius( int radius ) = 0;
  virtual int getEnvelopeVertexRadius( ) const = 0;

public Q_SLOTS:
  virtual void showLowerEnvelope( bool b ) = 0;
  virtual void showUpperEnvelope( bool b ) = 0;

protected:
  EnvelopeCallbackBase( QObject* parent );
}; // class EnvelopeCallbackBase

/**
   Updates and draws the lower and upper envelopes of an observed arrangement.
*/
template < typename Arr_, typename Traits = typename Arr_::Geometry_traits_2 >
class EnvelopeCallback : public EnvelopeCallbackBase
{
public:
  typedef Arr_                                          Arrangement;
  typedef typename Arrangement::Edge_iterator           Edge_iterator;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename Kernel::Point_2                      Kernel_point_2;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Kernel::Segment_2                    Segment_2;
  typedef typename Kernel::Ray_2                        Ray_2;
  typedef typename Kernel::Line_2                       Line_2;

  /**
     Construct an envelope callback observing the given arrangement.
  */
  EnvelopeCallback( Arrangement* arr_, QObject* parent );

  /**
     Enable/disable drawing the lower envelope.
  */
  void showLowerEnvelope( bool show );

  /**
     Enable/disable drawing the lower envelope.
  */
  void showUpperEnvelope( bool show );

  /**
     Slot: Update and redraw the envelopes.
  */
  void slotModelChanged( ) override;

  void setEnvelopeEdgeColor( const QColor& color )
  {
    this->lowerEnvelope->setEdgeColor( color );
    this->upperEnvelope->setEdgeColor( color );
  }

  const QColor& getEnvelopeEdgeColor( ) const
  {
    return this->lowerEnvelope->edgeColor( );
  }

  void setEnvelopeEdgeWidth( int width )
  {
    this->lowerEnvelope->setEdgeWidth( width );
    this->upperEnvelope->setEdgeWidth( width );
  }

  int getEnvelopeEdgeWidth( ) const
  {
    return this->lowerEnvelope->edgeWidth( );
  }

  void setEnvelopeVertexColor( const QColor& color )
  {
    this->lowerEnvelope->setVertexColor( color );
    this->upperEnvelope->setVertexColor( color );
  }

  const QColor& getEnvelopeVertexColor( ) const
  {
    return this->lowerEnvelope->vertexColor( );
  }

  void setEnvelopeVertexRadius( int radius )
  {
    this->lowerEnvelope->setVertexRadius( radius );
    this->upperEnvelope->setVertexRadius( radius );
  }

  int getEnvelopeVertexRadius( ) const
  {
    return this->lowerEnvelope->vertexRadius( );
  }

protected:
  void setScene( QGraphicsScene* scene_ ) override;

  /**
     Helper method to update the upper/lower envelope.
  */
  void updateEnvelope( bool lower );

  template < typename TTraits >
  void updateEnvelope( bool lower, TTraits traits );

  template < typename Kernel_ >
  void updateEnvelope( bool lower, 
                      CGAL::Arr_linear_traits_2< Kernel_ > traits );

#if 0
  template < typename CircularKernel >
  void updateEnvelope(bool lower,
                      CGAL::Arr_circular_arc_traits_2<CircularKernel> traits);
#endif
  template < typename Coefficient_ >
  void updateEnvelope(bool lower,
                      CGAL::Arr_algebraic_segment_traits_2<Coefficient_> traits);
 
  Construct_x_monotone_subcurve_2< Traits > construct_x_monotone_subcurve_2;
  Arrangement* arr;
  CGAL::Qt::CurveGraphicsItem< Traits >* lowerEnvelope;
  CGAL::Qt::CurveGraphicsItem< Traits >* upperEnvelope;
  bool showLower;
  bool showUpper;
}; // class EnvelopeCallback

#endif // ENVELOPE_CALLBACK_H
