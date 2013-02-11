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
// $URL: $
// $Id: $
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#ifndef ENVELOPE_CALLBACK_H
#define ENVELOPE_CALLBACK_H

#include "Callback.h"
#include "CurveGraphicsItem.h"

#include <CGAL/envelope_2.h>
#include <CGAL/Envelope_diagram_1.h>

#include <list>
#include "Utils.h"

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

public slots:
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
  // typedef typename Arrangement::Geometry_traits_2    Traits;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  // typedef typename Traits::Construct_x_monotone_curve_2
  //   Construct_x_monotone_curve_2;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename Kernel::Point_2                      Kernel_point_2;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Kernel::Segment_2                    Segment_2;
  typedef typename Kernel::Ray_2                        Ray_2;
  typedef typename Kernel::Line_2                       Line_2;
  typedef CGAL::Envelope_diagram_1< Traits >            Diagram_1;

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
  void slotModelChanged( );

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

  //  shouldn't need this here, since it is in the base class Callback
  //    void setScene( QGraphicsScene* scene_ );

protected:
  /**
     Helper method to update the upper/lower envelope.
  */
  void updateEnvelope( bool lower );

  template < typename TTraits >
  void updateEnvelope( bool lower, TTraits traits );

  template < typename CircularKernel >
  void updateEnvelope(bool lower,
                      CGAL::Arr_circular_arc_traits_2<CircularKernel> traits);

  template < typename Coefficient_ >
  void updateEnvelope(bool lower,
                      CGAL::Arr_algebraic_segment_traits_2<Coefficient_>
                      traits);
  
  Construct_x_monotone_subcurve_2< Traits > construct_x_monotone_subcurve_2;
  Arrangement* arr;
  CGAL::Qt::CurveGraphicsItem< Traits >* lowerEnvelope;
  CGAL::Qt::CurveGraphicsItem< Traits >* upperEnvelope;
  using CGAL::Qt::Callback::scene;
}; // class EnvelopeCallback

template < typename Arr_, typename Traits >
EnvelopeCallback<Arr_, Traits>::EnvelopeCallback(Arrangement* arr_,
                                                 QObject* parent) :
  EnvelopeCallbackBase( parent ),
  arr( arr_ ),
  lowerEnvelope( new CGAL::Qt::CurveGraphicsItem< Traits >( ) ),
  upperEnvelope( new CGAL::Qt::CurveGraphicsItem< Traits >( ) )
{
  this->lowerEnvelope->hide( );
  this->upperEnvelope->hide( );
}

template < typename Arr_, typename Traits >
void EnvelopeCallback< Arr_, Traits >::slotModelChanged( )
{
  this->updateEnvelope( true );
  this->updateEnvelope( false );
}

template < typename Arr_, typename Traits >
void EnvelopeCallback< Arr_, Traits >::updateEnvelope( bool lower )
{
  this->updateEnvelope( lower, Traits( ) );
}

template < typename Arr_, typename Traits >
template < typename TTraits >
void EnvelopeCallback< Arr_, Traits >::updateEnvelope(bool lower,
                                                      TTraits /* traits */)
{
  CGAL::Qt::CurveGraphicsItem< Traits >* envelopeToUpdate;
  if ( lower )
  {
    envelopeToUpdate = this->lowerEnvelope;
  }
  else
  {
    envelopeToUpdate = this->upperEnvelope;
  }
  envelopeToUpdate->clear( );

  std::list< X_monotone_curve_2 > curves;
  Edge_iterator eit;
  for (eit = this->arr->edges_begin( ); eit != this->arr->edges_end( ); ++eit)
  {
    curves.push_back( eit->curve( ) );
  }
  Diagram_1 diagram;
  if ( lower )
  {
    CGAL::lower_envelope_x_monotone_2(curves.begin(), curves.end(), diagram);
  }
  else
  {
    CGAL::upper_envelope_x_monotone_2(curves.begin(), curves.end(), diagram);
  }

  typename Diagram_1::Edge_const_handle e = diagram.leftmost( );
  typename Diagram_1::Vertex_const_handle v;
  QRectF clipRect = this->viewportRect( );
  CGAL::Qt::Converter< Kernel > convert( clipRect );
  while ( e != diagram.rightmost( ) )
  {
    if ( ! e->is_empty( ) )
    {
      // The edge is not empty: draw a representative curve.
      // Note that the we only draw the portion of the curve
      // that overlaps the x-range defined by the two vertices
      // that are incident to this edge.

      // TODO: generate a subcurve instead of just making a segment

      Point_2 leftPoint, rightPoint;
      if ( e->left( ) != NULL )
      {
        leftPoint = e->left( )->point( );
      }
      else
      {
        // std::cout << "handle unbounded curve" << std::endl;
        v = e->right( );
        e = v->right( );
        continue;
      }

      if ( e->right( ) != NULL )
      {
        rightPoint = e->right( )->point( );
      }
      else
      {
        // std::cout << "pRight is null; should never get here..."
        //           << std::endl;
      }
      X_monotone_curve_2 curve =
        this->construct_x_monotone_subcurve_2(e->curve(),
                                              leftPoint, rightPoint);
      envelopeToUpdate->insert( curve );
      envelopeToUpdate->insert( leftPoint );
      envelopeToUpdate->insert( rightPoint );
    }
    v = e->right( );

    // TODO: Draw the point associated with the current vertex.
    e = v->right( );
  }
  envelopeToUpdate->modelChanged( );
}

template < typename Arr_, typename Traits >
template < typename CircularKernel >
void EnvelopeCallback< Arr_, Traits >::
updateEnvelope(bool lower,
               CGAL::Arr_circular_arc_traits_2< CircularKernel > /* traits */)
{
  typedef typename Traits::Point_2 Arc_point_2;
  CGAL::Qt::CurveGraphicsItem< Traits >* envelopeToUpdate;
  if ( lower )
  {
    envelopeToUpdate = this->lowerEnvelope;
  }
  else
  {
    envelopeToUpdate = this->upperEnvelope;
  }
  envelopeToUpdate->clear( );

  std::list< X_monotone_curve_2 > curves;
  for ( Edge_iterator eit =
          this->arr->edges_begin( ); eit != this->arr->edges_end( ); ++eit )
  {
    curves.push_back( eit->curve( ) );
  }
  Diagram_1 diagram;
  if ( lower )
  {
    CGAL::lower_envelope_x_monotone_2( curves.begin(), curves.end(), diagram);
  }
  else
  {
    CGAL::upper_envelope_x_monotone_2(curves.begin( ), curves.end(), diagram);
  }

  typename Diagram_1::Edge_const_handle e = diagram.leftmost( );
  typename Diagram_1::Vertex_const_handle v;
  QRectF clipRect = this->viewportRect( );
  CGAL::Qt::Converter< Kernel > convert( clipRect );
  while ( e != diagram.rightmost( ) )
  {
    if ( ! e->is_empty( ) )
    {
      // The edge is not empty: draw a representative curve.
      // Note that the we only draw the portion of the curve
      // that overlaps the x-range defined by the two vertices
      // that are incident to this edge.

      // TODO: generate a subcurve instead of just making a segment

      Arc_point_2 leftPoint, rightPoint;
      if ( e->left( ) != NULL )
      {
        leftPoint = e->left( )->point( );
      }
      else
      {
        // std::cout << "handle unbounded curve" << std::endl;
        v = e->right( );
        e = v->right( );
        continue;
      }

      if ( e->right( ) != NULL )
      {
        rightPoint = e->right( )->point( );
      }
      else
      {
        // std::cout << "pRight is null; should never get here..."
        //           << std::endl;
      }
      X_monotone_curve_2 curve =
        this->construct_x_monotone_subcurve_2( e->curve( ),
                                               leftPoint, rightPoint );
      envelopeToUpdate->insert( curve );
    }
    v = e->right( );

    // TODO: Draw the point associated with the current vertex.
    e = v->right( );
  }
  envelopeToUpdate->modelChanged( );
}

template < typename Arr_, typename Traits >
template < typename Coefficient_ >
void EnvelopeCallback< Arr_, Traits >::
updateEnvelope(bool /* lower */,
               CGAL::Arr_algebraic_segment_traits_2<Coefficient_> /* traits */)
{
  // std::cout << "alg seg envelope stub" << std::endl;
}

template < typename Arr_, typename Traits >
void EnvelopeCallback< Arr_, Traits >::showLowerEnvelope( bool show )
{
  if ( show )
  {
    // std::cout << "Show lower envelope" << std::endl;
    this->scene->addItem( this->lowerEnvelope );
  }
  else
  {
    // std::cout << "Hide lower envelope" << std::endl;
    this->scene->removeItem( this->lowerEnvelope );
  }
}

template < typename Arr_, typename Traits >
void EnvelopeCallback< Arr_, Traits >::showUpperEnvelope( bool show )
{
  if ( show )
  {
    // std::cout << "Show upper envelope" << std::endl;
    this->scene->addItem( this->upperEnvelope );
  }
  else
  {
    // std::cout << "Hide upper envelope" << std::endl;
    this->scene->removeItem( this->upperEnvelope );
  }
}

#endif // ENVELOPE_CALLBACK_H
