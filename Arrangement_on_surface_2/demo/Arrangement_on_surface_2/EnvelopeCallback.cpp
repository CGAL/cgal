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

#include "EnvelopeCallback.h"
#include "Utils.h"

#include <CGAL/envelope_2.h>
#include <CGAL/Envelope_diagram_1.h>

#include <list>

EnvelopeCallbackBase::EnvelopeCallbackBase( QObject* parent ) :
  CGAL::Qt::Callback( parent )
{ }

template < typename Arr_, typename Traits >
EnvelopeCallback<Arr_, Traits>::EnvelopeCallback(Arrangement* arr_,
                                                 QObject* parent) :
  EnvelopeCallbackBase( parent ),
  arr( arr_ ),
  lowerEnvelope( new CGAL::Qt::CurveGraphicsItem< Traits >( ) ),
  upperEnvelope( new CGAL::Qt::CurveGraphicsItem< Traits >( ) ),
  showLower( false ),
  showUpper( false )
{
  this->lowerEnvelope->hide( );
  this->upperEnvelope->hide( );
}

template < typename Arr_, typename Traits >
void EnvelopeCallback<Arr_, Traits>::setScene( QGraphicsScene* scene_ )
{
  EnvelopeCallbackBase::setScene(scene_);
  lowerEnvelope->setScene(scene_);
  upperEnvelope->setScene(scene_);
}

template < typename Arr_, typename Traits >
void EnvelopeCallback< Arr_, Traits >::slotModelChanged( )
{
  if (showLower)
    this->updateEnvelope( true );
  if (showUpper)
    this->updateEnvelope( false );
}

template < typename Arr_, typename Traits >
void EnvelopeCallback< Arr_, Traits >::updateEnvelope( bool lower )
{
  this->updateEnvelope( lower, Traits( ) );
}

template < typename Arr_, typename Traits >
template < typename Kernel_ >
void EnvelopeCallback< Arr_, Traits >::
updateEnvelope( bool lower,
                    CGAL::Arr_linear_traits_2< Kernel_ > traits )
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

  typedef CGAL::Envelope_diagram_1< Traits > Diagram_1;
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

  typedef CGAL::Arr_linear_traits_2< Kernel_ > Trait;
  Arr_compute_y_at_x_2< Trait > compute_y_at_x_2;

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
        // v = e->right( );
        // e = v->right( );
        // continue;
        double leftPoint_y = compute_y_at_x_2.approx(e->curve(), clipRect.left());
        leftPoint = Point_2(clipRect.left(), leftPoint_y);
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

  if (e == diagram.rightmost( ))
  {
    if ( ! e->is_empty( ) )
    {
      Point_2 leftPoint, rightPoint;
      if ( e->left( ) != NULL )
      {
        leftPoint = e->left( )->point( );
      }
      else
      {
        double leftPoint_y = compute_y_at_x_2.approx(e->curve(), clipRect.left());
        leftPoint = Point_2(clipRect.left(), leftPoint_y);
      }

      if ( e->right( ) != NULL )
      {
        rightPoint = e->right( )->point( );
      }
      else
      {
        // std::cout << "pRight is null; should never get here..."
        //           << std::endl;
        double rightPoint_y = compute_y_at_x_2.approx(e->curve(), clipRect.right());
        rightPoint = Point_2(clipRect.right(), rightPoint_y);
      }

      X_monotone_curve_2 curve =
        this->construct_x_monotone_subcurve_2(e->curve(),
                                              leftPoint, rightPoint);
      envelopeToUpdate->insert( curve );

    }
    else
    {
    }
  }

  envelopeToUpdate->modelChanged( );
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

  typedef CGAL::Envelope_diagram_1< Traits > Diagram_1;
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
template < typename Coefficient_ >
void EnvelopeCallback< Arr_, Traits >::
updateEnvelope(bool lower,
                CGAL::Arr_algebraic_segment_traits_2<Coefficient_ > )
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
  for ( Edge_iterator eit =
          this->arr->edges_begin( ); eit != this->arr->edges_end( ); ++eit )
  {
    curves.push_back( eit->curve( ) );
  }

  typedef CGAL::Envelope_diagram_1< Traits > Diagram_1;
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

  Arr_compute_y_at_x_2< Traits > compute_y_at_x_2;

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
        // v = e->right( );
        // e = v->right( );
        // continue;
        double leftPoint_y = compute_y_at_x_2.approx(e->curve(), clipRect.left());
        leftPoint = Point_2(clipRect.left(), leftPoint_y);
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

  if (e == diagram.rightmost( ))
  {
    if ( ! e->is_empty( ) )
    {
      Point_2 leftPoint, rightPoint;
      if ( e->left( ) != NULL )
      {
        leftPoint = e->left( )->point( );
      }
      else
      {
        double leftPoint_y = compute_y_at_x_2.approx(e->curve(), clipRect.left());
        leftPoint = Point_2(clipRect.left(), leftPoint_y);
      }

      if ( e->right( ) != NULL )
      {
        rightPoint = e->right( )->point( );
      }
      else
      {
        // std::cout << "pRight is null; should never get here..."
        //           << std::endl;
        double rightPoint_y = compute_y_at_x_2.approx(e->curve(), clipRect.right());
        rightPoint = Point_2(clipRect.right(), rightPoint_y);
      }

      X_monotone_curve_2 curve =
        this->construct_x_monotone_subcurve_2(e->curve(),
                                              leftPoint, rightPoint);
      envelopeToUpdate->insert( curve );

    }
  }

}

#if 0
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
#endif
// template < typename Arr_, typename Traits >
// template < typename Coefficient_ >
// void EnvelopeCallback< Arr_, Traits >::
// updateEnvelope(bool  lower ,
//                CGAL::Arr_algebraic_segment_traits_2<Coefficient_> /* traits */)
// {
//   // std::cout << "alg seg envelope stub" << std::endl;
// }

template < typename Arr_, typename Traits >
void EnvelopeCallback< Arr_, Traits >::showLowerEnvelope( bool show )
{
  this->showLower = show;
  if ( showLower )
  {
    this->scene->addItem( this->lowerEnvelope );
    this->updateEnvelope( true );
  }
  else
  {
    this->scene->removeItem( this->lowerEnvelope );
  }
}

template < typename Arr_, typename Traits >
void EnvelopeCallback< Arr_, Traits >::showUpperEnvelope( bool show )
{
  this->showUpper = show;
  if ( showUpper )
  {
    this->scene->addItem( this->upperEnvelope );
    this->updateEnvelope( false );
  }
  else
  {
    this->scene->removeItem( this->upperEnvelope );
  }
}

template class EnvelopeCallback<Seg_arr, Seg_traits>;
template class EnvelopeCallback<Pol_arr, Pol_traits>;
template class EnvelopeCallback<Conic_arr, Conic_traits>;
template class EnvelopeCallback<Lin_arr, Lin_traits>;
template class EnvelopeCallback<Arc_arr, Arc_traits>;
template class EnvelopeCallback<Alg_seg_arr, Alg_seg_traits>;
