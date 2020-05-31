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

#include <vector>

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
  this->lowerEnvelope->hide();
  this->upperEnvelope->hide();
}

template < typename Arr_, typename Traits >
void EnvelopeCallback<Arr_, Traits>::setEnvelopeEdgeColor( const QColor& color )
{
  this->lowerEnvelope->setEdgeColor( color );
  this->upperEnvelope->setEdgeColor( color );
}

template < typename Arr_, typename Traits >
const QColor& EnvelopeCallback<Arr_, Traits>::getEnvelopeEdgeColor( ) const
{
  return this->lowerEnvelope->edgeColor( );
}

template < typename Arr_, typename Traits >
void EnvelopeCallback<Arr_, Traits>::setEnvelopeEdgeWidth( int width )
{
  this->lowerEnvelope->setEdgeWidth( width );
  this->upperEnvelope->setEdgeWidth( width );
}

template < typename Arr_, typename Traits >
int EnvelopeCallback<Arr_, Traits>::getEnvelopeEdgeWidth( ) const
{
  return this->lowerEnvelope->edgeWidth( );
}

template < typename Arr_, typename Traits >
void EnvelopeCallback<Arr_, Traits>::setEnvelopeVertexColor( const QColor& color )
{
  this->lowerEnvelope->setVertexColor( color );
  this->upperEnvelope->setVertexColor( color );
}

template < typename Arr_, typename Traits >
const QColor& EnvelopeCallback<Arr_, Traits>::getEnvelopeVertexColor( ) const
{
  return this->lowerEnvelope->vertexColor( );
}

template < typename Arr_, typename Traits >
void EnvelopeCallback<Arr_, Traits>::setEnvelopeVertexRadius( int radius )
{
  this->lowerEnvelope->setVertexRadius( radius );
  this->upperEnvelope->setVertexRadius( radius );
}

template < typename Arr_, typename Traits >
int EnvelopeCallback<Arr_, Traits>::getEnvelopeVertexRadius( ) const
{
  return this->lowerEnvelope->vertexRadius( );
}

template < typename Arr_, typename Traits >
void EnvelopeCallback<Arr_, Traits>::setScene( QGraphicsScene* scene_ )
{
  EnvelopeCallbackBase::setScene(scene_);
  lowerEnvelope->setScene(scene_);
  upperEnvelope->setScene(scene_);
  this->scene->addItem( this->lowerEnvelope );
  this->scene->addItem( this->upperEnvelope );
}

template < typename Arr_, typename Traits >
void EnvelopeCallback< Arr_, Traits >::slotModelChanged( )
{
  if (showLower) this->updateEnvelope(true);
  if (showUpper) this->updateEnvelope(false);
}

template < typename Arr_, typename Traits >
void EnvelopeCallback< Arr_, Traits >::updateEnvelope( bool lower )
{
  std::vector<X_monotone_curve_2> curves;
  for (auto it = this->arr->edges_begin(); it != this->arr->edges_end(); ++it)
  {
    curves.push_back(it->curve());
  }

  typedef CGAL::Envelope_diagram_1<Traits> Diagram_1;
  typedef typename Diagram_1::Edge_const_handle Edge_const_handle;

  Diagram_1 diagram;
  CGAL::Qt::CurveGraphicsItem<Traits>* envelopeToUpdate;

  if (lower)
  {
    envelopeToUpdate = this->lowerEnvelope;
    CGAL::lower_envelope_x_monotone_2(curves.begin(), curves.end(), diagram);
  }
  else
  {
    envelopeToUpdate = this->upperEnvelope;
    CGAL::upper_envelope_x_monotone_2(curves.begin(), curves.end(), diagram);
  }
  envelopeToUpdate->clear( );

  auto next_edge = [](const auto& e) -> Edge_const_handle {
    auto&& v = e->right();
    if (v) return v->right();
    else   return nullptr;
  };

  for (Edge_const_handle e = diagram.leftmost(); e; e = next_edge(e))
  {
    if (!e->is_empty())
    {
      boost::optional<Point_2> leftPoint, rightPoint;
      if (e->left())
      {
        leftPoint = e->left()->point();
      }

      if (e->right())
      {
        rightPoint = e->right()->point();
      }

      X_monotone_curve_2 curve =
        this->construct_x_monotone_subcurve_2(e->curve(),
                                              leftPoint, rightPoint);
      envelopeToUpdate->insert(curve);
      if (leftPoint)  envelopeToUpdate->insert(*leftPoint);
      if (rightPoint) envelopeToUpdate->insert(*rightPoint);
    }
  }

  envelopeToUpdate->modelChanged();
}

template < typename Arr_, typename Traits >
void EnvelopeCallback< Arr_, Traits >::showLowerEnvelope( bool show )
{
  this->showLower = show;
  if (this->showLower)
  {
    this->updateEnvelope(true);
    this->lowerEnvelope->show();
  }
  else
  {
    this->lowerEnvelope->hide();
  }
}

template < typename Arr_, typename Traits >
void EnvelopeCallback< Arr_, Traits >::showUpperEnvelope( bool show )
{
  this->showUpper = show;
  if (this->showUpper)
  {
    this->updateEnvelope(false);
    this->upperEnvelope->show();
  }
  else
  {
    this->upperEnvelope->hide();
  }
}

template class EnvelopeCallback<Seg_arr, Seg_traits>;
template class EnvelopeCallback<Pol_arr, Pol_traits>;
template class EnvelopeCallback<Conic_arr, Conic_traits>;
template class EnvelopeCallback<Lin_arr, Lin_traits>;
template class EnvelopeCallback<Alg_seg_arr, Alg_seg_traits>;
