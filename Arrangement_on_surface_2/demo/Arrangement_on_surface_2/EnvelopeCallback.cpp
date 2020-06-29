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
#include "CurveGraphicsItem.h"
#include "Utils.h"

#include <CGAL/envelope_2.h>
#include <CGAL/Envelope_diagram_1.h>

#include <vector>

EnvelopeCallbackBase::EnvelopeCallbackBase( QObject* parent ) :
  CGAL::Qt::Callback( parent )
{ }

template < typename Arr_>
EnvelopeCallback<Arr_>::EnvelopeCallback(Arrangement* arr_,
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

template < typename Arr_ >
void EnvelopeCallback<Arr_>::setEnvelopeEdgeColor( const QColor& color )
{
  this->lowerEnvelope->setEdgeColor( color );
  this->upperEnvelope->setEdgeColor( color );
}

template < typename Arr_ >
const QColor& EnvelopeCallback<Arr_>::getEnvelopeEdgeColor( ) const
{
  return this->lowerEnvelope->edgeColor( );
}

template < typename Arr_ >
void EnvelopeCallback<Arr_>::setEnvelopeEdgeWidth( int width )
{
  this->lowerEnvelope->setEdgeWidth( width );
  this->upperEnvelope->setEdgeWidth( width );
}

template < typename Arr_ >
int EnvelopeCallback<Arr_>::getEnvelopeEdgeWidth( ) const
{
  return this->lowerEnvelope->edgeWidth( );
}

template < typename Arr_ >
void EnvelopeCallback<Arr_>::setEnvelopeVertexColor( const QColor& color )
{
  this->lowerEnvelope->setVertexColor( color );
  this->upperEnvelope->setVertexColor( color );
}

template < typename Arr_ >
const QColor& EnvelopeCallback<Arr_>::getEnvelopeVertexColor( ) const
{
  return this->lowerEnvelope->vertexColor( );
}

template < typename Arr_ >
void EnvelopeCallback<Arr_>::setEnvelopeVertexRadius( int radius )
{
  this->lowerEnvelope->setVertexRadius( radius );
  this->upperEnvelope->setVertexRadius( radius );
}

template < typename Arr_ >
int EnvelopeCallback<Arr_>::getEnvelopeVertexRadius( ) const
{
  return this->lowerEnvelope->vertexRadius( );
}

template < typename Arr_ >
void EnvelopeCallback<Arr_>::setScene( QGraphicsScene* scene_ )
{
  EnvelopeCallbackBase::setScene(scene_);
  lowerEnvelope->setScene(scene_);
  upperEnvelope->setScene(scene_);
  this->scene->addItem( this->lowerEnvelope );
  this->scene->addItem( this->upperEnvelope );
}

template < typename Arr_ >
void EnvelopeCallback< Arr_>::slotModelChanged( )
{
  if (showLower) this->updateEnvelope(true);
  if (showUpper) this->updateEnvelope(false);
}

template < typename Arr_ >
void EnvelopeCallback< Arr_>::updateEnvelope( bool lower )
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
        leftPoint = e->left()->point();

      if (e->right())
        rightPoint = e->right()->point();

      Construct_x_monotone_subcurve_2<Traits> construct_x_monotone_subcurve_2;
      X_monotone_curve_2 curve =
        construct_x_monotone_subcurve_2(e->curve(), leftPoint, rightPoint);

      envelopeToUpdate->insert(curve);
      // TODO: visually show leftPoint and rightPoint
      // note: leftPoint and rightPoint are not the actual curve points
    }
  }

  envelopeToUpdate->modelChanged();
}

template < typename Arr_ >
void EnvelopeCallback< Arr_>::showLowerEnvelope( bool show )
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

template < typename Arr_ >
void EnvelopeCallback< Arr_>::showUpperEnvelope( bool show )
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

template class EnvelopeCallback<Seg_arr>;
template class EnvelopeCallback<Pol_arr>;
template class EnvelopeCallback<Conic_arr>;
template class EnvelopeCallback<Lin_arr>;
template class EnvelopeCallback<Alg_seg_arr>;
template class EnvelopeCallback<Bezier_arr>;
