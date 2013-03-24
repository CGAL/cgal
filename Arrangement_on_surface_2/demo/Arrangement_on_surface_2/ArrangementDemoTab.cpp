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

#include "ArrangementDemoTab.h"

#include <QGridLayout>

ArrangementDemoTabBase::ArrangementDemoTabBase( QWidget* parent ) :
  QWidget( parent ),
  graphicsView( new ArrangementDemoGraphicsView( this ) ),
  scene( new QGraphicsScene( -100, -100, 100, 100 ) ),
  layout( new QGridLayout( this ) ),
  arrangementGraphicsItem( NULL ),
  curveInputCallback( NULL ),
  deleteCurveCallback( NULL ),
  pointLocationCallback( NULL ),
  verticalRayShootCallback( NULL ),
  mergeEdgeCallback( NULL ),
  splitEdgeCallback( NULL ),
  envelopeCallback( NULL ),
  fillFaceCallback( NULL )
{
  this->setupUi( );
}

ArrangementDemoTabBase::~ArrangementDemoTabBase( )
{ 
}

void ArrangementDemoTabBase::setupUi( )
{
  this->layout->addWidget( this->graphicsView, 0, 0, 1, 1 );
  this->graphicsView->setScene( this->scene );
  this->graphicsView->setMouseTracking( true );
}

QGraphicsScene* ArrangementDemoTabBase::getScene( ) const
{
  return this->scene;
}

ArrangementDemoGraphicsView* ArrangementDemoTabBase::getView( ) const
{
  return this->graphicsView;
}

CGAL::Qt::ArrangementGraphicsItemBase*
ArrangementDemoTabBase::getArrangementGraphicsItem( ) const
{
  return this->arrangementGraphicsItem;
}

CGAL::Qt::GraphicsViewCurveInputBase*
ArrangementDemoTabBase::getCurveInputCallback( ) const
{
  return this->curveInputCallback;
}

CGAL::Qt::Callback* ArrangementDemoTabBase::getDeleteCurveCallback( ) const
{
  return this->deleteCurveCallback;
}

CGAL::Qt::Callback* ArrangementDemoTabBase::getPointLocationCallback( ) const
{
  return this->pointLocationCallback;
}

VerticalRayShootCallbackBase*
ArrangementDemoTabBase::getVerticalRayShootCallback( ) const
{
  return this->verticalRayShootCallback;
}

CGAL::Qt::Callback* ArrangementDemoTabBase::getMergeEdgeCallback( ) const
{
  return this->mergeEdgeCallback;
}

SplitEdgeCallbackBase* ArrangementDemoTabBase::getSplitEdgeCallback( ) const
{
  return this->splitEdgeCallback;
}

EnvelopeCallbackBase* ArrangementDemoTabBase::getEnvelopeCallback( ) const
{
  return this->envelopeCallback;
}

FillFaceCallbackBase* ArrangementDemoTabBase::getFillFaceCallback( ) const
{
  return this->fillFaceCallback;
}
