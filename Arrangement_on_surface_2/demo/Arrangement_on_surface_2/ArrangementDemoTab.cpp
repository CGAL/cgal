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

#include "ArrangementDemoTab.h"

#include <QGridLayout>

ArrangementDemoTabBase::ArrangementDemoTabBase( QWidget* parent ) :
  QWidget( parent ),
  graphicsView( new ArrangementDemoGraphicsView( this ) ),
  // scene( new QGraphicsScene( -100, -100, 100, 100 ) ),
  // scene( new QGraphicsScene( 0, 0, 100, 100 ) ),
  scene( new QGraphicsScene( ) ),

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
  // int viewWidth = this->graphicsView->width();
  // int viewHeight = this->graphicsView->height();

  double viewWidth = 0.0001;
  double viewHeight = 0.0001;

  // this->scene->setSceneRect(-viewWidth/2, -viewHeight/2, viewWidth, viewHeight);
  this->scene->setSceneRect(0, 0, viewWidth, viewHeight);
  this->layout->addWidget( this->graphicsView, 0, 0 );
  this->graphicsView->setScene( this->scene );
  this->graphicsView->setMouseTracking( true );
}

QGraphicsScene* ArrangementDemoTabBase::getScene( ) const
{
  return this->scene;
}

//! Feedback from setting up a new tab based on various options
/*!
	\return call from the ArrangementDemoTab
*/
ArrangementDemoGraphicsView* ArrangementDemoTabBase::getView( ) const
{
  return this->graphicsView;
}

CGAL::Qt::ArrangementGraphicsItemBase*
ArrangementDemoTabBase::getArrangementGraphicsItem( ) const
{
  return this->arrangementGraphicsItem;
}

//! Determining the points of the arrangement
/*!
	\return call to the ArrangementCurveInputCallback
*/
CGAL::Qt::GraphicsViewCurveInputBase*
ArrangementDemoTabBase::getCurveInputCallback( ) const
{
  return this->curveInputCallback;
}

//! eraser option i.e. to delete the selected curve.
/*!
	\return the drawing after the selected curve has been removed
*/
CGAL::Qt::Callback* ArrangementDemoTabBase::getDeleteCurveCallback( ) const
{
  return this->deleteCurveCallback;
}

//! Returns the point where the mouse is selected.
/*!
	\return the point from the curve originates
*/
CGAL::Qt::Callback* ArrangementDemoTabBase::getPointLocationCallback( ) const
{
  return this->pointLocationCallback;
}

//! Vertical ray offshoot feedback
/*!
	\return the ray in the direction closest to the edge of the screen
*/
VerticalRayShootCallbackBase*
ArrangementDemoTabBase::getVerticalRayShootCallback( ) const
{
  return this->verticalRayShootCallback;
}

//! Merging the segments
/*!
	\return the curves after merging them back together
*/
CGAL::Qt::Callback* ArrangementDemoTabBase::getMergeEdgeCallback( ) const
{
  return this->mergeEdgeCallback;
}

//! Splitting the curves drawn in the screen with points.
/*!
	\return the points of splitting
*/
SplitEdgeCallbackBase* ArrangementDemoTabBase::getSplitEdgeCallback( ) const
{
  return this->splitEdgeCallback;
}

//! feedback after the envelope call.
/*!
	\return result of the envelope call
*/
EnvelopeCallbackBase* ArrangementDemoTabBase::getEnvelopeCallback( ) const
{
  return this->envelopeCallback;
}

//! member function to fill the viewport
/*!
	\return result after calling the fill color option
*/
FillFaceCallbackBase* ArrangementDemoTabBase::getFillFaceCallback( ) const
{
  return this->fillFaceCallback;
}
