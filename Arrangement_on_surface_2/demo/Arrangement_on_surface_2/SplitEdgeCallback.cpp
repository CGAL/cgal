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

#include "SplitEdgeCallback.h"

void SplitEdgeCallbackBase::setSnappingEnabled( bool b )
{
  this->snappingEnabled = b;
}

void SplitEdgeCallbackBase::setSnapToGridEnabled( bool b )
{
  this->snapToGridEnabled = b;
}

void SplitEdgeCallbackBase::setColor( QColor c )
{
  this->color = c;
}

QColor SplitEdgeCallbackBase::getColor( ) const
{
  return this->color;
}

SplitEdgeCallbackBase::SplitEdgeCallbackBase( QObject* parent ) :
  CGAL::Qt::Callback( parent ),
  snappingEnabled( false ),
  snapToGridEnabled( false ),
  color( ::Qt::blue )
{ }
