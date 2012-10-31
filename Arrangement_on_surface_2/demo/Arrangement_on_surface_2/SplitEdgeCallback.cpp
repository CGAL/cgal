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
