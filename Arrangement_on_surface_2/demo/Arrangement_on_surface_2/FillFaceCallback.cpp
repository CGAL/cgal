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

#include "FillFaceCallback.h"

FillFaceCallbackBase::FillFaceCallbackBase( QObject* parent ) :
  CGAL::Qt::Callback( parent ),
  fillColor( ::Qt::black )
{ }

void FillFaceCallbackBase::setColor( QColor c )
{
  this->fillColor = c;
  Q_EMIT modelChanged( );
}

QColor FillFaceCallbackBase::getColor( ) const
{
  return this->fillColor;
}
