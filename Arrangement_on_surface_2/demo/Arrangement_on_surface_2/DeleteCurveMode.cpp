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

#include "DeleteCurveMode.h"

#include <QString>

DeleteCurveMode::DeleteCurveMode( ) : m_mode( DELETE_CURVE ) { }

DeleteCurveMode::DeleteCurveMode( const DeleteCurveMode& dcm ) :
  m_mode( dcm.mode( ) )
{ }

DeleteCurveMode::DeleteCurveMode( Mode mode ) : m_mode( mode ) { }

DeleteCurveMode::~DeleteCurveMode( ) { }

DeleteCurveMode::Mode DeleteCurveMode::mode( ) const
{
  return this->m_mode;
}

void DeleteCurveMode::setMode( Mode mode )
{
  this->m_mode = mode;
}

QString DeleteCurveMode::ToString( const DeleteCurveMode& mode )
{
  return ( mode.mode( ) == DELETE_CURVE ) ?
    QString("Delete Curve") : QString("Delete Edge");
}
