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
