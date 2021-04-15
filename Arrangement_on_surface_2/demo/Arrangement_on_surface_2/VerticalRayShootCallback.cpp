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

#include "VerticalRayShootCallback.h"

VerticalRayShootCallbackBase::VerticalRayShootCallbackBase(QObject* parent_) :
  CGAL::Qt::Callback( parent_ ),
  shootingUp( true )
{ }

void VerticalRayShootCallbackBase::setShootingUp( bool isShootingUp )
{
  this->shootingUp = isShootingUp;
}

