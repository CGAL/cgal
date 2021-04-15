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

#include "DeleteCurveModeItemEditor.h"

DeleteCurveModeItemEditor::DeleteCurveModeItemEditor( QWidget* parent ) :
  QComboBox( parent )
{
  this->setFrame( false );

  QVariant deleteCurveOption = QVariant::fromValue( DeleteCurveMode( ) );
  QVariant deleteEdgeOption =
    QVariant::fromValue( DeleteCurveMode( DeleteCurveMode::DELETE_EDGE ) );
  this->insertItem( 0, "Delete Curve", deleteCurveOption );
  this->insertItem( 1, "Delete Edge", deleteEdgeOption );
}

DeleteCurveMode DeleteCurveModeItemEditor::mode( ) const
{
  return this->itemData(this->currentIndex( ), Qt::UserRole ).value< DeleteCurveMode >();
}

void DeleteCurveModeItemEditor::setMode( DeleteCurveMode m )
{
  if ( m.mode( ) == DeleteCurveMode::DELETE_CURVE )
  {
    this->setCurrentIndex( 0 );
  }
  else
  {
    this->setCurrentIndex( 1 );
  }
}
