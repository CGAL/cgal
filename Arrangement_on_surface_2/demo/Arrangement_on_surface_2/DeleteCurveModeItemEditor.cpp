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
  return qVariantValue<DeleteCurveMode >(this->itemData(this->currentIndex( ),
                                                        Qt::UserRole ) );
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
