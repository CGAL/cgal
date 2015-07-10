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

#include "ColorItemEditor.h"
#include "DeleteCurveModeItemEditor.h"
#include "PropertyValueDelegate.h"
#include "DeleteCurveMode.h"

#include <iostream>
#include <QItemEditorFactory>


PropertyValueDelegate::PropertyValueDelegate( QObject* parent ):
  QItemDelegate( parent )
{ 
  QItemEditorFactory* factory = new QItemEditorFactory;
  QItemEditorCreatorBase* creator =
    new QStandardItemEditorCreator< PositiveSpinBox >( );
  factory->registerEditor( QVariant::UInt, creator );
  this->setItemEditorFactory( factory );
}

QWidget* PropertyValueDelegate::
createEditor( QWidget* parent, const QStyleOptionViewItem& option,
              const QModelIndex& index ) const
{
  QWidget* editor;
  QVariant myData = index.data( Qt::UserRole );

  // check for data types we need to handle ourselves
 /*
  if ( qVariantCanConvert< QColor >( myData ) )
  {
    ColorItemEditor* colorEditor = new ColorItemEditor( parent );
    editor = colorEditor;

    QObject::connect( colorEditor, SIGNAL(confirmed()), this, SLOT(commit()));
  }
  else if ( qVariantCanConvert< DeleteCurveMode >( myData ) )
  {
    DeleteCurveModeItemEditor* modeEditor =
      new DeleteCurveModeItemEditor( parent );
    modeEditor->setMode( qVariantValue< DeleteCurveMode >( myData ) );
    editor = modeEditor;

    QObject::connect( modeEditor, SIGNAL( currentIndexChanged( int ) ), this,
                      SLOT( commit( ) ) );
  }
  else
  { // default handler
    editor = QItemDelegate::createEditor( parent, option, index );
  }*/

 if ( myData.canConvert< QColor >() )
  {
    ColorItemEditor* colorEditor = new ColorItemEditor( parent );
    editor = colorEditor;

    QObject::connect( colorEditor, SIGNAL(confirmed()), this, SLOT(commit()));
  }
  else if ( myData.canConvert< DeleteCurveMode >() )
  {
    DeleteCurveModeItemEditor* modeEditor =
      new DeleteCurveModeItemEditor( parent );
    modeEditor->setMode( myData.value< DeleteCurveMode >() );
    editor = modeEditor;

    QObject::connect( modeEditor, SIGNAL( currentIndexChanged( int ) ), this,
                      SLOT( commit( ) ) );
  }
  else
  { // default handler
    editor = QItemDelegate::createEditor( parent, option, index );
  }

  return editor;
}

void PropertyValueDelegate::setModelData( QWidget* editor,
                                          QAbstractItemModel* model,
                                          const QModelIndex& index ) const
{
  ColorItemEditor* colorEditor = qobject_cast<ColorItemEditor*>(editor);
  if (colorEditor)
  {
    // std::cout << "set color model data" << std::endl;
    model->setData(index, colorEditor->color(), Qt::DisplayRole);
    model->setData(index, colorEditor->color(), Qt::DecorationRole);
    model->setData(index, QVariant::fromValue(colorEditor->color()),
                    Qt::UserRole);
    return;
  }
  DeleteCurveModeItemEditor* modeEditor =
    qobject_cast<DeleteCurveModeItemEditor*>(editor);
  if (modeEditor) {
    model->setData(index, DeleteCurveMode::ToString(modeEditor->mode()),
                   Qt::DisplayRole);
    model->setData(index, QVariant::fromValue(modeEditor->mode()),
                   Qt::UserRole);
    return;
  }
  QItemDelegate::setModelData(editor, model, index);
}

bool PropertyValueDelegate::eventFilter( QObject* object, QEvent* event )
{
  QWidget* editor = qobject_cast<QWidget*>(object);
  if ((event->type() == QEvent::FocusOut) ||
      (event->type() == QEvent::Hide && editor->isWindow()))
  {
    ColorItemEditor* colorEditor = qobject_cast<ColorItemEditor*>(editor);
    if (colorEditor)
      return false;
    DeleteCurveModeItemEditor* modeEditor =
      qobject_cast<DeleteCurveModeItemEditor*>(editor);
    if (modeEditor)
      return false;
  }
  return QItemDelegate::eventFilter( object, event );
}

void PropertyValueDelegate::commit( )
{
  // std::cout << "commit selection" << std::endl;
  QWidget* editor = qobject_cast< QWidget* >( sender( ) );
  if ( editor )
  {
    emit( commitData( editor ) );
    emit( closeEditor( editor ) );
  }
}

PositiveSpinBox::PositiveSpinBox( QWidget* parent ) :
  QSpinBox( parent )
{
  this->setMinimum( 1 );
}

void PositiveSpinBox::setValue( unsigned int val )
{
  QSpinBox::setValue( val );
}

unsigned int PositiveSpinBox::value( ) const
{
  return QSpinBox::value( );
}
