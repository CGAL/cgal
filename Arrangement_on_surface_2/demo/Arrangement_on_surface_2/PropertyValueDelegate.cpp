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

#include "ColorItemEditor.h"
#include "PropertyValueDelegate.h"

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

 if ( myData.canConvert< QColor >() )
  {
    ColorItemEditor* colorEditor = new ColorItemEditor( parent );
    editor = colorEditor;

    QObject::connect( colorEditor, SIGNAL(confirmed()), this, SLOT(commit()));
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
  }
  return QItemDelegate::eventFilter( object, event );
}

void PropertyValueDelegate::commit( )
{
  QWidget* editor = qobject_cast< QWidget* >( sender( ) );
  if ( editor )
  {
    Q_EMIT( commitData( editor ) );
    Q_EMIT( closeEditor( editor ) );
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
