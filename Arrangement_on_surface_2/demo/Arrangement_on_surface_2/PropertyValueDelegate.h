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

#ifndef PROPERTY_VALUE_DELEGATE_H
#define PROPERTY_VALUE_DELEGATE_H

#include <QtGui>
#include <QItemDelegate>
#include <QSpinBox>


class PropertyValueDelegate : public QItemDelegate
{
  Q_OBJECT

  public:
  PropertyValueDelegate( QObject* parent = nullptr );

public:
  QWidget* createEditor( QWidget* parent, const QStyleOptionViewItem& option,
                         const QModelIndex& index ) const;
  void setModelData( QWidget* editor, QAbstractItemModel* model,
                     const QModelIndex& index ) const;
  bool eventFilter( QObject* object, QEvent* event );

public Q_SLOTS:
  void commit( );

};

class PositiveSpinBox : public QSpinBox
{
  Q_OBJECT
  Q_PROPERTY( unsigned int value READ value WRITE setValue USER true )

    public:
    PositiveSpinBox( QWidget* parent );
  void setValue( unsigned int );
  unsigned int value( ) const;
};

#endif // PROPERTY_VALUE_DELEGATE_H
