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

#ifndef PROPERTY_VALUE_DELEGATE_H
#define PROPERTY_VALUE_DELEGATE_H

#include <QtGui>

class PropertyValueDelegate : public QItemDelegate
{
  Q_OBJECT

  public:
  PropertyValueDelegate( QObject* parent = 0 );

public:
  QWidget* createEditor( QWidget* parent, const QStyleOptionViewItem& option,
                         const QModelIndex& index ) const;
  void setModelData( QWidget* editor, QAbstractItemModel* model,
                     const QModelIndex& index ) const;
  bool eventFilter( QObject* object, QEvent* event );

public slots:
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
