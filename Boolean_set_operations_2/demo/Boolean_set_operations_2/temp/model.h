/****************************************************************************
**
** Copyright (C) 2009 Nokia Corporation and/or its subsidiary(-ies).
** All rights reserved.
** Contact: Nokia Corporation (qt-info@nokia.com)
**
** This file is part of the documentation of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial Usage
** Licensees holding valid Qt Commercial licenses may use this file in
** accordance with the Qt Commercial License Agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and Nokia.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 2.1 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL included in the
** packaging of this file.  Please review the following information to
** ensure the GNU Lesser General Public License version 2.1 requirements
** will be met: http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html.
**
** In addition, as a special exception, Nokia gives you certain additional
** rights.  These rights are described in the Nokia Qt LGPL Exception
** version 1.1, included in the file LGPL_EXCEPTION.txt in this package.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3.0 as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL included in the
** packaging of this file.  Please review the following information to
** ensure the GNU General Public License version 3.0 requirements will be
** met: http://www.gnu.org/copyleft/gpl.html.
**
** If you have questions regarding the use of this file, please contact
** Nokia at qt-info@nokia.com.
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef MODEL_H
#define MODEL_H

#include <QAbstractListModel>
#include <QObject>
#include <QStringList>


class StringListModel : public QAbstractListModel
{
  Q_OBJECT

public:
  StringListModel(const QStringList &strings, QObject *parent = 0)
    : QAbstractListModel(parent), stringList(strings) {}

  int rowCount(const QModelIndex &parent = QModelIndex()) const;
  QVariant data(const QModelIndex &index, int role) const;
  QVariant headerData(int section, Qt::Orientation orientation,

    int role = Qt::DisplayRole) const;



  Qt::ItemFlags flags(const QModelIndex &index) const;
  bool setData(const QModelIndex &index, const QVariant &value,

    int role = Qt::EditRole);



  bool insertRows(int position, int rows, const QModelIndex &index = QModelIndex());
  bool removeRows(int position, int rows, const QModelIndex &index = QModelIndex());



private:
  QStringList stringList;
};


#endif

