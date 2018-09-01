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
// Author(s) : Saar Katz <kats.saar@gmail.com>

#ifndef CGAL_POLYGONTABLEMODEL_H
#define CGAL_POLYGONTABLEMODEL_H

#include <QAbstractTableModel>

#include "PolygonListItem.h"

class PolygonTableModel : public QAbstractTableModel {
  Q_OBJECT
public:
  PolygonTableModel(QObject *parent = 0);

  int rowCount(const QModelIndex &parent = QModelIndex()) const;
  int columnCount(const QModelIndex &parent = QModelIndex()) const;

  QVariant data(const QModelIndex &index, int role) const;
  QVariant headerData(int section, Qt::Orientation orientation, 
    int role = Qt::DisplayRole) const;

  Qt::ItemFlags flags(const QModelIndex &index) const;
  bool setData(const QModelIndex &index, const QVariant &value,
    int role = Qt::EditRole);

  bool insertRows(int position, int rows, 
                  const QModelIndex &index = QModelIndex());
  bool removeRows(int position, int rows, 
                  const QModelIndex &index = QModelIndex());
  bool insertItem(int position, PolygonListItem* item);
  bool appendItem(PolygonListItem* item);

  PolygonListItem* getPolygonItem(int index);

private:
  QList<PolygonListItem*> m_polygonList;
};

#endif // CGAL_POLYGONTABLEMODEL_H
