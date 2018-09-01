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

#ifndef CGAL_EVENTFILTERMANAGER_H
#define CGAL_EVENTFILTERMANAGER_H

#include <QWidget>
#include <QMap>
#include <QDebug>

class EventFilterManager : public QObject {
  Q_OBJECT

public:
  EventFilterManager(QObject* parent = 0);
  ~EventFilterManager();

  void addFilterWidget(const QString name, QObject* filter);
  QObject* getFilterWidget(const QString name);
  QObject* removeFilterWidget(const QString name);
  QString const getCurrentFilterName() const;
  void setCurrentFilterName(const QString name);
  QList<QString> getFilterNames();
  void clear();

protected:
  bool eventFilter(QObject* target, QEvent* event);

private:
  QMap<QString, QObject*> m_filterMap;
  QString m_currentFilter;
};

#endif // CGAL_EVENTFILTERMANAGER_H
