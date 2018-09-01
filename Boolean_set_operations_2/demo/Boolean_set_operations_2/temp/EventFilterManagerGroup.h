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

#ifndef CGAL_EVENTFILTERMANAGERGROUP_H
#define CGAL_EVENTFILTERMANAGERGROUP_H

#include <QWidget>
#include <QMap>
#include <QDebug>

#include "EventFilterManager.h"

class EventFilterManagerGroup {
public:
  EventFilterManagerGroup();
  ~EventFilterManagerGroup();

  void addObjectToWatch(const QString name, QObject* object);
  void removeObjectFromWatch(const QString name);
  void addFilterWidget(const QString targetName, const QString name,
                       QObject* filter);
  QObject* removeFilterWidget(const QString name);
  QObject* getFilterWidget(const QString name);
  QString const getCurrentFilterName() const;
  void setCurrentFilterName(const QString name);
  QList<QString> getFilterNames();
  void clear();

private:
  QMap<QString, EventFilterManager*> m_objectsMap;
  QMap<QString, QString> m_filtersMap;
  QString m_currentFilter;
  QString m_currentWatchedObject;
};

#endif // CGAL_EVENTFILTERMANAGERGROUP_H
