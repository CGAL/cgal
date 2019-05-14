// Copyright (c) 2018 GeometryFactory Sarl (France)
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
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Maxime GIMENO

#ifndef THREE_H
#define THREE_H

#include <CGAL/license/Three.h>

#include <QString>
#include <QObject>
#include <QDockWidget>
#include <CGAL/Three/Scene_interface.h>
#include <QMainWindow>

#ifdef three_EXPORTS
#  define THREE_EXPORT Q_DECL_EXPORT
#else
#  define THREE_EXPORT Q_DECL_IMPORT
#endif

namespace CGAL{
namespace Three{
class Polyhedron_demo_plugin_interface;
class THREE_EXPORT Three{
public:

  Three();
  virtual ~Three(){}
  static QMainWindow* mainWindow();
  static Scene_interface* scene();
  static QObject* connectableScene();
  static Three* messages();
  /*! \brief Adds a dock widget to the interface
   *
   * Adds a dock widget in the left section of the MainWindow. If the slot is already taken, the dock widgets will be tabified.
   */
  void addDockWidget(QDockWidget* dock_widget);

  /*! \brief Gets an item of the templated type.
   * \returns the first `SceneType` item found in the scene's list of currently selected items;
   * \returns NULL if there is no `SceneType` in the list.
   */
  template<class SceneType>
  static SceneType* getSelectedItem();

  /*! \brief Automatically connects each action of `plugin` to the corresponding slot.
   *
   * \attention Each action named `ActionName` in the plugin's `actions()` list must have a corresponding slot named `on_ActionsName_triggered()`
   * in the plugin.
   */
  static void autoConnectActions(CGAL::Three::Polyhedron_demo_plugin_interface* plugin);
protected:
  static QMainWindow* s_mainwindow;
  static Scene_interface* s_scene;
  static QObject* s_connectable_scene;
  static Three* s_three;

};
}
}

#endif // THREE_H
