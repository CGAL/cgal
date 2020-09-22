// Copyright (c) 2018 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
#include <CGAL/Three/Viewer_interface.h>
#include <QMainWindow>
#include <QApplication>

#ifdef three_EXPORTS
#  define THREE_EXPORT Q_DECL_EXPORT
#else
#  define THREE_EXPORT Q_DECL_IMPORT
#endif

#if QT_VERSION < QT_VERSION_CHECK(5, 15, 0)
#define CGAL_QT_SKIP_EMPTY_PARTS QString::SkipEmptyParts
#else
#define CGAL_QT_SKIP_EMPTY_PARTS ::Qt::SkipEmptyParts
#endif

namespace CGAL{
namespace Three{
//define enum depending on Qt version
class Polyhedron_demo_plugin_interface;
class THREE_EXPORT Three{
public:

  Three();
  virtual ~Three(){}
  static QMainWindow* mainWindow();
  static Viewer_interface* mainViewer();
  static Viewer_interface* currentViewer();
  static void setCurrentViewer(CGAL::Three::Viewer_interface* viewer);
  static Viewer_interface* activeViewer();
  static Scene_interface* scene();
  static QObject* connectableScene();
  static RenderingMode defaultSurfaceMeshRenderingMode();
  static RenderingMode defaultPointSetRenderingMode();
  static QString modeName(RenderingMode mode);
  static RenderingMode modeFromName(QString name);
  static int getDefaultPointSize();
  static int getDefaultNormalLength();
  static int getDefaultLinesWidth();
  /*! \brief Adds a dock widget to the interface
   *
   * Adds a dock widget in the left section of the MainWindow. If the slot is already
   * taken, the dock widgets will be tabified.
   */
  void addDockWidget(QDockWidget* dock_widget);

  /*! \brief Gets an item of the templated type.
   * \returns the first `SceneType` item found in the scene's list of currently selected
   * items;
   * \returns nullptr if there is no `SceneType` in the list.
   */
  template<class SceneType>
  static SceneType* getSelectedItem();

  /*! \brief Automatically connects each action of `plugin` to the corresponding slot.
   *
   * \attention Each action named `ActionName` in the plugin's `actions()` list must have
   *  a corresponding slot named `on_ActionsName_triggered()`
   * in the plugin.
   */
  static void autoConnectActions(CGAL::Three::Polyhedron_demo_plugin_interface* plugin);
  /*!
   * Displays in the console a blue text preceded by the mention
   * "INFO: ".
   */
  static void information(QString);
  /*!
   * Displays in the console an orange text preceded by the mention "WARNING: ".
   */
  static void warning(QString);
  /*!
   * Displays in the console a red text preceded by the mention "ERROR: ".
   */
  static void error(QString);
  /*!
   * Displays an information popup.
   */
  static void information(QString title, QString message);
  /*!
   * Displays a warning popup.
   */
  static void warning(QString title, QString message);
  /*!
   * Displays an error popup.
   */
  static void error(QString title, QString message);
protected:
  static QMainWindow* s_mainwindow;
  static Viewer_interface* s_mainviewer;
  static Viewer_interface* s_currentviewer;
  static Scene_interface* s_scene;
  static QObject* s_connectable_scene;
  static Three* s_three;
  static RenderingMode s_defaultSMRM;
  static RenderingMode s_defaultPSRM;
  static int default_point_size;
  static int default_normal_length;
  static int default_lines_width;

public:
  struct CursorScopeGuard
  {
    CursorScopeGuard(QCursor cursor)
    {
      QApplication::setOverrideCursor(cursor);
    }
    ~CursorScopeGuard()
    {
      QApplication::restoreOverrideCursor();
    }
  };
};
}
}

#endif // THREE_H
