// Copyright (c) 2009-2015  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent RINEAU

//! \file Polyhedron_demo_plugin_interface.h
#ifndef POLYHEDRON_DEMO_PLUGIN_INTERFACE_H
#define POLYHEDRON_DEMO_PLUGIN_INTERFACE_H

#include <CGAL/license/Three.h>


#include <QString>
#include <QList>
#include <QtPlugin>
#include <QDebug>

class QAction;
class QMainWindow;
class Messages_interface;
namespace CGAL {
namespace Three {
class Scene_interface;
  /*!
   * This virtual class provides the basic functions used for making a plugin.
   */
class Polyhedron_demo_plugin_interface
{
public:
  //! \brief initializes the plugin
  //! This function acts like a constructor. This is where the attributes must be initialized.
  //! The Message_interface allows to print warnings or errors on the screen and the `Console` widget.
  virtual void init(QMainWindow* , CGAL::Three::Scene_interface* , Messages_interface*) = 0;

  //! \brief indicates if an action is usable or not.
  //! This function usually tests the type of the selected item to determine if `action` can be applied to it,
  //! but not necessarly.
  //! @returns \c true if `action` can be called in the current state, \c false
  //! otherwise
  virtual bool applicable(QAction* action) const = 0;
  //!Contains all the plugin's actions.
  virtual QList<QAction*> actions() const = 0;
  //!\brief is called when the application is closed.
  //! Override this function if you need to perform a specific action
  //! when the application is closed, like hide the widgets if you don't want
  //! their visibility to be saved.
  virtual void closure() {
 }
protected :
};
}
}
Q_DECLARE_INTERFACE(CGAL::Three::Polyhedron_demo_plugin_interface,
                    "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

#endif // POLYHEDRON_DEMO_PLUGIN_INTERFACE_H
