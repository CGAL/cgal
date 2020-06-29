// Copyright (c) 2009,2010,2012  GeometryFactory Sarl (France)
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
//! \file Polyhedron_demo_io_plugin_interface.h
#ifndef POLYHEDRON_DEMO_IO_PLUGIN_INTERFACE_H
#define POLYHEDRON_DEMO_IO_PLUGIN_INTERFACE_H

#include <CGAL/license/Three.h>


#include <QFileInfo>
#include <QStringList>
#include <QtPlugin>
class QMainWindow;
class Messages_interface;
namespace CGAL{
namespace Three {
class Scene_item;
class Scene_interface;
  /*!
   * This class provides a base for creating a new IO plugin.
   */
class Polyhedron_demo_io_plugin_interface
{
public:
  //! \brief Initializes the plugin
  //! This function is called in the constructor of the MainWindow.
  //! Whatever initialization the plugin needs can be done here. Default
  //! behavior is to do nothing.
  virtual void init(){}
  //!Returns the name of the plugin
  //!It is used by the loading system.
  virtual QString name() const = 0;
  virtual ~Polyhedron_demo_io_plugin_interface() {}
  /*! The filters for the names of the files that can be used
   * by the plugin.
   * Example : to filter OFF files : return "OFF files (*.off)"
*/
  virtual QString nameFilters() const = 0;
  //!Returns only the filters used for saving. The default is `nameFilters()`.
  //! You must override this function to change its behavior.
  virtual QString saveNameFilters() const {return nameFilters();}
  //!Returns only the filters used for loading. The default is `nameFilters()`.
  //! You must override this function to change its behavior.
  //! If multiple plugins have the same load filters, only once will be kept,
  //! so be careful not to use one that already exists.
  virtual QString loadNameFilters() const {return nameFilters();}

  //! Specifies if the io_plugin is able to load an item or not.
  //! This must be overriden.
  virtual bool canLoad(QFileInfo fileinfo) const = 0;
  //! Loads one or more item(s) from a file. `ok` is `true` if the loading
  //! was successful, `false` otherwise.
  //! New items will be added to the scene if `add_to_scene` is `true`.
  //! You don't want that when you reload an item, for example,
  //! as it will be added at some other point of the process.
  //! This must be overriden.
  virtual QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true) = 0;
  //!Specifies if the io_plugin can save the item or not.
  //!This must be overriden.
  virtual bool canSave(const Scene_item*) = 0;
  //!Saves one or more items in the file corresponding to the path
  //!contained in fileinfo. Returns false if error.
  //! This must be overriden.
  //! @attention When a file is successfully saved, it must be removed from the
  //! list.
  virtual bool save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& ) = 0;

  //! If this returns `true`, then the loader will be chosen as default in the
  //! list of available loaders when saving a file, which means it will be the
  //! first in the list.
  virtual bool isDefaultLoader(const Scene_item*) const { return false; }
  //! If this returns `true`, then the loader will be chosen as default in the
  //! list of available loaders when loading a file, which means it will be the
  //! first in the list.
  //! @param name is the extension without the dot (e.g. "off" for a .off file)
  virtual bool isDefaultLoader(const QString& name) const { Q_UNUSED(name); return false; }
};
}
}
Q_DECLARE_INTERFACE(CGAL::Three::Polyhedron_demo_io_plugin_interface,
                    "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.90")

#endif // POLYHEDRON_DEMO_IO_PLUGIN_INTERFACE_H
