// Copyright (c) 2009,2010,2012  GeometryFactory Sarl (France)
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
//
//
// Author(s)     : Laurent RINEAU
//! \file Polyhedron_demo_io_plugin_interface.h
#ifndef POLYHEDRON_DEMO_IO_PLUGIN_INTERFACE_H
#define POLYHEDRON_DEMO_IO_PLUGIN_INTERFACE_H

#include <CGAL/license/Three.h>


#include <QFileInfo>
#include <QStringList>

namespace CGAL{
namespace Three {
class Scene_item;
  /*!
   * This class provides a base for creating a new IO plugin.
   */
class Polyhedron_demo_io_plugin_interface 
{
public:
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
  virtual bool canLoad() const = 0;
  //!  Loads an item from a file.
  //! This must be overriden.
  virtual Scene_item* load(QFileInfo fileinfo) = 0;
  //!Specifies if the io_plugin can save the item or not.
  //!This must be overriden.
  virtual bool canSave(const Scene_item*) = 0;
  //!Saves the item in the file corresponding to the path
  //!contained in fileinfo. Returns false if error.
  //! This must be overriden.
  virtual bool save(const Scene_item*, QFileInfo fileinfo) = 0;
};
}
}
Q_DECLARE_INTERFACE(CGAL::Three::Polyhedron_demo_io_plugin_interface,
                    "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0")

#endif // POLYHEDRON_DEMO_IO_PLUGIN_INTERFACE_H
