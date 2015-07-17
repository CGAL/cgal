// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : 
//******************************************************************************

//#include <CGAL_demo/Plugin_interface.h>
//#include <CGAL_demo/Plugin_helper.h>
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_io_plugin_interface.h"

#include "implicit_functions/Implicit_function_interface.h"
#include "Scene_implicit_function_item.h"
#include "ui_Function_dialog.h"

#include <iostream>

#include <QAction>
#include <QMainWindow>
#include <QPluginLoader>
#include <QDir>
#include <QApplication>
#include <QMenu>
#include <QList>
#include <QLibrary>

class Io_implicit_function_plugin :
  public QObject, 
  // public Polyhedron_demo_plugin_interface,
  protected Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  Io_implicit_function_plugin();
  virtual ~Io_implicit_function_plugin() {}

  bool applicable(QAction*) const { return true; }

  QString name() const { return "implicit functions"; }
  // QString nameFilters() const { return ""; }
  // bool canLoad() const { return false; }
  

  typedef Polyhedron_demo_plugin_helper Plugin_helper;
  
  using Plugin_helper::init;
  virtual void init(QMainWindow* mainWindow, Scene_interface* scene_interface);
  
  QList<QAction*> actions() const
  {
    return QList<QAction*>();
  }
  
public Q_SLOTS:
  void load_function() const;
  
private:
  void load_function_plugins();
  
private:
  QList<Implicit_function_interface*> functions_;
};



Io_implicit_function_plugin::
Io_implicit_function_plugin()
{
  load_function_plugins();
}


void
Io_implicit_function_plugin::
init(QMainWindow* mainWindow, Scene_interface* scene_interface)
{
  this->scene = scene_interface;
  this->mw = mainWindow;
  
  QAction* actionLoadFunction = new QAction("Load &implicit function", mw);
  if( NULL != actionLoadFunction )
  {
    connect(actionLoadFunction, SIGNAL(triggered()), this, SLOT(load_function()));
  }
  
  QMenu* menuFile = mw->findChild<QMenu*>("menuFile");
  if ( NULL != menuFile )
  {
    QList<QAction*> menuFileActions = menuFile->actions();
    
    // Look for action just after "Load..." action
    QAction* actionAfterLoad = NULL;
    for ( QList<QAction*>::iterator it_action = menuFileActions.begin(), 
         end = menuFileActions.end() ; it_action != end ; ++ it_action ) //Q_FOREACH( QAction* action, menuFileActions)
    {
      if ( NULL != *it_action && (*it_action)->text().contains("Load") )
      {
        ++it_action;
        if ( it_action != end && NULL != *it_action )
        {
          actionAfterLoad = *it_action;
        }
      }
    }
    
    // Insert "Load implicit function" action
    if ( NULL != actionAfterLoad )
    {
      menuFile->insertAction(actionAfterLoad,actionLoadFunction);      
    }
    
  }
}


void
Io_implicit_function_plugin::
load_function() const
{
  QDialog dialog(mw);
  Ui::FunctionDialog ui;
  ui.setupUi(&dialog);
  
  connect(ui.buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
  connect(ui.buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));
  
  // Add loaded functions to the dialog
  int i=0;
  Q_FOREACH( Implicit_function_interface* f, functions_ )
  {
    ui.functionList->insertItem(i++,f->name());
  }

  // Open window
  int return_code = dialog.exec();
  if(return_code == QDialog::Rejected) { return; }
  
  // Get selected function
  i = ui.functionList->currentIndex();
  Implicit_function_interface* function = functions_[i];
  
  // Create Scene_implicit_function object and add it to the framework
  Scene_implicit_function_item* item =
    new Scene_implicit_function_item(function);
  
  item->setName(tr("%1").arg(function->name()));
  item->setRenderingMode(FlatPlusEdges);

  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  scene->itemChanged(index);
    
  Scene_interface::Item_id new_item_id = scene->addItem(item);
  scene->setSelectedItem(new_item_id);
}

void
Io_implicit_function_plugin::
load_function_plugins()
{
  QDir pluginsDir(qApp->applicationDirPath());
  QString dirname = pluginsDir.dirName();
  if ( !pluginsDir.cd("implicit_functions") ) { 
    // In that case, dirname may be "Debug" or "Release" and one has to
    // search in ../implicit_functions/Debug or
    // ../implicit_functions/Release
    QString newDir = QString("../implicit_functions/") + dirname;
    if( !pluginsDir.cd(newDir) ) return; 
  }
  
  Q_FOREACH (QString fileName, pluginsDir.entryList(QDir::Files))
  {
    if ( fileName.contains("plugin") && QLibrary::isLibrary(fileName) )
    {
      qDebug("    + Loading Function \"%s\"...", fileName.toUtf8().data());
      QPluginLoader loader;
      loader.setFileName(pluginsDir.absoluteFilePath(fileName));
      QObject *function_plugin = loader.instance();
      if ( NULL != function_plugin )
      {
        Implicit_function_interface* function =
          qobject_cast<Implicit_function_interface*>(function_plugin);
        
        if ( NULL != function )
        {
          functions_ << function;
        }
      }
    }
  }
}

#include "Io_implicit_function_plugin.moc"
