// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//******************************************************************************

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>

#include "implicit_functions/Implicit_function_interface.h"
#include "Scene_implicit_function_item.h"
#include "ui_Function_dialog.h"

#include <iostream>

#include <QAction>
#include <QMainWindow>
#include <QPluginLoader>
#include <QDir>
#include <QMenu>
#include <QList>
#include <QLibrary>
using namespace CGAL::Three;
class Io_implicit_function_plugin :
  public QObject,
  protected Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "implicit_function_io_plugin.json")

public:
  Io_implicit_function_plugin();
  virtual ~Io_implicit_function_plugin() {}

  bool applicable(QAction*) const { return true; }

  QString name() const { return "implicit functions"; }
  // QString nameFilters() const { return ""; }
  // bool canLoad() const { return false; }


  typedef Polyhedron_demo_plugin_helper Plugin_helper;

  using Plugin_helper::init;
  virtual void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*);

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
init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface *)
{
  this->scene = scene_interface;
  this->mw = mainWindow;

  QAction* actionLoadFunction = new QAction("Generate &Implicit Function", mw);
  if( NULL != actionLoadFunction )
  {
    connect(actionLoadFunction, SIGNAL(triggered()), this, SLOT(load_function()));
  }

  QMenu* menuFile = mw->findChild<QMenu*>("menuFile");
  QMenu* menu = menuFile->findChild<QMenu*>("menuGenerateObject");
  if(!menu){
    QAction* actionLoad = mw->findChild<QAction*>("actionLoadPlugin");
    menu = new QMenu(tr("Generate &Object"), menuFile);
    menu->setObjectName("menuGenerateObject");
    menuFile->insertMenu(actionLoad, menu);
  }
  // Insert "Generate Implicit Function" action
    menu->addAction(actionLoadFunction);
}


void
Io_implicit_function_plugin::
load_function() const
{
  QDialog dialog(mw);
  Ui::FunctionDialog ui;
  ui.setupUi(&dialog);
  dialog.setWindowFlags(Qt::Dialog|Qt::CustomizeWindowHint|Qt::WindowCloseButtonHint);
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

  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  scene->itemChanged(index);

  CGAL::Three::Scene_interface::Item_id new_item_id = scene->addItem(item);
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

#include "Implicit_function_io_plugin.moc"
