#include "config.h"
#include "MainWindow.h"
#include "Scene.h"
#include "Scene_item.h"
#include <CGAL/Qt/debug.h>

#include <QtDebug>
#include <QUrl>
#include <QFileDialog>
#include <QFileInfo>
#include <QSettings>
#include <QHeaderView>
#include <QMenu>
#include <QAction>
#include <QLibrary>
#include <QPluginLoader>
#include <QMessageBox>
#include <QScrollBar>

#include "Polyhedron_demo_plugin_interface.h"
#include "Polyhedron_demo_io_plugin_interface.h"

#include "ui_MainWindow.h"

MainWindow::~MainWindow()
{
  delete ui;
}

MainWindow::MainWindow(QWidget* parent)
  : CGAL::Qt::DemosMainWindow(parent)
{
  ui = new Ui::MainWindow;
  ui->setupUi(this);

  // Saves some pointers from ui, for latter use.
  treeView = ui->treeView;
  viewer = ui->viewer;

  // Setup the submenu of the View menu that can toggle the dockwidgets
  Q_FOREACH(QDockWidget* widget, findChildren<QDockWidget*>()) {
    ui->menuDockWindows->addAction(widget->toggleViewAction());
  }
  ui->menuDockWindows->removeAction(ui->dummyAction);

  // do not save the state of the viewer (anoying)
  viewer->setStateFileName(QString::null);

  // accept drop events
  setAcceptDrops(true);

  // setup scene
  scene = new Scene(this);
  viewer->setScene(scene);
  treeView->setModel(scene);

  // setup the treeview: delegation and columns sizing...
  treeView->setItemDelegate(new SceneDelegate(this));

  treeView->header()->setStretchLastSection(false);
  treeView->header()->setResizeMode(Scene::NameColumn, QHeaderView::Stretch);
  treeView->header()->setResizeMode(Scene::NameColumn, QHeaderView::Stretch);
  treeView->header()->setResizeMode(Scene::ColorColumn, QHeaderView::ResizeToContents);
  treeView->header()->setResizeMode(Scene::RenderingModeColumn, QHeaderView::Fixed);
  treeView->header()->setResizeMode(Scene::VisibleColumn, QHeaderView::Fixed);

  treeView->resizeColumnToContents(Scene::ColorColumn);
  treeView->resizeColumnToContents(Scene::RenderingModeColumn);
  treeView->resizeColumnToContents(Scene::VisibleColumn);

  // setup connections
  connect(scene, SIGNAL(dataChanged(const QModelIndex &, const QModelIndex & )),
          this, SLOT(updateInfo()));

  connect(scene, SIGNAL(dataChanged(const QModelIndex &, const QModelIndex & )),
          viewer, SLOT(updateGL()));

  connect(scene, SIGNAL(updated()),
          viewer, SLOT(update()));

  connect(scene, SIGNAL(itemAboutToBeDestroyed(Scene_item*)),
          this, SLOT(removeManipulatedFrame(Scene_item*)));

  connect(scene, SIGNAL(updated_bbox()),
          this, SLOT(updateViewerBBox()));

  connect(treeView->selectionModel(),
          SIGNAL(selectionChanged ( const QItemSelection & , const QItemSelection & ) ),
          this, SLOT(updateInfo()));

  connect(treeView->selectionModel(),
          SIGNAL(selectionChanged ( const QItemSelection & , const QItemSelection & ) ),
          this, SLOT(selectionChanged()));

  connect(viewer, SIGNAL(selected(int)),
          this, SLOT(selectSceneItem(int)));

  connect(ui->actionAntiAliasing, SIGNAL(toggled(bool)),
          viewer, SLOT(setAntiAliasing(bool)));

  connect(ui->actionDraw_two_sides, SIGNAL(toggled(bool)),
          viewer, SLOT(setTwoSides(bool)));

  // enable anti-aliasing by default
  ui->actionAntiAliasing->setChecked(true);

  // add the "About CGAL..." and "About demo..." entries
  this->addAboutCGAL();
  this->addAboutDemo(":/cgal/Point_set_demo/about.html");

  // Connect the button "addButton" with actionFileOpen
  ui->addButton->setDefaultAction(ui->actionFileOpen);
  // Same with "removeButton" and "duplicateButton"
  ui->removeButton->setDefaultAction(ui->actionFileClose);
  ui->duplicateButton->setDefaultAction(ui->actionDuplicate);

  // Connect actionQuit (Ctrl+Q) and qApp->quit()
  connect(ui->actionQuit, SIGNAL(triggered()),
          this, SLOT(quit()));

  // Recent files menu
  this->addRecentFiles(ui->menuFile, ui->actionQuit);
  connect(this, SIGNAL(openRecentFile(QString)),
	  this, SLOT(open(QString)));

  // Empty the menus implemented by plugins:
  // "Analysis", "Processing", "Reconstruction".
  clearMenu(ui->menuAnalysis);
  clearMenu(ui->menuProcessing);
  clearMenu(ui->menuReconstruction);

  // Loads plugins, and re-enable actions that need it.
  loadPlugins();

  readSettings(); // Among other things, the column widths are stored.
}

void MainWindow::loadPlugins()
{
  Q_FOREACH(QObject *obj, QPluginLoader::staticInstances())
  {
    initPlugin(obj);
    initIOPlugin(obj);
  }

  QDir pluginsDir(qApp->applicationDirPath());
  Q_FOREACH (QString fileName, pluginsDir.entryList(QDir::Files)) {
    if(fileName.contains("plugin") && QLibrary::isLibrary(fileName)) {
      qDebug("### Loading \"%s\"...", fileName.toUtf8().data());
      QPluginLoader loader;
      loader.setFileName(pluginsDir.absoluteFilePath(fileName));
      QObject *obj = loader.instance();
      if(obj) {
        initPlugin(obj);
        initIOPlugin(obj);
      }
      else {
        qDebug("Error loading \"%s\": %s",
               qPrintable(fileName),
               qPrintable(loader.errorString()));
      }
    }
  }
}

bool MainWindow::initPlugin(QObject* obj)
{
  QObjectList childs = this->children();
  Polyhedron_demo_plugin_interface* plugin =
    qobject_cast<Polyhedron_demo_plugin_interface*>(obj);
  if(plugin) {
    // Calls plugin's init() method
    plugin->init(this, this->scene, this);

    Q_FOREACH(QAction* action, plugin->actions()) {
      // If action does not belong to the menus, add it to "Edit" menu.
      // TODO: implement something less naive.
      if(!childs.contains(action)) {
        std::cerr << "Add " << action->text().toStdString() << " menu item to the Edit menu\n";
        ui->menuEdit->addAction(action);
      }
      // Show and enable menu item
      addAction(action);
    }
    return true;
  }
  else
    return false;
}

bool MainWindow::initIOPlugin(QObject* obj)
{
  Polyhedron_demo_io_plugin_interface* plugin =
    qobject_cast<Polyhedron_demo_io_plugin_interface*>(obj);
  if(plugin) {
//     std::cerr << "I/O plugin\n";
    io_plugins << plugin;
    return true;
  }
  else
    return false;
}

void MainWindow::clearMenu(QMenu* menu)
{
  Q_FOREACH(QAction* action, menu->actions())
  {
    QMenu* menu = action->menu();
    if(menu) {
      clearMenu(menu);
    }
    action->setVisible(false);
  }
  menu->menuAction()->setEnabled(false);
}

void MainWindow::addAction(QAction* action)
{
  if(!action) return;

  action->setVisible(true);
  action->setEnabled(true);
  Q_FOREACH(QWidget* widget, action->associatedWidgets())
  {
//     qDebug() << QString("%1 (%2)\n")
//       .arg(widget->objectName())
//       .arg(widget->metaObject()->className());
    QMenu* menu = qobject_cast<QMenu*>(widget);
    if(menu)
    {
      addAction(menu->menuAction());
    }
  }
}

void MainWindow::message(QString message, QString colorName, QString /*font*/) {
  if (message.endsWith('\n')) {
    message.remove(message.length()-1, 1);
  }
  statusBar()->showMessage(message, 5000);
  message = "<font color=\"" + colorName + "\" >" + message + "</font><br>";
  message = "[" + QTime::currentTime().toString() + "] " + message;
  ui->consoleTextEdit->insertHtml(message);
  ui->consoleTextEdit->verticalScrollBar()->setValue(ui->consoleTextEdit->verticalScrollBar()->maximum());
}

void MainWindow::information(QString text) {
  this->message("INFO: " + text, "");
}

void MainWindow::warning(QString text) {
  this->message("WARNING: " + text, "blue");
}

void MainWindow::error(QString text) {
  this->message("ERROR: " + text, "red");
}


void MainWindow::updateViewerBBox()
{
  const Scene::Bbox bbox = scene->bbox();
  const double xmin = bbox.xmin;
  const double ymin = bbox.ymin;
  const double zmin = bbox.zmin;
  const double xmax = bbox.xmax;
  const double ymax = bbox.ymax;
  const double zmax = bbox.zmax;
  // qDebug() << QString("Bounding box: (%1, %2, %3) - (%4, %5, %6)\n")
  // .arg(xmin).arg(ymin).arg(zmin).arg(xmax).arg(ymax).arg(zmax);
  qglviewer::Vec
    vec_min(xmin, ymin, zmin),
    vec_max(xmax, ymax, zmax);
  viewer->setSceneBoundingBox(vec_min,
                              vec_max);
  viewer->camera()->showEntireScene();
}

void MainWindow::open(QString filename)
{
  QFileInfo fileinfo(filename);
  Scene_item* item = 0;
  if(fileinfo.isFile() && fileinfo.isReadable()) {
    Q_FOREACH(Polyhedron_demo_io_plugin_interface* plugin,
              io_plugins)
    {
      if(plugin->canLoad()) {
        item = plugin->load(fileinfo);
        if(item) break; // go out of the loop
      }
    }
    if(item) {
      Scene::Item_id index = scene->addItem(item);
      QSettings settings;
      settings.setValue("Point set open directory",
                        fileinfo.absoluteDir().absolutePath());
      this->addToRecentFiles(filename);
      selectSceneItem(index);
    }
    else {
      QMessageBox::critical(this,
                            tr("Cannot open file"),
                            tr("File %1 has not a known file format.")
                            .arg(filename));
    }
  }
  else {
    QMessageBox::critical(this,
                          tr("Cannot open file"),
                          tr("File %1 is not a readable file.")
                          .arg(filename));
  }
}

void MainWindow::selectSceneItem(int i)
{
  if(i < 0) return;
  if((unsigned int)i >= scene->numberOfEntries()) return;

  treeView->selectionModel()->select(scene->createSelection(i),
                                     QItemSelectionModel::ClearAndSelect);
}

int MainWindow::getSelectedSceneItemIndex() const
{
  QModelIndexList selectedRows = treeView->selectionModel()->selectedRows();
  if(selectedRows.empty())
    return -1;
  else
    return selectedRows.first().row();
}

void MainWindow::selectionChanged()
{
  scene->setSelectedItem(getSelectedSceneItemIndex());
  Scene_item* item = scene->item(getSelectedSceneItemIndex());
  if(item != NULL && item->manipulatable()) {
    viewer->setManipulatedFrame(item->manipulatedFrame());
    connect(viewer->manipulatedFrame(), SIGNAL(modified()),
            this, SLOT(updateInfo()));
  }

  viewer->updateGL();
}

void MainWindow::removeManipulatedFrame(Scene_item* item)
{
  if(item->manipulatable() &&
     item->manipulatedFrame() == viewer->manipulatedFrame()) {
    viewer->setManipulatedFrame(0);
  }
}

void MainWindow::updateInfo() {
  Scene_item* item = scene->item(getSelectedSceneItemIndex());
  if(item)
    ui->infoLabel->setText(item->toolTip());
  else
    ui->infoLabel->clear();
}

void MainWindow::readSettings()
{
  this->readState("MainWindow", Size|State);
}

void MainWindow::writeSettings()
{
  this->writeState("MainWindow");
  std::cerr << "Write setting... done.\n";
}

void MainWindow::quit()
{
  close();
}

void MainWindow::closeEvent(QCloseEvent *event)
{
  writeSettings();
  event->accept();
}

void MainWindow::on_actionFileOpen_triggered()
{
  QStringList filters;
  Q_FOREACH(Polyhedron_demo_io_plugin_interface* plugin, io_plugins) {
    if(plugin->canLoad()) {
      filters += plugin->nameFilters();
    }
  }
  filters << tr("All files (*)");

  QSettings settings;
  QString directory = settings.value("Point set open directory",
                                     QDir::current().dirName()).toString();
  QStringList filenames =
    QFileDialog::getOpenFileNames(this,
                                  tr("Open File..."),
                                  directory,
                                  filters.join(";;"));
  if(!filenames.isEmpty()) {
    Q_FOREACH(QString filename, filenames) {
      open(filename);
    }
  }
}

void MainWindow::on_actionSaveAs_triggered()
{
  QModelIndexList selectedRows = treeView->selectionModel()->selectedRows();
  if(selectedRows.size() != 1)
    return;
  Scene_item* item = scene->item(getSelectedSceneItemIndex());
  if(!item)
    return;

  QVector<Polyhedron_demo_io_plugin_interface*> canSavePlugins;
  QStringList filters;
  Q_FOREACH(Polyhedron_demo_io_plugin_interface* plugin, io_plugins) {
    if(plugin->canSave(item)) {
      canSavePlugins << plugin;
      filters += plugin->nameFilters();
    }
  }
  filters << tr("All files (*)");

  if(canSavePlugins.isEmpty()) {
    QMessageBox::warning(this,
                         tr("Cannot save"),
                         tr("The selected object %1 cannot be saved.")
                         .arg(item->name()));
    return;
  }

  QSettings settings;
  QString directory = settings.value("Point set save directory",
                                     QDir::current().dirName()).toString();
  QString filename =
    QFileDialog::getSaveFileName(this,
                                 tr("Save As..."),
                                 directory,
                                 filters.join(";;"));
  if (filename.isEmpty())
    return;
    
  QFileInfo fileinfo(filename);

  bool saved = false;
  Q_FOREACH(Polyhedron_demo_io_plugin_interface* plugin, canSavePlugins) {
    if(plugin->save(item, fileinfo)) {
      saved = true;
      break;
    }
  }
  if (saved) {
    settings.setValue("Point set save directory",
                      fileinfo.absoluteDir().absolutePath());
  }
  else {
    QMessageBox::warning(this,
                         tr("Cannot save"),
                         tr("Error while saving object %1 as %2.")
                         .arg(item->name())
                         .arg(filename));
  }
}

bool MainWindow::on_actionFileClose_triggered()
{
  int index = scene->erase(getSelectedSceneItemIndex());
  selectSceneItem(index);
  return index >= 0;
}

void MainWindow::on_actionFileCloseAll_triggered()
{
  while(on_actionFileClose_triggered()) {
  }
}

void MainWindow::on_actionDuplicate_triggered()
{
  int index = scene->duplicate(getSelectedSceneItemIndex());
  selectSceneItem(index);
}

void MainWindow::on_actionConvertToPointSet_triggered()
{
  int index = scene->convertToPointSet(getSelectedSceneItemIndex());
  selectSceneItem(index);
}

void MainWindow::on_actionDeleteSelection_triggered()
{
  scene->deleteSelection(getSelectedSceneItemIndex());
}

void MainWindow::on_actionResetSelection_triggered()
{
  scene->resetSelection(getSelectedSceneItemIndex());
}

void MainWindow::on_actionShowHide_triggered()
{
  Q_FOREACH(QModelIndex index, treeView->selectionModel()->selectedRows())
  {
    int i = index.row();
    Scene_item* item = scene->item(i);
    item->setVisible(!item->visible());
    scene->itemChanged(i);
  }
}

void MainWindow::setAddKeyFrameKeyboardModifiers(::Qt::KeyboardModifiers m)
{
  viewer->setAddKeyFrameKeyboardModifiers(m);
}
