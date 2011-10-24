#include "config.h"
#include "MainWindow.h"
#include <CGAL_demo/Scene.h>
#include <CGAL_demo/Scene_item.h>
#include <CGAL/Qt/debug.h>

#include <QDragEnterEvent>
#include <QDropEvent>
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
#include <QClipboard>

#include <CGAL_demo/Plugin_interface.h>
#include <CGAL_demo/Io_plugin_interface.h>

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

  // Save some pointers from ui, for latter use.
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
  treeView->header()->setResizeMode(Scene::ABColumn, QHeaderView::Fixed);
  treeView->header()->setResizeMode(Scene::VisibleColumn, QHeaderView::Fixed);

  treeView->resizeColumnToContents(Scene::ColorColumn);
  treeView->resizeColumnToContents(Scene::RenderingModeColumn);
  treeView->resizeColumnToContents(Scene::ABColumn);
  treeView->resizeColumnToContents(Scene::VisibleColumn);

  // setup connections
  connect(scene, SIGNAL(dataChanged(const QModelIndex &, const QModelIndex & )),
          this, SLOT(updateInfo()));
  
  connect(scene, SIGNAL(dataChanged(const QModelIndex &, const QModelIndex & )),
          this, SLOT(updateDisplayInfo()));

  connect(scene, SIGNAL(dataChanged(const QModelIndex &, const QModelIndex & )),
          viewer, SLOT(updateGL()));

  connect(scene, SIGNAL(updated()),
          viewer, SLOT(update()));
  
  connect(scene, SIGNAL(selectionChanged()),
          this, SLOT(selectSceneItem()));

  connect(scene, SIGNAL(itemAboutToBeDestroyed(Scene_item*)),
          this, SLOT(removeManipulatedFrame(Scene_item*)));

  connect(scene, SIGNAL(updated_bbox()),
          this, SLOT(updateViewerBBox()));

  connect(treeView->selectionModel(), 
          SIGNAL(selectionChanged ( const QItemSelection & , const QItemSelection & ) ),
          this, SLOT(updateInfo()));
  
  connect(treeView->selectionModel(), 
          SIGNAL(selectionChanged ( const QItemSelection & , const QItemSelection & ) ),
          this, SLOT(updateDisplayInfo()));

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
  this->addAboutDemo(":/cgal/Mesh_3/about.html");

  // Connect the button "addButton" with actionLoad
  ui->addButton->setDefaultAction(ui->actionLoad);
  // Same with "removeButton" and "duplicateButton"
  ui->removeButton->setDefaultAction(ui->actionErase);
  ui->duplicateButton->setDefaultAction(ui->actionDuplicate);

  // Connect actionQuit (Ctrl+Q) and qApp->quit()
  connect(ui->actionQuit, SIGNAL(triggered()),
          this, SLOT(quit()));

  // Recent files menu
  this->addRecentFiles(ui->menuFile, ui->actionQuit);
  connect(this, SIGNAL(openRecentFile(QString)),
	  this, SLOT(open(QString)));

  // Reset the "Operation menu"
  clearMenu(ui->menuOperations);

  // Load plugins, and re-enable actions that need it.
  loadPlugins();

  readSettings(); // Among other things, the column widths are stored.

  this->dumpObjectTree();
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
  Plugin_interface* plugin =
    qobject_cast<Plugin_interface*>(obj);
  if(plugin) {
    // Call plugin's init() method
    plugin->init(this, this->scene, this);

    Q_FOREACH(QAction* action, plugin->actions()) {
      // If action does not belong to the menus, add it to "Operations" menu
      if(!childs.contains(action)) {
        ui->menuOperations->addAction(action);
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
  Io_plugin_interface* plugin =
    qobject_cast<Io_plugin_interface*>(obj);
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
  
  // Moves cursor to the end of the block
  ui->consoleTextEdit->moveCursor(QTextCursor::EndOfBlock);
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

void MainWindow::dragEnterEvent(QDragEnterEvent *event)
{
  if (event->mimeData()->hasFormat("text/uri-list"))
    event->acceptProposedAction();
}

void MainWindow::dropEvent(QDropEvent *event)
{
  Q_FOREACH(QUrl url, event->mimeData()->urls()) {
    QString filename = url.toLocalFile();
    if(!filename.isEmpty()) {
      qDebug() << QString("dropEvent(\"%1\")\n").arg(filename);
      open(filename);
    }
  }
  event->acceptProposedAction();
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
    Q_FOREACH(Io_plugin_interface* plugin, 
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
      settings.setValue("OFF open directory",
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

void
MainWindow::selectSceneItem()
{
  selectSceneItem(scene->mainSelectionIndex());
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
//    connect(viewer->manipulatedFrame(), SIGNAL(modified()),
//            this, SLOT(updateInfo()));
//    connect(viewer->manipulatedFrame(), SIGNAL(modified()),
//            this, SLOT(updateDisplayInfo()));
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

void MainWindow::updateDisplayInfo() {
  Scene_item* item = scene->item(getSelectedSceneItemIndex());
  if(item)
    ui->displayLabel->setPixmap(item->graphicalToolTip());
  else 
    ui->displayLabel->clear();
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

void MainWindow::on_actionLoad_triggered()
{
  QStringList filters;
  Q_FOREACH(Io_plugin_interface* plugin, io_plugins) {
    if(plugin->canLoad()) {
      filters += plugin->nameFilters();
    }
  }
  filters << tr("All files (*)");

  QSettings settings;
  QString directory = settings.value("OFF open directory",
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

  QVector<Io_plugin_interface*> canSavePlugins;
  QStringList filters;
  Q_FOREACH(Io_plugin_interface* plugin, io_plugins) {
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
  
  QString filename = 
    QFileDialog::getSaveFileName(this,
                                 tr("Save to File..."),
                                 QString(),
                                 filters.join(";;"));
  
  QFileInfo fileinfo(filename);
  if(!fileinfo.isFile() ||
     QMessageBox::warning(this,
                          tr("File exists"),
                          tr("The file %1 already exists! Continue?")
                          .arg(filename),
                          QMessageBox::Yes|QMessageBox::No) == 
     QMessageBox::Yes)
  {

    Q_FOREACH(Io_plugin_interface* plugin, canSavePlugins) {
      if(plugin->save(item, fileinfo))
        break;
    }
  }
}

bool MainWindow::on_actionErase_triggered()
{
  int index = scene->erase(getSelectedSceneItemIndex());
  selectSceneItem(index);
  return index >= 0;
}

void MainWindow::on_actionEraseAll_triggered()
{
  while(on_actionErase_triggered()) {
  }
}

void MainWindow::on_actionDuplicate_triggered()
{
  int index = scene->duplicate(getSelectedSceneItemIndex());
  selectSceneItem(index);
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

void MainWindow::on_actionSetPolyhedronA_triggered()
{
  int i = getSelectedSceneItemIndex();
  scene->setItemA(i);
}

void MainWindow::on_actionSetPolyhedronB_triggered()
{
  int i = getSelectedSceneItemIndex();
  scene->setItemB(i);
}

void MainWindow::setAddKeyFrameKeyboardModifiers(::Qt::KeyboardModifiers m)
{
  viewer->setAddKeyFrameKeyboardModifiers(m);
}

void MainWindow::on_actionCopy_snapshot_triggered()
{
  // copy snapshot to clipboard
	QApplication::setOverrideCursor(Qt::WaitCursor);
  QClipboard *qb = QApplication::clipboard();
  viewer->makeCurrent();
  viewer->raise();
  QImage snapshot = viewer->grabFrameBuffer(true);
  qb->setImage(snapshot);
	QApplication::restoreOverrideCursor();
}

void MainWindow::on_actionSave_snapshot_triggered()
{
	// save snapshot to file
	QApplication::setOverrideCursor(Qt::WaitCursor);
        viewer->saveSnapshot(false, false);
	QApplication::restoreOverrideCursor();
}



