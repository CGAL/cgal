#include "MainWindow.h"
#include "Scene.h"
#include <CGAL/Qt/debug.h>

#include <QDragEnterEvent>
#include <QDropEvent>
#include <QTextStream>
#include <QUrl>
#include <QFileDialog>
#include <QSettings>
#include <QHeaderView>

MainWindow::MainWindow(QWidget* parent)
  : CGAL::Qt::DemosMainWindow(parent)
{
  setupUi(this);

  addDockWidget(::Qt::LeftDockWidgetArea, polyhedraDockWidget);
  menuDockWindows->addAction(polyhedraDockWidget->toggleViewAction());
  menuDockWindows->removeAction(dummyAction);

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
  treeView->resizeColumnToContents(Scene::ColorColumn);
  treeView->resizeColumnToContents(Scene::RenderingModeColumn);
  treeView->resizeColumnToContents(Scene::ABColumn);
  treeView->resizeColumnToContents(Scene::ActivatedColumn);

  treeView->header()->setStretchLastSection(false);
  treeView->header()->setResizeMode(Scene::NameColumn, QHeaderView::Stretch);
  treeView->header()->setResizeMode(Scene::ColorColumn, QHeaderView::ResizeToContents);
  treeView->header()->setResizeMode(Scene::RenderingModeColumn, QHeaderView::Fixed);
  treeView->header()->setResizeMode(Scene::ABColumn, QHeaderView::Fixed);
  treeView->header()->setResizeMode(Scene::ActivatedColumn, QHeaderView::Fixed);

  connect(scene, SIGNAL(dataChanged(const QModelIndex &, const QModelIndex & )),
          viewer, SLOT(updateGL()));

  connect(scene, SIGNAL(updated()),
          viewer, SLOT(update()));

  connect(scene, SIGNAL(updated_bbox()),
          this, SLOT(updateViewerBBox()));

  connect(treeView->selectionModel(), 
          SIGNAL(selectionChanged ( const QItemSelection & , const QItemSelection & ) ),
          this, SLOT(selectionChanged()));

  connect(actionAntiAliasing, SIGNAL(toggled(bool)),
          viewer, SLOT(setAntiAliasing(bool)));

  actionAntiAliasing->setChecked(true);

  // add the "About CGAL..." and "About demo..." entries
  this->addAboutCGAL();
  this->addAboutDemo(":/cgal/Polyhedron_3/about.html");

  // Connect the button "addButton" with actionLoadPolyhedron
  addButton->setDefaultAction(actionLoadPolyhedron);
  // Same with "removeButton" and "duplicateButton"
  removeButton->setDefaultAction(actionErasePolyhedron);
  duplicateButton->setDefaultAction(actionDuplicatePolyhedron);

  // Connect actionQuit (Ctrl+Q) and qApp->quit()
  connect(actionQuit, SIGNAL(triggered()),
          this, SLOT(quit()));

  // recent files...
  for (int i = 0; i < MaxRecentFiles; ++i) {
    recentFileActs[i] = new QAction(this);
    recentFileActs[i]->setVisible(false);
    connect(recentFileActs[i], SIGNAL(triggered()),
            this, SLOT(openRecentFile()));
    menuFile->insertAction(actionQuit, recentFileActs[i]);
  }
  recentFilesSeparator = menuFile->insertSeparator(actionQuit);
  recentFilesSeparator->setVisible(false);
  updateRecentFileActions();

  readSettings(); // Among other things, the column widths are stored.
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
      QTextStream(stderr) << QString("dropEvent(\"%1\")\n").arg(filename);
      open(filename);
    }
  }
  event->acceptProposedAction();
}

void MainWindow::updateViewerBBox()
{
  const CGAL::Bbox_3 bbox = scene->bbox();
  const double xmin = bbox.xmin();
  const double ymin = bbox.ymin();
  const double zmin = bbox.zmin();
  const double xmax = bbox.xmax();
  const double ymax = bbox.ymax();
  const double zmax = bbox.zmax();
  QTextStream(stderr)
    << QString("Bounding box: (%1, %2, %3) - (%4, %5, %6)\n")
    .arg(xmin).arg(ymin).arg(zmin).arg(xmax).arg(ymax).arg(zmax);
  qglviewer::Vec 
    vec_min(xmin, ymin, zmin),
    vec_max(xmax, ymax, zmax);
  viewer->setSceneBoundingBox(vec_min,
                              vec_max);
  viewer->camera()->showEntireScene();
}

void MainWindow::open(QString filename)
{
  int index = scene->open(filename);
  if(index >= 0) {
    setCurrentFile(filename);
    selectPolyhedron(index);
  }
}



void MainWindow::selectPolyhedron(int i)
{
  treeView->selectionModel()->select(scene->createSelection(i),
                                     QItemSelectionModel::ClearAndSelect);
}

bool MainWindow::onePolygonIsSelected() const
{
  return treeView->selectionModel()->selectedRows().size() == 1;
}

int MainWindow::getSelectedPolygonIndex() const
{
  QModelIndexList selectedRows = treeView->selectionModel()->selectedRows();
  if(selectedRows.empty())
    return -1;
  else
    return selectedRows.first().row();
}

Polyhedron* MainWindow::getSelectedPolygon()
{
  // scene->getPolyhedron(...) returns 0 iff the index is not valid
  return scene->getPolyhedron(getSelectedPolygonIndex());
}

void MainWindow::selectionChanged()
{
  if(onePolygonIsSelected()) {
    scene->setSelectedItem(getSelectedPolygonIndex());
  }
  else {
    scene->setSelectedItem(-1);
  }
  viewer->updateGL();
}

void MainWindow::openRecentFile()
{
  QAction *action = qobject_cast<QAction *>(sender());
  if (action)
    open(action->data().toString());
}

void MainWindow::setCurrentFile(const QString &fileName)
{
  QSettings settings;
  QStringList files = settings.value("recentFileList").toStringList();
  files.removeAll(fileName);
  files.prepend(fileName);
  while (files.size() > MaxRecentFiles)
    files.removeLast();

  settings.setValue("recentFileList", files);

  updateRecentFileActions();
}

QString MainWindow::strippedName(const QString &fullFileName)
{
  return QFileInfo(fullFileName).fileName();
}

void MainWindow::updateRecentFileActions()
{
  QSettings settings;
  QStringList files = settings.value("recentFileList").toStringList();

  int numRecentFiles = qMin(files.size(), (int)MaxRecentFiles);

  for (int i = 0; i < numRecentFiles; ++i) {
    QString text = tr("&%1 %2").arg(i).arg(strippedName(files[i]));
    recentFileActs[i]->setText(text);
    recentFileActs[i]->setData(files[i]);
    recentFileActs[i]->setVisible(true);
  }
  for (int j = numRecentFiles; j < MaxRecentFiles; ++j)
    recentFileActs[j]->setVisible(false);

  recentFilesSeparator->setVisible(numRecentFiles > 0);
}

void MainWindow::quit()
{
  writeSettings();
  close();
}

void MainWindow::readSettings()
{
  QSettings settings;
  
  settings.beginGroup("MainWindow");
  resize(settings.value("size", QSize(400, 400)).toSize());
  move(settings.value("pos", QPoint(0, 0)).toPoint());
  settings.endGroup();

  settings.beginGroup("SceneView");
  treeView->setColumnWidth(Scene::NameColumn,
                           settings.value("NameColumn width").toInt());
  settings.endGroup();
}

void MainWindow::writeSettings()
{
  QSettings settings;

  settings.beginGroup("MainWindow");
  settings.setValue("size", size());
  settings.setValue("pos", pos());
  settings.endGroup();

  settings.beginGroup("SceneView");
  settings.setValue("NameColumn width", treeView->columnWidth(Scene::NameColumn));
  settings.endGroup();
  std::cerr << "Write setting... done.\n";
}

void MainWindow::closeEvent(QCloseEvent *event)
{
  writeSettings();
  event->accept();
}

void MainWindow::on_actionLoadPolyhedron_triggered()
{
  QStringList filenames = 
    QFileDialog::getOpenFileNames(this,
                                  tr("Load polyhedron..."),
                                  QString(),
                                  tr("OFF files (*.off)\n"
                                     "All files (*)"));
  if(!filenames.isEmpty()) {
    Q_FOREACH(QString filename, filenames) {
      open(filename);
    }
  }
}

void MainWindow::on_actionSaveAs_triggered()
{
  if(!onePolygonIsSelected())
	  return;

  QString filename = 
    QFileDialog::getSaveFileName(this,
                                 tr("Save polyhedron..."),
                                 QString(),
                                 tr("OFF files (*.off)\n"
                                    "All files (*)"));
  if(!filename.isEmpty())
        scene->save(getSelectedPolygonIndex(),filename);
}




void MainWindow::on_actionErasePolyhedron_triggered()
{
  if(onePolygonIsSelected()) {
    int index = scene->erase(getSelectedPolygonIndex());
    selectPolyhedron(index);
  }
}

void MainWindow::on_actionDuplicatePolyhedron_triggered()
{
  if(onePolygonIsSelected()) {
    int index = scene->duplicate(getSelectedPolygonIndex());
    selectPolyhedron(index);
  }
}

