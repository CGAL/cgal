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

#include "ui_MainWindow.h"

MainWindow::MainWindow(QWidget* parent)
  : CGAL::Qt::DemosMainWindow(parent)
{
  ui = new Ui::MainWindow;
  ui->setupUi(this);

  // Save some pointers from ui, for latter use.
  treeView = ui->treeView;
  viewer = ui->viewer;

  addDockWidget(::Qt::LeftDockWidgetArea, ui->polyhedraDockWidget);
  ui->menuDockWindows->addAction(ui->polyhedraDockWidget->toggleViewAction());
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
  treeView->header()->setResizeMode(Scene::ActivatedColumn, QHeaderView::Fixed);

  treeView->resizeColumnToContents(Scene::ColorColumn);
  treeView->resizeColumnToContents(Scene::RenderingModeColumn);
  treeView->resizeColumnToContents(Scene::ABColumn);
  treeView->resizeColumnToContents(Scene::ActivatedColumn);


  connect(scene, SIGNAL(dataChanged(const QModelIndex &, const QModelIndex & )),
          viewer, SLOT(updateGL()));

  connect(scene, SIGNAL(updated()),
          viewer, SLOT(update()));

  connect(scene, SIGNAL(updated_bbox()),
          this, SLOT(updateViewerBBox()));

  connect(treeView->selectionModel(), 
          SIGNAL(selectionChanged ( const QItemSelection & , const QItemSelection & ) ),
          this, SLOT(selectionChanged()));

  connect(viewer, SIGNAL(selected(int)),
          this, SLOT(selectPolyhedron(int)));

  connect(ui->actionAntiAliasing, SIGNAL(toggled(bool)),
          viewer, SLOT(setAntiAliasing(bool)));

  ui->actionAntiAliasing->setChecked(true);

  connect(ui->actionViewEdges, SIGNAL(toggled(bool)),
          scene, SLOT(setViewEdges(bool)));

  ui->actionViewEdges->setChecked(true);

  // add the "About CGAL..." and "About demo..." entries
  this->addAboutCGAL();
  this->addAboutDemo(":/cgal/Polyhedron_3/about.html");

  // Connect the button "addButton" with actionLoadPolyhedron
  ui->addButton->setDefaultAction(ui->actionLoadPolyhedron);
  // Same with "removeButton" and "duplicateButton"
  ui->removeButton->setDefaultAction(ui->actionErasePolyhedron);
  ui->duplicateButton->setDefaultAction(ui->actionDuplicatePolyhedron);

  // Connect actionQuit (Ctrl+Q) and qApp->quit()
  connect(ui->actionQuit, SIGNAL(triggered()),
          this, SLOT(quit()));

  this->addRecentFiles(ui->menuFile, ui->actionQuit);
  connect(this, SIGNAL(openRecentFile(QString)),
	  this, SLOT(open(QString)));

  readSettings(); // Among other things, the column widths are stored.
}

MainWindow::~MainWindow()
{
  delete ui;
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
  const Scene::Bbox bbox = scene->bbox();
  const double xmin = bbox.xmin;
  const double ymin = bbox.ymin;
  const double zmin = bbox.zmin;
  const double xmax = bbox.xmax;
  const double ymax = bbox.ymax;
  const double zmax = bbox.zmax;
  // QTextStream(stderr) << QString("Bounding box: (%1, %2, %3) - (%4, %5, %6)\n")
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
  if(fileinfo.isFile() && fileinfo.isReadable()) {
    int index = scene->open(filename);
    if(index >= 0) {
      QSettings settings;
      settings.setValue("OFF open directory",
			fileinfo.absoluteDir().absolutePath());
	this->addToRecentFiles(filename);
      selectPolyhedron(index);
    }
  }
}

void MainWindow::selectPolyhedron(int i)
{
  if(i < 0) return;
  if(i >= scene->numberOfPolyhedra()) return;

  treeView->selectionModel()->select(scene->createSelection(i),
                                     QItemSelectionModel::ClearAndSelect);
}

bool MainWindow::onePolygonIsSelected() const
{
  QModelIndexList selectedRows = treeView->selectionModel()->selectedRows();
  if(selectedRows.size() == 1)
  {
    int i = selectedRows.first().row();
    if(scene->polyhedronType(i) == Scene::POLYHEDRON_ENTRY)
      return true;
  }
  return false;
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
  writeSettings();
  close();
}

void MainWindow::closeEvent(QCloseEvent *event)
{
  writeSettings();
  event->accept();
}

void MainWindow::on_actionLoadPolyhedron_triggered()
{
  QSettings settings;
  QString directory = settings.value("OFF open directory",
				     QDir::current().dirName()).toString();
  QStringList filenames = 
    QFileDialog::getOpenFileNames(this,
                                  tr("Load polyhedron..."),
                                  directory,
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




bool MainWindow::on_actionErasePolyhedron_triggered()
{
  int index = scene->erase(getSelectedPolygonIndex());
  selectPolyhedron(index);
  return index >= 0;
}

void MainWindow::on_actionEraseAll_triggered()
{
  while(on_actionErasePolyhedron_triggered()) {
  }
}

void MainWindow::on_actionDuplicatePolyhedron_triggered()
{
  int index = scene->duplicate(getSelectedPolygonIndex());
  selectPolyhedron(index);
}

void MainWindow::on_actionActivatePolyhedron_triggered()
{
  Q_FOREACH(QModelIndex index, treeView->selectionModel()->selectedRows())
  {
    int i = index.row();
    scene->setPolyhedronActivated(i,
                                  !scene->isPolyhedronActivated(i));
  }
}

void MainWindow::on_actionSetPolyhedronA_triggered()
{
  int i = getSelectedPolygonIndex();
  scene->setPolyhedronA(i);
}

void MainWindow::on_actionSetPolyhedronB_triggered()
{
  int i = getSelectedPolygonIndex();
  scene->setPolyhedronB(i);
}

void MainWindow::setAddKeyFrameKeyboardModifiers(::Qt::KeyboardModifiers m)
{
  viewer->setAddKeyFrameKeyboardModifiers(m);
}
