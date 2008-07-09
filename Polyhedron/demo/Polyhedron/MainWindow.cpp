#include "MainWindow.h"
#include "Scene.h"
#include <CGAL/Qt/debug.h>

#include <QDragEnterEvent>
#include <QDropEvent>
#include <QTextStream>
#include <QUrl>
#include <QFileDialog>

MainWindow::MainWindow(QWidget* parent)
  : CGAL::Qt::DemosMainWindow(parent)
{
  setupUi(this);

  addDockWidget(Qt::LeftDockWidgetArea, polyhedraDockWidget);
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

  connect(scene, SIGNAL(dataChanged(const QModelIndex &, const QModelIndex & )),
          viewer, SLOT(updateGL()));

  connect(scene, SIGNAL(updated()),
          viewer, SLOT(update()));

  connect(scene, SIGNAL(updated_bbox()),
          this, SLOT(updateViewerBBox()));

  connect(treeView->selectionModel(), 
          SIGNAL(selectionChanged ( const QItemSelection & , const QItemSelection & ) ),
          this, SLOT(on_treeView_itemSelectionChanged()));

  connect(actionAntiAliasing, SIGNAL(toggled(bool)),
          viewer, SLOT(setAntiAliasing(bool)));

  // add the "About CGAL..." entry
  this->addAboutCGAL();

  // Connect the button "addButton" with actionLoadPolyhedron
  addButton->setDefaultAction(actionLoadPolyhedron);
  // Same with "removeButton" and "duplicateButton"
  removeButton->setDefaultAction(actionErasePolyhedron);
  duplicateButton->setDefaultAction(actionDuplicatePolyhedron);

  // Connect actionQuit (Ctrl+Q) and qApp->quit()
  connect(actionQuit, SIGNAL(triggered()),
          qApp, SLOT(quit()));
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
      scene->open(filename);
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
  scene->open(filename);
}


Polyhedron* MainWindow::getSelectedPolygon()
{
  // scene->getPolyhedron(...) returns 0 iff the index is not valid
  return scene->getPolyhedron(getSelectedPolygonIndex());
}

void MainWindow::on_treeView_itemSelectionChanged()
{
  if(onePolygonIsSelected()) {
    scene->setSelectedItem(getSelectedPolygonIndex());
  }
  else {
    scene->setSelectedItem(-1);
  }
  viewer->updateGL();
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
      scene->open(filename);
    }
  }
}

void MainWindow::on_actionErasePolyhedron_triggered()
{
  if(onePolygonIsSelected())
    scene->erase(getSelectedPolygonIndex());
}

void MainWindow::on_actionDuplicatePolyhedron_triggered()
{
  if(onePolygonIsSelected())
    scene->duplicate(getSelectedPolygonIndex());
}

