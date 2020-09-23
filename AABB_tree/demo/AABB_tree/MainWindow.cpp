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
#include <QClipboard>

#include "ui_MainWindow.h"

MainWindow::MainWindow(QWidget* parent)
: CGAL::Qt::DemosMainWindow(parent)
{
  ui = new Ui::MainWindow;
  ui->setupUi(this);

  // saves some pointers from ui, for latter use.
  m_pViewer = ui->viewer;

  // does not save the state of the viewer
  m_pViewer->setStateFileName(QString());

  // accepts drop events
  setAcceptDrops(true);
  // setups scene
        m_pScene = new Scene();
  m_pViewer->setScene(m_pScene);
        m_pViewer->setManipulatedFrame(m_pScene->manipulatedFrame());

  // connects actionQuit (Ctrl+Q) and qApp->quit()
  connect(ui->actionQuit, SIGNAL(triggered()),
    this, SLOT(quit()));

  this->addRecentFiles(ui->menuFile, ui->actionQuit);
  connect(this, SIGNAL(openRecentFile(QString)),
    this, SLOT(open(QString)));

  readSettings();
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
  m_pScene->update_bbox();
  const Scene::Bbox bbox = m_pScene->bbox();
  const double xmin = bbox.xmin();
  const double ymin = bbox.ymin();
  const double zmin = bbox.zmin();
  const double xmax = bbox.xmax();
  const double ymax = bbox.ymax();
  const double zmax = bbox.zmax();
  CGAL::qglviewer::Vec
    vec_min(xmin, ymin, zmin),
    vec_max(xmax, ymax, zmax);
  m_pViewer->setSceneBoundingBox(vec_min,vec_max);
  m_pViewer->camera()->showEntireScene();
}

void MainWindow::open(QString filename)
{
  QFileInfo fileinfo(filename);
  if(fileinfo.isFile() && fileinfo.isReadable())
  {
    int index = m_pScene->open(filename);
    if(index >= 0)
    {
      QSettings settings;
      settings.setValue("OFF open directory",
        fileinfo.absoluteDir().absolutePath());
      this->addToRecentFiles(filename);

      // update bbox
      updateViewerBBox();
      m_pViewer->update();
    }
  }
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


void MainWindow::setAddKeyFrameKeyboardModifiers(::Qt::KeyboardModifiers m)
{
  m_pViewer->setAddKeyFrameKeyboardModifiers(m);
}

void MainWindow::on_actionInside_points_triggered()
{
  bool ok;

  const unsigned int nb_points = (unsigned)
    QInputDialog::getInt(NULL, "#Points",
    "#Points:",10000,1,100000000,9,&ok);

  if(!ok)
    return;

  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->generate_inside_points(nb_points);
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
}

void MainWindow::on_actionPoints_in_interval_triggered()
{
  bool ok;

  const unsigned int nb_points = (unsigned)
                QInputDialog::getInt(NULL, "#Points",
    "#Points:",10000,1,100000000,9,&ok);

  if(!ok)
    return;

  const double min =
    QInputDialog::getDouble(NULL, "min",
    "Min:",-0.1,-1000.0,1000.0,9,&ok);
  if(!ok)
    return;
  const double max =
    QInputDialog::getDouble(NULL, "max",
    "Max:",0.1,-1000.0,1000.0,9,&ok);
  if(!ok)
    return;

  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->generate_points_in(nb_points,min,max);
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
}

void MainWindow::on_actionBoundary_segments_triggered()
{
  bool ok;

  const unsigned int nb_slices = (unsigned)
    QInputDialog::getInt(NULL, "#Slices",
    "Slices:",100,1,1000000,8,&ok);

  if(!ok)
    return;

  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->generate_boundary_segments(nb_slices);
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
}

void MainWindow::on_actionBoundary_points_triggered()
{
  bool ok;

  const unsigned int nb_points = (unsigned)
    QInputDialog::getInt(NULL, "#Points",
    "Points:",1000,1,10000000,8,&ok);

  if(!ok)
    return;

  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->generate_boundary_points(nb_points);
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
}

void MainWindow::on_actionEdge_points_triggered()
{
  bool ok;

  const unsigned int nb_points = (unsigned)
    QInputDialog::getInt(NULL, "#Points",
    "Points:",1000,1,10000000,8,&ok);

  if(!ok)
    return;

  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->generate_edge_points(nb_points);
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
}

void MainWindow::on_actionBench_distances_triggered()
{
  bool ok;
  const double duration = QInputDialog::getDouble(NULL, "Duration",
    "Duration (s):",1.0,0.01,1000,8,&ok);
  if(!ok)
    return;

  QApplication::setOverrideCursor(Qt::WaitCursor);
  std::cout << std::endl << "Benchmark distances" << std::endl;
  m_pScene->benchmark_distances(duration);
  QApplication::restoreOverrideCursor();
}

void MainWindow::on_actionBench_intersections_triggered()
{
  bool ok;
  const double duration = QInputDialog::getDouble(NULL, "Duration",
    "Duration (s):",1.0,0.01,1000.0,8,&ok);
  if(!ok)
    return;

  QApplication::setOverrideCursor(Qt::WaitCursor);
  std::cout << std::endl << "Benchmark intersections" << std::endl;
  m_pScene->benchmark_intersections(duration);
  QApplication::restoreOverrideCursor();
}

void MainWindow::on_actionUnsigned_distance_function_to_facets_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
    m_pScene->activate_cutting_plane();
  m_pScene->unsigned_distance_function();
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
}

void MainWindow::on_actionUnsigned_distance_function_to_edges_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
    m_pScene->activate_cutting_plane();
  m_pScene->unsigned_distance_function_to_edges();
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
}

void MainWindow::on_actionSigned_distance_function_to_facets_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
    m_pScene->activate_cutting_plane();
  m_pScene->signed_distance_function();
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
}

void MainWindow::on_actionIntersection_cutting_plane_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->activate_cutting_plane();
  m_pScene->cut_segment_plane();
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
}

void MainWindow::on_actionCutting_plane_none_triggered()
{
  m_pScene->clear_cutting_plane();
  m_pViewer->update();
}

void MainWindow::on_actionView_polyhedron_triggered()
{
  m_pScene->toggle_view_poyhedron();
  m_pViewer->update();
}

void MainWindow::on_actionView_points_triggered()
{
  m_pScene->toggle_view_points();
  m_pViewer->update();
}

void MainWindow::on_actionView_segments_triggered()
{
  m_pScene->toggle_view_segments();
  m_pViewer->update();
}

void MainWindow::on_actionView_cutting_plane_triggered()
{
  m_pScene->toggle_view_plane();
  m_pViewer->update();
}

void MainWindow::on_actionClear_points_triggered()
{
  m_pScene->clear_points();
  m_pViewer->update();
}

void MainWindow::on_actionClear_segments_triggered()
{
  m_pScene->clear_segments();
  m_pViewer->update();
}

void MainWindow::on_actionClear_cutting_plane_triggered()
{
  m_pScene->clear_cutting_plane();
  m_pViewer->update();
}

void MainWindow::on_actionRefine_bisection_triggered()
{
  bool ok;
  const double max_len =
    QInputDialog::getDouble(NULL, "Max edge len",
    "Max edge len:",0.1,0.001,100.0,9,&ok);
  if(!ok)
    return;

  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->refine_bisection(max_len * max_len);
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
}

void MainWindow::on_actionBench_memory_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->bench_memory();
  QApplication::restoreOverrideCursor();
}

void MainWindow::on_actionBench_construction_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->bench_construction();
  QApplication::restoreOverrideCursor();
}

void MainWindow::on_actionBench_intersections_vs_nbt_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->bench_intersections_vs_nbt();
  QApplication::restoreOverrideCursor();
}

void MainWindow::on_actionBench_distances_vs_nbt_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->bench_distances_vs_nbt();
  QApplication::restoreOverrideCursor();
}

void MainWindow::on_actionRefine_loop_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->refine_loop();
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
}

void MainWindow::on_actionSave_snapshot_triggered()
{
  return;
}
void MainWindow::on_actionCopy_snapshot_triggered()
{
  // copy snapshot to clipboard
  QApplication::setOverrideCursor(Qt::WaitCursor);
  QClipboard *qb = QApplication::clipboard();
  m_pViewer->makeCurrent();
  m_pViewer->raise();
  QImage snapshot = m_pViewer->grabFramebuffer();
  qb->setImage(snapshot);
  QApplication::restoreOverrideCursor();
}




