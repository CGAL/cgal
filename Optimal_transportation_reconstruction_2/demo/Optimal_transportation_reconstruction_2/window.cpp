// STL
#include <fstream>

// Qt
#include <QtGui>
#include <QDialog>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <QClipboard>

// local
#include "window.h"
#include "ui_options.h"
#include "dialog_options.h"

MainWindow::MainWindow() :
QMainWindow(), Ui_MainWindow(),
maxNumRecentFiles(15), recentFileActs(15)
{
  setupUi(this);

  // init scene
  m_scene = new Scene;
  viewer->set_scene(m_scene);

  // options
  m_verbose = 0;
  m_mchoice = 0;
  m_relevance = 1.0;
  m_use_flip = true;
  m_percent  = 100.0;
  m_relocation = 2;

  // accepts drop events
  setAcceptDrops(true);

  // Handling actions
  addRecentFiles(menuFile, actionQuit);
  connect(actionQuit, SIGNAL(triggered()), this, SLOT(close()));
  connect(this, SIGNAL(openRecentFile(QString)), this, SLOT(open(QString)));
  Q_EMIT update();
  viewer->update();
}

MainWindow::~MainWindow()
{
  delete m_scene;
}

void MainWindow::addToRecentFiles(QString fileName)
{
  QSettings settings;
  QStringList files = settings.value("recentFileList").toStringList();
  files.removeAll(fileName);
  files.prepend(fileName);
  while (files.size() > (int)maxNumRecentFiles)
    files.removeLast();
  settings.setValue("recentFileList", files);
  updateRecentFileActions();
}

void MainWindow::dropEvent(QDropEvent *event)
{
  Q_FOREACH(QUrl url, event->mimeData()->urls())
  {
    QString filename = url.toLocalFile();
    if (!filename.isEmpty())
    {
      QTextStream(stderr) << QString("dropEvent(\"%1\")\n").arg(filename);
      open(filename);
    }
  }
  event->acceptProposedAction();
}

void MainWindow::closeEvent(QCloseEvent *event)
{
  event->accept();
}

void MainWindow::dragEnterEvent(QDragEnterEvent *event)
{
  if (event->mimeData()->hasFormat("text/uri-list"))
    event->acceptProposedAction();
}

void MainWindow::openRecentFile_aux()
{
  QAction* action = qobject_cast<QAction*>(sender());
  if (action)
    Q_EMIT openRecentFile(action->data().toString());
}

void MainWindow::updateRecentFileActions()
{
  QSettings settings;
  QStringList files = settings.value("recentFileList").toStringList();

  int numRecentFiles = qMin(files.size(), (int)maxNumberOfRecentFiles());
  for (int i = 0; i < numRecentFiles; ++i) {
    QString strippedName = QFileInfo(files[i]).fileName();
    QString text = tr("&%1 %2").arg(i).arg(strippedName);
    recentFileActs[i]->setText(text);
    recentFileActs[i]->setData(files[i]);
    recentFileActs[i]->setVisible(true);
  }
  for (unsigned int j = numRecentFiles; j < maxNumberOfRecentFiles(); ++j)
    recentFileActs[j]->setVisible(false);

  recentFilesSeparator->setVisible(numRecentFiles > 0);
}

void MainWindow::addRecentFiles(QMenu* menu, QAction* insertBeforeAction)
{
  if (insertBeforeAction)
    recentFilesSeparator = menu->insertSeparator(insertBeforeAction);
  else
    recentFilesSeparator = menu->addSeparator();
  recentFilesSeparator->setVisible(false);

  for (unsigned int i = 0; i < maxNumberOfRecentFiles(); ++i) {
    recentFileActs[i] = new QAction(this);
    recentFileActs[i]->setVisible(false);
    connect(recentFileActs[i], SIGNAL(triggered()), this, SLOT(openRecentFile_aux()));
    if (insertBeforeAction)
      menu->insertAction(insertBeforeAction, recentFileActs[i]);
    else
      menu->addAction(recentFileActs[i]);
  }
  updateRecentFileActions();
}

void MainWindow::open(const QString& filename)
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_scene->load(filename, this);
  QApplication::restoreOverrideCursor();
  addToRecentFiles(filename);
  on_actionRecenter_triggered();
  update();
}

void MainWindow::save(const QString& filename)
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_scene->save(filename);
  QApplication::restoreOverrideCursor();
}

void MainWindow::update()
{
  viewer->repaint();
}

void MainWindow::on_actionClear_triggered()
{
  m_scene->clear();
  update();
}

void MainWindow::on_actionLoadPoints_triggered()
{
  QString fileName =
      QFileDialog::getOpenFileName(this, tr("Open point set"), ".");
  if (fileName.isEmpty()) return;
  open(fileName);
}

void MainWindow::on_actionSave_triggered()
{
  QString filename =
    QFileDialog::getSaveFileName(this, tr("Save point set"), ".", tr("*.xy"));
  if (filename.isEmpty()) return;
  save(filename);
}


void MainWindow::on_actionInsertPoint_toggled()
{
  viewer->toggle_insert_points();
  update();
}

void MainWindow::on_actionRecenter_triggered()
{
  double center_x, center_y, scale;
  m_scene->compute_bbox(center_x, center_y, scale);
  viewer->set_camera(center_x, center_y, 1. / scale);
  update();
}


///////////////////////////
// PREDEFINED POINT SETS //
///////////////////////////

void MainWindow::on_actionCircle_triggered()
{
  bool ok;
  int density = QInputDialog::getInt(
      this, tr("Density"), tr("Density:"), 300, 10, 10000, 1, &ok);
  if (!ok) return;

  float x = QInputDialog::getDouble(
      this, tr("Center x"), tr("Center x:"), 0.0, 0.0, 10, 1, &ok);
  if (!ok) return;

  float y = QInputDialog::getDouble(
      this, tr("Center y"), tr("Center y:"), 0.0, 0.0, 10, 1, &ok);
  if (!ok) return;

  float radius = QInputDialog::getDouble(
      this, tr("Radius"), tr("Radius:"), 0.5, 0.1, 10, 1, &ok);
  if (!ok) return;

  m_scene->append_predefined_circle(density, x, y, radius);
  update();
}

void MainWindow::on_actionHalf_circle_triggered()
{
  bool ok;
  unsigned int density = QInputDialog::getInt(
      this, tr("Density"), tr("Density:"), 150, 10, 10000, 1, &ok);
  if (!ok) return;

  m_scene->append_predefined_half_circle(density);
  update();
}

void MainWindow::on_actionSpiral_triggered()
{
  bool ok;

  int loops = QInputDialog::getInt(
      this, tr("Loops"), tr("Number:"), 3, 1, 10000, 1, &ok);
  if (!ok) return;

  int density = QInputDialog::getInt(
      this, tr("Density"), tr("Density:"), 10, 1, 10000, 1, &ok);
  if (!ok) return;

  m_scene->append_predefined_spiral(loops,density);
  update();
}


void MainWindow::on_actionLine_triggered()
{
  bool ok;
  unsigned int density = QInputDialog::getInt(
      this, tr("Density"), tr("Density:"), 50, 1, 10000, 1, &ok);
  if (!ok) return;

  m_scene->append_predefined_line(density);
  update();
}

void MainWindow::on_actionBox_triggered()
{
  bool ok;
  int density = QInputDialog::getInt(
      this, tr("Density"), tr("Density:"), 100, 1, 10000, 1, &ok);
  if (!ok) return;

  float x = QInputDialog::getDouble(
      this, tr("x"), tr("x:"), 0.0, 0.0, 10.0, 1, &ok);
  if (!ok) return;

  float y = QInputDialog::getDouble(
      this, tr("y"), tr("y:"), 0.0, 0.0, 10.0, 1, &ok);
  if (!ok) return;

  float sx = QInputDialog::getDouble(
      this, tr("size x"), tr("size x:"), 1.0, 0.1, 10.0, 1, &ok);
  if (!ok) return;

  float sy = QInputDialog::getDouble(
      this, tr("size y"), tr("size y:"), 1.0, 0.1, 10.0, 1, &ok);
  if (!ok) return;

  m_scene->append_predefined_box(density, x, y, sx, sy);
  update();
}

void MainWindow::on_actionBoxes_triggered()
{
  bool ok;
  unsigned int density = QInputDialog::getInt(
      this, tr("Density"), tr("Density:"), 10, 1, 10000, 1, &ok);
  if (!ok) return;

  m_scene->append_predefined_boxes(density);
  update();
}

void MainWindow::on_actionParallel_lines_triggered()
{
  bool ok;

  int lines = QInputDialog::getInt(
      this, tr("Lines"), tr("Number:"), 3, 1, 10000, 1, &ok);
  if(!ok) return;

  float space = QInputDialog::getDouble(
      this, tr("Space"), tr("Space:"), 0.2, 0.1, 10, 1, &ok);
  if (!ok) return;

  int density = QInputDialog::getInt(
      this, tr("Density"), tr("Density:"), 10, 1, 10000, 1, &ok);
  if(!ok) return;

  m_scene->append_predefined_parallel_lines(lines, space, density);
  update();
}

void MainWindow::on_actionBox_with_boundaries_triggered()
{
  bool ok;
  unsigned int density = QInputDialog::getInt(
      this, tr("Density"), tr("Density:"), 10, 1, 10000, 1, &ok);
  if (!ok) return;

  m_scene->append_predefined_box_with_boundaries(density);
  update();
}

void MainWindow::on_actionBox_with_missing_corners_triggered()
{
  bool ok;
  unsigned int density = QInputDialog::getInt(
      this, tr("Density"), tr("Density:"), 10, 1, 10000, 1, &ok);
  if (!ok) return;

  m_scene->append_predefined_box_with_missing_corners(density);
  update();
}

void MainWindow::on_actionStar_triggered()
{
  bool ok;
  int nb_branches = QInputDialog::getInt(
      this, tr("Branches"), tr("Branches:"), 20, 2, 10000, 1, &ok);
  if(!ok) return;

  int density = QInputDialog::getInt(
      this, tr("Density"), tr("Density:"), 100, 10, 10000, 1, &ok);
  if(!ok) return;

  m_scene->append_star(nb_branches,density);
  update();
}


void MainWindow::on_actionStair_triggered()
{
  bool ok;
  unsigned int density = QInputDialog::getInt(
      this, tr("Density"), tr("Density:"), 30, 2, 10000, 1, &ok);
  if (!ok) return;

  m_scene->append_predefined_stair(density);
  on_actionRecenter_triggered();
  update();
}

void MainWindow::on_actionSkyline_triggered()
{
  bool ok;
  unsigned int density = QInputDialog::getInt(
      this, tr("Density"), tr("Density:"), 10, 1, 10000, 1, &ok);
  if (!ok) return;

  m_scene->append_predefined_skyline(density);
  on_actionRecenter_triggered();
  update();
}

void MainWindow::on_actionIncreasingly_sharp_angles_triggered()
{
  bool ok;
  unsigned int density = QInputDialog::getInt(
      this, tr("Density"), tr("Density:"), 100, 1, 10000, 1, &ok);
  if (!ok) return;

  double min_angle = QInputDialog::getDouble(
      this, tr("Min Angle"), tr("Min Angle:"), 1.0, 1.0, 360.0, 1.0, &ok);
  if (!ok) return;

  m_scene->append_predefined_increasingly_sharp_angles(density, min_angle);
  on_actionRecenter_triggered();
  update();
}

void MainWindow::on_actionWidely_variable_sampling_triggered()
{
  bool ok;
  float d1 = QInputDialog::getDouble(
      this, tr("Delta-angle"), tr("Delta-angle:"), 1, 0.01, 20.0, 2, &ok);
  if (!ok) return;

  float d2 = QInputDialog::getDouble(
      this, tr("Delta-angle"), tr("Delta-angle:"), 10, 0.01, 30.0, 2, &ok);
  if (!ok) return;

  m_scene->append_widely_variable_sampling(d1, d2);
  update();
}

void MainWindow::on_actionNoise_triggered()
{
  bool ok;
  double noise = QInputDialog::getDouble(
      this, tr("Noise"), tr("Amount:"), 0.003, 0.0, 1000.0, 8, &ok);
  if (!ok) return;

  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_scene->noise(noise);
  QApplication::restoreOverrideCursor();
  update();
}

void MainWindow::on_actionAdd_outliers_triggered()
{
  bool ok;
  unsigned int nb = (unsigned)
        QInputDialog::getInt(
            this, tr("Outliers"), tr("How many:"), 10, 1, 100000, 1, &ok);
  if (!ok) return;

  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_scene->add_outliers(nb);
  QApplication::restoreOverrideCursor();
  update();
}

void MainWindow::on_actionSubdivide_triggered()
{
  m_scene->subdivide();
  update();
}

void MainWindow::on_actionDecimate_triggered()
{
  bool ok;
  int percent = QInputDialog::getInt(
      this, tr("Decimate"), tr("Percentage:"), 50, 1, 100, 1, &ok);
  if (!ok) return;

  QApplication::setOverrideCursor(Qt::WaitCursor);
  const double percentage = double(percent) / 100.0;
  m_scene->decimate(percentage);
  QApplication::restoreOverrideCursor();
  update();
}


////////////////////
// RECONSTRUCTION //
////////////////////

void MainWindow::set_scene_options()
{
  m_scene->set_options(m_verbose, m_mchoice, m_use_flip,
      m_relocation, m_relevance, percentage());
}

void MainWindow::on_actionReconstruction_reinit_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  set_scene_options();
  m_scene->init_reconstruction(percentage());
  QApplication::restoreOverrideCursor();
  update();
}

void MainWindow::on_actionReconstruction_one_step_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  set_scene_options();
  m_scene->reconstruct(1);
  QApplication::restoreOverrideCursor();
  update();
}

void MainWindow::on_actionReconstruction_10_steps_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  set_scene_options();
  m_scene->reconstruct(10);
  QApplication::restoreOverrideCursor();
  update();
}

void MainWindow::on_actionReconstruction_100_steps_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  set_scene_options();
  m_scene->reconstruct(100);
  QApplication::restoreOverrideCursor();
  update();
}

void MainWindow::on_actionReconstruction_1000_steps_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  set_scene_options();
  m_scene->reconstruct(1000);
  QApplication::restoreOverrideCursor();
  update();
}

void MainWindow::on_actionReconstruction_until_triggered()
{
  bool ok;
  int nb_points = QInputDialog::getInt(
      this, tr("Number of Points"), tr("Nb:"), 4, 1, 1000000, 1, &ok);
  if (!ok) return;

  set_scene_options();

  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_scene->reconstruct_until(nb_points);
  QApplication::restoreOverrideCursor();

  update();
}

void MainWindow::on_actionReconstruction_Wasserstein_tolerance_triggered()
{
  bool ok;
  double tolerance = QInputDialog::getDouble(
      this, tr("Wasserstein tolerance"), tr("Tolerance:"), 0.01, 0., 1000000., 5, &ok);
  if (!ok) return;

  set_scene_options();

  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_scene->reconstruct_wasserstein_tolerance (tolerance);
  QApplication::restoreOverrideCursor();

  update();
}

void MainWindow::on_actionRelocate_vertices_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_scene->relocate_all_points();
  QApplication::restoreOverrideCursor();
  update();
}

void MainWindow::on_actionOutput_console_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_scene->output_console();
  QApplication::restoreOverrideCursor();
  update();
}



//////////
// VIEW //
//////////

void MainWindow::on_actionView_points_toggled()
{
  viewer->toggle_view_points();
  update();
}

void MainWindow::on_actionView_tolerance_toggled()
{
  viewer->toggle_view_tolerance();
  update();
}

void MainWindow::on_actionView_vertices_toggled()
{
  viewer->toggle_view_vertices();
  update();
}

void MainWindow::on_actionView_edges_toggled()
{
  viewer->toggle_view_edges();
  update();
}

void MainWindow::on_actionView_ghost_toggled()
{
  viewer->toggle_view_ghost_edges();
  update();
}

void MainWindow::on_actionView_edge_cost_toggled()
{
  viewer->toggle_view_edge_cost();
  update();
}

void MainWindow::on_actionView_edge_priority_toggled()
{
  viewer->toggle_view_edge_priority();
  update();
}

void MainWindow::on_actionView_bins_toggled()
{
  viewer->toggle_view_bins();
  update();
}

void MainWindow::on_actionView_foot_points_toggled()
{
  viewer->toggle_view_foot_points();
  update();
}

void MainWindow::on_actionView_relocation_toggled()
{
  viewer->toggle_view_relocation();
  update();
}

void MainWindow::on_actionView_relevance_toggled()
{
  viewer->toggle_view_edge_relevance();
  update();
}

void MainWindow::on_actionSet_options_triggered()
{
  Dialog_options dlg;
  dlg.set_all_ranges();
  dlg.set_verbose(m_verbose);
  dlg.set_random_sample_size(m_mchoice);
  dlg.set_percent(m_percent);
  dlg.set_relocation(m_relocation);
  dlg.set_relevance(m_relevance);
  dlg.set_use_flip(m_use_flip);
  dlg.set_multiple_choice_checkbox(m_mchoice != 0);
  dlg.set_line_thickness(viewer->line_thickness());
  dlg.set_point_size(viewer->point_size());
  dlg.set_vertex_size(viewer->vertex_size());

  if (dlg.exec() == QDialog::Accepted)
  {
    m_verbose = dlg.get_verbose();
    m_mchoice = dlg.get_mchoice();
    m_percent = dlg.get_percent();
    m_relocation = dlg.get_relocation();
    m_relevance  = dlg.get_relevance();
    m_use_flip   = dlg.get_use_flip();

    set_scene_options();
    viewer->line_thickness() = dlg.get_line_thickness();
    viewer->point_size()  = dlg.get_point_size();
    viewer->vertex_size() = dlg.get_vertex_size();
    update();
  }
}
