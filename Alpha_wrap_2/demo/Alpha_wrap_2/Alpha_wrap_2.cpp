// #define CGAL_AW2_DEBUG_PP

#include "Alpha_wrap_2_options.h"
#include "Screenshot_options.h"
#include "scene.h"
#include "ui_Alpha_wrap_2.h"

#include <CGAL/Real_timer.h>

#include <CGAL/Qt/utility.h>
#include <CGAL/Qt/DemosMainWindow.h>

#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <QWidget>

#include <string>

class MainWindow
  : public CGAL::Qt::DemosMainWindow,
    public Ui_MainWindow
{
  Q_OBJECT

private:
  Scene* m_scene;

  unsigned int maxNumRecentFiles;
  QAction* recentFilesSeparator;
  QVector<QAction*> recentFileActs;

  QString currentScreenshotDir = "./";
  QString currentScreenshotFilename = "screenshot";
  int currentScreenshotNumber = 0;

  CGAL::Real_timer timer;

public:
  MainWindow();
  ~MainWindow();

  void update();
  void update_viewer_camera();

protected:
  void addRecentFiles(QMenu* menu, QAction* insertBefore = 0);
  unsigned int maxNumberOfRecentFiles() const;

protected Q_SLOTS:
  // drag & drop
  void dropEvent(QDropEvent *event);
  void closeEvent(QCloseEvent *event);
  void dragEnterEvent(QDragEnterEvent *event);

  // recent files
  void openRecentFile_aux();
  void updateRecentFileActions();
  void addToRecentFiles(QString fileName);

  // io
  void open(QString file);

  void setStepByStepAction(bool b);
  void update_stats(const Alpha_wrapper& wrapper);
  void resetAction();

public Q_SLOTS:
  void on_actionLoad_triggered();
  void on_actionSaveInput_triggered();
  void on_actionClear_triggered();
  void on_actionScreenshot_triggered();
  void on_actionScreenshotSettings_triggered();

  void on_actionExplodeInput_triggered();
  void on_actionSmoothInput_triggered();
  void on_actionAlphaWrap_triggered();
  void on_actionNextStep_triggered();
  void on_actionNext10Steps_triggered();
  void on_actionTerminate_triggered();

  void on_actionViewInput_triggered();
  void on_actionViewAlphaWrap_triggered();
  void on_actionViewInsideOutside_triggered();
  void on_actionViewTriangulationEdge_triggered();
  void on_actionViewVoronoiDiagram_triggered();
  void on_actionViewEmptyCircles_triggered();
  void on_actionViewSteinerPoint_triggered();
  void on_actionViewNextGate_triggered();
  void on_actionViewNextGateEmptyCircle_triggered();
  void on_actionViewStats_triggered();

Q_SIGNALS:
  void openRecentFile(QString filename);
};

MainWindow::MainWindow()
  : DemosMainWindow(),
    Ui_MainWindow(),
    maxNumRecentFiles(15),
    recentFileActs(15)
{
  setupUi(this);

  // init scene
  m_scene = new Scene;
  viewer->set_scene(m_scene);

  // accepts drop events
  setAcceptDrops(true);
  connect(actionQuit, SIGNAL(triggered()), this, SLOT(close()));
  addRecentFiles(menuInput, actionQuit);
  connect(this, SIGNAL(openRecentFile(QString)), this, SLOT(open(QString)));

  stats->setVisible(false);
}

MainWindow::~MainWindow()
{
  delete m_scene;
}

void MainWindow::update()
{
  viewer->repaint();
}

void MainWindow::update_viewer_camera()
{
  CGAL::Bbox_2 bbox = m_scene->get_bbox();
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()));
  viewer->set_camera(0.5 * (bbox.xmin() + bbox.xmax()),
                     0.5 * (bbox.ymin() + bbox.ymax()),
                     1.5 / diag_length);
}

void MainWindow::addRecentFiles(QMenu* menu, QAction* insertBeforeAction)
{
  if(insertBeforeAction) {
    recentFilesSeparator = menu->insertSeparator(insertBeforeAction);
  } else {
    recentFilesSeparator = menu->addSeparator();
  }
  recentFilesSeparator->setVisible(false);

  for(unsigned int i=0; i<maxNumberOfRecentFiles(); ++i) {
    recentFileActs[i] = new QAction(this);
    recentFileActs[i]->setVisible(false);
    connect(recentFileActs[i], SIGNAL(triggered()), this, SLOT(openRecentFile_aux()));
    if(insertBeforeAction)
      menu->insertAction(insertBeforeAction, recentFileActs[i]);
    else
      menu->addAction(recentFileActs[i]);
  }
  updateRecentFileActions();
}

unsigned int MainWindow::maxNumberOfRecentFiles() const { return maxNumRecentFiles; }

// drag & drop
void MainWindow::dropEvent(QDropEvent *event)
{
  Q_FOREACH(QUrl url, event->mimeData()->urls())
  {
    QString filename = url.toLocalFile();
    if(!filename.isEmpty())
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
  if(event->mimeData()->hasFormat("text/uri-list"))
    event->acceptProposedAction();
}

// recent files
void MainWindow::openRecentFile_aux()
{
  QAction* action = qobject_cast<QAction*>(sender());
  if(action)
    Q_EMIT openRecentFile(action->data().toString());
}

void MainWindow::updateRecentFileActions()
{
  QSettings settings;
  QStringList files = settings.value("recentFileList").toStringList();

  int numRecentFiles = qMin(files.size(), static_cast<int>(maxNumberOfRecentFiles()));
  for(int i=0; i<numRecentFiles; ++i) {
    QString strippedName = QFileInfo(files[i]).fileName();
    QString text = tr("&%1 %2").arg(i).arg(strippedName);
    recentFileActs[i]->setText(text);
    recentFileActs[i]->setData(files[i]);
    recentFileActs[i]->setVisible(true);
  }
  for(unsigned int j = numRecentFiles; j < maxNumberOfRecentFiles(); ++j)
    recentFileActs[j]->setVisible(false);

  recentFilesSeparator->setVisible(numRecentFiles > 0);
}

void MainWindow::addToRecentFiles(QString fileName)
{
    QSettings settings;
    QStringList files = settings.value("recentFileList").toStringList();
    files.removeAll(fileName);
    files.prepend(fileName);
    while (files.size() > static_cast<int>(maxNumRecentFiles))
      files.removeLast();
    settings.setValue("recentFileList", files);
    updateRecentFileActions();
}

// io
void MainWindow::open(QString filename)
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  bool success = m_scene->load(filename);
  QApplication::restoreOverrideCursor();
  if(!success)
    return;
  this->actionExplodeInput->setEnabled(true);
  this->actionSmoothInput->setEnabled(true);
  this->actionAlphaWrap->setEnabled(true);
  addToRecentFiles(filename);
  update_viewer_camera();
  update();
  resetAction();
}

void MainWindow::setStepByStepAction(bool b)
{
  this->actionNextStep->setEnabled(b);
  this->actionNext10Steps->setEnabled(b);
  this->actionTerminate->setEnabled(b);
}

void MainWindow::update_stats(const Alpha_wrapper& wrapper)
{
  if(!wrapper.queue().empty())
    return;

  actionViewStats->setEnabled(true);
  stats->setVisible(actionViewStats->isChecked());
  // output # of cmp is # of edges for now
  complexity->setText(QString::fromStdString("Alpha Wrap complexity : "+ std::to_string(m_scene->get_alpha_wrap().size()) +" edges"));
  time->setText(QString::fromStdString("Total execution time : "+ std::to_string(timer.time()) +" s"));
  double percentage_nb_proj = (m_scene->get_nb_proj()*100.)/m_scene->get_nb_steiner();
  double percentage_nb_inter = 100. - percentage_nb_proj;
  steiner->setText(QString::fromStdString("Steiner points : "+ std::to_string(percentage_nb_proj) +"% projection - "+ std::to_string(percentage_nb_inter) + "% intersection"));
  gate->setText(QString::fromStdString("Number of gates traversed : "+ std::to_string(m_scene->get_nb_gate_traversed())));
}

void MainWindow::resetAction()
{
  timer.reset();
  setStepByStepAction(false);
  actionViewStats->setEnabled(false);
  stats->setVisible(false);
  m_scene->reset_stats();
}

void MainWindow::on_actionLoad_triggered()
{
  QSettings settings;
  QString directory =
      settings.value("Open PSLG", QDir::current().dirName())
          .toString();
  QString fileName = QFileDialog::getOpenFileName(this, tr("Open PSLG"), directory, tr("DAT files (*.dat)"));
  if(!fileName.isEmpty()) {
    open(fileName);
  }
}

void MainWindow::on_actionSaveInput_triggered()
{
  QSettings settings;
  QString directory =
      settings.value("Save PSLG", QDir::current().dirName())
          .toString();
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save PSLG"), directory, tr("DAT files (*.dat)"));
  if(!fileName.isEmpty()) {
    AW2::IO::write_output_polylines_file(fileName.toStdString(), m_scene->get_input());
  }
}

void MainWindow::on_actionClear_triggered()
{
  m_scene->clear();
  this->actionExplodeInput->setEnabled(false);
  this->actionSmoothInput->setEnabled(false);
  this->actionAlphaWrap->setEnabled(false);
  resetAction();
  update();
}

void MainWindow::on_actionScreenshot_triggered()
{
  QPixmap pixmap(viewer->size());
  viewer->render(&pixmap, QPoint(), QRegion(viewer->rect()));
  QString filepath = currentScreenshotDir + QDir::separator() + currentScreenshotFilename + "_" +  QString::number(currentScreenshotNumber) + ".png";
  pixmap.save(filepath);
  currentScreenshotNumber++;
}

void MainWindow::on_actionScreenshotSettings_triggered()
{
  ScreenshotOptions dialog_box_screenshot(this, m_scene->get_screenshot_folder(),
                                          m_scene->get_screenshot_filename(),
                                          m_scene->get_screenshot_number());
  if(dialog_box_screenshot.exec() != QDialog::Accepted) return;
  currentScreenshotDir = dialog_box_screenshot.get_current_screenshot_dir();
  currentScreenshotFilename = dialog_box_screenshot.get_current_screenshot_filename();
  currentScreenshotNumber = dialog_box_screenshot.get_current_screenshot_numbering_start();
  m_scene->set_screenshot_folder(currentScreenshotDir);
  m_scene->set_screenshot_filename(currentScreenshotFilename);
  m_scene->set_screenshot_number(currentScreenshotNumber);
}


void MainWindow::on_actionExplodeInput_triggered()
{
  m_scene->get_alpha_wrap().clear();
  m_scene->get_wrapper().clear();
  m_scene->explode_input();
  resetAction();
  update_viewer_camera();
  update();
}

void MainWindow::on_actionSmoothInput_triggered()
{
  m_scene->get_alpha_wrap().clear();
  m_scene->get_wrapper().clear();
  m_scene->smooth_input();
  resetAction();
  update();
}

void MainWindow::on_actionAlphaWrap_triggered()
{
  AlphaWrapOptions dialog_box_alpha_wrap(this, m_scene->get_alpha(), m_scene->get_offset(),
                                         m_scene->get_is_alpha_relative(), m_scene->get_is_offset_relative(),
                                         m_scene->get_is_step_by_step());
  if(dialog_box_alpha_wrap.exec() != QDialog::Accepted)
    return;

  resetAction();
  m_scene->set_alpha(dialog_box_alpha_wrap.get_alpha_value());
  m_scene->set_offset(dialog_box_alpha_wrap.get_offset_value());
  m_scene->set_is_alpha_relative(dialog_box_alpha_wrap.is_alpha_relative());
  m_scene->set_is_offset_relative(dialog_box_alpha_wrap.is_offset_relative());
  m_scene->set_is_step_by_step(dialog_box_alpha_wrap.is_step_by_step());

  // Get the alpha wrap options
  double alpha_value = dialog_box_alpha_wrap.get_alpha_value();
  double offset_value = dialog_box_alpha_wrap.get_offset_value();
  if(dialog_box_alpha_wrap.is_alpha_relative() ||
     dialog_box_alpha_wrap.is_offset_relative()) {
    double longest_diag_length = m_scene->get_diagonal_bbox();
    if(dialog_box_alpha_wrap.is_alpha_relative())
      alpha_value = longest_diag_length / alpha_value;
    if(dialog_box_alpha_wrap.is_offset_relative())
      offset_value = longest_diag_length / offset_value;
  }

  timer.start();
  m_scene->init_alpha_data_structure(alpha_value, offset_value);
  if(!dialog_box_alpha_wrap.is_step_by_step()) {
    m_scene->alpha_flood_fill();
  }
  timer.stop();
  m_scene->extract_pslg_2_soup();
  setStepByStepAction(dialog_box_alpha_wrap.is_step_by_step());
  update_viewer_camera();
  update();
  update_stats(m_scene->get_wrapper());
}

void MainWindow::on_actionNextStep_triggered()
{
  Alpha_wrapper& wrapper = m_scene->get_wrapper();
  if(wrapper.queue().empty())
    return;
  timer.start();
  m_scene->alpha_flood_fill(1);
  timer.stop();
  m_scene->extract_pslg_2_soup();
  update();
  update_stats(wrapper);
}

void MainWindow::on_actionNext10Steps_triggered()
{
  Alpha_wrapper& wrapper = m_scene->get_wrapper();
  if(wrapper.queue().empty())
    return;
  timer.start();
  m_scene->alpha_flood_fill(10);
  timer.stop();
  m_scene->extract_pslg_2_soup();
  update();
  update_stats(wrapper);
}

void MainWindow::on_actionTerminate_triggered()
{
  Alpha_wrapper& wrapper = m_scene->get_wrapper();
  if(wrapper.queue().empty())
    return;
  timer.start();
  m_scene->alpha_flood_fill();
  timer.stop();
  m_scene->extract_pslg_2_soup();
  update();
  update_stats(wrapper);
}

void MainWindow::on_actionViewInput_triggered() {
  m_scene->toggle_view_input();
  update();
}

void MainWindow::on_actionViewAlphaWrap_triggered() {
  m_scene->toggle_view_alpha_wrap();
  update();
}

void MainWindow::on_actionViewInsideOutside_triggered() {
  m_scene->toggle_view_dt2_inside_outside();
  update();
}

void MainWindow::on_actionViewTriangulationEdge_triggered() {
  m_scene->toggle_view_dt2_edge();
  update();
}

void MainWindow::on_actionViewVoronoiDiagram_triggered() {
  m_scene->toggle_view_voronoi();
  update();
}

void MainWindow::on_actionViewEmptyCircles_triggered() {
  m_scene->toggle_view_empty_alpha_pencils();
  update();
}

void MainWindow::on_actionViewSteinerPoint_triggered() {
  m_scene->toggle_view_steiner_point();
  update();
}

void MainWindow::on_actionViewNextGate_triggered() {
  m_scene->toggle_view_next_gate();
  update();
}

void MainWindow::on_actionViewNextGateEmptyCircle_triggered() {
  m_scene->toggle_view_next_gate_pencil();
  update();
}

void MainWindow::on_actionViewStats_triggered() {
  stats->setVisible(!stats->isVisible());
}

#include "Alpha_wrap_2.moc"
#include <CGAL/Qt/resources.h>

int main(int argc, char** argv)
{
  QApplication app(argc, argv);

  app.setApplicationName("2D Alpha Wrapping");

  // Import resources from libCGAL (Qt6).
  CGAL_QT_INIT_RESOURCES;

  MainWindow window;
  window.show();
  return app.exec();
}
