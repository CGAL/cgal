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
  m_pViewer = ui->viewer;

  // do not save the state of the viewer (anoying)
  m_pViewer->setStateFileName(QString::null);

  // accept drop events
  setAcceptDrops(true);

  // setup scene
  m_pScene = new Scene;
  m_pViewer->setScene(m_pScene);

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
  const Scene::Bbox bbox = m_pScene->bbox();
  const double xmin = bbox.xmin;
  const double ymin = bbox.ymin;
  const double zmin = bbox.zmin;
  const double xmax = bbox.xmax;
  const double ymax = bbox.ymax;
  const double zmax = bbox.zmax;
  qglviewer::Vec 
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
