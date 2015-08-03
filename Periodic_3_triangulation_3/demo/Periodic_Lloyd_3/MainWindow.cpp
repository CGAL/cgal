
#include "MainWindow.h"
#include <ctime>

MainWindow::MainWindow(QWidget* parent): CGAL::Qt::DemosMainWindow(parent)
{
  setupUi(this);
  this->viewer->setScene(&scene);

  process = new QProcess(this);

  connectActions();
  this->addAboutDemo(":/cgal/help/about_Periodic_Lloyd_3.html");
  this->addAboutCGAL();

  scene.eight_copies=false;
  scene.two_dimensional=false;

  qtimer = new QTimer(this);
  connect(qtimer, SIGNAL(timeout()), this, SLOT(lloydStep()));
}

void
MainWindow::connectActions()
{
  QObject::connect(this->actionNew_Point_Set, SIGNAL(triggered()), 
		   this, SLOT(newPointSet()));

  QObject::connect(this->actionLoad_points, SIGNAL(triggered()), 
		   this, SLOT(loadPoints()));

  QObject::connect(this->actionSave_points, SIGNAL(triggered()), 
		   this, SLOT(savePoints()));

  QObject::connect(this->speedSlider, SIGNAL(valueChanged(int)), 
		   this, SLOT(speedChanged(int)));

  QObject::connect(this->viewer, SIGNAL(valueChanged(int)),
                   this, SLOT(speedChanged(int)));

  QObject::connect(this, SIGNAL(sceneChanged()), 
		   this->viewer, SLOT(sceneChanged()));
  
  QObject::connect(this->actionStep, SIGNAL(triggered()),
                   this, SLOT(lloydStep()));

  QObject::connect(this->actionPlay, SIGNAL(toggled(bool)),
                   this, SLOT(togglePause(bool)));

  QObject::connect(this->actionShow_8_Copies, SIGNAL(toggled(bool)),
                   this, SLOT(toggle8Copies(bool)));

  QObject::connect(this->action2D_version, SIGNAL(toggled(bool)),
                   this, SLOT(toggle2D(bool)));

  QObject::connect(this->actionQuit, SIGNAL(triggered()), 
		   qApp, SLOT(quit()));

  QObject::connect(this->actionDemo_Help, SIGNAL(triggered()),
                   this, SLOT(help()));
}

void
MainWindow::togglePause(bool p)
{
  if (p) {
    int speed = (100-(speedSlider->value()))*100;
    qtimer->start(speed);
  }
  else qtimer->stop();
}

void
MainWindow::toggle8Copies(bool on)
{
  scene.eight_copies = on;
  emit(sceneChanged());
}

void
MainWindow::toggle2D(bool on)
{
  scene.two_dimensional = on;
  if (on) newPointSet();
  else emit(sceneChanged());
}

void
MainWindow::newPointSet()
{
  bool ok;

  int numberOfPoints = QInputDialog::getInt(this,
      "Periodic Lloyd", "Number of points: ", 100, 0, 2147483647, 1, &ok );
  
  if (ok) newPoints(numberOfPoints);
}

void
MainWindow::loadPoints()
{
  QString fileName = QFileDialog
    ::getOpenFileName(this, tr("Open point set"),
	".", tr("All files (*)"));
  if(fileName.isEmpty()) return;
  
  std::ifstream ifs(fileName.toLatin1().data() );
  scene.points.clear();
  Iso_cuboid_3 dom;
  ifs >> dom;
  std::copy(std::istream_iterator<Point_3>(ifs), 
      std::istream_iterator<Point_3>(),
      std::back_inserter(scene.points));

  scene.periodic_triangulation.set_domain(dom);
  scene.periodic_triangulation.insert(scene.points.begin(), scene.points.end());

  FT cx(0),cy(0),cz(0);
  for (int i=0 ; i<8 ; i++) {
    cx += dom[i].x();
    cy += dom[i].y();
    cy += dom[i].y();
  }
  qglviewer::Vec center(cx/8.,cy/8.,cz/8.);
  viewer->setSceneCenter(center);
  viewer->setSceneRadius(std::sqrt(
	  ((dom.xmax()-dom.xmin())*(dom.xmax()-dom.xmin()))
	  + ((dom.xmax()-dom.xmin())*(dom.xmax()-dom.xmin()))
	  + ((dom.xmax()-dom.xmin())*(dom.xmax()-dom.xmin()))));

  speedSlider->setRange(0,100);
  speedSlider->setSliderPosition(100);

  emit (sceneChanged()); 
}

void
MainWindow::savePoints()
{
  QString fileName = QFileDialog
    ::getSaveFileName(this, tr("Save point set"),
	".", tr("*.pts"));
  if(fileName.isEmpty()) return;
  
  std::ofstream ofs(fileName.toLatin1().data() );
  ofs << scene.periodic_triangulation.domain() << '\n';
  for (std::list<Point_3>::iterator pit = scene.points.begin() ;
       pit != scene.points.end() ; ++pit) ofs << *pit << '\n';
}

void MainWindow::lloydStep() {
  scene.lloyd_step();
  viewer->updateGL();
  viewer->changed();
  }

void 
MainWindow::speedChanged(int i)
{
  int speed = (100-i)*100;
  if (qtimer->isActive()) {
     qtimer->stop();
     qtimer->start(speed);
  }
}

void
MainWindow::newPoints(int n)
{
  scene.periodic_triangulation.clear();
  scene.points.clear();

  CGAL::Random rnd(std::time(NULL));
  CGAL::Random_points_in_cube_3<Point_3, Creator> in_cube(1,rnd);

  for (int i=0 ; i<n ; i++) 
    if (scene.two_dimensional) {
      Point_3 rdpt = *in_cube++;
      scene.points.push_back(Point_3(rdpt.x(),rdpt.y(),0.));
    } else 
      scene.points.push_back(*in_cube++);

  Iso_cuboid_3 dom(-1,-1,-1,1,1,1);
  scene.periodic_triangulation.set_domain(dom);
  scene.periodic_triangulation.insert(scene.points.begin(), scene.points.end());

  FT cx(0),cy(0),cz(0);
  for (int i=0 ; i<8 ; i++) {
    cx += dom[i].x();
    cy += dom[i].y();
    cy += dom[i].y();
  }
  qglviewer::Vec center(cx/8.,cy/8.,cz/8.);
  viewer->setSceneCenter(center);
  viewer->setSceneRadius(std::sqrt(
	  ((dom.xmax()-dom.xmin())*(dom.xmax()-dom.xmin()))
	  + ((dom.xmax()-dom.xmin())*(dom.xmax()-dom.xmin()))
	  + ((dom.xmax()-dom.xmin())*(dom.xmax()-dom.xmin()))));

  speedSlider->setRange(0,100);
  speedSlider->setSliderPosition(100);

  emit (sceneChanged()); 
}

void MainWindow::help() {
  QString app = QLibraryInfo::location(QLibraryInfo::BinariesPath)
    + QDir::separator();
#if !defined(Q_OS_MAC)
  app += QString("assistant");
#else
  app += QString("Assistant.app/Contents/MacOS/Assistant");
#endif

  QStringList args;
  QString help_path = QCoreApplication::applicationDirPath() 
    + QDir::separator()
    + QString("./Periodic_Lloyd_3.qhc");
  args << QString("-collectionFile") << help_path;
  process->start(app, args);
  if (!process->waitForStarted()) {
    QMessageBox::critical(this, tr("Remote Control"),
      tr("Could not start Qt Assistant from %1.").arg(app));
  }
}


