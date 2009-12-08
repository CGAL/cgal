
#include "MainWindow.h"

MainWindow::MainWindow(QWidget* parent): CGAL::Qt::DemosMainWindow(parent)
{
  setupUi(this);
  this->viewer->setScene(&scene);
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

  QObject::connect(this->speedSlider, SIGNAL(valueChanged(int)), 
		   this, SLOT(speedChanged(int)));

  QObject::connect(this->viewer, SIGNAL(valueChanged(int)),
                   this, SLOT(speedChanged(int)));

  QObject::connect(this, SIGNAL(sceneChanged()), 
		   this->viewer, SLOT(sceneChanged()));
  
  QObject::connect(this->actionStep, SIGNAL(triggered()),
                   this, SLOT(lloydStep()));

  QObject::connect(this->actionStop, SIGNAL(toggled(bool)),
                   this, SLOT(togglePause(bool)));

  QObject::connect(this->actionShow_8_Copies, SIGNAL(toggled(bool)),
                   this, SLOT(toggle8Copies(bool)));

  QObject::connect(this->action2D_version, SIGNAL(toggled(bool)),
                   this, SLOT(toggle2D(bool)));

  QObject::connect(this->actionQuit, SIGNAL(triggered()), 
		   qApp, SLOT(quit()));
}

void
MainWindow::togglePause(bool p)
{
  if (p)
    qtimer->stop();
  else {
    int speed = (100-(speedSlider->value()))*100;
    qtimer->start(speed);
  }
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
  int numberOfPoints = QInputDialog::getInteger(this,
      "Periodic Lloyd", "Number of points: ", 100, 0, 2147483647, 1, &ok );
  
  if (ok) newPoints(numberOfPoints);
}

void MainWindow::lloydStep() {
  scene.lloyd_step();
  viewer->updateGL();
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

  scene.periodic_triangulation.set_domain(Iso_cuboid_3(-1,-1,-1,1,1,1));
  scene.periodic_triangulation.insert(scene.points.begin(), scene.points.end());

  speedSlider->setRange(0,100);
  speedSlider->setSliderPosition(100);

  emit (sceneChanged()); 

}


#include "MainWindow.moc"

