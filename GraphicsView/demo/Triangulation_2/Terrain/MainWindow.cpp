
#include "MainWindow.h"

MainWindow::MainWindow(QWidget* parent): CGAL::Qt::DemosMainWindow(parent)
{
  setupUi(this);
  this->viewer->setScene(&scene);
  connectActions();
  this->addAboutDemo(":/cgal/help/about_Terrain.html");
  this->addAboutCGAL();

  this->addRecentFiles(this->menuFile, this->actionQuit);
  connect(this, SIGNAL(openRecentFile(QString)),
	  this, SLOT(open(QString)));
}


void
MainWindow::connectActions()
{
  QObject::connect(this->actionLoad_New_File, SIGNAL(triggered()), 
		   this, SLOT(open_file()));


  QObject::connect(this, SIGNAL(sceneChanged()), 
		   this->viewer, SLOT(sceneChanged()));

  QObject::connect(this->actionQuit, SIGNAL(triggered()), 
		   qApp, SLOT(quit()));
}

void
MainWindow::open_file()
{

  QString fileName = QFileDialog::getOpenFileName(this,
						  tr("Open Points and polylines file"),
						  "./data",
						  tr("pap files (*.pap)"));

  if(! fileName.isEmpty()){
    open(fileName);
  }
}


/* Yet another file format 

number of points
x y z

number of polyline vertices
x y z

number of polyline vertices
x y z

*/

void
MainWindow::open(const QString& fileName)
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  scene.terrain.clear();
  std::ifstream ifs(qPrintable(fileName));
  int n;
  ifs >> n;
  Point_3 p;
  std::vector<Point_3> points;
  points.resize(n);
  for(int i=0; i<n; i++){
    ifs >> p;
    points[i] = p;
  }
  scene.bbox = points[0].bbox();
  
  for(int i=0; i<n; i++){
    scene.bbox = scene.bbox + points[i].bbox();
  }
  scene.terrain.insert(points.begin(), points.end());

  while(ifs >> n) {
    Vertex_handle vh, wh;
    ifs >> p;
    vh = scene.terrain.insert(p);
    // read a polyline with n points
    for(int i=1; i < n; i++){
      ifs >> p;
      wh = scene.terrain.insert(p, vh->face());
      scene.terrain.insert_constraint(vh, wh);
      vh = wh;
    }
  }

  this->addToRecentFiles(fileName);
  QApplication::restoreOverrideCursor();
  emit (sceneChanged());
}


#include "MainWindow.moc"

