#include "MainWindow.h"
#include <QFileDialog>

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
  scene.points.resize(n);
  for(int i=0; i<n; i++){
    ifs >> p;
    scene.points[i] = p;
  }
  scene.bbox = scene.points[0].bbox();
  
  for(int i=0; i<n; i++){
    scene.bbox = scene.bbox + scene.points[i].bbox();
  }
  scene.polylines.clear();
  while(ifs >> n) {
    std::vector<Point_3> line;
    // read a polyline with n points
    for(int i=0; i < n; i++){
      ifs >> p;
      line.push_back(p);
    }
    scene.polylines.push_back(line);
  }

  scene.refresh();

  this->addToRecentFiles(fileName);
  QApplication::restoreOverrideCursor();
  emit (sceneChanged());
}


#include "MainWindow.moc"

