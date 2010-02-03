
#include "MainWindow.h"

MainWindow::MainWindow(QWidget* parent): CGAL::Qt::DemosMainWindow(parent), nbcube(0)
{
  setupUi(this);
  this->viewer->setScene(&scene);
  connectActions();
  this->addAboutDemo(":/cgal/help/about_Combinatorial_map_3.html");
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

  QObject::connect(this->actionSubdivide, SIGNAL(triggered()), 
		   this, SLOT(subdivide()));

  QObject::connect(this->actionCreateCube, SIGNAL(triggered()), 
		   this, SLOT(createCube()));

  QObject::connect(this, SIGNAL(sceneChanged()), 
		   this->viewer, SLOT(sceneChanged()));


  QObject::connect(this->actionQuit, SIGNAL(triggered()), 
		   qApp, SLOT(quit()));
}

void
MainWindow::open_file()
{

  QString fileName = QFileDialog::getOpenFileName(this,
						  tr("Open File"),
						  "./off",
						  tr("off files (*.off)"));

  if(! fileName.isEmpty()){
    open(fileName);
  }
}


void
MainWindow::createCube()
{  
  make_cube(scene.map, Point_3(nbcube, nbcube, nbcube), 1);
  ++nbcube;
  
  emit (sceneChanged());
}

void
MainWindow::subdivide()
{  
  // do the subdivision
  emit (sceneChanged());
}

void
MainWindow::open(const QString& fileName)
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  //scene.map.clear();

  std::ifstream ifs(qPrintable(fileName));

  CGAL::import_from_polyhedron_flux<Map>(scene.map,ifs);


  this->addToRecentFiles(fileName);
  QApplication::restoreOverrideCursor();
  emit (sceneChanged());
}


#include "MainWindow.moc"

