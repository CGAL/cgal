
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
	  this, SLOT(load_off(QString)));
}


void
MainWindow::connectActions()
{
  QObject::connect(this->actionImportOFF, SIGNAL(triggered()), 
		   this, SLOT(import_off()));

  QObject::connect(this->actionAddOFF, SIGNAL(triggered()), 
		   this, SLOT(add_off()));

  QObject::connect(this->actionSubdivide, SIGNAL(triggered()), 
		   this, SLOT(subdivide()));

  QObject::connect(this->actionCreateCube, SIGNAL(triggered()), 
		   this, SLOT(create_cube()));

  QObject::connect(this, SIGNAL(sceneChanged()), 
		   this->viewer, SLOT(sceneChanged()));


  QObject::connect(this->actionQuit, SIGNAL(triggered()), 
		   qApp, SLOT(quit()));
}

void
MainWindow::import_off()
{
  QString fileName = QFileDialog::getOpenFileName(this,
						  tr("Import OFF"),
						  "./off",
						  tr("off files (*.off)"));

  if(! fileName.isEmpty())
    {
      scene.map.clear();
      load_off(fileName);
    }
}

void
MainWindow::add_off()
{
  QString fileName = QFileDialog::getOpenFileName(this,
						  tr("Add OFF"),
						  "./off",
						  tr("off files (*.off)"));

  if(! fileName.isEmpty())
    {
      load_off(fileName);
    }
}


void
MainWindow::create_cube()
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
MainWindow::load_off(const QString& fileName)
{
  QApplication::setOverrideCursor(Qt::WaitCursor);

  std::ifstream ifs(qPrintable(fileName));

  CGAL::import_from_polyhedron_flux<Map>(scene.map,ifs);

  this->addToRecentFiles(fileName);
  QApplication::restoreOverrideCursor();
  emit (sceneChanged());
}


#include "MainWindow.moc"

