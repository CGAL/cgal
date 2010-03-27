
#include "MainWindow.h"
#include <CGAL/Delaunay_triangulation_3.h>

// Function defined in map_3_subivision.cpp
void subdivide_map_3(Map& m);
  
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
  
  QObject::connect(this->actionImport3DTDS, SIGNAL(triggered()), 
		   this, SLOT(import_3DTDS()));

  QObject::connect(this->actionQuit, SIGNAL(triggered()), 
		   qApp, SLOT(quit()));
  
  QObject::connect(this->actionSubdivide, SIGNAL(triggered()), 
		   this, SLOT(subdivide()));

  QObject::connect(this->actionCreate3Cubes, SIGNAL(triggered()), 
		   this, SLOT(create_3cubes()));

  QObject::connect(this->actionCreate2Volumes, SIGNAL(triggered()), 
		   this, SLOT(create_2volumes()));

  QObject::connect(this, SIGNAL(sceneChanged()), 
		   this, SLOT(onSceneChanged()));

  QObject::connect(this->actionClear, SIGNAL(triggered()), 
		   this, SLOT(clear()));

}

void
MainWindow::onSceneChanged()
{
  int mark = scene.map.get_new_mark();
  scene.map.negate_mask_mark(mark);
  
  unsigned int nb0, nb1, nb2, nb3, nb4;

  scene.map.count_cells(mark, &nb0, &nb1, &nb2, &nb3, &nb4, NULL); 

  std::ostringstream os;
  os<<"Darts: "<< scene.map.size_of_darts()
    <<",  Vertices:"<<nb0
    <<",  Edges:"<<nb1
    <<",  Faces:"<<nb2
    <<",  Volumes:"<<nb3
    <<",  Connected components:"<<nb4;
  
  scene.map.free_mark(mark);

  viewer->sceneChanged();

  statusBar()->showMessage(os.str().c_str());
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
      load_off(fileName,true);
    }
}

void
MainWindow::import_3DTDS()
{
  QString fileName = QFileDialog::getOpenFileName(this,
						  tr("Import 3DTDS"),
						  ".",
						  tr("Data file (*)"));

  if(! fileName.isEmpty())
    {
      load_3DTDS(fileName,true);
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
      load_off(fileName,false);
    }
}

void
MainWindow::load_off(const QString& fileName, bool clear)
{
  QApplication::setOverrideCursor(Qt::WaitCursor);

  if (clear) scene.map.clear();

  std::ifstream ifs(qPrintable(fileName));

  CGAL::import_from_polyhedron_flux<Map>(scene.map,ifs);

  this->addToRecentFiles(fileName);
  QApplication::restoreOverrideCursor();
  emit (sceneChanged());
}

void
MainWindow::load_3DTDS(const QString& fileName, bool clear)
{
  QApplication::setOverrideCursor(Qt::WaitCursor);

  if (clear) scene.map.clear();

  typedef CGAL::Delaunay_triangulation_3<Kernel>  Triangulation;
  Triangulation T;
    
  std::ifstream ifs(qPrintable(fileName));
  std::istream_iterator<Point_3> begin(ifs), end;
  T.insert(begin, end);

  CGAL::import_from_triangulation_3<Map, Triangulation>(scene.map, T);

  QApplication::restoreOverrideCursor();
  emit (sceneChanged());
}

void
MainWindow::create_3cubes()
{  
  Dart_handle d1=make_cube(scene.map, Point_3(nbcube, nbcube, nbcube), 1);
  Dart_handle d2=make_cube(scene.map, Point_3(nbcube+1, nbcube, nbcube), 1);
  Dart_handle d3=make_cube(scene.map, Point_3(nbcube, nbcube+1, nbcube), 1);

  scene.map.sew3(d1->beta(1,1)->beta(2),d2->beta(2));
  scene.map.sew3(d1->beta(2)->beta(1,1)->beta(2),d3);

  ++nbcube;
  
  emit (sceneChanged());
}

void
MainWindow::create_2volumes()
{
  Dart_handle d1=make_cube(scene.map, Point_3(nbcube, nbcube, nbcube), 1);
  Dart_handle d2=make_cube(scene.map, Point_3(nbcube+1, nbcube, nbcube), 1);
  Dart_handle d3=make_cube(scene.map, Point_3(nbcube, nbcube+1, nbcube), 1);
  Dart_handle d4=make_cube(scene.map, Point_3(nbcube+1, nbcube+1, nbcube), 1);

  scene.map.sew3(d1->beta(1,1)->beta(2),d2->beta(2));
  scene.map.sew3(d1->beta(2)->beta(1,1)->beta(2),d3);

  scene.map.sew3(d3->beta(1,1)->beta(2),d4->beta(2));
  scene.map.sew3(d2->beta(2)->beta(1,1)->beta(2),d4);

  remove_face_3(scene.map,d3);
  remove_face_3(scene.map,d2->beta(2));
  
  ++nbcube;
  
  emit (sceneChanged());
}

void
MainWindow::subdivide()
{
  subdivide_map_3(scene.map);
  emit (sceneChanged());
}

void
MainWindow::clear()
{  
  scene.map.clear();
  emit (sceneChanged());
}

#include "MainWindow.moc"

