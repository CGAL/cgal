
#include "MainWindow.h"


#include <QMainWindow>
#include <QActionGroup>
#include <QFileDialog>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/point_generators_2.h>

#include "TriangulationVerticesGraphicsItem_2.h"

  
MainWindow::MainWindow()
{
  setupUi(this);
  setupStatusBar();

  // The navigation adds zooming and translation functionality to the QGraphicsView
  navigation = new CGAL::Navigation(this->graphicsView);
  this->graphicsView->viewport()->installEventFilter(navigation);

  // Add GraphicItems for the Delaunay triangulation, the input points and the Voronoi diagram
#ifdef DELAUNAY_VORONOI
  sdt = new CGAL::QTriangulation_2<Delaunay> (&dt);
  dgi = new CGAL::TriangulationGraphicsItem_2<Delaunay>(&dt);
#else 
  sdt = new CGAL::QConstrainedTriangulation_2<Delaunay> (&dt);
  dgi = new CGAL::ConstrainedTriangulationGraphicsItem_2<Delaunay>(&dt);
#endif    
    
  QObject::connect(sdt, SIGNAL(changed()),
		   dgi, SLOT(modelChanged()));

  CGAL::TriangulationVerticesGraphicsItem_2<Delaunay> * dvgi;
  dvgi = new  CGAL::TriangulationVerticesGraphicsItem_2<Delaunay>(&dt);
  dvgi->setPen(QPen(Qt::red, 2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(dvgi);

  QObject::connect(sdt, SIGNAL(changed()),
		   dvgi, SLOT(modelChanged()));

#ifdef DELAUNAY_VORONOI
  vgi = new CGAL::VoronoiGraphicsItem_2<Delaunay>(&dt);
  vgi->setPen(QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  QObject::connect(sdt, SIGNAL(changed()),
  	   vgi, SLOT(modelChanged()));
#endif    
    

  // Setup input handlers. They get events before the scene gets them
  // and the input they generate is passed to the triangulation with 
  // the signal/slot mechanism
    
  pi = new CGAL::PolylineInput_2<K>(&scene, 0, false); // inputs polylines which are not closed

#ifdef DELAUNAY_VORONOI
  pi->setNumberOfVertices(1);  // In this case we only want to insert points
#endif

  QObject::connect(pi, SIGNAL(produce(CGAL::Object)),
		   sdt, SLOT(consume(CGAL::Object)));
    
  mp = new CGAL::TriangulationMovingPoint_2<Delaunay>(&dt);
  QObject::connect(mp, SIGNAL(produce(CGAL::Object)),
		   sdt, SLOT(consume(CGAL::Object)));
  

  connectActions();


  QObject::connect(navigation, SIGNAL(mouseCoordinates(QString)),
		   this, SLOT(updateMouseCoordinates(QString)));

  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  scene.setSceneRect(-200, -20, 400, 400);
    
  this->graphicsView->setCacheMode(QGraphicsView::CacheBackground);
  this->graphicsView->setRenderHint(QPainter::Antialiasing);
  this->graphicsView->setScene(&scene);
  this->graphicsView->setMinimumSize(400, 400);
}


void
MainWindow::connectActions()
{
  // Connect File Menu actions 
    
  QObject::connect(this->actionClear, SIGNAL(triggered()), 
		   this, SLOT(clear()));

  QObject::connect(this->actionLoadConstraints, SIGNAL(triggered()), 
		   this, SLOT(loadConstraints()));

  QObject::connect(this->actionSaveConstraints, SIGNAL(triggered()), 
		   this, SLOT(saveConstraints()));

  QObject::connect(this->actionExit, SIGNAL(triggered()), 
		   this, SLOT(close()));

  // Connect Edit Menu actions

  QObject::connect(this->actionInsertRandomPoints, SIGNAL(triggered()),
		   this, SLOT(insertRandomPoints()));
    

  // Connect Help Menu actions
  QObject::connect(this->actionAbout, SIGNAL(triggered()),
		   this, SLOT(about()));

  QObject::connect(this->actionAboutCGAL, SIGNAL(triggered()),
		   this, SLOT(aboutCGAL()));


  // Connect Tool Menu actions
QObject::connect(this->actionShowDelaunay, SIGNAL(toggled(bool)),
		   this, SLOT(showDelaunay(bool)));

QObject::connect(this->actionShowVoronoi, SIGNAL(toggled(bool)),
		   this, SLOT(showVoronoi(bool)));
 this->actionShowDelaunay->setChecked(true);
  

  QObject::connect(this->actionInsertPolyline, SIGNAL(toggled(bool)),
		   this, SLOT(insertPolyline(bool)));

  QObject::connect(this->actionMovingPoint, SIGNAL(toggled(bool)),
		   this, SLOT(movingPoint(bool)));

  // We put mutually exclusive actions in an QActionGroup
  QActionGroup* ag = new QActionGroup(this);
  ag->addAction(this->actionInsertPolyline);
  ag->addAction(this->actionMovingPoint);

  // The following function call also emits a toggled signal
  this->actionInsertPolyline->setChecked(true);
}  

void
MainWindow::insertPolyline(bool checked)
{
  if(checked){
    scene.installEventFilter(pi);
  } else {
    scene.removeEventFilter(pi);
  }
}


void
MainWindow::movingPoint(bool checked)
{

  if(checked){
    scene.installEventFilter(mp);
  } else {
    scene.removeEventFilter(mp);
  }
}


void
MainWindow::showDelaunay(bool checked)
{
  if(checked){
    scene.addItem(dgi);
  } else {  
    scene.removeItem(dgi);
  }
}


void
MainWindow::showVoronoi(bool checked)
{
#ifdef DELAUNAY_VORONOI
  if(checked){
    scene.addItem(vgi);
  } else {  
    scene.removeItem(vgi);
  }
#endif
}


void
MainWindow::clear()
{
  sdt->clear();
}


void
MainWindow::loadConstraints()
{
  QString fileName = QFileDialog::getOpenFileName(this,
						  tr("Open Constraint File"),
						  ".",
						  tr("Poly files (*.poly)\n"
						     "Edge files (*.edg)"));
  if(! fileName.isEmpty()){
    loadConstraints(fileName);
  }
}


void
MainWindow::loadConstraints(QString fileName)
{
}

void
MainWindow::saveConstraints()
{
  QString fileName = QFileDialog::getSaveFileName(this,
						  tr("Save Constraints"),
						  ".",
						  tr("Poly files (*.poly)\n"
						     "Edge files (*.edg)"));
  if(! fileName.isEmpty()){
    saveConstraints(fileName);
  }
}


void
MainWindow::saveConstraints(QString fileName)
{
}

void
MainWindow::insertRandomPoints()
{
  typedef CGAL::Creator_uniform_2<double,Point_2>  Creator;
  CGAL::Random_points_in_disc_2<Point_2,Creator> g( 100.0);
  
  std::vector<Point_2> points;
  points.reserve(10);
  for(int i = 0; i < 10; i++){
    points.push_back(*g++);
  }
  sdt->insert(points.begin(), points.end());
}

void
MainWindow::about()
{
  QMessageBox::about(this, tr(" About the Demo"),
		     tr("<h2>Constrained Delaunay Triangulation</h2>"
			"<p>Copyright &copy; 2008 GeometryFactory"
			"<p>This application illustrates the 2D Constrained Delaunay triangulation "
			"of CGAL. <p>See also <a href='http://www.cgal.org/Pkg/Triangulation2'>the online manual</a>"));
}

void
MainWindow::aboutCGAL()
{
  QMessageBox::about(this, tr(" About CGAL"),
		     tr("<h2>CGAL - Computational Geometry Algorithms Library</h2>"
			"<p>CGAL provides efficient and reliable geometric algorithms in the form of a C++ library."
                        "<p>For more information visit <a href='http://www.cgal.org/'>www.cgal.org</a>"));
}
void
MainWindow::setupStatusBar()
{
  xycoord = new QLabel(" -0.00000 , -0.00000 ");
  xycoord->setAlignment(Qt::AlignHCenter);
  xycoord->setMinimumSize(xycoord->sizeHint());
  this->statusbar->addWidget(xycoord);
}
void
MainWindow::updateMouseCoordinates(QString s)
{
  xycoord->setText(s);
}

#include "MainWindow.rcc"
 
#include "MainWindow.moc"
