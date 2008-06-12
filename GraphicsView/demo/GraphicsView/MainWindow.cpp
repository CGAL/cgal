#include <fstream>
#include "MainWindow.h"


#include <QMainWindow>
#include <QActionGroup>
#include <QFileDialog>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/point_generators_2.h>

#include "QTriangulationVerticesGraphicsItem_2.h"
#include "QTriangulationCircumcenter_2.h"

  
MainWindow::MainWindow()
{
  setupUi(this);
  setupStatusBar();

  // The navigation adds zooming and translation functionality to the QGraphicsView
  navigation = new CGAL::QNavigation(this->graphicsView);
  this->graphicsView->viewport()->installEventFilter(navigation);

  //  navigation2 = new CGAL::Navigation2(this->graphicsView);
  //this->graphicsView->installEventFilter(navigation2);

  // Add GraphicItems for the Delaunay triangulation, the input points and the Voronoi diagram
#ifdef DELAUNAY_VORONOI
  dgi = new CGAL::QTriangulationGraphicsItem_2<Delaunay>(&dt);
#else 
  dgi = new CGAL::QConstrainedTriangulationGraphicsItem_2<Delaunay>(&dt);
#endif    
    
  QObject::connect(this, SIGNAL(changed()),
		   dgi, SLOT(modelChanged()));

  CGAL::QTriangulationVerticesGraphicsItem_2<Delaunay> * dvgi;
  dvgi = new  CGAL::QTriangulationVerticesGraphicsItem_2<Delaunay>(&dt);
  dvgi->setPen(QPen(Qt::red, 1, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(dvgi);

  QObject::connect(this, SIGNAL(changed()),
  	   dvgi, SLOT(modelChanged()));

#ifdef DELAUNAY_VORONOI
  vgi = new CGAL::QVoronoiGraphicsItem_2<Delaunay>(&dt);
  vgi->setPen(QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  QObject::connect(this, SIGNAL(changed()),
		   vgi, SLOT(modelChanged()));
#endif    
    

  // Setup input handlers. They get events before the scene gets them
  // and the input they generate is passed to the triangulation with 
  // the signal/slot mechanism
    
  pi = new CGAL::QPolylineInput_2<K>(&scene, 0, false); // inputs polylines which are not closed
  tcc = new CGAL::QTriangulationCircumcenter_2<Delaunay>(&scene, &dt);
  tcc->setPen(QPen(Qt::red, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

#ifdef DELAUNAY_VORONOI
  pi->setNumberOfVertices(1);  // In this case we only want to insert points
#endif

  QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
		   this, SLOT(process(CGAL::Object)));
    
  mp = new CGAL::QTriangulationMovingPoint_2<Delaunay>(&dt);
  QObject::connect(mp, SIGNAL(generate(CGAL::Object)),
		   this, SLOT(process(CGAL::Object)));
  

  connectActions();

  actionShowDelaunay->setChecked(true);


  QObject::connect(navigation, SIGNAL(mouseCoordinates(QString)),
		   this, SLOT(updateMouseCoordinates(QString)));

  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  scene.setSceneRect(0,0, 100, 100);

  //  QMatrix m(1.0, 0.0, 0.0, -1.0, 0.0, 0.0);  // identity matrix with reversed Y
  //this->graphicsView->setMatrix(m); 
  this->graphicsView->setTransformationAnchor(QGraphicsView::NoAnchor);
  this->graphicsView->setCacheMode(QGraphicsView::CacheBackground);
  this->graphicsView->setRenderHint(QPainter::Antialiasing);
  this->graphicsView->setScene(&scene);
  this->graphicsView->setMinimumSize(200, 200); // this is the size in pixels on the screen
}


void
MainWindow::connectActions()
{
  QObject::connect(this->actionExit, SIGNAL(triggered()), 
		   this, SLOT(close()));

  // We put mutually exclusive actions in an QActionGroup
  QActionGroup* ag = new QActionGroup(this);
  ag->addAction(this->actionInsertPolyline);
  ag->addAction(this->actionMovingPoint);

  // The following function call also emits a toggled signal
  this->actionInsertPolyline->setChecked(true);
}  

void
MainWindow::process(CGAL::Object o)
{
  std::list<Point_2> points;
  if(CGAL::assign(points, o)){
    dt.insert(points.begin(), points.end());
  }
  emit(changed());
}

void
MainWindow::on_actionInsertPolyline_toggled(bool checked)
{
  if(checked){
    scene.installEventFilter(pi);
  } else {
    scene.removeEventFilter(pi);
  }
}


void
MainWindow::on_actionMovingPoint_toggled(bool checked)
{

  if(checked){
    scene.installEventFilter(mp);
  } else {
    scene.removeEventFilter(mp);
  }
}


void
MainWindow::on_actionShowDelaunay_toggled(bool checked)
{
  if(checked){
    scene.addItem(dgi);
  } else {  
    scene.removeItem(dgi);
  }
}


void
MainWindow::on_actionShowVoronoi_toggled(bool checked)
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
MainWindow::on_actionCircumcenter_toggled(bool checked)
{
#ifdef DELAUNAY_VORONOI
  if(checked){
    scene.installEventFilter(tcc);
    this->graphicsView->setMouseTracking(true);
    tcc->show();
  } else {  
    scene.removeEventFilter(vgi);
    this->graphicsView->setMouseTracking(false);
    tcc->hide();
  }
#endif
}


void
MainWindow::on_actionClear_triggered()
{
  dt.clear();
  emit(changed());
  scene.update(); // do we need that?
}


void
MainWindow::on_actionLoadConstraints_triggered()
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
  std::ifstream ifs(qPrintable(fileName));

  std::list<K::Point_2> points;
  std::list< std::pair<K::Point_2,K::Point_2> > segments;
  K::Point_2 p,q;

  std::vector<K::Point_2> P;
  P.resize(2);
  char c;
  while(ifs >> c){
    ifs >>p >> q;
    if(p != q){
      segments.push_back(std::make_pair(p,q));
    }
  }
  // dt.insert_constraints(segments.begin(), segments.end());
  /*
  while(ifs >> p){
    points.push_back(p);
  }
  //std::cout << "Read " << points.size() << " points" << std::endl;
  dt.insert_polyline(points.begin(), points.end());
  emit(changed());
  */
     this->graphicsView->ensureVisible(dgi);
}

void
MainWindow::on_actionSaveConstraints_triggered()
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
MainWindow::on_actionInsertRandomPoints_triggered()
{
  typedef CGAL::Creator_uniform_2<double,Point_2>  Creator;
  CGAL::Random_points_in_disc_2<Point_2,Creator> g( 100.0);
  
  std::vector<Point_2> points;
  points.reserve(10);
  for(int i = 0; i < 10; i++){
    points.push_back(*g++);
  }
  dt.insert(points.begin(), points.end());
  emit(changed());
}

void
MainWindow::on_actionAbout_triggered()
{
  QMessageBox::about(this, tr(" About the Demo"),
		     tr("<h2>Constrained Delaunay Triangulation</h2>"
			"<p>Copyright &copy; 2008 GeometryFactory"
			"<p>This application illustrates the 2D Constrained Delaunay triangulation "
			"of CGAL. <p>See also <a href='http://www.cgal.org/Pkg/Triangulation2'>the online manual</a>"));
}

void
MainWindow::on_actionAboutCGAL_triggered()
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
