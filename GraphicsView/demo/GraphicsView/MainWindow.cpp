#include <fstream>
#include "MainWindow.h"


#include <QMainWindow>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <QLabel>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/point_generators_2.h>

#include "QTriangulationCircumcenter_2.h"
#include "QTriangulationMovingPoint_2.h"
#include <CGAL/IO/QtPolylineInput.h>
#include <CGAL/IO/QtTriangulationGraphicsItem.h>
#include <CGAL/IO/QtConstrainedTriangulationGraphicsItem.h>
#include <CGAL/IO/QtVoronoiGraphicsItem.h>

#include <CGAL/IO/QtNavigation.h>

#include <QGLWidget>


  
MainWindow::MainWindow()
{
  setupUi(this);
  setupStatusBar();

  // Add GraphicItems for the Delaunay triangulation, the input points and the Voronoi diagram
#ifdef DELAUNAY_VORONOI
  dgi = new CGAL::QtTriangulationGraphicsItem<Delaunay>(&dt);
#else 
  dgi = new CGAL::QtConstrainedTriangulationGraphicsItem<Delaunay>(&dt);
#endif    
    
  QObject::connect(this, SIGNAL(changed()),
		   dgi, SLOT(modelChanged()));

  dgi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

#ifdef DELAUNAY_VORONOI
  vgi = new CGAL::QtVoronoiGraphicsItem<Delaunay>(&dt);
  vgi->setPen(QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  QObject::connect(this, SIGNAL(changed()),
		   vgi, SLOT(modelChanged()));
#endif    
    

  // Setup input handlers. They get events before the scene gets them
  // and the input they generate is passed to the triangulation with 
  // the signal/slot mechanism
    
  pi = new CGAL::QtPolylineInput<K>(&scene, 0, false); // inputs polylines which are not closed
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


  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  scene.setSceneRect(0,0, 100, 100);

  //  QMatrix m(1.0, 0.0, 0.0, -1.0, 0.0, 0.0);  // identity matrix with reversed Y
  //this->graphicsView->setMatrix(m); 
  this->graphicsView->setTransformationAnchor(QGraphicsView::NoAnchor);
  this->graphicsView->setCacheMode(QGraphicsView::CacheBackground);
  this->graphicsView->setRenderHint(QPainter::Antialiasing);
  this->graphicsView->setScene(&scene);
  this->graphicsView->setMinimumSize(200, 200); // this is the size in pixels on the screen

  // The navigation adds zooming and translation functionality to the QGraphicsView
  navigation = new CGAL::QtNavigation(this->graphicsView);
  this->graphicsView->viewport()->installEventFilter(navigation);
  this->graphicsView->installEventFilter(navigation);

  QObject::connect(navigation, SIGNAL(mouseCoordinates(QString)),
		   xycoord, SLOT(setText(QString)));
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
    if(points.size() == 1) {
      dt.insert(points.front());
    }
#ifndef DELAUNAY_VORONOI
    else {
      insert_polyline(points.begin(), points.end());
    }
#endif // ifndef DELAUNAY_VORONOI
  }
  emit(changed());
}

void
MainWindow::on_actionUse_OpenGL_toggled(bool checked)
{ 
  if(checked) {
    QGLWidget* new_viewport = new QGLWidget;

    // Setup the format to allow antialiasing with OpenGL:
    // one need to activate the SampleBuffers, if the graphic driver allows
    // this.
    QGLFormat glformat = new_viewport->format();
    glformat.setOption(QGL::SampleBuffers);
    new_viewport->setFormat(glformat);

    this->graphicsView->setViewport(new_viewport);
    statusBar()->showMessage(tr("OpenGL activated"), 1000);
  }
  else {
    this->graphicsView->setViewport(new QWidget);
    statusBar()->showMessage(tr("OpenGL deactivated"), 1000);
  }
  this->graphicsView->viewport()->installEventFilter(navigation);
  this->graphicsView->setFocus();
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
    dgi->setDrawEdges(true);
  } else {
    dgi->setDrawEdges(false);
    if(!dgi->drawVertices())
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
  typedef std::list< std::pair<Point_2, Point_2> > Segments;
  Segments segments;
  K::Point_2 p,q;
  while(ifs >> p) {
    points.push_back(p);
  }
#ifndef DELAUNAY_VORONOI
  insert_polyline(points.begin(), points.end());
#endif

  std::vector<K::Point_2> P;
//   P.resize(2);
//   char c;
//   while(ifs >> c){
//     ifs >>p >> q;
//     if(p != q){
//       segments.push_back(std::make_pair(p,q));
//     }
//   }
// #ifndef DELAUNAY_VORONOI
//   for(Segments::const_iterator it = segments.begin();
//       it != segments.end(); ++it)
//   {
//     dt.insert_constraint(it->first, it->second);
//   }
// #endif
  /*
  while(ifs >> p){
    points.push_back(p);
  }
  //std::cout << "Read " << points.size() << " points" << std::endl;
  dt.insert_polyline(points.begin(), points.end());
  */
  // recenter
  actionRecenter->trigger();
  emit(changed());
}

void
MainWindow::on_actionRecenter_triggered()
{
  this->graphicsView->setSceneRect(dgi->boundingRect());
  this->graphicsView->fitInView(dgi->boundingRect(), Qt::KeepAspectRatio);  
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
  this->graphicsView->fitInView(dgi->boundingRect(), Qt::KeepAspectRatio);
}

void
MainWindow::on_actionInsertRandomPoints_triggered()
{
  typedef CGAL::Creator_uniform_2<double,Point_2>  Creator;
  CGAL::Random_points_in_disc_2<Point_2,Creator> g( 100.0);
  
  const int number_of_points = QInputDialog::getInteger(this, 
                                                        tr("Number of random points"),
                                                        tr("Enter number of random points"));
  std::vector<Point_2> points;
  points.reserve(number_of_points);
  for(int i = 0; i < number_of_points; i++){
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
  xycoord = new QLabel(" -0.00000 , -0.00000 ", this);
  xycoord->setAlignment(Qt::AlignHCenter);
  xycoord->setMinimumSize(xycoord->sizeHint());
  this->statusbar->addWidget(xycoord);
}

#include "MainWindow.moc"
