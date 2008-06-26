#include <fstream>
#include "Constrained_Delaunay_triangulation_MainWindow.h"

#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/point_generators_2.h>

#include "TriangulationCircumcircle.h"
#include "TriangulationMovingPoint.h"
#include <CGAL/Qt/GraphicsViewPolylineInput.h>
#include <CGAL/Qt/TriangulationGraphicsItem.h>
#include <CGAL/Qt/ConstrainedTriangulationGraphicsItem.h>
  
Constrained_Delaunay_triangulation_MainWindow::Constrained_Delaunay_triangulation_MainWindow()
  : DemosMainWindow()
{
  setupUi(this);

  // Add a GraphicItem for the Delaunay triangulation
  dgi = new CGAL::Qt::ConstrainedTriangulationGraphicsItem<Delaunay>(&dt);
    
  QObject::connect(this, SIGNAL(changed()),
		   dgi, SLOT(modelChanged()));

  dgi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(dgi);

  // Setup input handlers. They get events before the scene gets them
  // and the input they generate is passed to the triangulation with 
  // the signal/slot mechanism    
  pi = new CGAL::Qt::GraphicsViewPolylineInput<K>(this, &scene, 0, false); // inputs polylines which are not closed

  QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
		   this, SLOT(processInput(CGAL::Object)));
    
  mp = new CGAL::Qt::TriangulationMovingPoint<Delaunay>(&dt, this);
  // TriangulationMovingPoint<Delaunay> generates an empty Object() each
  // time the moving point moves.
  // The following connection is for the purpose of emitting changed().
  QObject::connect(mp, SIGNAL(generate(CGAL::Object)),
		   this, SLOT(processInput(CGAL::Object)));

  tcc = new CGAL::Qt::TriangulationCircumcircle<Delaunay>(&scene, &dt, this);
  tcc->setPen(QPen(Qt::red, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  

  // 
  // Manual handling of actions
  //
  QObject::connect(this->actionExit, SIGNAL(triggered()), 
		   this, SLOT(close()));

  // We put mutually exclusive actions in an QActionGroup
  QActionGroup* ag = new QActionGroup(this);
  ag->addAction(this->actionInsertPolyline);
  ag->addAction(this->actionMovingPoint);

  // Check two actions 
  this->actionInsertPolyline->setChecked(true);
  this->actionShowDelaunay->setChecked(true);

  //
  // Setup the scene and the view
  //
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  scene.setSceneRect(0,0, 100, 100);
  this->graphicsView->setScene(&scene);

  // Uncomment the following line to get antialiasing by default.
//   actionUse_Antialiasing->setChecked(true);

  // Turn the vertical axe upside down
  this->graphicsView->matrix().scale(1, -1);
                                                      
  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_constrained_Delaunay_triangulation.html");
  this->addAboutCGAL();
}

void
Constrained_Delaunay_triangulation_MainWindow::processInput(CGAL::Object o)
{
  std::list<Point_2> points;
  if(CGAL::assign(points, o)){
    if(points.size() == 1) {
      dt.insert(points.front());
    }
    else {
      insert_polyline(points.begin(), points.end());
    }
  }
  emit(changed());
}

void
Constrained_Delaunay_triangulation_MainWindow::on_actionInsertPolyline_toggled(bool checked)
{
  if(checked){
    scene.installEventFilter(pi);
  } else {
    scene.removeEventFilter(pi);
  }
}


void
Constrained_Delaunay_triangulation_MainWindow::on_actionMovingPoint_toggled(bool checked)
{

  if(checked){
    scene.installEventFilter(mp);
  } else {
    scene.removeEventFilter(mp);
  }
}


void
Constrained_Delaunay_triangulation_MainWindow::on_actionShowDelaunay_toggled(bool checked)
{
  dgi->setDrawEdges(checked);
}


void
Constrained_Delaunay_triangulation_MainWindow::on_actionCircumcenter_toggled(bool checked)
{
  if(checked){
    scene.installEventFilter(tcc);
    this->graphicsView->setMouseTracking(true);
    tcc->show();
  } else {  
    scene.removeEventFilter(tcc);
    this->graphicsView->setMouseTracking(false);
    tcc->hide();
  }
}


void
Constrained_Delaunay_triangulation_MainWindow::on_actionClear_triggered()
{
  dt.clear();
  emit(changed());
}


void
Constrained_Delaunay_triangulation_MainWindow::on_actionLoadConstraints_triggered()
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
Constrained_Delaunay_triangulation_MainWindow::loadConstraints(QString fileName)
{
  std::ifstream ifs(qPrintable(fileName));

  std::list<K::Point_2> points;
  typedef std::list< std::pair<Point_2, Point_2> > Segments;
  Segments segments;
  K::Point_2 p,q;
  while(ifs >> p) {
    points.push_back(p);
  }
  insert_polyline(points.begin(), points.end());

  std::vector<K::Point_2> P;

  actionRecenter->trigger();
  emit(changed());
}

void
Constrained_Delaunay_triangulation_MainWindow::on_actionRecenter_triggered()
{
  this->graphicsView->setSceneRect(dgi->boundingRect());
  this->graphicsView->fitInView(dgi->boundingRect(), Qt::KeepAspectRatio);  
}

void
Constrained_Delaunay_triangulation_MainWindow::on_actionSaveConstraints_triggered()
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
Constrained_Delaunay_triangulation_MainWindow::saveConstraints(QString fileName)
{
  QMessageBox::warning(this,
                       tr("saveConstraints"),
                       tr("Not implemented!"));
}

void
Constrained_Delaunay_triangulation_MainWindow::on_actionInsertRandomPoints_triggered()
{
  typedef CGAL::Creator_uniform_2<double,Point_2>  Creator;
  CGAL::Random_points_in_disc_2<Point_2,Creator> g( 100.0);
  
  const int number_of_points = 
    QInputDialog::getInteger(this, 
                             tr("Number of random points"),
                             tr("Enter number of random points"));

  std::vector<Point_2> points;
  points.reserve(number_of_points);
  for(int i = 0; i < number_of_points; ++i){
    points.push_back(*g++);
  }
  dt.insert(points.begin(), points.end());
  emit(changed());
}

#include "Constrained_Delaunay_triangulation_MainWindow.moc"
