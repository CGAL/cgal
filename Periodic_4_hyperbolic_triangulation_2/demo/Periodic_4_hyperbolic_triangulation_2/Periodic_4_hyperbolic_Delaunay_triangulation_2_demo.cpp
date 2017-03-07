#include <fstream>

// CGAL headers
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h> 
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h> 
#include <CGAL/point_generators_2.h>
#include <CGAL/Hyperbolic_octagon_translation_matrix.h>
#include <CGAL/Hyperbolic_random_points_in_disc_2.h>

#include <CGAL/Qt/HyperbolicPainterOstream.h>
// for viewportsBbox
#include <CGAL/Qt/utility.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <QGraphicsEllipseItem>

// for filtering
#include <set>
#include <string>

// GraphicsView items and event filters (input classes)
#include <CGAL/Qt/TriangulationCircumcircle.h>
#include <CGAL/Qt/TriangulationPointInputAndConflictZone.h>
#include <CGAL/Qt/TriangulationGraphicsItem.h>    
#include <CGAL/Qt/VoronoiGraphicsItem.h>
#include <CGAL/Qt/DemosMainWindow.h>

#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>

// the two base classes
#include "ui_Periodic_4_hyperbolic_Delaunay_triangulation_2.h"

#include <CGAL/Timer.h>


typedef CORE::Expr                                                              NT;
typedef CGAL::Cartesian<NT>                                                     Kernel;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel>     Traits;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>            Triangulation;
typedef Hyperbolic_octagon_translation_matrix<NT>                               Octagon_matrix;
typedef Kernel::Point_2                                                         Point;
typedef Triangulation::Vertex_handle                                            Vertex_handle;
typedef Traits::Side_of_fundamental_octagon                                     Side_of_fundamental_octagon;


class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Periodic_4_hyperbolic_Delaunay_triangulation_2
{
  Q_OBJECT
  
private:    

  Triangulation                                                      dt;
  QGraphicsEllipseItem                                             * disk;
  QGraphicsScene                                                     scene;  

  CGAL::Qt::TriangulationGraphicsItem<Triangulation>               * dgi;
  CGAL::Qt::VoronoiGraphicsItem<Triangulation>                     * vgi;

  CGAL::Qt::TriangulationPointInputAndConflictZone<Triangulation>  * pi;
  CGAL::Qt::TriangulationCircumcircle<Triangulation>               * tcc;
public:
  MainWindow();

public slots:

  void processInput(CGAL::Object o);
  
  void on_actionCircumcenter_toggled(bool checked);

  void on_actionShowTriangulation_toggled(bool checked);

  void on_actionShowVoronoi_toggled(bool checked);

  void on_actionInsertPoint_toggled(bool checked);
  
  void on_actionInsertRandomPoints_triggered();

  void on_actionLoadPoints_triggered();

  void on_actionClear_triggered();

  void on_actionRecenter_triggered();


signals:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow(), dt(Traits())
{

  dt.insert_dummy_points(true);
  
  setupUi(this);

  this->graphicsView->setAcceptDrops(false);
  
  // Add Poincaré disk
  qreal origin_x = 0, origin_y = 0, radius = 1, diameter = 2*radius;
  qreal left_top_corner_x = origin_x - radius;
  qreal left_top_corner_y = origin_y - radius;
  qreal width = diameter, height = diameter;
  
  // set background
  qreal eps = 0.01;
  QGraphicsRectItem* rect = new QGraphicsRectItem(left_top_corner_x - eps, left_top_corner_y - eps, width + 2*eps, height + 2*eps);
  rect->setPen(Qt::NoPen);
  rect->setBrush(Qt::white);
  scene.addItem(rect);
  
  // Add a GraphicItem for the Triangulation triangulation
  dgi = new CGAL::Qt::TriangulationGraphicsItem<Triangulation>(&dt);

  QObject::connect(this, SIGNAL(changed()), dgi, SLOT(modelChanged()));

  dgi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  dgi->setEdgesPen(QPen(QColor(200, 200, 0), 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(dgi);

  // Add a GraphicItem for the Voronoi diagram
  vgi = new CGAL::Qt::VoronoiGraphicsItem<Triangulation>(&dt);

  QObject::connect(this, SIGNAL(changed()), vgi, SLOT(modelChanged()));

  vgi->setEdgesPen(QPen(Qt::blue, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(vgi);
  vgi->hide();


  // Setup input handlers. They get events before the scene gets them
  // and the input they generate is passed to the triangulation with 
  // the signal/slot mechanism    
  pi = new CGAL::Qt::TriangulationPointInputAndConflictZone<Triangulation>(&scene, &dt, this );

  QObject::connect(pi, SIGNAL(generate(CGAL::Object)), this, SLOT(processInput(CGAL::Object)));

  tcc = new CGAL::Qt::TriangulationCircumcircle<Triangulation>(&scene, &dt, this);
  tcc->setPen(QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

 
  // 
  // Manual handling of actions
  //

  QObject::connect(this->actionQuit, SIGNAL(triggered()), this, SLOT(close()));

  // We put mutually exclusive actions in an QActionGroup
  QActionGroup* ag = new QActionGroup(this);
  ag->addAction(this->actionInsertPoint);
  ag->addAction(this->actionCircumcenter);

  // Check two actions 
  this->actionInsertPoint->setChecked(true);
  this->actionShowTriangulation->setChecked(true);

  // //
  // // Setup the scene and the view
  // //
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  scene.setSceneRect(left_top_corner_x, left_top_corner_y, width, height);
  this->graphicsView->setScene(&scene);
  this->graphicsView->setMouseTracking(true);
  this->graphicsView->shear(230, 230);
  this->graphicsView->rotate(90);

  // // The navigation adds zooming and translation functionality to the
  // // QGraphicsView
  this->addNavigation(this->graphicsView);
  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_Triangulation_triangulation_2.html");
  this->addAboutCGAL();

}


void
MainWindow::processInput(CGAL::Object o)
{

  Point p;
  if(CGAL::assign(p, o)){
    Vertex_handle v = dt.insert(p);
  }
  emit(changed());

}




/* 
 *  Qt Automatic Connections
 *  http://doc.trolltech.com/4.4/designer-using-a-component.html#automatic-connections
 * 
 *  setupUi(this) generates connections to the slots named
 *  "on_<action_name>_<signal_name>"
 */
void
MainWindow::on_actionInsertPoint_toggled(bool checked)
{
  if(checked){
    scene.installEventFilter(pi);
  } else {
    scene.removeEventFilter(pi);
  }
}



void
MainWindow::on_actionCircumcenter_toggled(bool checked)
{
  if(checked){
    scene.installEventFilter(tcc);
    tcc->show();
  } else {  
    scene.removeEventFilter(tcc);
    tcc->hide();
  }
}



void
MainWindow::on_actionShowTriangulation_toggled(bool checked)
{
  dgi->setVisibleEdges(checked);
}


void
MainWindow::on_actionShowVoronoi_toggled(bool checked)
{
  vgi->setVisible(checked);
}


void
MainWindow::on_actionClear_triggered()
{
  dt.clear();
  emit(changed());
}


void
MainWindow::on_actionInsertRandomPoints_triggered()
{
  bool ok = false;
  const int number_of_points = 
    QInputDialog::getInt(this, 
                        tr("Number of random points"),
                        tr("Enter number of random points"),
			     100,
			     0,
			     std::numeric_limits<int>::max(),
			     1,
			     &ok);

  if(!ok) {
    return;
  }

  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);


  // typedef CGAL::Creator_uniform_2<double, Point> Creator;
  // CGAL::Random_points_in_disc_2<Point, Creator> g( 1.0 );
  typedef CGAL::Cartesian<double>::Point_2  Point_d;
  vector<Point_d> v;
  Hyperbolic_random_points_in_disc_2_double(v, 5*number_of_points, -1, 0.159);

  Traits::Side_of_fundamental_octagon pred;

  vector<Point> pts;
  int cnt = 0;
  for (int i = 0; cnt < number_of_points && i < v.size(); i++) {
    if (pred(v[i]) != CGAL::ON_UNBOUNDED_SIDE) {
      pts.push_back(Point(v[i].x(), v[i].y()));
      cnt++;
    }
  } 
  
  if (pts.size() < number_of_points) {
    cout << "Creation of random points failed! Please try again..." << endl;
    return;
  }

  CGAL::Timer tt;
  tt.start();
  for (int i = 0; i < pts.size(); i++) {
    dt.insert(pts[i]);
  }
  tt.stop();

  cout << "Time elapsed for the insertion of " << number_of_points << " points: " << tt.time() << " secs." << endl;
  cout << "Number of vertices in the triangulation: " << dt.number_of_vertices() << endl;

  QApplication::restoreOverrideCursor();
  emit(changed());
}


void
MainWindow::on_actionLoadPoints_triggered()
{
  QString fileName = QFileDialog::getOpenFileName(this,
						  tr("Open Points file"),
						  ".");
  if(! fileName.isEmpty()){
    //open(fileName);
  }
}


void
MainWindow::on_actionRecenter_triggered()
{
  this->graphicsView->setSceneRect(dgi->boundingRect());
  this->graphicsView->fitInView(dgi->boundingRect(), Qt::KeepAspectRatio);  
}


#include "Periodic_4_hyperbolic_Delaunay_triangulation_2_demo.moc"

int main(int argc, char **argv)
{

  QApplication app(argc, argv);
  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Periodic_4_hyperbolic_Delaunay_triangulation_2 testing name");
  MainWindow mainWindow;
  mainWindow.show();
  QStringList args = app.arguments();
  args.removeAt(0);

  return app.exec();
}
