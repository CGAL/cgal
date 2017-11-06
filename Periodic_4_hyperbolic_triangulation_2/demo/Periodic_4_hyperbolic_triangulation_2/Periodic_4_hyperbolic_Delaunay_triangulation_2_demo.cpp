#include <fstream>

// CGAL headers 

#include <CGAL/point_generators_2.h>
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
#include <CGAL/Qt/TriangulationGraphicsItem.h>    
#include <CGAL/Qt/VoronoiGraphicsItem.h>
#include <CGAL/Qt/TriangulationConflictZone.h>
#include <CGAL/Qt/DemosMainWindow.h>
#include <CGAL/Qt/TriangulationPointInput.h>

#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h> 
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>

#include <CGAL/Hyperbolic_octagon_translation.h>

// the two base classes
#include "ui_Periodic_4_hyperbolic_Delaunay_triangulation_2.h"

#include <CGAL/Timer.h>


typedef CORE::Expr                                                              NT;
typedef CGAL::Cartesian<NT>                                                     Kernel;

typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel>     Traits; 
typedef Traits::Segment_2                                                       Segment_2;
typedef Traits::Circle_2                                                        Circle_2;
typedef Traits::Circular_arc_2                                                  Circular_arc_2;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>            Triangulation;
typedef Kernel::Point_2                                                         Point;
typedef Triangulation::Vertex_handle                                            Vertex_handle;
typedef Traits::Side_of_original_octagon                                        Side_of_original_octagon;
typedef Traits::Construct_hyperbolic_segment_2                                  Construct_hyperbolic_segment_2;
typedef Traits::Construct_hyperbolic_line_2                                     Construct_hyperbolic_line_2;
typedef Triangulation::Hyperbolic_translation                                   Hyperbolic_translation;
typedef Triangulation::Locate_type                                              Locate_type;
typedef Triangulation::Face_handle                                              Face_handle;
typedef Triangulation::Face_iterator                                            Face_iterator;
typedef Traits::Construct_point_2                                               Make_point;


using CGAL::to_double;
using std::cout;
using std::endl;

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Periodic_4_hyperbolic_Delaunay_triangulation_2
{
  Q_OBJECT
  
private:    

  Point                                                              movingPoint;

  Triangulation                                                      dt;
  QGraphicsEllipseItem                                             * disk;
  QGraphicsScene                                                     scene;  

  CGAL::Qt::TriangulationGraphicsItem<Triangulation>               * dgi;
  CGAL::Qt::VoronoiGraphicsItem<Triangulation>                     * vgi;
  CGAL::Qt::TriangulationConflictZone<Triangulation>               * cz;
  CGAL::Qt::TriangulationPointInput<Triangulation>                 * pi;
  CGAL::Qt::TriangulationCircumcircle<Triangulation>               * tcc;

  Point     source;
  Point     target;
  double    time;
  bool      go;

  void        animate();
  double      updateTime();
  Face_handle last_location;
  Point       get_image(Point, Point, double);
  double      timestep;
  Hyperbolic_translation      last_Hyperbolic_translation;
  int         idx;

public:
  MainWindow();

public slots:

  void processInput(CGAL::Object o);

  void updateMovingPoint(CGAL::Object o);
  
  void on_actionCircumcenter_toggled(bool checked);

  void on_actionPlayDemo_toggled(bool checked);

  void on_actionShowTriangulation_toggled(bool checked);

  void on_actionShowConflictZone_toggled(bool checked);

  void on_actionShowVoronoi_toggled(bool checked);

  void on_actionShowOctagon_toggled(bool checked);

  void on_actionInsertPoint_toggled(bool checked);
  
  void on_actionInsertRandomPoints_triggered();

  void on_actionLoadPoints_triggered();

  void on_actionClear_triggered();

  void on_actionRecenter_triggered();

  virtual void open(QString fileName);

signals:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow(), dt(Traits())
{
  
  idx = 1;
    //dt.insert_dummy_points(true);

  std::srand(std::time(0));
  double rx1, rx2, ry1, ry2;
  rx1 = ((double)std::rand()/(double)RAND_MAX)*0.6;
  rx2 = ((double)std::rand()/(double)RAND_MAX)*0.6;
  ry1 = ((double)std::rand()/(double)RAND_MAX)*0.6;
  ry2 = ((double)std::rand()/(double)RAND_MAX)*0.6;
  source = Point(rx1, ry1);
  target = Point(rx2, ry2);
  Segment_2 seg = Construct_hyperbolic_line_2()(source, target);
  Circular_arc_2* carc = boost::get<Circular_arc_2>(&seg);
  source = carc->source();
  target = carc->target();

  timestep = 0.005;
  time = updateTime();
  last_location = dt.periodic_locate(get_image(source, target, time), last_Hyperbolic_translation); 
  //cout << "last location: face " << last_location->get_number() << " with Hyperbolic_translation " << last_Hyperbolic_translation << endl; 


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
  //QObject::connect(dgi,   SIGNAL(generate(CGAL::Object)), this, SLOT(processInput(CGAL::Object)));

  dgi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  dgi->setEdgesPen(QPen(QColor(200, 200, 0), 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  dgi->setSource(source);
  dgi->setTarget(target);
  dgi->setMovingPoint(source);
  dgi->setVisibleDemo(false);
  dgi->setVisibleConflictZone(false);
  dgi->setVisibleEdges(false);
  dgi->setVisibleCopies(false);
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
  pi = new CGAL::Qt::TriangulationPointInput<Triangulation>(&scene, &dt, this );

  QObject::connect(pi, SIGNAL(generate(CGAL::Object)), this, SLOT(processInput(CGAL::Object)));

  tcc = new CGAL::Qt::TriangulationCircumcircle<Triangulation>(&scene, &dt, this);
  tcc->setPen(QPen(Qt::blue, 0.01, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));


  cz = new CGAL::Qt::TriangulationConflictZone<Triangulation>(&scene, &dt, this);
  QObject::connect(cz, SIGNAL(generate(CGAL::Object)), this, SLOT(updateMovingPoint(CGAL::Object)));


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
  this->actionShowOctagon->setChecked(true);

  // //
  // // Setup the scene and the view
  // //
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  scene.setSceneRect(left_top_corner_x, left_top_corner_y, width, height);
  this->graphicsView->setScene(&scene);
  this->graphicsView->setMouseTracking(true);
  //this->graphicsView->scale(240, 240);

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
    if (last_location != Face_handle()) {
      for (Triangulation::Face_iterator fit = dt.faces_begin(); fit != dt.faces_end(); fit++)
        fit->tds_data().clear();
    }

    // make sure that the vertices and faces of the triangulation are clean
    for (Triangulation::Face_iterator fi = dt.faces_begin(); fi != dt.faces_end(); fi++) {
      fi->tds_data().clear();
    }
    for (Triangulation::Vertex_iterator vi = dt.vertices_begin(); vi != dt.vertices_end(); vi++) {
      vi->remove_translation();
    }

    Vertex_handle v = dt.insert(p);
  }
  emit(changed());
  qApp->processEvents();

}

void
MainWindow::updateMovingPoint(CGAL::Object o) {
  CGAL::assign(movingPoint, o);
  dgi->setMovingPoint(movingPoint);
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


void MainWindow::on_actionPlayDemo_toggled(bool checked) {
  on_actionInsertPoint_toggled(!checked);
  if (checked) {
    dgi->setVisibleDemo(true);
    go = true;
    animate();
  } else {
    dgi->setVisibleDemo(false);
    go = false;
  }
}


Point
MainWindow::get_image(Point src, Point tgt, double time) {
  Segment_2 seg = Construct_hyperbolic_segment_2()(src, tgt);
  Circular_arc_2* carc = boost::get<Circular_arc_2>(&seg);
  Circle_2 crc = carc->circle();
  
  double sx = to_double(((src.x()) - crc.center().x())/sqrt(crc.squared_radius()));
  double sy = to_double(((src.y()) - crc.center().y())/sqrt(crc.squared_radius()));
  double tx = to_double(((tgt.x()) - crc.center().x())/sqrt(crc.squared_radius()));
  double ty = to_double(((tgt.y()) - crc.center().y())/sqrt(crc.squared_radius()));
  
  double dot = to_double(sx*tx + sy*ty);
  double n1 = sqrt(sx*sx+sy*sy);
  double n2 = sqrt(tx*tx+ty*ty);
  double theta = acos(dot/n1/n2);

  double x = sin((1.0-time)*theta)/sin(theta)*sx + sin(time*theta)/sin(theta)*tx;
  double y = sin((1.0-time)*theta)/sin(theta)*sy + sin(time*theta)/sin(theta)*ty;
  
  x = to_double(x*sqrt(crc.squared_radius()) + crc.center().x());
  y = to_double(y*sqrt(crc.squared_radius()) + crc.center().y());

  Point p(x, y);

  return p;
}


double 
MainWindow::updateTime() {
  double t = 0.1;
  Side_of_original_octagon check;
  while(check(get_image(source, target, t)) == CGAL::ON_UNBOUNDED_SIDE) {
    t += timestep;
  }
  return t;
}


void
MainWindow::animate() {
  
  Point p = get_image(source, target, time);

  Side_of_original_octagon check;
  if (check(p) != CGAL::ON_UNBOUNDED_SIDE) {
    Locate_type lt;
    int li;
    last_location = dt.periodic_locate(p, lt, li, last_Hyperbolic_translation, last_location);
    dgi->setMovingPoint(p);
    dgi->setSource(source);
    dgi->setTarget(target);
    
    time += timestep;

  } else {

    Hyperbolic_translation o;
    for (int i = 0; i < 3; i++) {
      if (last_Hyperbolic_translation.is_identity()) {
        o = last_location->translation(i).inverse();  
        
      } else {
        o = last_location->translation(i);
      }
      
      if (check(Make_point()(p,o)) == CGAL::ON_BOUNDED_SIDE) {
        break;
      }
    }

    assert(check(Make_point()(p,o)) == CGAL::ON_BOUNDED_SIDE);

    source = Make_point()(source,o);
    target = Make_point()(target,o);
    
    Segment_2 seg = Construct_hyperbolic_line_2()(source, target);
    Circular_arc_2* carc = boost::get<Circular_arc_2>(&seg);
    if (sqrt(squared_distance(source, carc->source())) < sqrt(squared_distance(source, carc->target()))) {
      source = carc->source();
      target = carc->target();
    } else {
      source = carc->target();
      target = carc->source();
    }

    time = updateTime();
  }

  emit(changed());
  qApp->processEvents();

  if (go) {
    //boost::this_thread::sleep(boost::posix_time::milliseconds(50));

    // Uncomment here to generate snapshots!
    // QPixmap sshot = this->grab(QRect(QPoint(560, 210), QSize(780, 770)));
    // //QPixmap sshot = this->grab(QRect(QPoint(320, 55), QSize(625, 615)));
    // std::stringstream ss;
    // ss << "/Users/iordanov/Desktop/shots/sshot" << (idx++) << ".png";
    // sshot.save(QString(ss.str().c_str()));
    
    animate();
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
MainWindow::on_actionShowConflictZone_toggled(bool checked)
{
  dgi->setVisibleConflictZone(checked);
  if(checked){
    scene.installEventFilter(cz);
  } else {  
    scene.removeEventFilter(cz);
  }
}


void
MainWindow::on_actionShowVoronoi_toggled(bool checked)
{
  vgi->setVisible(checked);
}

void
MainWindow::on_actionShowOctagon_toggled(bool checked)
{
  dgi->setVisibleOctagon(checked);
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
  std::vector<Point_d> v;
  Hyperbolic_random_points_in_disc_2_double(v, 5*number_of_points, -1, 0.159);

  Traits::Side_of_original_octagon pred;

  std::vector<Point> pts;
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
  dt.insert(pts.begin(), pts.end());
  // for (int i = 0; i < pts.size(); i++) {
  //   dt.insert(pts[i]);
  // }
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
    open(fileName);
  }
}


void
MainWindow::open(QString fileName)
{
  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);
  std::ifstream ifs(qPrintable(fileName));

  Point p;
  std::vector<Point> points;
  while(ifs >> p) {
    points.push_back(p);
  }
  cout << "Read " << points.size() << " input points." << endl << "Inserting!" << endl;
  dt.insert(points.begin(), points.end());
  cout << "Done inserting!" << endl;

  // default cursor
  QApplication::restoreOverrideCursor();
  emit(changed());    
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
  app.setApplicationName("P4HDT2 Demo");
  MainWindow mainWindow;
  mainWindow.show();
  QStringList args = app.arguments();
  args.removeAt(0);

  return app.exec();
}
