#include <fstream>

// CGAL headers

#include <CGAL/point_generators_2.h>
#include <internal/Qt/HyperbolicPainterOstream.h>
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
#include <CGAL/Qt/DemosMainWindow.h>
#include <internal/Qt/TriangulationCircumcircle.h>
#include <internal/Qt/TriangulationGraphicsItem.h>
#include <internal/Qt/VoronoiGraphicsItem.h>
#include <internal/Qt/TriangulationConflictZone.h>
#include <internal/Qt/TriangulationPointInput.h>

#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>

#include <CGAL/Hyperbolic_octagon_translation.h>

// the two base classes
#include "ui_P4HDT2.h"

#include <CGAL/Timer.h>

typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<>           Traits;
typedef Traits::Segment_2                                                       Segment_2;
typedef Traits::Circle_2                                                        Circle_2;
typedef Traits::Circular_arc_2                                                  Circular_arc_2;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<>                  Triangulation;
typedef Triangulation::Point                                                    Point;
typedef Triangulation::Vertex_handle                                            Vertex_handle;
typedef Traits::Side_of_original_octagon                                        Side_of_original_octagon;
typedef Traits::Construct_hyperbolic_segment_2                                  Construct_hyperbolic_segment_2;
typedef Traits::Construct_inexact_intersection_2                                Construct_inexact_intersection_2;
typedef Triangulation::Hyperbolic_translation                                   Hyperbolic_translation;
typedef Triangulation::Locate_type                                              Locate_type;
typedef Triangulation::Face_handle                                              Face_handle;
typedef Triangulation::Face_iterator                                            Face_iterator;
typedef Traits::Construct_point_2                                               Make_point;

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Periodic_4_hyperbolic_Delaunay_triangulation_2
{
  Q_OBJECT

private:

  Point                                                              movingPoint;

  Triangulation                                                      dt;
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

  void        initialize_animation_parameters();
  void        animate();
  double      updateTime();
  Point       get_image(Point, Point, double);
  Circle_2    poincare;
  Face_handle last_location;
  double      timestep;
  Hyperbolic_translation      last_loc_translation;

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


#include <internal/hyperbolic_free_motion_animation.h>


MainWindow::MainWindow()
  : DemosMainWindow(), dt(Traits())
{

  initialize_animation_parameters();

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
  this->addAboutDemo(":/cgal/help/icons/about_P4HD_Triangulation_2.html");
  this->addAboutCGAL();

}


void
MainWindow::processInput(CGAL::Object o)
{

  Point p;
  if(CGAL::assign(p, o)){
    for(Triangulation::Vertex_iterator vi = dt.vertices_begin(); vi != dt.vertices_end(); vi++) {
      vi->clear_translation();
    }

    dt.insert(p);
    dt.try_to_remove_dummy_vertices();
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
  if(checked) {
    dgi->setVisibleDemo(true);
    go = true;
    animate();
  } else {
    dgi->setVisibleDemo(false);
    go = false;
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
           (std::numeric_limits<int>::max)(),
           1,
           &ok);

  if(!ok) {
    return;
  }

  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);


  typedef CGAL::Cartesian<double>::Point_2                   Point_double;
  typedef CGAL::Creator_uniform_2<double, Point_double >     Creator;

  CGAL::Random_points_in_disc_2<Point_double, Creator> g(0.85);

  Traits::Side_of_original_octagon pred;
  std::vector<Point> pts;
  int cnt = 0;
  do {
      Point_double pd = *(++g);
      if(pred(pd) != CGAL::ON_UNBOUNDED_SIDE) {
          Point pt = Point(pd.x(), pd.y());
          pts.push_back(pt);
          cnt++;
      }
  } while(cnt < number_of_points);


  CGAL::Timer tt;
  tt.start();
  dt.insert(pts.begin(), pts.end());
  tt.stop();

  std::cout << "Time elapsed for the insertion of " << number_of_points << " points: " << tt.time() << " secs." << std::endl;
  std::cout << "Number of vertices in the triangulation: " << dt.number_of_vertices() << std::endl;

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
  std::cout << "Read " << points.size() << " input points." << std::endl << "Inserting!" << std::endl;
  dt.insert(points.begin(), points.end());
  std::cout << "Done inserting!" << std::endl;

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

#include "P4HDT2.moc"
