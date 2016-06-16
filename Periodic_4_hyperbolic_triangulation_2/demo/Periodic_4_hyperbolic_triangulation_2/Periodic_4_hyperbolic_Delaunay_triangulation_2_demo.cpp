#include <fstream>

// CGAL headers
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h> 
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h> 
#include <CGAL/point_generators_2.h>
#include <CGAL/Hyperbolic_octagon_translation_matrix.h>

// unique words
#include <CGAL/Square_root_2_field.h>
// to be deleted (iiordano: why?)
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
#include <CGAL/Qt/TriangulationMovingPoint.h>
#include <CGAL/Qt/TriangulationConflictZone.h>
#include <CGAL/Qt/TriangulationRemoveVertex.h>
#include <CGAL/Qt/TriangulationPointInputAndConflictZone.h>
#include <CGAL/Qt/VoronoiGraphicsItem.h>
#include <CGAL/Qt/TriangulationGraphicsItemWithColorInfo.h>     // Visualise color
//#include <CGAL/TranslationInfo.h>                               // Store color
#include <CGAL/Qt/DemosMainWindow.h>


#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>

#define INITIAL_RECURSION_DEPTH 4

// dummy points
#include <CGAL/Periodic_4_hyperbolic_triangulation_dummy_14.h>

// the two base classes
#include "ui_Periodic_4_hyperbolic_Delaunay_triangulation_2.h"



typedef CORE::Expr                                                              NT;         
typedef CGAL::Cartesian<NT>                                                     Kernel;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel>     Traits;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>            Triangulation;
typedef Hyperbolic_octagon_translation_matrix<Traits>                           Octagon_matrix;
typedef Kernel::Point_2                                                         Point;
typedef Triangulation::Vertex_handle                                            Vertex_handle;
typedef Traits::Side_of_fundamental_octagon                                     Side_of_fundamental_octagon;


struct PointsComparator {
  static double eps; 
  
  bool operator() (const Point& l, const Point& r) const
  {
    if(l.x() < r.x() - eps) {
      return true;
    }
    if(l.x() < r.x() + eps) {
      if(l.y() < r.y() - eps) {
        return true;
      }
    }
    return false;
  }
};


double PointsComparator::eps = 0.0001;
string glabels[] = { "a", "\\bar{b}", "c", "\\bar{d}", "\\bar{a}", "b", "\\bar{c}", "d" };


void recurr(vector<Octagon_matrix>& v, vector<Octagon_matrix> g, int depth = 1) {
  if (depth > 1) {
    
    recurr(v, g, depth-1);

    vector<Octagon_matrix> tmp;
    vector<string> tmpw;
    for (int i = 0; i < v.size(); i++) {
      tmp.push_back(v[i]);
    }

    for (int i = 0; i < tmp.size(); i++) {
      for (int j = 0; j < g.size(); j++) {
        v.push_back(tmp[i]*g[j]);
      }
    }
  } else if (depth == 1) {
    for (int i = 0; i < g.size(); i++) {
      v.push_back(g[i]);
    }
  } 
}


void my_unique_words(std::vector<Point>& p, Point input, int depth) {
  std::vector<Octagon_matrix> g;
  get_generators(g);
  std::vector<Octagon_matrix> v;
  recurr(v, g, depth);
  std::set<Octagon_matrix> s;

  for (int i = 0; i < v.size(); i++) {
    s.insert( v[i] );
  }

  cout << "Translating... " << endl;
  for (set<Octagon_matrix>::iterator it = s.begin(); it != s.end(); it++) {
    Octagon_matrix m = *it;
    Point res;
    res = m.apply(input);
    p.push_back( res );
  }
  cout << "Done! Now I need to draw " << p.size() << " points..." << endl;

}




class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Periodic_4_hyperbolic_Delaunay_triangulation_2
{
  Q_OBJECT
  
private:  


  int                                                           cidx;
  std::vector<int>                                              ccol;
  bool                                                          dummy_mode;

  Triangulation                                                 dt;
  QGraphicsEllipseItem                                        * disk;
  QGraphicsScene                                                scene;  

  CGAL::Qt::TriangulationGraphicsItem<Triangulation>          * dgi;
  
  
  // for drawing Voronoi diagram of the orbit of the origin

  CGAL::Qt::TriangulationMovingPoint<Triangulation>                * mp;
  CGAL::Qt::TriangulationRemoveVertex<Triangulation>               * trv;
  CGAL::Qt::TriangulationPointInputAndConflictZone<Triangulation>  * pi;
  CGAL::Qt::TriangulationCircumcircle<Triangulation>               * tcc;
public:
  MainWindow();

public slots:

  void processInput(CGAL::Object o);
  
  void on_actionCircumcenter_toggled(bool checked);

  void on_actionShowTriangulation_toggled(bool checked);

  void on_actionInsertPoint_toggled(bool checked);
  
  void on_actionInsertRandomPoints_triggered();

  void on_actionLoadPoints_triggered();

  void on_actionClear_triggered();

  void on_actionRecenter_triggered();


signals:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow(), dt(Traits(1))
{

  dt.insert_dummy_points();

  cidx = 0;
  for (int i = 0; i < 14; i++)
    ccol.push_back(i);
  
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
  
  // set disk
  disk = new QGraphicsEllipseItem(left_top_corner_x, left_top_corner_y, width, height);
  QPen pen;  // creates a default pen
  pen.setWidth(0);
  //pen.setBrush(Qt::black);
  pen.setBrush(Qt::black);
  disk->setPen(pen);

  scene.addItem(disk);
  
  // Add a GraphicItem for the Triangulation triangulation
  dgi = new CGAL::Qt::TriangulationGraphicsItem<Triangulation>(&dt);

  QObject::connect(this, SIGNAL(changed()),
		   dgi, SLOT(modelChanged()));

  dgi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  dgi->setEdgesPen(QPen(QColor(200, 200, 0), 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(dgi);

  // Setup input handlers. They get events before the scene gets them
  // and the input they generate is passed to the triangulation with 
  // the signal/slot mechanism    
  pi = new CGAL::Qt::TriangulationPointInputAndConflictZone<Triangulation>(&scene, &dt, this );

  QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
		   this, SLOT(processInput(CGAL::Object)));

  mp = new CGAL::Qt::TriangulationMovingPoint<Triangulation>(&dt, this);
  // TriangulationMovingPoint<Triangulation> emits a modelChanged() signal each
  // time the moving point moves.
  // The following connection is for the purpose of emitting changed().
  QObject::connect(mp, SIGNAL(modelChanged()),
		   this, SIGNAL(changed()));

  trv = new CGAL::Qt::TriangulationRemoveVertex<Triangulation>(&dt, this);
  QObject::connect(trv, SIGNAL(modelChanged()),
		   this, SIGNAL(changed()));

  tcc = new CGAL::Qt::TriangulationCircumcircle<Triangulation>(&scene, &dt, this);
  tcc->setPen(QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

  //cz = new CGAL::Qt::TriangulationConflictZone<Triangulation>(&scene, &dt, this);

  // 
  // Manual handling of actions
  //

  QObject::connect(this->actionQuit, SIGNAL(triggered()), 
		   this, SLOT(close()));

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

  //this->addRecentFiles(this->menuFile, this->actionQuit);
  //connect(this, SIGNAL(openRecentFile(QString)),
	//        this, SLOT(open(QString)));


}


void
MainWindow::processInput(CGAL::Object o)
{

  Point p;
  if(CGAL::assign(p, o)){
    Vertex_handle v = dt.insert(p);
  }
  emit(changed());
  assert(dt.is_valid());

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
    scene.installEventFilter(trv);
  } else {
    scene.removeEventFilter(pi);
    scene.removeEventFilter(trv);
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
MainWindow::on_actionClear_triggered()
{
  dt.clear();
  dt.insert_dummy_points();
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

  vector<Point> pts;
  pts.reserve(number_of_points);

  typedef CGAL::Creator_uniform_2<double, Point> Creator;

  CGAL::Random_points_in_disc_2<Point, Creator> g( 1.0 );
  CGAL::cpp11::copy_n( g, number_of_points, std::back_inserter(pts));

  for(int i = 0; i < number_of_points; ++i){
    processInput(make_object(pts[i]));
  }
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

  Q_INIT_RESOURCE(res);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Periodic_4_hyperbolic_Delaunay_triangulation_2 demo");

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  //Q_INIT_RESOURCE(File);
  //Q_INIT_RESOURCE(Triangulation_2);
  //Q_INIT_RESOURCE(Input);
  //Q_INIT_RESOURCE(CGAL);

  MainWindow mainWindow;
  mainWindow.show();

  QStringList args = app.arguments();
  args.removeAt(0);
  //Q_FOREACH(QString filename, args) {
    //mainWindow.open(filename);
  //}

  return app.exec();
}
