#include <fstream>

// CGAL headers
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h> 
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/point_generators_2.h>
// unique words
#include <CGAL/Square_root_2_field.h>
#include <CGAL/Hyperbolic_octagon_group.h>
#include <CGAL/Hyperbolic_random_points_in_disc_2.h>
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
#include "CGAL/Qt/TriangulationCircumcircle.h"
#include "CGAL/Qt/TriangulationMovingPoint.h"
#include "CGAL/Qt/TriangulationConflictZone.h"
#include "CGAL/Qt/TriangulationRemoveVertex.h"
#include "CGAL/Qt/TriangulationPointInputAndConflictZone.h"
#include <CGAL/Qt/VoronoiGraphicsItem.h>
#include <CGAL/Qt/TriangulationGraphicsItemWithColorInfo.h>     // Visualise color
#include <CGAL/TranslationInfo.h>                               // Store color
#include <CGAL/Qt/DemosMainWindow.h>


#define OPTION_INSERT_DUMMY_POINTS 0


// dummy points
#if OPTION_INSERT_DUMMY_POINTS == 1
  #include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_dummy.h>
#endif
  
// the two base classes
#include "ui_Periodic_4_hyperbolic_Delaunay_triangulation_2.h"




typedef CGAL::Exact_predicates_inexact_constructions_kernel   InR;
typedef CGAL::Exact_predicates_exact_constructions_kernel     R;

typedef CGAL::Triangulation_hyperbolic_traits_2<R>            THT2;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<R> 
                                                              K;
typedef K::Point_2                                            Point;
// keep color
typedef TranslationInfo<std::string>                          Vb_info;
typedef CGAL::Triangulation_vertex_base_with_info_2< Vb_info, K > 
                                                              Vb;
typedef CGAL::Triangulation_face_base_with_info_2 <CGAL::Hyperbolic_face_info_2, K > 
                                                              Fb;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2< K, CGAL::Triangulation_data_structure_2<Vb, Fb> > 
                                                              Delaunay;
typedef Delaunay::Vertex_handle                               Vertex_handle;


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


void recurr(vector<Octagon_group>& v, vector<Octagon_group> g, int depth = 1) {
  if (depth > 1) {
    
    recurr(v, g, depth-1);

    vector<Octagon_group> tmp;
    vector<string> tmpw;
    for (int i = 0; i < v.size(); i++) {
      tmp.push_back(v[i]);
    }

    for (int i = 0; i < tmp.size(); i++) {
      for (int j = 0; j < g.size(); j++) {
        v.push_back(tmp[i]*g[j]);
      }
    }
  } else {
    for (int i = 0; i < g.size(); i++) {
      v.push_back(g[i]);
    }
  }
}


void my_unique_words(std::vector<Point>& p, Point input, int depth) {
  std::vector<Octagon_group> g;
  get_generators(g);
  std::vector<Octagon_group> v;
  recurr(v, g, depth);
  std::set<Octagon_group> s;

  for (int i = 0; i < v.size(); i++) {
    s.insert( v[i] );
  }

  //cout << "Original point and images: " << endl;
  //cout << input.x() << ", " << input.y() << endl;
  for (set<Octagon_group>::iterator it = s.begin(); it != s.end(); it++) {
    Octagon_group m = *it;
    pair<double, double> res;
    res = m.apply(to_double(input.x()), to_double(input.y()));
    //cout << res.first << ", " << res.second << endl;
    p.push_back( Point(res.first, res.second) );
  }

}



void apply_unique_words(std::vector<Point>& points, Point input = Point(0, 0), double threshold = 10, int word_length = 6, double d = .998)
{

  cout << "apply_unique_words called with threshold = " << threshold << ", word_length = " << word_length << ", d = " << d << endl;

  static vector<Octagon_group> unique_words;
  static bool generated = false;
  if(generated == false) {
    generate_unique_words(unique_words, threshold, word_length);
    generated = true;
  }
  
  //points.resize(unique_words.size());
  pair<double, double> res;
  for(size_t i = 0; i < unique_words.size(); i++) {
    pair<double, double> res;
    res = unique_words[i].apply(to_double(input.x()), to_double(input.y()));

    double dist = res.first*res.first + res.second*res.second;
    if(dist < d) {
      points.push_back( Point(res.first, res.second) );
    }
  }
}

void apply_unique_words_G(std::vector<Point>& points, Point input = Point(0, 0), double threshold = 6/*13.5*/, int word_length = 6/*20*/)
{
  static vector<Octagon_group> unique_words;
  static bool generated = false;
  if(generated == false) {
    generate_unique_words(unique_words, threshold, word_length);
    generated = true;
    
    // to generate all words
    /*
    ofstream fwords("m_w");
    fwords << "local words;\nPrint(\"words!\\n\");\nwords := [\n";
    for(size_t i = 0; i < unique_words.size() - 1; i++) {
      fwords << unique_words[i].label << "," <<endl;
    }
    fwords << unique_words[unique_words.size() - 1].label << endl;
    fwords << "];\nreturn words;";
    fwords.close();
    */
    //
  }
  
  ifstream findices("indicesG8_1");
  long index;
  vector<long> indices;
  while(findices >> index) {
    indices.push_back(index);
  }
  cout << "indices size " << indices.size() << endl; 

  pair<double, double> res;
  for(size_t i = 0; i < indices.size(); i++) {
    pair<double, double> res;
    if(indices[i] < unique_words.size() /*&& unique_words[indices[i]].length() <  13.*/) {
      res = unique_words[indices[i]].apply(to_double(input.x()), to_double(input.y()));
      points.push_back(Point(res.first, res.second));
    }
  }
  cout << "nb of points to insert " << points.size() << endl;
}


class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Periodic_4_hyperbolic_Delaunay_triangulation_2
{
  Q_OBJECT
  
private:  

  int              recursion_depth;
  int              cidx;
  std::vector<int> ccol;

  Delaunay                                                      dt;
  QGraphicsEllipseItem                                        * disk;
  QGraphicsScene                                                scene;  

  CGAL::Qt::TriangulationGraphicsItem<Delaunay>               * dgi;
  CGAL::Qt::VoronoiGraphicsItem<Delaunay>                     * vgi;
  
  // for drawing Voronoi diagram of the orbit of the origin
  CGAL::Qt::VoronoiGraphicsItem<Delaunay>                     * origin_vgi;

  CGAL::Qt::TriangulationMovingPoint<Delaunay>                * mp;
  CGAL::Qt::TriangulationConflictZone<Delaunay>               * cz;
  CGAL::Qt::TriangulationRemoveVertex<Delaunay>               * trv;
  CGAL::Qt::TriangulationPointInputAndConflictZone<Delaunay>  * pi;
  CGAL::Qt::TriangulationCircumcircle<Delaunay>               * tcc;
public:
  MainWindow();

public slots:

  void processInput(CGAL::Object o);

  void on_actionMovingPoint_toggled(bool checked);

  void on_actionShowConflictZone_toggled(bool checked);

  void on_actionCircumcenter_toggled(bool checked);

  void on_actionShowDelaunay_toggled(bool checked);

  void on_actionShowVoronoi_toggled(bool checked);

  void on_actionInsertPoint_toggled(bool checked);
  
  void on_actionInsertRandomPoints_triggered();

  void on_actionModifyDepth_triggered();

  void on_actionLoadPoints_triggered();

  void on_actionSavePoints_triggered();

  void on_actionClear_triggered();

  void on_actionRecenter_triggered();

  virtual void open(QString fileName);

signals:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow(), dt(K(1))
{
  recursion_depth = 1;
  cidx = 0;
  for (int i = 0; i < 10; i++)
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
  pen.setBrush(QColor(200, 200, 0));
  disk->setPen(pen);

  scene.addItem(disk);
  
  // another input point, instead of the origin
  
  double phi = CGAL_PI / 8.;
  double psi = CGAL_PI / 3.;
  double rho = std::sqrt(cos(psi)*cos(psi) - sin(phi)*sin(phi));
  
  Point origin = Point(0, 0);
  const Point a(cos(phi)*cos(phi + psi)/rho, sin(phi)*cos(phi + psi)/rho);
  
  
  // dt to form the octagon tessellation
  vector<Point> origin_orbit;
  apply_unique_words(origin_orbit, origin, 10, 6);
  origin_orbit.push_back(Point(0, 0));
  // for(long i = 0; i < origin_orbit.size(); i++) {
  //   cout << origin_orbit[i] << endl;
  // }
  cout << "nb of points on the orbit of the origin: " << origin_orbit.size() << endl;

  Delaunay* dtO = new Delaunay(K(1));
  dtO->insert(origin_orbit.begin(), origin_orbit.end());
  
  origin_vgi = new CGAL::Qt::VoronoiGraphicsItem<Delaunay>(dtO);
  origin_vgi->setVisible(true);
  
  QObject::connect(this, SIGNAL(changed()),
                   origin_vgi, SLOT(modelChanged()));
  
  QColor br(149, 179, 179);
  origin_vgi->setEdgesPen(QPen(br, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(origin_vgi);
  
  
  // Add a GraphicItem for the Delaunay triangulation
  dgi = new CGAL::Qt::TriangulationGraphicsItem<Delaunay>(&dt);

  QObject::connect(this, SIGNAL(changed()),
		   dgi, SLOT(modelChanged()));

  dgi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  dgi->setEdgesPen(QPen(QColor(200, 200, 0), 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(dgi);

  // Add a GraphicItem for the Voronoi diagram
  vgi = new CGAL::Qt::VoronoiGraphicsItem<Delaunay>(&dt);

  QObject::connect(this, SIGNAL(changed()),
		   vgi, SLOT(modelChanged()));

  QColor brown(139, 69, 19);
  vgi->setEdgesPen(QPen(brown, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(vgi);
  vgi->hide();

  // Setup input handlers. They get events before the scene gets them
  // and the input they generate is passed to the triangulation with 
  // the signal/slot mechanism    
  pi = new CGAL::Qt::TriangulationPointInputAndConflictZone<Delaunay>(&scene, &dt, this );

  QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
		   this, SLOT(processInput(CGAL::Object)));

  mp = new CGAL::Qt::TriangulationMovingPoint<Delaunay>(&dt, this);
  // TriangulationMovingPoint<Delaunay> emits a modelChanged() signal each
  // time the moving point moves.
  // The following connection is for the purpose of emitting changed().
  QObject::connect(mp, SIGNAL(modelChanged()),
		   this, SIGNAL(changed()));

  trv = new CGAL::Qt::TriangulationRemoveVertex<Delaunay>(&dt, this);
  QObject::connect(trv, SIGNAL(modelChanged()),
		   this, SIGNAL(changed()));

  tcc = new CGAL::Qt::TriangulationCircumcircle<Delaunay>(&scene, &dt, this);
  tcc->setPen(QPen(Qt::red, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

  cz = new CGAL::Qt::TriangulationConflictZone<Delaunay>(&scene, &dt, this);

  // 
  // Manual handling of actions
  //

  QObject::connect(this->actionQuit, SIGNAL(triggered()), 
		   this, SLOT(close()));

  // We put mutually exclusive actions in an QActionGroup
  QActionGroup* ag = new QActionGroup(this);
  ag->addAction(this->actionInsertPoint);
  ag->addAction(this->actionMovingPoint);
  ag->addAction(this->actionCircumcenter);
  ag->addAction(this->actionShowConflictZone);

  // Check two actions 
  this->actionInsertPoint->setChecked(true);
  this->actionShowDelaunay->setChecked(true);

  // //
  // // Setup the scene and the view
  // //
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  scene.setSceneRect(left_top_corner_x, left_top_corner_y, width, height);
  this->graphicsView->setScene(&scene);
  this->graphicsView->setMouseTracking(true);
  this->graphicsView->shear(230, 230);

  // // The navigation adds zooming and translation functionality to the
  // // QGraphicsView
  this->addNavigation(this->graphicsView);
  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_Delaunay_triangulation_2.html");
  this->addAboutCGAL();

  this->addRecentFiles(this->menuFile, this->actionQuit);
  connect(this, SIGNAL(openRecentFile(QString)),
	        this, SLOT(open(QString)));


#if OPTION_INSERT_DUMMY_POINTS == 1

  std::vector<Point> pts;
  cout << "Inserting dummy points now! " << endl;
  dt.insert_dummy_points(pts);
  for (int i = 0; i < pts.size(); i++) {
    processInput(make_object(pts[i]));
  }
  cout << "Dummy points inserted! " << endl;
  emit(changed());

#endif

}


void
MainWindow::processInput(CGAL::Object o)
{

  typedef CGAL::Exact_predicates_inexact_constructions_kernel GT;
  typedef GT::Point_2 Point_2;
  Point_2 pp;
  if (CGAL::assign(pp, o)) {
  	Point tmp(pp.x(), pp.y());
  	o = make_object(tmp);
  }

  Point p;
  if(CGAL::assign(p, o)){
    QPointF qp(CGAL::to_double(p.x()), CGAL::to_double(p.y()));
    
    // note that if the point is on the boundary then the disk contains the point
    if(disk->contains(qp)){
      //dt.insert(p);
      
      //delete
      vector<Point> points;
      //apply_unique_words(points, p, 4, 1, .998);
      my_unique_words(points, p, recursion_depth);

      points.push_back(p);
      Vertex_handle v;
      for(size_t j = 0; j < points.size() ; j++) {
        v = dt.insert(points[j]);
        v->info().setColor(ccol[cidx]);
      }
      cidx = (cidx + 1) % ccol.size();
      //
    }
    // delete
    else {
      static double phi = CGAL_PI / 8.;
      static double psi = CGAL_PI / 3.;
      static double rho = std::sqrt(cos(psi)*cos(psi) - sin(phi)*sin(phi));
      
      static Point origin = Point(0, 0);
      static Point a(cos(phi)*cos(phi + psi)/rho, sin(phi)*cos(phi + psi)/rho);
      static Point b((cos(psi)-sin(phi))/rho, 0.);
      
      static Point current_point_a = origin;
      static Point current_point_b = origin;
      static double dT = 0.05;
      
      static double dx_a = dT * CGAL::to_double(a.x());
      static double dy_a = dT * CGAL::to_double(a.y());
      current_point_a = Point(current_point_a.x() + dx_a, current_point_a.y() + dy_a);
      if(current_point_a.x() > a.x()) {
        current_point_a = origin;
      }
      
      static double dx_b = dT * CGAL::to_double(b.x());
      static double dy_b = dT * CGAL::to_double(b.y());
      current_point_b = Point(current_point_b.x() + dx_b, current_point_b.y() + dy_b);
      std::cout << current_point_b.x() << " : " << current_point_b.y() << std::endl;
      if(current_point_b.x() > b.x()) {
        current_point_b = origin;
      }
      
      vector<Point> points;
    
      // current_point_b <-> current_point_b
      Point current_point = current_point_a;
      apply_unique_words(points, current_point);
      points.push_back(current_point);
      dt.clear();
      dt.insert(points.begin(), points.end());
      
    }
    
  }
  emit(changed());

  cout << "v = " << dt.number_of_vertices() << ", f = " << dt.number_of_faces() << endl;

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
MainWindow::on_actionMovingPoint_toggled(bool checked)
{

  if(checked){
    scene.installEventFilter(mp);
  } else {
    scene.removeEventFilter(mp);
  }
}


void
MainWindow::on_actionShowConflictZone_toggled(bool checked)
{

  if(checked){
    scene.installEventFilter(cz);
  } else {
    scene.removeEventFilter(cz);
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
MainWindow::on_actionShowDelaunay_toggled(bool checked)
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

  typedef CGAL::Exact_predicates_inexact_constructions_kernel GT;
  typedef GT::Point_2 Point_2;
  typedef GT::FT FT;

  vector<Point_2> pts;
  Hyperbolic_random_points_in_disc_2<GT>(pts, number_of_points);

  for(int i = 0; i < number_of_points; ++i){
    processInput(make_object(pts[i]));
  }
  QApplication::restoreOverrideCursor();
  emit(changed());
}


void
MainWindow::on_actionModifyDepth_triggered()
{
  bool ok = false;
  const int result = 
    QInputDialog::getInt(this, 
                        tr("Modify recursion depth"),
                        tr("Enter new recursion depth"),
           recursion_depth,
           0,
           10,
           1,
           &ok);

  if(!ok) {
    return;
  }

  recursion_depth = result;
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
  
  K::Point p;
  std::vector<K::Point> points;
  while(ifs >> p) {
    points.push_back(p);
  }
  dt.insert(points.begin(), points.end());

  // default cursor
  QApplication::restoreOverrideCursor();
  this->addToRecentFiles(fileName);
  actionRecenter->trigger();
  emit(changed());
    
}

void
MainWindow::on_actionSavePoints_triggered()
{
/*  // dump points
  QString fileName = QFileDialog::getSaveFileName(this,
						  tr("Save points"),
						  ".");
  if(! fileName.isEmpty()){
    std::ofstream ofs(qPrintable(fileName));
    for(Delaunay::Finite_vertices_iterator 
          vit = dt.finite_vertices_begin(),
          end = dt.finite_vertices_end();
        vit!= end; ++vit)
    {
      ofs << vit->point() << std::endl;
    }
  }
*/  
  
  // take a snapshot
  std::cout << "snapshot...";
  
  const QRect viewerRect = graphicsView->mapFromScene(scene.sceneRect()).boundingRect();
  const QRect imageRect = QRect(QPoint(0, 0), viewerRect.size()); 
  
  QImage snapshot(imageRect.size(), QImage::Format_ARGB32);//QImage::Format_ARGB32_Premultiplied
  
  QPainter painter(&snapshot);
  painter.setRenderHint(QPainter::Antialiasing);
  //painter.setRenderHint(QPainter::SmoothPixmapTransform);
  
  graphicsView->render(&painter, imageRect, viewerRect);
  bool saved = snapshot.save("mysnap.png", "PNG", 100);
  assert(saved == true);
  
  std::cout << "done" << std::endl;
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
  app.setApplicationName("Periodic_4_hyperbolic_Delaunay_triangulation_2 demo");

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  Q_INIT_RESOURCE(File);
  Q_INIT_RESOURCE(Triangulation_2);
  Q_INIT_RESOURCE(Input);
  Q_INIT_RESOURCE(CGAL);

  MainWindow mainWindow;
  mainWindow.show();

  QStringList args = app.arguments();
  args.removeAt(0);
  Q_FOREACH(QString filename, args) {
    mainWindow.open(filename);
  }

  return app.exec();
}
