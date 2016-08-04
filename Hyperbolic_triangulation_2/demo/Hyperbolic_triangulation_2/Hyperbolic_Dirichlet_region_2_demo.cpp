#include <fstream>

// CGAL headers
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_triangulation_traits_2.h>

// to be deleted
#include <CGAL/Qt/HyperbolicPainterOstream.h>
//

#include <CGAL/point_generators_2.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <QGraphicsEllipseItem>

#include <string>

// GraphicsView items and event filters (input classes)
#include "TriangulationCircumcircle.h"

#include "TriangulationMovingPoint.h"
#include "TriangulationConflictZone.h"
#include "TriangulationRemoveVertex.h"
#include "TriangulationPointInputAndConflictZone.h"
//#include <CGAL/Qt/TriangulationGraphicsItem.h>
#include <CGAL/Qt/VoronoiGraphicsItem.h>

// store color
#include <TranslationInfo.h>
// visualise color
#include <CGAL/Qt/TriangulationGraphicsItemWithColorInfo.h>

// unique words
#include <temp.h>

// dummy points
#include <CGAL/Periodic_2_hyperbolic_triangulation_dummy.h>

// for filtering
#include <set>

// for viewportsBbox
#include <CGAL/Qt/utility.h>
  
// the two base classes
#include "ui_Hyperbolic_translations_2.h"
#include <CGAL/Qt/DemosMainWindow.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel InR;
typedef CGAL::Exact_predicates_exact_constructions_kernel R;
typedef CGAL::Hyperbolic_triangulation_traits_2<R> K;

typedef K::Point_2 Point_2;
typedef Point_2 Point;
typedef K::Iso_rectangle_2 Iso_rectangle_2;

// keep color
typedef TranslationInfo<std::string> Vb_info;

typedef CGAL::Triangulation_vertex_base_with_info_2< Vb_info, K > Vb;
typedef CGAL::Hyperbolic_triangulation_face_base_2<K> Fb;

typedef CGAL::Hyperbolic_Delaunay_triangulation_2< K, CGAL::Triangulation_data_structure_2<Vb, Fb> > Delaunay;

typedef Delaunay::Vertex_handle Vertex_handle;
 
//typedef CGAL::Hyperbolic_Delaunay_triangulation_2<K> Delaunay;

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

void apply_unique_words(std::vector<Point>& points, Point input = Point(0, 0), double threshold = 10, int word_length = 4/*13*/, double d = 1.)
{
  static vector<OctagonMatrix> unique_words;
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

void apply_unique_words_G(std::vector<Point>& points, Point input = Point(0, 0), double threshold = 13.5, int word_length = 20)
{
  static vector<OctagonMatrix> unique_words;
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
  public Ui::Delaunay_triangulation_2
{
  Q_OBJECT
  
private:  
  Delaunay dt;
  QGraphicsEllipseItem* disk;
  QGraphicsScene scene;  

  CGAL::Qt::TriangulationGraphicsItem<Delaunay> * dgi;
  CGAL::Qt::VoronoiGraphicsItem<Delaunay> * vgi;
  
  // for drawing Voronoi diagram of the orbit of the origin
  CGAL::Qt::VoronoiGraphicsItem<Delaunay> * origin_vgi;

  CGAL::Qt::TriangulationMovingPoint<Delaunay> * mp;
  CGAL::Qt::TriangulationConflictZone<Delaunay> * cz;
  CGAL::Qt::TriangulationRemoveVertex<Delaunay> * trv;
  CGAL::Qt::TriangulationPointInputAndConflictZone<Delaunay> * pi;
  CGAL::Qt::TriangulationCircumcircle<Delaunay> *tcc;
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

  void on_actionLoadPoints_triggered();

  void on_actionSavePoints_triggered();

  void on_actionClear_triggered();

  void on_actionRecenter_triggered();

  virtual void open(QString fileName);

signals:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow()
{
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
  
  //QPen diskPen = disk->pen();
  //diskPen.setWidth(1);
  //std::cout << "width "<< diskPen.width() << std::endl;
  //disk->setPen(diskPen);
  
  scene.addItem(disk);
  
  // another input point, instead of the origin
  
  double phi = CGAL_PI / 8.;
  double psi = CGAL_PI / 3.;
  double rho = std::sqrt(cos(psi)*cos(psi) - sin(phi)*sin(phi));
  
  Point origin = Point(0, 0);
  const Point_2 a(cos(phi)*cos(phi + psi)/rho, sin(phi)*cos(phi + psi)/rho);
  
  
  // dt to form the octagon tessellation
  vector<Point> origin_orbit;
  //apply_unique_words(origin_orbit, origin, 8, 8, 0.998);
  origin_orbit.push_back(Point(0, 0));
  cout << "nb of points on the orbit of the origin: " << origin_orbit.size() << endl;
  for(long i = 0; i < origin_orbit.size(); i++) {
    cout << origin_orbit[i] << endl;
  }
  
  Delaunay* dtO = new Delaunay();
  dtO->insert(origin_orbit.begin(), origin_orbit.end());
  
  origin_vgi = new CGAL::Qt::VoronoiGraphicsItem<Delaunay>(dtO);
  origin_vgi->setVisible(true);
  
  QObject::connect(this, SIGNAL(changed()),
                   origin_vgi, SLOT(modelChanged()));
  
  QColor br(139, 69, 19);
  origin_vgi->setEdgesPen(QPen(br, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(origin_vgi);
  
  // extra points
  vector<Point> inner;
  /*
  inner.push_back(Point(0, 0));
  inner.push_back(Point(0.3, 0.3));
  inner.push_back(Point(-0.4, 0.2));
  inner.push_back(Point(-0.1, -0.25));
  inner.push_back(Point(0.2, -0.4));
  inner.push_back(Point(0.02, -0.1572));
  inner.push_back(Point(-0.098, 0.137));
  inner.push_back(Point(-0.1867, 0.005));
  inner.push_back(Point(-0.24, -0.188));*/
  
  vector<Point> points;
  Vertex_handle v;
  Point p;
  for(size_t i = 0, col = 0; i < inner.size(); i++, col++) {
    points.resize(0);
    //apply_unique_words(points, inner[i], 8, 8, /*0.992*/0.99);
    points.push_back(inner[i]);
    for(size_t j = 0; j < points.size(); j++) {
      v = dt.insert(points[j]);
      if (col == 5) {
        col = 6;
      }
      v->info().setColor(col);
    }
  }
  
  // insert many points of GXX O
  
  /*
  vector<Point> points;
  apply_unique_words_G(points, origin);
  
  Delaunay temp_dt = Delaunay(K(1));
  temp_dt.insert(points.begin(), points.end());
  
  points.push_back(origin);
  Delaunay::Vertex_handle origin_v = temp_dt.insert(origin);
  
  vector<Point> incident_to_origin;
  Delaunay::Vertex_circulator start = temp_dt.incident_vertices(origin_v);
  Delaunay::Vertex_circulator next = start;
  do {
    incident_to_origin.push_back(next->point());
    next++;
  } while(next != start);
  incident_to_origin.push_back(origin);
  
  dt.insert(incident_to_origin.begin(), incident_to_origin.end());
  */
   
  /*
  vector<Point> points;
  apply_unique_words(points, origin, 7, 4);
  points.push_back(origin);
  dt.insert(points.begin(), points.end());
  */
  
  // insert many points
/*
  vector<Point> points;
  apply_unique_words(points, Point(0, 0), 5);
  points.push_back(Point(0, 0));
  
  Vertex_handle v;
  for(size_t j = 0; j < points.size(); j++) {
    v = dt.insert(points[j]);
    v->info().setColor(6);
  }
  
  // dummy points
  vector<Point> inner, on_boundary, on_vertex;
  CGAL::compute_redundant_dummy_points<K>(inner, on_boundary, on_vertex);
  
  // filter
  set<Point, PointsComparator> filter;
  set<Point, PointsComparator>::iterator it;
  
  // inner
  for(size_t i = 0; i < inner.size() -1; i++) {
    points.resize(0);
    apply_unique_words(points, inner[i], 5);
    points.push_back(inner[i]);
    for(size_t j = 0; j < points.size(); j++) {
      v = dt.insert(points[j]);
      v->info().setColor(1);
    }
  }
  
  // on_boundary
  filter.clear();
  long nb = 0;
  for(size_t i = 0; i < on_boundary.size(); i++) {
    points.resize(0);
    apply_unique_words(points, on_boundary[i], 5);
    points.push_back(on_boundary[i]);
    nb += points.size();
    filter.insert(points.begin(), points.end());
  }
  cout << "before " << nb;
  cout << "after " << filter.size();
  
  for(it = filter.begin(); it != filter.end(); it++) {
    v = dt.insert(*it);
    v->info().setColor(7);
  }
  
  // on_vertex
  filter.clear();
  nb = 0;
  for(size_t i = 0; i < on_vertex.size(); i++) {
    points.resize(0);
    apply_unique_words(points, on_vertex[i], 5);
    points.push_back(on_vertex[i]);
    nb += points.size();
    filter.insert(points.begin(), points.end());
  }
  cout << "before " << nb;
  cout << "after " << filter.size();
  
  for(it = filter.begin(); it != filter.end(); it++) {
    v = dt.insert(*it);
    v->info().setColor(3);
  }
*/    
  // Add a GraphicItem for the Delaunay triangulation
  dgi = new CGAL::Qt::TriangulationGraphicsItem<Delaunay>(&dt);

  QObject::connect(this, SIGNAL(changed()),
		   dgi, SLOT(modelChanged()));

  dgi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
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

  //
  // Setup the scene and the view
  //
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  scene.setSceneRect(left_top_corner_x, left_top_corner_y, width, height);
  this->graphicsView->setScene(&scene);
  this->graphicsView->setMouseTracking(true);
  
  // we want to adjust the coordinates of QGraphicsView to the coordinates of QGraphicsScene
  // the following line must do this:
  //   this->graphicsView->fitInView( scene.sceneRect(), Qt::KeepAspectRatio);
  // It does not do this sufficiently well.
  // Current solution:
  
  QRect viewerRect = graphicsView->mapFromScene(scene.sceneRect()).boundingRect();
  std::cout << "resolution before " << viewerRect.width() << " " << viewerRect.height() << std::endl;
  
  // shear 230 230
  this->graphicsView->shear(407, 407);
  //this->graphicsView->shear(307, 307);
  
  viewerRect = graphicsView->mapFromScene(scene.sceneRect()).boundingRect();
  std::cout << "resolution after " << viewerRect.width() << " " << viewerRect.height() << std::endl;
  
  // Turn the vertical axis upside down
  this->graphicsView->matrix().scale(1, -1);
                                                      
  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_Delaunay_triangulation_2.html");
  this->addAboutCGAL();

  this->addRecentFiles(this->menuFile, this->actionQuit);
  connect(this, SIGNAL(openRecentFile(QString)),
	  this, SLOT(open(QString)));
}


void
MainWindow::processInput(CGAL::Object o)
{
  Point_2 p;
  if(CGAL::assign(p, o)){
    QPointF qp(CGAL::to_double(p.x()), CGAL::to_double(p.y()));
    
    // note that if the point is on the boundary then the disk contains the point
    if(disk->contains(qp)){
      //dt.insert(p);
      
      //delete
      vector<Point> points;
      apply_unique_words/*_G*/(points, p, 6, 4);;
      points.push_back(p);
      Vertex_handle v;
      for(size_t j = 0; j < points.size(); j++) {
        v = dt.insert(points[j]);
        v->info().setColor(1);
      }
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
      apply_unique_words/*_G*/(points, current_point, 6, 4);
      points.push_back(current_point);
      dt.clear();
      dt.insert(points.begin(), points.end());
      
      /*
      static vector<Vertex_handle> old_vertices;
      if(dt.number_of_vertices() != 0) {
        for(size_t i = 0; i < old_vertices.size(); i++) {
          dt.remove(old_vertices[i]);
        }
        // we call it because dt.remove() is not enough in hyperbolic case
        dt.mark_faces();
      }
      old_vertices.resize(0);
      
      Vertex_handle new_vertex;
      for(size_t i = 0; i < points.size(); i++) {
        new_vertex = dt.insert(points[i]);
        old_vertices.push_back(new_vertex);
      }
      std::cout << "old_vertices.size()" << old_vertices.size() << std::endl;
      */
      /*
      Vertex_handle v;
      for(size_t j = 0; j < points.size(); j++) {
        v = dt.insert(points[j]);
        v->info().setColor(6);
      }
      
      static vector<Point> origin_images;
      static bool filled = false;
      if(filled == false) {
        apply_unique_words(origin_images, Point(0, 0), 6, 4);
        origin_images.push_back(Point(0, 0));
        filled = true;
      }
      dt.insert(origin_images.begin(), origin_images.end());
      */
    }
    
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
  QRectF rect = CGAL::Qt::viewportsBbox(&scene);
  CGAL::Qt::Converter<K> convert;  
  Iso_rectangle_2 isor = convert(rect);
  CGAL::Random_points_in_iso_rectangle_2<Point_2> pg(isor.min(), isor.max());
  bool ok = false;
  const int number_of_points = 
    QInputDialog::getInteger(this, 
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
  std::vector<Point_2> points;
  points.reserve(number_of_points);
  for(int i = 0; i < number_of_points; ++i){
    points.push_back(*pg++);
  }
  dt.insert(points.begin(), points.end());
  // default cursor
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
  
  K::Point_2 p;
  std::vector<K::Point_2> points;
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


#include "Hyperbolic_Dirichlet_region_2_demo.moc"

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Delaunay_triangulation_2 demo");

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
