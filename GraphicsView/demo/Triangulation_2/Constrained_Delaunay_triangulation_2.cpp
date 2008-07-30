#include <fstream>
#include <vector>

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Timer.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>

// GraphicsView items and event filters (input classes)
#include "TriangulationCircumcircle.h"
#include <CGAL/Qt/GraphicsViewPolylineInput.h>
#include <CGAL/Qt/ConstrainedTriangulationGraphicsItem.h>
  
// the two base classes
#include "ui_Constrained_Delaunay_triangulation_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef CGAL::Triangulation_vertex_base_2<K>  Vertex_base;
//typedef CGAL::Constrained_triangulation_face_base_2<K> Face_base;
typedef CGAL::Delaunay_mesh_face_base_2<K> Face_base;
typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base>  TDS;
typedef CGAL::Exact_predicates_tag              Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;


typedef std::pair<std::vector<Point_2>*, int> Point_iterator;



template <typename Kernel, typename Iterator>
struct Sort_traits_2 {

  Kernel k;

  Sort_traits_2 (const Kernel &kernel = Kernel())
      : k (kernel)
  {}

  typedef Iterator Point_2;

  struct Less_x_2 {
    Kernel k;
    Less_x_2 (const Kernel &kernel = Kernel())
        : k (kernel)
    {}
    bool operator() (const Point_2 &p, const Point_2 &q) const
    {
      return k.less_x_2_object() ((*(p.first))[p.second], (*(q.first))[q.second]);
    }
  };

  Less_x_2
  less_x_2_object() const
  {
    return Less_x_2(k);
  }

  struct Less_y_2 {
    Kernel k;
    Less_y_2 (const Kernel &kernel = Kernel())
        : k (kernel)
    {}
    bool operator() (const Point_2 &p, const Point_2 &q) const
    {
      return k.less_y_2_object() ((*(p.first))[p.second], (*(q.first))[q.second]);
    }
  };


  Less_y_2
  less_y_2_object() const
  {
    return Less_y_2(k);
  }
};


typedef CGAL::Hilbert_sort_2<Sort_traits_2<K, Point_iterator> > Hilbert_sort_2;
typedef CGAL::Multiscale_sort<Hilbert_sort_2> Spatial_sort_2;




class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Constrained_Delaunay_triangulation_2
{
  Q_OBJECT
  
private:  
  CDT cdt; 
  QGraphicsScene scene;  

  CGAL::Qt::ConstrainedTriangulationGraphicsItem<CDT> * dgi;

  CGAL::Qt::GraphicsViewPolylineInput<K> * pi;
  CGAL::Qt::TriangulationCircumcircle<CDT> *tcc;
public:
  MainWindow();

private:
  template <typename Iterator> 
  void insert_polyline(Iterator b, Iterator e)
  {
    Point_2 p, q;
    CDT::Vertex_handle vh, wh;
    Iterator it = b;
    vh = cdt.insert(*it);
    p = *it;
    ++it;
    for(; it != e; ++it){
      q = *it;
      if(p != q){
        wh = cdt.insert(*it);
        cdt.insert_constraint(vh,wh);
        vh = wh;
        p = q;
      } else {
        std::cout << "duplicate point: " << p << std::endl; 
      }
    }
    emit(changed());
  }

protected slots:
void open(QString);

public slots:

  void processInput(CGAL::Object o);

  void on_actionShowDelaunay_toggled(bool checked);

  void on_actionInsertPolyline_toggled(bool checked);
  
  void on_actionCircumcenter_toggled(bool checked);

  void on_actionClear_triggered();

  void on_actionRecenter_triggered();

  void on_actionLoadConstraints_triggered();

  void loadPolyConstraints(QString);

  void loadEdgConstraints(QString);

  void on_actionSaveConstraints_triggered();

  void saveConstraints(QString);

  void on_actionMakeGabrielConform_triggered();

  void on_actionMakeDelaunayConform_triggered();

  void on_actionMakeDelaunayMesh_triggered();

  void on_actionInsertRandomPoints_triggered();

signals:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow()
{
  setupUi(this);

  // Add a GraphicItem for the CDT triangulation
  dgi = new CGAL::Qt::ConstrainedTriangulationGraphicsItem<CDT>(&cdt);
    
  QObject::connect(this, SIGNAL(changed()),
		   dgi, SLOT(modelChanged()));

  dgi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(dgi);

  // Setup input handlers. They get events before the scene gets them
  // and the input they generate is passed to the triangulation with 
  // the signal/slot mechanism    
  pi = new CGAL::Qt::GraphicsViewPolylineInput<K>(this, &scene, 0, true); // inputs polylines which are not closed

  QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
		   this, SLOT(processInput(CGAL::Object)));
    

  tcc = new CGAL::Qt::TriangulationCircumcircle<CDT>(&scene, &cdt, this);
  tcc->setPen(QPen(Qt::red, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  

  // 
  // Manual handling of actions
  //
  QObject::connect(this->actionExit, SIGNAL(triggered()), 
		   this, SLOT(close()));

  // We put mutually exclusive actions in an QActionGroup
  QActionGroup* ag = new QActionGroup(this);
  ag->addAction(this->actionInsertPolyline);

  // Check two actions 
  this->actionInsertPolyline->setChecked(true);
  this->actionShowDelaunay->setChecked(true);

  //
  // Setup the scene and the view
  //
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  scene.setSceneRect(-100, -100, 100, 100);
  this->graphicsView->setScene(&scene);
  this->graphicsView->setMouseTracking(true);

  // Turn the vertical axis upside down
  this->graphicsView->matrix().scale(1, -1);
                                                      
  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_Constrained_Delaunay_triangulation_2.html");
  this->addAboutCGAL();

  this->addRecentFiles(this->menuFile, this->actionExit);
  connect(this, SIGNAL(openRecentFile(QString)),
	  this, SLOT(open(QString)));
}


void
MainWindow::processInput(CGAL::Object o)
{
  std::list<Point_2> points;
  if(CGAL::assign(points, o)){
    if(points.size() == 1) {
      cdt.insert(points.front());
    }
    else {
      insert_polyline(points.begin(), points.end());
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
MainWindow::on_actionInsertPolyline_toggled(bool checked)
{
  if(checked){
    scene.installEventFilter(pi);
  } else {
    scene.removeEventFilter(pi);
  }
}


void
MainWindow::on_actionShowDelaunay_toggled(bool checked)
{
  dgi->setDrawEdges(checked);
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
MainWindow::on_actionClear_triggered()
{
  cdt.clear();
  emit(changed());
}


void 
MainWindow::open(QString fileName)
{
  if(! fileName.isEmpty()){
    if(fileName.endsWith(".poly")){
      loadPolyConstraints(fileName);
      this->addToRecentFiles(fileName);
    } else if(fileName.endsWith(".edg")){
      loadEdgConstraints(fileName);
      this->addToRecentFiles(fileName);
    }
  }
}

void
MainWindow::on_actionLoadConstraints_triggered()
{
  QString fileName = QFileDialog::getOpenFileName(this,
						  tr("Open Constraint File"),
						  ".",
						  tr("Edge files (*.edg)\n"
						     "Poly files (*.poly)"));
  open(fileName);
}

void
MainWindow::loadPolyConstraints(QString fileName)
{
  std::ifstream ifs(qPrintable(fileName));
  bool first=true;
  int n;
  ifs >> n;
  
  K::Point_2 p,q, qold;
  CDT::Vertex_handle vp, vq, vqold;
  while(ifs >> p) {
    ifs >> q;
    if((!first) && (p == qold)){
      vp = vqold;
    } else {
      vp = cdt.insert(p);
    }
    vq = cdt.insert(q, vp->face());
    cdt.insert_constraint(vp,vq);
    qold = q;
    vqold = vq;
    first = false;
  }

  emit(changed());
  actionRecenter->trigger();
}


void
MainWindow::loadEdgConstraints(QString fileName)
{
  CGAL::Timer tim;
  tim.start();
  std::ifstream ifs(qPrintable(fileName));
  bool first=true;
  int n;
  ifs >> n;
  
  K::Point_2 p,q, qold;
#if 1
  CDT::Vertex_handle vp, vq, vqold;
  while(ifs >> p) {
    ifs >> q;
    if(p == q){
      std::cout << "Ignore zero length segment" << std::endl;
      continue;
    }
    if((!first) && (p == qold)){
      vp = vqold;
    } else {
      vp = cdt.insert(p);
    }
    vq = cdt.insert(q, vp->face());
    cdt.insert_constraint(vp,vq);
    qold = q;
    vqold = vq;
    first = false;
  }

#else 

  Spatial_sort_2 sort_2;

  std::vector<Point_2> points;
  std::vector<bool> bop; // beginning of polyline
  std::vector<std::pair<std::vector<Point_2>*,int> > iterators;
  std::vector<CDT::Vertex_handle> vertices;
  CDT::Vertex_handle vh;

  points.reserve(n); // As n is the number of segments the vectors might become twice as big
  bop.reserve(n);

  while(ifs >> p){
    ifs >> q;
    if(p == q){
      std::cout << "Ignore zero length segment" << std::endl;
      continue;
    }
    if(first || (p != qold)){
      // start a new polyline
      bop.push_back(true);
      points.push_back(p);
    } 
    bop.push_back(false);
    points.push_back(q);
    first = false;
  }
  iterators.reserve(points.size());
  for(int i=0; i < points.size(); i++){
    iterators.push_back(std::make_pair(&points,i));
  }

  sort_2(iterators.begin(), iterators.end());

  // insert the points in the spatial sort order
  first = true;
  vertices.resize(points.size());
  for(std::vector<std::pair<std::vector<Point_2>*,int> >::iterator it = iterators.begin();
      it != iterators.end(); 
      it++) {
    if(first){
      vh = vertices[it->second] = cdt.insert(points[it->second]);
      first = false;
    } else {
      vh = vertices[it->second] = cdt.insert(points[it->second], vh->face());
    }
  }

  // insert the constraints
  CDT::Vertex_handle vp, vq;
  for(int i = 0; i < vertices.size(); i++){
    if(!bop[i]){
      cdt.insert_constraint(vertices[i-1], vertices[i]);
    }
  }
#endif
  tim.stop();
  statusBar()->showMessage(QString("Insertion took %1 seconds").arg(tim.time()), 2000);
  emit(changed());
  actionRecenter->trigger();
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
  QMessageBox::warning(this,
                       tr("saveConstraints"),
                       tr("Not implemented!"));
}


void
MainWindow::on_actionMakeGabrielConform_triggered()
{
  int nv = cdt.number_of_vertices();
  CGAL::make_conforming_Gabriel_2(cdt);
  nv = cdt.number_of_vertices() - nv;
  statusBar()->showMessage(QString("Added %1 vertices").arg(nv), 2000);
  emit(changed());
}


void
MainWindow::on_actionMakeDelaunayConform_triggered()
{
  int nv = cdt.number_of_vertices();
  CGAL::make_conforming_Delaunay_2(cdt);
  nv = cdt.number_of_vertices() - nv;
  statusBar()->showMessage(QString("Added %1 vertices").arg(nv), 2000);
  emit(changed());
}


void
MainWindow::on_actionMakeDelaunayMesh_triggered()
{
  double len = 0;
  int cc = 0;
  for(CDT::Finite_edges_iterator it = cdt.finite_edges_begin();
      it != cdt.finite_edges_end();
      ++it){
    if(cdt.is_constrained(*it)){
      ++cc;
      len+= sqrt(cdt.segment(*it).squared_length());
    }
  }
  double al = len/cc;
  al /= 2.0;
  int nv = cdt.number_of_vertices();
  CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, al));
  nv = cdt.number_of_vertices() - nv;
  statusBar()->showMessage(QString("Added %1 vertices").arg(nv), 2000);
  emit(changed());
}


void
MainWindow::on_actionInsertRandomPoints_triggered()
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
  cdt.insert(points.begin(), points.end());
  emit(changed());
}

#include "Constrained_Delaunay_triangulation_2.moc"

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Constrained_Delaunay_triangulation_2 demo");

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  Q_INIT_RESOURCE(File);
  Q_INIT_RESOURCE(Triangulation_2);
  Q_INIT_RESOURCE(Input);
  Q_INIT_RESOURCE(CGAL);

  MainWindow mainWindow;
  mainWindow.show();
  return app.exec();
}
