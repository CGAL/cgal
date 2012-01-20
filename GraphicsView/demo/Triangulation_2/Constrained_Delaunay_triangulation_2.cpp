//#define CGAL_USE_BOOST_BIMAP

#include <fstream>
#include <vector>

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Lipschitz_sizing_field_2.h>
#include <CGAL/Lipschitz_sizing_field_criteria_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Random.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Timer.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <QDragEnterEvent>
#include <QDropEvent>

// GraphicsView items and event filters (input classes)
#include "TriangulationCircumcircle.h"
#include <CGAL/Qt/GraphicsViewPolylineInput.h>
#include <CGAL/Qt/DelaunayMeshTriangulationGraphicsItem.h>
#include <CGAL/Qt/Converter.h>

// the two base classes
#include "ui_Constrained_Delaunay_triangulation_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

// for viewportsBbox(QGraphicsScene*)
#include <CGAL/Qt/utility.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef CGAL::Triangulation_vertex_base_2<K>  Vertex_base;
typedef CGAL::Constrained_triangulation_face_base_2<K> Face_base;

template <class Gt,
          class Fb >
class Enriched_face_base_2 : public Fb {
public:
  typedef Gt Geom_traits;
  typedef typename Fb::Vertex_handle Vertex_handle;
  typedef typename Fb::Face_handle Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other Fb2;
    typedef Enriched_face_base_2<Gt,Fb2> Other;
  };

protected:
  int status;

public:
  Enriched_face_base_2(): Fb(), status(-1) {};

  Enriched_face_base_2(Vertex_handle v0, 
		       Vertex_handle v1, 
		       Vertex_handle v2)
    : Fb(v0,v1,v2), status(-1) {};

  Enriched_face_base_2(Vertex_handle v0, 
		       Vertex_handle v1, 
		       Vertex_handle v2,
		       Face_handle n0, 
		       Face_handle n1, 
		       Face_handle n2)
    : Fb(v0,v1,v2,n0,n1,n2), status(-1) {};

  inline
  bool is_in_domain() const { return (status%2 == 1); };

  inline
  void set_in_domain(const bool b) { status = (b ? 1 : 0); };

  inline 
  void set_counter(int i) { status = i; };

  inline 
  int counter() const { return status; };

  inline 
  int& counter() { return status; };
}; // end class Enriched_face_base_2

typedef Enriched_face_base_2<K, Face_base> Fb;
typedef CGAL::Triangulation_data_structure_2<Vertex_base, Fb>  TDS;
typedef CGAL::Exact_predicates_tag              Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;

typedef CGAL::Lipschitz_sizing_field_2<K> Lipschitz_sizing_field;
typedef CGAL::Lipschitz_sizing_field_criteria_2<CDT, Lipschitz_sizing_field> Lipschitz_criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Lipschitz_criteria> Lipschitz_mesher;

typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Face_handle Face_handle;
typedef CDT::All_faces_iterator All_faces_iterator;


void
initializeID(const CDT& ct)
{
  for(All_faces_iterator it = ct.all_faces_begin(); it != ct.all_faces_end(); ++it){
    it->set_counter(-1);
  }
}


void 
discoverComponent(const CDT & ct, 
		  Face_handle start, 
		  int index, 
		  std::list<CDT::Edge>& border )
{
  if(start->counter() != -1){
    return;
  }
  std::list<Face_handle> queue;
  queue.push_back(start);

  while(! queue.empty()){
    Face_handle fh = queue.front();
    queue.pop_front();
    if(fh->counter() == -1){
      fh->counter() = index;
      fh->set_in_domain(index%2 == 1);
      for(int i = 0; i < 3; i++){
	CDT::Edge e(fh,i);
	Face_handle n = fh->neighbor(i);
	if(n->counter() == -1){
	  if(ct.is_constrained(e)){
	    border.push_back(e);
	  } else {
	    queue.push_back(n);
	  }
	}
	
      }
    }
  }
}

void 
discoverComponents(const CDT & ct)
{
  if (ct.dimension()!=2) return;
  int index = 0;
  std::list<CDT::Edge> border;
  discoverComponent(ct, ct.infinite_face(), index++, border);
  while(! border.empty()){
    CDT::Edge e = border.front();
    border.pop_front();
    Face_handle n = e.first->neighbor(e.second);
    if(n->counter() == -1){
      discoverComponent(ct, n, e.first->counter()+1, border);
    }
  }
} 



class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Constrained_Delaunay_triangulation_2
{
  Q_OBJECT
  
private:  
  CDT cdt; 
  QGraphicsScene scene;
  std::list<Point_2> seeds;

  CGAL::Qt::DelaunayMeshTriangulationGraphicsItem<CDT> * dgi;

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

public slots:
  void open(QString);

  void processInput(CGAL::Object o);

  void on_actionShowDelaunay_toggled(bool checked);

  void on_actionShow_constrained_edges_toggled(bool checked);

  void on_actionShow_faces_in_domain_toggled(bool checked);

  void on_actionInsertPolyline_toggled(bool checked);
  
  void on_actionCircumcenter_toggled(bool checked);

  void on_actionClear_triggered();

  void on_actionRecenter_triggered();

  void on_actionLoadConstraints_triggered();

  void loadFile(QString);

  void loadPolygonConstraints(QString);

  void loadEdgConstraints(QString);

  void on_actionSaveConstraints_triggered();

  void saveConstraints(QString);

  void on_actionMakeGabrielConform_triggered();

  void on_actionMakeDelaunayConform_triggered();

  void on_actionMakeDelaunayMesh_triggered();

  void on_actionMakeLipschitzDelaunayMesh_triggered();

  void on_actionInsertRandomPoints_triggered();

signals:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow()
{
  setupUi(this);

  this->graphicsView->setAcceptDrops(false);

  // Add a GraphicItem for the CDT triangulation
  dgi = new CGAL::Qt::DelaunayMeshTriangulationGraphicsItem<CDT>(&cdt);
  QColor facesColor(::Qt::blue);
  facesColor.setAlpha(150);
  dgi->setFacesInDomainBrush(facesColor);
    
  QObject::connect(this, SIGNAL(changed()),
		   dgi, SLOT(modelChanged()));

  dgi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  dgi->setZValue(-1);
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
  QObject::connect(this->actionQuit, SIGNAL(triggered()), 
		   this, SLOT(close()));

  // We put mutually exclusive actions in an QActionGroup
  QActionGroup* ag = new QActionGroup(this);
  ag->addAction(this->actionInsertPolyline);

  // Check two actions 
  this->actionInsertPolyline->setChecked(true);
  this->actionShowDelaunay->setChecked(true);
  this->actionShow_faces_in_domain->setChecked(true);
  this->actionShow_constrained_edges->setChecked(true);

  //
  // Setup the scene and the view
  //
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  scene.setSceneRect(-100, -100, 100, 100);
  this->graphicsView->setScene(&scene);
  this->graphicsView->setMouseTracking(true);

  // Turn the vertical axis upside down
  this->graphicsView->scale(1, -1);
                                                      
  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_Constrained_Delaunay_triangulation_2.html");
  this->addAboutCGAL();

  this->addRecentFiles(this->menuFile, this->actionQuit);
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
      /*
      std::cout.precision(12);
      std::cout << points.size() << std::endl;
      for( std::list<Point_2>::iterator it =  points.begin(); it != points.end(); ++it){
	std::cout << *it << std::endl;
      }
      */
      insert_polyline(points.begin(), points.end());
    }
  }


  initializeID(cdt);
  discoverComponents(cdt);
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
  dgi->setVisibleEdges(checked);
  update();
}

void
MainWindow::on_actionShow_constrained_edges_toggled(bool checked)
{
  dgi->setVisibleConstraints(checked);
  update();
}

void
MainWindow::on_actionShow_faces_in_domain_toggled(bool checked)
{
  dgi->setVisibleFacesInDomain(checked);
  update();
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
    if(fileName.endsWith(".cgal")){
      loadFile(fileName);
      this->addToRecentFiles(fileName);
    } else if(fileName.endsWith(".plg")){
      loadPolygonConstraints(fileName);
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
						  tr("Edge files (*.edg);;"
						     "Poly files (*.plg)"));
  open(fileName);
}

void
MainWindow::loadFile(QString fileName)
{
  std::ifstream ifs(qPrintable(fileName));
  ifs >> cdt;
  if(!ifs) abort();
  initializeID(cdt);
  discoverComponents(cdt);
  emit(changed());
  actionRecenter->trigger();
}

void
MainWindow::loadPolygonConstraints(QString fileName)
{
  K::Point_2 p,q, first;
  CDT::Vertex_handle vp, vq, vfirst;
  std::ifstream ifs(qPrintable(fileName));
  int n;
  // int counter = 0;
  while(ifs >> n){
    int poly_size = n;
    ifs >> first;
    p = first;
    vfirst = vp = cdt.insert(p);
    n--;
    while(n--){
      ifs >> q;
      vq = cdt.insert(q, vp->face());
      if(vp != vq) {
        cdt.insert_constraint(vp,vq);
        // std::cerr << "inserted constraint #" << counter++ << std::endl;
      }
      p = q;
      vp = vq;
    }
    if(poly_size != 2 && vp != vfirst) {
      cdt.insert_constraint(vp, vfirst);
    }
  }
  
  
  initializeID(cdt);
  discoverComponents(cdt);
  emit(changed());
  actionRecenter->trigger();
}


void
MainWindow::loadEdgConstraints(QString fileName)
{
  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);
  CGAL::Timer tim;
  tim.start();
  std::ifstream ifs(qPrintable(fileName));
  bool first=true;
  int n;
  ifs >> n;
  
  K::Point_2 p,q, qold;

  CDT::Vertex_handle vp, vq, vqold;
  while(ifs >> p) {
    ifs >> q;
    if(p == q){
      continue;
    }
    if((!first) && (p == qold)){
      vp = vqold;
    } else {
      vp = cdt.insert(p);
    }
    vq = cdt.insert(q, vp->face());
    if(vp != vq) {
      cdt.insert_constraint(vp,vq);
    }
    qold = q;
    vqold = vq;
    first = false;
  }


  tim.stop();
  statusBar()->showMessage(QString("Insertion took %1 seconds").arg(tim.time()), 2000);
  // default cursor
  QApplication::restoreOverrideCursor();
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
  std::ofstream output(qPrintable(fileName));
  if (output) output << cdt;
}


void
MainWindow::on_actionMakeGabrielConform_triggered()
{
  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);
  std::size_t nv = cdt.number_of_vertices();
  CGAL::make_conforming_Gabriel_2(cdt);
  nv = cdt.number_of_vertices() - nv;
  initializeID(cdt);
  discoverComponents(cdt);
  statusBar()->showMessage(QString("Added %1 vertices").arg(nv), 2000);
  // default cursor
  QApplication::restoreOverrideCursor();
  emit(changed());
}


void
MainWindow::on_actionMakeDelaunayConform_triggered()
{
  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);
  std::size_t nv = cdt.number_of_vertices();
  CGAL::make_conforming_Delaunay_2(cdt);
  initializeID(cdt);
  discoverComponents(cdt);
  nv = cdt.number_of_vertices() - nv;
  statusBar()->showMessage(QString("Added %1 vertices").arg(nv), 2000);
   // default cursor
  QApplication::restoreOverrideCursor();
  emit(changed());
}


void
MainWindow::on_actionMakeDelaunayMesh_triggered()
{
  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);
  double edge_length = 0;
  CGAL::Timer timer;
  timer.start();
  initializeID(cdt);
  discoverComponents(cdt);

  std::size_t nv = cdt.number_of_vertices();
  CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, edge_length), true);
  timer.stop();
  nv = cdt.number_of_vertices() - nv;
  initializeID(cdt);
  discoverComponents(cdt);
  statusBar()->showMessage(QString("Added %1 vertices in %2 seconds").arg(nv).arg(timer.time()), 2000);
  // default cursor
  QApplication::restoreOverrideCursor();
  emit(changed());

}

void
MainWindow::on_actionMakeLipschitzDelaunayMesh_triggered()
{
  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);
  std::set<Point_2> points;
  for(CDT::Finite_edges_iterator it = cdt.finite_edges_begin();
      it != cdt.finite_edges_end();
      ++it){
    if(cdt.is_constrained(*it)){
      Segment_2 s = cdt.segment(*it);
      points.insert(s.source());
      points.insert(s.target());
    }
  }

  initializeID(cdt);
  discoverComponents(cdt);

  std::size_t nv = cdt.number_of_vertices();
  Lipschitz_sizing_field field(points.begin(), points.end(), 0.7 ); // k-lipschitz with k=1
  Lipschitz_criteria criteria(0.125, &field);
  Lipschitz_mesher mesher(cdt);
  mesher.set_criteria(criteria);
  mesher.init(true);
  //  mesher.set_seeds(m_seeds.begin(),m_seeds.end(),false);
  mesher.refine_mesh();
  nv = cdt.number_of_vertices() - nv;
  statusBar()->showMessage(QString("Added %1 vertices").arg(nv), 2000);
  // default cursor
  QApplication::restoreOverrideCursor();
  emit(changed());
}



void
MainWindow::on_actionInsertRandomPoints_triggered()
{
  QRectF rect = CGAL::Qt::viewportsBbox(&scene);
  CGAL::Qt::Converter<K> convert;
  Iso_rectangle_2 isor = convert(rect);
  CGAL::Random_points_in_iso_rectangle_2<Point_2> pg((isor.min)(), (isor.max)());
  bool ok = false;
  const int number_of_points = 
    QInputDialog::getInteger(this, 
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
  std::vector<Point_2> points;
  points.reserve(number_of_points);
  for(int i = 0; i < number_of_points; ++i){
    points.push_back(*pg++);
  }
  cdt.insert(points.begin(), points.end());
  // default cursor
  QApplication::restoreOverrideCursor();
  emit(changed());
}

#include "Constrained_Delaunay_triangulation_2.moc"
#include <CGAL/Qt/resources.h>

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Constrained_Delaunay_triangulation_2 demo");

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  CGAL_QT4_INIT_RESOURCES;

  MainWindow mainWindow;
  mainWindow.show();

  QStringList args = app.arguments();
  args.removeAt(0);
  Q_FOREACH(QString filename, args) {
    mainWindow.open(filename);
  }

  return app.exec();
}
