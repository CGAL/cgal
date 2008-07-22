#include <fstream>

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>

#include <CGAL/point_generators_2.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>

// GraphicsView items and event filters (input classes)
#include "TriangulationCircumcircle.h"
#include "TriangulationMovingPoint.h"
#include <CGAL/Qt/GraphicsViewPolylineInput.h>
#include <CGAL/Qt/ConstrainedTriangulationGraphicsItem.h>
  
// the two base classes
#include "ui_Constrained_Delaunay_triangulation_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Triangulation_vertex_base_2<K>  Vertex_base;
typedef CGAL::Constrained_triangulation_face_base_2<K> Face_base;
typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base>  TDS;
typedef CGAL::Exact_predicates_tag              Itag;


typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Constrained_Delaunay_triangulation_2
{
  Q_OBJECT
  
private:  
  CDT cdt; 
  QGraphicsScene scene;  

  CGAL::Qt::ConstrainedTriangulationGraphicsItem<CDT> * dgi;

  CGAL::Qt::TriangulationMovingPoint<CDT> * mp;
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

  void processInput(CGAL::Object o);

  void on_actionMovingPoint_toggled(bool checked);

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
  pi = new CGAL::Qt::GraphicsViewPolylineInput<K>(this, &scene, 0, false); // inputs polylines which are not closed

  QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
		   this, SLOT(processInput(CGAL::Object)));
    
  mp = new CGAL::Qt::TriangulationMovingPoint<CDT>(&cdt, this);
  // TriangulationMovingPoint<CDT> generates an empty Object() each
  // time the moving point moves.
  // The following connection is for the purpose of emitting changed().
  QObject::connect(mp, SIGNAL(generate(CGAL::Object)),
		   this, SIGNAL(changed()));

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
  ag->addAction(this->actionMovingPoint);

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
MainWindow::on_actionLoadConstraints_triggered()
{
  QString fileName = QFileDialog::getOpenFileName(this,
						  tr("Open Constraint File"),
						  ".",
						  tr("Edge files (*.edg)\n"
						     "Poly files (*.poly)"));
  if(! fileName.isEmpty()){
    if(fileName.endsWith(".poly")){
      loadPolyConstraints(fileName);
    } else if(fileName.endsWith(".edg")){
      loadEdgConstraints(fileName);
    }
  }
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

  actionRecenter->trigger();
  emit(changed());
}


void
MainWindow::loadEdgConstraints(QString fileName)
{
  std::ifstream ifs(qPrintable(fileName));
  bool first=true;
  int n;
  ifs >> n;
  
  K::Point_2 p,q, qold;
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
