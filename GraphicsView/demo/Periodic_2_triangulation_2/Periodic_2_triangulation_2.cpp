#include <fstream>

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_triangulation_traits_2.h>
#include <CGAL/point_generators_2.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>

// GraphicsView items and event filters (input classes)
#include "TriangulationPointInput.h"
#include "TriangulationMovingPoint.h"
#include "TriangulationRemoveVertex.h"
#include "PeriodicTriangulationLocate.h"
#include <CGAL/Qt/PeriodicTriangulationGraphicsItem.h>

// for viewportsBbox
#include <CGAL/Qt/utility.h>
  
// the two base classes
#include "ui_Periodic_2_triangulation_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

//typedef CGAL::Exact_predicates_inexact_constructions_kernel     EPIC;
struct EPIC : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Periodic_2_triangulation_traits_2<EPIC>           K;
typedef K::Point_2                                              Point_2;
typedef K::Iso_rectangle_2                                      Iso_rectangle_2;

#define NGHK_DELAUNAY
#ifdef NGHK_DELAUNAY
typedef CGAL::Periodic_2_Delaunay_triangulation_2<K>            Periodic_triangulation;
#else
typedef CGAL::Periodic_2_triangulation_2<K>                     Periodic_triangulation;
#endif

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Periodic_2_triangulation_2
{
  Q_OBJECT
  
private:  
  Periodic_triangulation triang; 
  QGraphicsScene scene;  

  CGAL::Qt::PeriodicTriangulationGraphicsItem<Periodic_triangulation> * pt_gi;

  CGAL::Qt::TriangulationMovingPoint<Periodic_triangulation> * pt_mp;
  CGAL::Qt::TriangulationRemoveVertex<Periodic_triangulation> * trv;
  CGAL::Qt::TriangulationPointInput<Periodic_triangulation> * pt_pi;
  CGAL::Qt::PeriodicTriangulationLocate<Periodic_triangulation> * pt_l;
public:
  MainWindow();

public slots:

  void processInput(CGAL::Object o);

  void on_actionMovingPoint_toggled(bool checked);

  void on_actionShowDelaunay_toggled(bool checked);

  void on_actionInsertPoint_toggled(bool checked);

  void on_actionShowConflictZone_toggled(bool checked);

  void on_actionInsertRandomPoints_triggered();

  void on_actionConvertTo9Cover_triggered();

  void on_actionConvertTo1Cover_triggered();

  void on_actionLoadPoints_triggered();

  void on_actionSavePoints_triggered();

  void on_actionClear_triggered();

  void on_actionRecenter_triggered();

  void open(const QString& fileName);

signals:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow()
{
  setupUi(this);

  // Add a GraphicItem for the Periodic triangulation
  pt_gi = new CGAL::Qt::PeriodicTriangulationGraphicsItem<Periodic_triangulation>(&triang);

  QObject::connect(this, SIGNAL(changed()),
		   pt_gi, SLOT(modelChanged()));

  pt_gi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(pt_gi);

  // Setup input handlers. They get events before the scene gets them
  // and the input they generate is passed to the triangulation with 
  // the signal/slot mechanism    
  pt_pi = new CGAL::Qt::TriangulationPointInput<Periodic_triangulation>(&scene, &triang, this );
  
  QObject::connect(pt_pi, SIGNAL(generate(CGAL::Object)),
		   this, SLOT(processInput(CGAL::Object)));
  
  pt_mp = new CGAL::Qt::TriangulationMovingPoint<Periodic_triangulation>(&triang, this);
  // TriangulationMovingPoint<Periodic_triangulation> emits a modelChanged() signal each
  // time the moving point moves.
  // The following connection is for the purpose of emitting changed().
  QObject::connect(pt_mp, SIGNAL(modelChanged()),
		   this, SIGNAL(changed()));

  trv = new CGAL::Qt::TriangulationRemoveVertex<Periodic_triangulation>(&triang, this);
  QObject::connect(trv, SIGNAL(modelChanged()),
		   this, SIGNAL(changed()));

  pt_l = new CGAL::Qt::PeriodicTriangulationLocate<Periodic_triangulation>(&scene, &triang, this);
  
  // 
  // Manual handling of actions
  //

  QObject::connect(this->actionQuit, SIGNAL(triggered()), 
		   this, SLOT(close()));

  // We put mutually exclusive actions in an QActionGroup
  QActionGroup* ag = new QActionGroup(this);
  ag->addAction(this->actionInsertPoint);
  ag->addAction(this->actionMovingPoint);
  ag->addAction(this->actionShowConflictZone);

  // Check two actions 
  this->actionInsertPoint->setChecked(true);
  this->actionShowDelaunay->setChecked(true);

  //
  // Setup the scene and the view
  //
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  scene.setSceneRect(0, 0, 1, 1);
  this->graphicsView->setScene(&scene);
  this->graphicsView->setMouseTracking(true);

  // Turn the vertical axis upside down
  this->graphicsView->matrix().scale(1, -1);
                                                      
  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_Periodic_2_triangulation_2.html");

  this->addRecentFiles(this->menuFile, this->actionQuit);
  connect(this, SIGNAL(openRecentFile(QString)),
	  this, SLOT(open(QString)));
  
  on_actionRecenter_triggered();
}


void
MainWindow::processInput(CGAL::Object o)
{
  bool was_empty = triang.empty();
  Point_2 p;
  if(CGAL::assign(p, o)) {
    double dx = triang.domain().xmax() - triang.domain().xmin();
    double dy = triang.domain().ymax() - triang.domain().ymin();
    p = Point_2(p.x()- std::floor(p.x()/dx),
                p.y()- std::floor(p.y()/dy));
    triang.insert(p);
  }
  emit(changed());

  if (was_empty)
    on_actionRecenter_triggered();
}


/* 
 *  Qt Automatic Connections
 *  http://doc.trolltech.com/4.4/designer-using-a-copt_mponent.html#automatic-connections
 * 
 *  setupUi(this) generates connections to the slots named
 *  "on_<action_name>_<signal_name>"
 */
void
MainWindow::on_actionInsertPoint_toggled(bool checked)
{
  if(checked){
    scene.installEventFilter(pt_pi);
    scene.installEventFilter(trv);
  } else {
    scene.removeEventFilter(pt_pi);
    scene.removeEventFilter(trv);
  }
}

void
MainWindow::on_actionShowConflictZone_toggled(bool checked)
{
  if(checked) {
    scene.installEventFilter(pt_l);
  } else {
    scene.removeEventFilter(pt_l);
  }
}


void
MainWindow::on_actionMovingPoint_toggled(bool checked)
{

  if(checked){
    scene.installEventFilter(pt_mp);
  } else {
    scene.removeEventFilter(pt_mp);
  }
}

void
MainWindow::on_actionShowDelaunay_toggled(bool checked)
{
  pt_gi->setVisibleEdges(checked);
}


void
MainWindow::on_actionClear_triggered()
{
  triang.clear();
  emit(changed());
}


void
MainWindow::on_actionInsertRandomPoints_triggered()
{
  CGAL::Random_points_in_iso_rectangle_2<Point_2> pg(triang.domain().min(), triang.domain().max());
  bool ok = false;
  const int number_of_points = 
    QInputDialog::getInteger(this, 
                             tr("Number of random points"),
                             tr("Enter number of random points"),
			     250,
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
  triang.insert(points.begin(), points.end(), true);

  // default cursor
  QApplication::restoreOverrideCursor();

  on_actionRecenter_triggered();
  emit(changed());
}


void
MainWindow::on_actionConvertTo9Cover_triggered() {
  if (triang.is_1_cover()) {
    triang.convert_to_9_sheeted_covering();
    emit(changed());
  }
}

void
MainWindow::on_actionConvertTo1Cover_triggered() {
  if (!triang.is_1_cover()) {
    triang.convert_to_1_sheeted_covering();
    emit(changed());
  }
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
MainWindow::open(const QString& fileName)
{
  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);
  std::ifstream ifs(qPrintable(fileName));
  
  K::Point_2 p;
  std::vector<K::Point_2> points;
  while(ifs >> p) {
    points.push_back(p);
  }
  triang.insert(points.begin(), points.end());

  // default cursor
  QApplication::restoreOverrideCursor();
  this->addToRecentFiles(fileName);
  actionRecenter->trigger();
  emit(changed());
    
}

void
MainWindow::on_actionSavePoints_triggered()
{
  QString fileName = QFileDialog::getSaveFileName(this,
						  tr("Save points"),
						  ".");
  if(! fileName.isEmpty()){
    std::ofstream ofs(qPrintable(fileName));
    for(Periodic_triangulation::Vertex_iterator 
          vit = triang.vertices_begin(),
          end = triang.vertices_end();
        vit!= end; ++vit)
    {
      ofs << vit->point() << std::endl;
    }
  }
}


void
MainWindow::on_actionRecenter_triggered()
{
  this->graphicsView->setSceneRect(pt_gi->boundingRect());
  this->graphicsView->fitInView(pt_gi->boundingRect(), Qt::KeepAspectRatio);  
}


#include "Periodic_2_triangulation_2.moc"

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("www.nghk.nl");
  app.setOrganizationName("Nico Kruithof");
  app.setApplicationName("Periodic_2_triangulation_2 demo");

  // Ipt_mport resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  Q_INIT_RESOURCE(File);
  Q_INIT_RESOURCE(Periodic_2_triangulation_2);
  Q_INIT_RESOURCE(Input);
  Q_INIT_RESOURCE(CGAL);

  MainWindow mainWindow;
  mainWindow.show();
  try {
    return app.exec();
  }
  catch (char const *str) {
    std::cerr << "EXCEPTION: " << str << std::endl;
    return -1;
  }
  catch (...) {
    std::cerr << "Unknown exception!" << std::endl;
    return -1;
  }
}
