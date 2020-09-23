#include <fstream>

#include <boost/config.hpp>
#include <boost/version.hpp>
// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/point_generators_2.h>
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
#include <CGAL/IO/WKT.h>
#endif

// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>

// GraphicsView items and event filters (input classes)
#include "TriangulationCircumcircle.h"
#include "TriangulationMovingPoint.h"
#include "TriangulationConflictZone.h"
#include "TriangulationRemoveVertex.h"
#include "TriangulationPointInputAndConflictZone.h"
#include <CGAL/Qt/PeriodicTriangulationGraphicsItem.h>
#include <CGAL/Qt/PeriodicVoronoiGraphicsItem.h>

// for viewportsBbox
#include <CGAL/Qt/utility.h>

// the two base classes
#include "ui_Periodic_2_triangulation_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

//typedef CGAL::Exact_predicates_inexact_constructions_kernel     EPIC;
struct EPIC : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<EPIC>  Gt;
typedef Gt::Point_2                                             Point_2;
typedef Gt::Iso_rectangle_2                                     Iso_rectangle_2;

typedef CGAL::Periodic_2_Delaunay_triangulation_2<Gt>           Periodic_DT;

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Periodic_2_triangulation_2
{
  Q_OBJECT

private:
  Periodic_DT triang;
  QGraphicsScene scene;

  typedef CGAL::Qt::PeriodicTriangulationGraphicsItem<Periodic_DT> PTGI;

  CGAL::Qt::PeriodicTriangulationGraphicsItem<Periodic_DT> * pt_gi;
  CGAL::Qt::PeriodicTriangulationVoronoiGraphicsItem<Periodic_DT> * vgi;

  CGAL::Qt::TriangulationMovingPoint<Periodic_DT> * pt_mp;
  CGAL::Qt::TriangulationConflictZone<Periodic_DT> * pt_cz;
  CGAL::Qt::TriangulationRemoveVertex<Periodic_DT> * pt_rv;
  CGAL::Qt::TriangulationPointInputAndConflictZone<Periodic_DT> * pt_pi;
  CGAL::Qt::TriangulationCircumcircle<Periodic_DT> *pt_cc;
public:
  MainWindow();

public Q_SLOTS:

  void processInput(CGAL::Object o);

  void on_actionClear_triggered();
  void on_actionInsertPoint_toggled(bool checked);
  void on_actionInsertRandomPoints_triggered();
  void on_actionConvertTo9Cover_triggered();
  void on_actionConvertTo1Cover_triggered();

  void on_actionMovingPoint_toggled(bool checked);
  void on_actionShowConflictZone_toggled(bool checked);
  void on_actionCircumcenter_toggled(bool checked);

  void on_actionShowDelaunay_toggled(bool checked);
  void on_actionShowVoronoi_toggled(bool checked);
  void on_actionNoneSimplicesEmphasized_triggered(bool checked);
  void on_actionUniqueSimplicesEmphasized_triggered(bool checked);
  void on_actionStoredSimplicesEmphasized_triggered(bool checked);
  void on_actionUniqueCoverDomainSimplicesEmphasized_triggered(bool checked);
  void on_actionStoredCoverDomainSimplicesEmphasized_triggered(bool checked);


  void on_actionLoadPoints_triggered();
  void on_actionSavePoints_triggered();

  void on_actionRecenter_triggered();

  virtual void open(QString fileName);

Q_SIGNALS:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow()
{
  setupUi(this);

  this->graphicsView->setAcceptDrops(false);

  // Add a GraphicItem for the Periodic triangulation
  pt_gi = new CGAL::Qt::PeriodicTriangulationGraphicsItem<Periodic_DT>(&triang);

  QObject::connect(this, SIGNAL(changed()),
                   pt_gi, SLOT(modelChanged()));

  pt_gi->setVerticesPen(QPen(Qt::red, 5, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(pt_gi);

  // Add a GraphicItem for the Voronoi diagram
  vgi = new CGAL::Qt::PeriodicTriangulationVoronoiGraphicsItem<Periodic_DT>(&triang);

  QObject::connect(this, SIGNAL(changed()),
                   vgi, SLOT(modelChanged()));

  vgi->setEdgesPen(QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(vgi);
  vgi->hide();

  // Setup input handlers. They get events before the scene gets them
  // and the input they generate is passed to the triangulation with
  // the signal/slot mechanism
  pt_pi = new CGAL::Qt::TriangulationPointInputAndConflictZone<Periodic_DT>(&scene, &triang, this );
  QObject::connect(pt_pi, SIGNAL(generate(CGAL::Object)),
                   this, SLOT(processInput(CGAL::Object)));

  pt_mp = new CGAL::Qt::TriangulationMovingPoint<Periodic_DT>(&triang, this);
  // TriangulationMovingPoint<Periodic_DT> emits a modelChanged() signal each
  // time the moving point moves.
  // The following connection is for the purpose of emitting changed().
  QObject::connect(pt_mp, SIGNAL(modelChanged()),
                   this, SIGNAL(changed()));

  pt_cz = new CGAL::Qt::TriangulationConflictZone<Periodic_DT>(&scene, &triang, this);
  QObject::connect(pt_cz, SIGNAL(modelChanged()),
                   this, SIGNAL(changed()));

  pt_rv = new CGAL::Qt::TriangulationRemoveVertex<Periodic_DT>(&triang, this);
  QObject::connect(pt_rv, SIGNAL(modelChanged()),
                   this, SIGNAL(changed()));

  pt_cc = new CGAL::Qt::TriangulationCircumcircle<Periodic_DT>(&scene, &triang, this);
  pt_cc ->setPen(QPen(::Qt::black, .01));
  QObject::connect(pt_cc, SIGNAL(modelChanged()),
                   this, SIGNAL(changed()));

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
  scene.setSceneRect(0, 0, 1, 1);
  this->graphicsView->setScene(&scene);
  this->graphicsView->setMouseTracking(true);

  // Turn the vertical axis upside down
  this->graphicsView->transform().scale(1, -1);

  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_Periodic_2_triangulation_2.html");
  this->addAboutCGAL();

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
  Q_EMIT( changed());

  if (was_empty)
    on_actionRecenter_triggered();
}


/*
 *  Qt Automatic Connections
 *  https://doc.qt.io/qt-5/designer-using-a-ui-file.html#automatic-connections
 *
 *  setupUi(this) generates connections to the slots named
 *  "on_<action_name>_<signal_name>"
 */
void
MainWindow::on_actionInsertPoint_toggled(bool checked)
{
  if(checked){
    scene.installEventFilter(pt_pi);
    scene.installEventFilter(pt_rv);
  } else {
    scene.removeEventFilter(pt_pi);
    scene.removeEventFilter(pt_rv);
  }
}

void
MainWindow::on_actionShowConflictZone_toggled(bool checked)
{
  if(checked) {
    scene.installEventFilter(pt_cz);
  } else {
    scene.removeEventFilter(pt_cz);
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
MainWindow::on_actionCircumcenter_toggled(bool checked)
{
  if(checked){
    scene.installEventFilter(pt_cc);
    pt_cc->show();
  } else {
    scene.removeEventFilter(pt_cc);
    pt_cc->hide();
  }
}


void
MainWindow::on_actionShowDelaunay_toggled(bool checked)
{
  pt_gi->setVisibleEdges(checked);
}

void
MainWindow::on_actionShowVoronoi_toggled(bool checked)
{
  vgi->setVisible(checked);
}

void
MainWindow::on_actionClear_triggered()
{
  triang.clear();
  Q_EMIT( changed());
}


void
MainWindow::on_actionInsertRandomPoints_triggered()
{
  CGAL::Random_points_in_iso_rectangle_2<Point_2> pg((triang.domain().min)(),
                                                     (triang.domain().max)());
  bool ok = false;

  const int number_of_points =
    QInputDialog::getInt(this,
                             tr("Number of random points"),
                             tr("Enter number of random points"),
                             250,
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
  triang.insert(points.begin(), points.end(), true);

  // default cursor
  QApplication::restoreOverrideCursor();

  on_actionRecenter_triggered();
  Q_EMIT( changed());
}


void
MainWindow::on_actionConvertTo9Cover_triggered() {
  if (triang.is_1_cover()) {
    triang.convert_to_9_sheeted_covering();
    Q_EMIT( changed());
  }
}

void
MainWindow::on_actionConvertTo1Cover_triggered() {
  if (!triang.is_1_cover()) {
    triang.convert_to_1_sheeted_covering();
    Q_EMIT( changed());
  }
}

void
MainWindow::on_actionLoadPoints_triggered()
{
  QString fileName = QFileDialog::getOpenFileName(this,
                                                  tr("Open Points file"),
                                                  ".",
                                                  tr("CGAL files (*.pts.cgal);;"
                                                   #if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
                                                     "WKT files (*.wkt *.WKT);;"
                                                   #endif
                                                     "All files (*)"));
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

  Point_2 p;
  std::vector<Point_2> points;
  if(fileName.endsWith(".wkt", Qt::CaseInsensitive))
  {
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
    CGAL::read_multi_point_WKT(ifs, points);
#endif
  }
  else
  {
    while(ifs >> p) {
      points.push_back(p);
    }
  }
  triang.clear();
  triang.insert(points.begin(), points.end());

  // default cursor
  QApplication::restoreOverrideCursor();
  this->addToRecentFiles(fileName);
  actionRecenter->trigger();
  Q_EMIT( changed());

}

void
MainWindow::on_actionSavePoints_triggered()
{
  QString fileName = QFileDialog::getSaveFileName(this,
                                                  tr("Save points"),
                                                  ".",
                                                  tr("CGAL files (*.pts.cgal);;"
                                                   #if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
                                                     "WKT files (*.wkt *.WKT);;"
                                                   #endif
                                                     "All files (*)"));
  if(! fileName.isEmpty()){
    std::ofstream ofs(qPrintable(fileName));
    if(fileName.endsWith(".wkt", Qt::CaseInsensitive))
    {
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
      std::vector<Point_2> points;
      points.reserve(std::distance(triang.unique_vertices_begin(),
                                   triang.unique_vertices_end()));
      for(Periodic_DT::Unique_vertex_iterator
          vit = triang.unique_vertices_begin(),
          end = triang.unique_vertices_end();
          vit!= end; ++vit)
      {
        points.push_back(vit->point());
      }
      CGAL::write_multi_point_WKT(ofs, points);
#endif
    }
    else
    {
      for(Periodic_DT::Unique_vertex_iterator
          vit = triang.unique_vertices_begin(),
          end = triang.unique_vertices_end();
          vit!= end; ++vit)
      {
        ofs << vit->point() << std::endl;
      }
    }
  }
}


void
MainWindow::on_actionRecenter_triggered()
{
  pt_gi->modelChanged();
  this->graphicsView->setSceneRect(pt_gi->boundingRect());
  this->graphicsView->fitInView(pt_gi->boundingRect(), Qt::KeepAspectRatio);
}

void MainWindow::on_actionNoneSimplicesEmphasized_triggered(bool)
{
  pt_gi->setEmphasizedSimplices(PTGI::NONE);
  Q_EMIT changed();
}
void MainWindow::on_actionUniqueSimplicesEmphasized_triggered(bool)
{
  pt_gi->setEmphasizedSimplices(PTGI::UNIQUE);
  Q_EMIT changed();
}
void MainWindow::on_actionStoredSimplicesEmphasized_triggered(bool)
{
  pt_gi->setEmphasizedSimplices(PTGI::STORED);
  Q_EMIT changed();
}
void MainWindow::on_actionUniqueCoverDomainSimplicesEmphasized_triggered(bool)
{
  pt_gi->setEmphasizedSimplices(PTGI::UNIQUE_COVER_DOMAIN);
  Q_EMIT changed();
}
void MainWindow::on_actionStoredCoverDomainSimplicesEmphasized_triggered(bool)
{
  pt_gi->setEmphasizedSimplices(PTGI::STORED_COVER_DOMAIN);
  Q_EMIT changed();
}


#include "Periodic_2_Delaunay_triangulation_2.moc"
#include <CGAL/Qt/resources.h>

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("www.nghk.nl");
  app.setOrganizationName("Nico Kruithof");
  app.setApplicationName("Periodic_2_Delaunay_triangulation_2 demo");

  // Import resources from libCGAL (Qt5).
  // See https://doc.qt.io/qt-5/qdir.html#Q_INIT_RESOURCE
  CGAL_QT_INIT_RESOURCES;

  MainWindow mainWindow;
  mainWindow.show();
  QStringList args = app.arguments();
  args.removeAt(0);
  Q_FOREACH(QString filename, args) {
    mainWindow.open(filename);
  }
  return app.exec();
}
