#include <fstream>
#include <boost/config.hpp>
#include <boost/version.hpp>

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
#include <CGAL/IO/WKT.h>
#endif

#include <CGAL/point_generators_2.h>
// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>

// GraphicsView items and event filters (input classes)

#include "RegularTriangulationRemoveVertex.h"
#include <CGAL/Qt/GraphicsViewCircleInput.h>
#include <CGAL/Qt/RegularTriangulationGraphicsItem.h>
#include <CGAL/Qt/PowerdiagramGraphicsItem.h>

// for viewportsBbox
#include <CGAL/Qt/utility.h>
// the two base classes
#include "ui_Regular_triangulation_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point_2;
typedef K::Weighted_point_2                                 Weighted_point_2;
typedef K::Point_2                                          Circle_2;
typedef K::Iso_rectangle_2                                  Iso_rectangle_2;

typedef CGAL::Regular_triangulation_2<K>                    Regular;

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Regular_triangulation_2
{
  Q_OBJECT

private:
  Regular dt;
  QGraphicsScene scene;

  CGAL::Qt::RegularTriangulationGraphicsItem<Regular> * dgi;
  CGAL::Qt::PowerdiagramGraphicsItem<Regular> * vgi;

  CGAL::Qt::RegularTriangulationRemoveVertex<Regular> * trv;
  CGAL::Qt::GraphicsViewCircleInput<K> * pi;
public:
  MainWindow();

public Q_SLOTS:

  void processInput(CGAL::Object o);

  void on_actionShowRegular_toggled(bool checked);

  void on_actionShowPowerdiagram_toggled(bool checked);

  void on_actionInsertPoint_toggled(bool checked);

  void on_actionInsertRandomPoints_triggered();

  void on_actionLoadPoints_triggered();

  void on_actionSavePoints_triggered();

  void on_actionClear_triggered();

  void on_actionRecenter_triggered();


Q_SIGNALS:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow()
{
  setupUi(this);

  // Add a GraphicItem for the regular triangulation
  dgi = new CGAL::Qt::RegularTriangulationGraphicsItem<Regular>(&dt);

  QObject::connect(this, SIGNAL(changed()),
                   dgi, SLOT(modelChanged()));

  dgi->setVerticesPen(QPen(Qt::red, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  dgi->setEdgesPen(QPen(Qt::black, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(dgi);

  // Add a GraphicItem for the Powerdiagram diagram
  vgi = new CGAL::Qt::PowerdiagramGraphicsItem<Regular>(&dt);

  QObject::connect(this, SIGNAL(changed()),
                   vgi, SLOT(modelChanged()));

  vgi->setEdgesPen(QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(vgi);
  vgi->hide();

  // Setup input handlers. They get events before the scene gets them
  // and the input they generate is passed to the triangulation with
  // the signal/slot mechanism
  pi = new CGAL::Qt::GraphicsViewCircleInput<K>(this, &scene, 1); // emits center/radius


  QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
                   this, SLOT(processInput(CGAL::Object)));

  trv = new CGAL::Qt::RegularTriangulationRemoveVertex<Regular>(&dt, this);
  QObject::connect(trv, SIGNAL(modelChanged()),
                   this, SIGNAL(changed()));

  //
  // Manual handling of actions
  //
  QObject::connect(this->actionQuit, SIGNAL(triggered()),
                   this, SLOT(close()));

  // We put mutually exclusive actions in an QActionGroup
  QActionGroup* ag = new QActionGroup(this);
  ag->addAction(this->actionInsertPoint);

  // Check two actions
  this->actionInsertPoint->setChecked(true);
  this->actionShowRegular->setChecked(true);

  //
  // Setup the scene and the view
  //
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  scene.setSceneRect(-100, -100, 100, 100);
  this->graphicsView->setScene(&scene);
  this->graphicsView->setMouseTracking(true);

  // Turn the vertical axis upside down
  this->graphicsView->transform().scale(1, -1);

  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_Regular_triangulation_2.html");
  this->addAboutCGAL();
}


void
MainWindow::processInput(CGAL::Object o)
{
  std::pair<Point_2, K::FT > center_sqr;
  if(CGAL::assign(center_sqr, o)){
    Regular::Point wp(center_sqr.first, center_sqr.second);
    dt.insert(wp);
  }

  Q_EMIT( changed());
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
    scene.installEventFilter(pi);
    scene.installEventFilter(trv);
  } else {
    scene.removeEventFilter(pi);
    scene.removeEventFilter(trv);
  }
}



void
MainWindow::on_actionShowRegular_toggled(bool checked)
{
  dgi->setVisibleEdges(checked);
}


void
MainWindow::on_actionShowPowerdiagram_toggled(bool checked)
{
  vgi->setVisible(checked);
}


void
MainWindow::on_actionClear_triggered()
{
  dt.clear();
  Q_EMIT( changed());
}


void
MainWindow::on_actionInsertRandomPoints_triggered()
{
  QRectF rect = CGAL::Qt::viewportsBbox(&scene);
  CGAL::Qt::Converter<K> convert;
  Iso_rectangle_2 isor = convert(rect);
  CGAL::Random_points_in_iso_rectangle_2<Point_2> pg((isor.min)(), (isor.max)());
  CGAL::Random rnd(CGAL::get_default_random());

  const int number_of_points =
    QInputDialog::getInt(this,
                             tr("Number of random points"),
                             tr("Enter number of random points"), 100, 0);

  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);
  std::vector<Weighted_point_2> points;
  points.reserve(number_of_points);
  for(int i = 0; i < number_of_points; ++i){
    Weighted_point_2 wp(*pg++, rnd.get_double(0, 500));
    points.push_back(wp);
  }
  dt.insert(points.begin(), points.end());
  // default cursor
  QApplication::setOverrideCursor(Qt::ArrowCursor);
  Q_EMIT( changed());
}


void
MainWindow::on_actionLoadPoints_triggered()
{
  QString fileName = QFileDialog::getOpenFileName(this,
                                                  tr("Open Points file"),
                                                  ".",
                                                  tr("Weighted Points (*.wpts.cgal);;"
                                                   #if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
                                                     "WKT files (*.wkt *.WKT);;"
                                                   #endif
                                                     "All (*)"));

  if(! fileName.isEmpty()){
    std::ifstream ifs(qPrintable(fileName));
    std::vector<Weighted_point_2> points;
    if(fileName.endsWith(".wkt",Qt::CaseInsensitive))
    {
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
      std::vector<K::Point_3> points_3;
      CGAL::read_multi_point_WKT(ifs, points_3);
      for(const K::Point_3& p : points_3)
      {
        points.push_back(Weighted_point_2(K::Point_2(p.x(), p.y()), p.z()));
      }
#endif
    }
    else
    {
      Weighted_point_2 p;
      while(ifs >> p) {
        points.push_back(p);
      }
    }
    dt.insert(points.begin(), points.end());

    actionRecenter->trigger();
    Q_EMIT( changed());
  }
}


void
MainWindow::on_actionSavePoints_triggered()
{
  QString fileName = QFileDialog::getSaveFileName(this,
                                                  tr("Save points"),
                                                  ".reg.cgal",
                                                  tr("Weighted Points (*.wpts.cgal);;"
                                                   #if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
                                                     "WKT files (*.wkt *.WKT);;"
                                                   #endif
                                                     "All (*)"));
  if(! fileName.isEmpty()){
    std::ofstream ofs(qPrintable(fileName));
    if(fileName.endsWith(".wkt",Qt::CaseInsensitive))
    {
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
      std::vector<K::Point_3> points_3;
      for(Regular::Finite_vertices_iterator
          vit = dt.finite_vertices_begin(),
          end = dt.finite_vertices_end();
          vit!= end; ++vit)
      {
        points_3.push_back(K::Point_3(vit->point().x(),
                                      vit->point().y(),
                                      vit->point().weight()));
      }
      CGAL::write_multi_point_WKT(ofs, points_3);
#endif
    }
    else
    {
      for(Regular::Finite_vertices_iterator
          vit = dt.finite_vertices_begin(),
          end = dt.finite_vertices_end();
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
  this->graphicsView->setSceneRect(dgi->boundingRect());
  this->graphicsView->fitInView(dgi->boundingRect(), Qt::KeepAspectRatio);
}


#include "Regular_triangulation_2.moc"
#include <CGAL/Qt/resources.h>

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Regular_triangulation_2 demo");

  // Import resources from libCGAL (Qt5).
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
