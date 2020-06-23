#include <fstream>
#include <cmath>
#include <boost/config.hpp>
#include <boost/version.hpp>

// CGAL headers
#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_circle_2_traits_2.h>
#include <CGAL/Min_ellipse_2.h>
#include <CGAL/Min_ellipse_2_traits_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/min_quadrilateral_2.h>
#include <CGAL/rectangular_p_center_2.h>
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
#include <CGAL/IO/WKT.h>
#endif

// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <QGraphicsEllipseItem>
#include <QGraphicsRectItem>

// GraphicsView items and event filters (input classes)

#include <CGAL/Qt/PointsGraphicsItem.h>
#include <CGAL/Qt/PolygonGraphicsItem.h>

// for viewportsBbox
#include <CGAL/Qt/utility.h>

#include <CGAL/Qt/GraphicsViewPolylineInput.h>

// the two base classes
#include "ui_Bounding_volumes.h"
#include <CGAL/Qt/DemosMainWindow.h>

#include "Ellipse.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Vector_2 Vector_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;

typedef CGAL::Polygon_2<K> Polygon_2;

typedef CGAL::Min_circle_2<CGAL::Min_circle_2_traits_2<K> > Min_circle;
typedef CGAL::Min_ellipse_2<CGAL::Min_ellipse_2_traits_2<K> > Min_ellipse;

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Bounding_volumes
{
  Q_OBJECT

private:
  Polygon_2 convex_hull, min_rectangle, min_parallelogram;
  Min_circle mc;
  Min_ellipse me;
  QGraphicsScene scene;

  std::vector<Point_2> points;
  CGAL::Qt::PointsGraphicsItem<std::vector<Point_2> > * pgi;
  CGAL::Qt::PolygonGraphicsItem<Polygon_2> * convex_hull_gi;
  CGAL::Qt::PolygonGraphicsItem<Polygon_2> * min_rectangle_gi;
  CGAL::Qt::PolygonGraphicsItem<Polygon_2> * min_parallelogram_gi;
  QGraphicsEllipseItem *cgi, *egi;

  const std::size_t P;
  QGraphicsRectItem *p_center[3];
  Iso_rectangle_2 p_center_iso_rectangle[3];
  CGAL::Qt::GraphicsViewPolylineInput<K> * pi;

public:
  MainWindow();

public Q_SLOTS:

  void update();

  void update_from_points();

  void processInput(CGAL::Object o);

  void on_actionShowMinCircle_toggled(bool checked);

  void on_actionShowMinEllipse_toggled(bool checked);

  void on_actionShowMinRectangle_toggled(bool checked);

  void on_actionShowMinParallelogram_toggled(bool checked);

  void on_actionShowConvexHull_toggled(bool checked);

  void on_actionShowPCenter_toggled(bool checked);

  void on_actionInsertPoint_toggled(bool checked);

  void on_actionInsertRandomPoints_triggered();

  void on_actionLoadPoints_triggered();

  void on_actionSavePoints_triggered();

  void on_actionClear_triggered();

  void on_actionRecenter_triggered();

  virtual void open(QString fileName);

Q_SIGNALS:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow(), P(3)
{
  setupUi(this);

  QObject::connect(this, SIGNAL(changed()), this, SLOT(update()));

  // Add a GraphicItem for the Min_circle
  cgi = new QGraphicsEllipseItem;
  cgi->setPen(QPen(Qt::red, 0, Qt::SolidLine));
  cgi->hide();
  scene.addItem(cgi);

  egi = new QGraphicsEllipseItem;
  egi->setPen(QPen(Qt::magenta, 0, Qt::SolidLine));
  egi->hide();
  scene.addItem(egi);

  for(std::size_t i =0; i < P; i++){
    p_center[i] = new QGraphicsRectItem;
    p_center[i]->setPen(QPen(Qt::cyan, 0, Qt::SolidLine));
    p_center[i]->hide();
    scene.addItem(p_center[i]);
  }

  // Graphics Item for the input point set
  pgi = new CGAL::Qt::PointsGraphicsItem<std::vector<Point_2> >(&points);

  QObject::connect(this, SIGNAL(changed()),
                   pgi, SLOT(modelChanged()));
  pgi->setVerticesPen(QPen(Qt::black, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(pgi);


  // Graphics Item for the convex hull
  convex_hull_gi = new CGAL::Qt::PolygonGraphicsItem<Polygon_2>(&convex_hull);

  QObject::connect(this, SIGNAL(changed()),
                   convex_hull_gi, SLOT(modelChanged()));
  convex_hull_gi->setEdgesPen(QPen(Qt::black, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(convex_hull_gi);


  // Graphics Item for the min rectangle
  min_rectangle_gi = new CGAL::Qt::PolygonGraphicsItem<Polygon_2>(&min_rectangle);

  QObject::connect(this, SIGNAL(changed()),
                   min_rectangle_gi, SLOT(modelChanged()));
  min_rectangle_gi->setEdgesPen(QPen(Qt::green, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(min_rectangle_gi);


  // Graphics Item for the min parallelogram
  min_parallelogram_gi = new CGAL::Qt::PolygonGraphicsItem<Polygon_2>(&min_parallelogram);

  QObject::connect(this, SIGNAL(changed()),
                   min_parallelogram_gi, SLOT(modelChanged()));
  min_parallelogram_gi->setEdgesPen(QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(min_parallelogram_gi);


  // Setup input handlers. They get events before the scene gets them
  // and the input they generate is passed to the triangulation with
  // the signal/slot mechanism
  pi = new CGAL::Qt::GraphicsViewPolylineInput<K>(this, &scene, 1);

  scene.installEventFilter(pi);

  QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
                   this, SLOT(processInput(CGAL::Object)));


  //
  // Manual handling of actions
  //

  QObject::connect(this->actionQuit, SIGNAL(triggered()),
                   this, SLOT(close()));

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
  this->addAboutDemo(":/cgal/help/about_Bounding_volumes.html");
  this->addAboutCGAL();

  this->addRecentFiles(this->menuFile, this->actionQuit);
  connect(this, SIGNAL(openRecentFile(QString)),
          this, SLOT(open(QString)));
}

void
MainWindow::update()
{
  if(this->actionShowConvexHull->isChecked()){
    convex_hull_gi->show();
  }else {
    convex_hull_gi->hide();
  }

  if(this->actionShowMinRectangle->isChecked()){
    min_rectangle_gi->show();
  }else {
    min_rectangle_gi->hide();
  }


  if(this->actionShowMinParallelogram->isChecked()){
    min_parallelogram_gi->show();
  }else {
    min_parallelogram_gi->hide();
  }

  CGAL::Qt::Converter<K> convert;

  if(this->actionShowPCenter->isChecked() && convex_hull.size()>=3){
    for(std::size_t i=0; i< P; i++){
      p_center[i]->setRect(convert(p_center_iso_rectangle[i]));
      p_center[i]->show();
    }
  }

  if (mc.is_degenerate() || (! this->actionShowMinCircle->isChecked())){
    cgi->hide();
  } else {
    K::Circle_2 c;
    if (mc.number_of_support_points() == 2)
      c = K::Circle_2(mc.support_point(0), mc.support_point(1));
    else
      c = K::Circle_2(mc.support_point(0), mc.support_point(1), mc.support_point(2));


    cgi->setRect(convert(c.bbox()));
    cgi->show();
  }

  if (me.is_degenerate()  || (! this->actionShowMinEllipse->isChecked()) ){
    egi->hide();
  } else {
    if (me.number_of_support_points() == 2) {
    } else {
      Ellipse_2<K> e(me);
      double half_width = sqrt(e.va() * e.va());
      double half_height = sqrt(e.vb() * e.vb());
      double angle = std::atan2( e.va().y(), e.va().x() ) * 180.0/CGAL_PI;
      Vector_2 wh(half_width, half_height);

      Iso_rectangle_2 isor(e.center()+ wh, e.center()-wh);
      egi->setRect(convert(isor));
      // Rotate an item 45 degrees around (x, y).
      double x = e.center().x();
      double y = e.center().y();
      egi->setTransform(QTransform().translate(x, y).rotate(angle).translate(-x, -y));
      egi->show();
    }
  }
}


void
MainWindow::update_from_points()
{
    convex_hull.clear();
    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(convex_hull));

    min_rectangle.clear();
    CGAL::min_rectangle_2(convex_hull.vertices_begin(), convex_hull.vertices_end(), std::back_inserter(min_rectangle));

    min_parallelogram.clear();
    CGAL::min_parallelogram_2(convex_hull.vertices_begin(), convex_hull.vertices_end(), std::back_inserter(min_parallelogram));

    std::vector<Point_2> center;
    double radius;

    CGAL::rectangular_p_center_2 (points.begin(), points.end(), std::back_inserter(center), radius, static_cast<int>(P));
    Vector_2 rvec(radius, radius);

    for(std::size_t i = 0; i < center.size(); i++){
      p_center_iso_rectangle[i] = Iso_rectangle_2(center[i]-rvec, center[i]+rvec);
    }
}


void
MainWindow::processInput(CGAL::Object o)
{
  std::list<Point_2> input;
  if(CGAL::assign(input, o)){
    Point_2 p = input.front();

    mc.insert(p);
    me.insert(p);
    points.push_back(p);

    convex_hull.push_back(p);
    Polygon_2 tmp;
    CGAL::convex_hull_2(convex_hull.vertices_begin(), convex_hull.vertices_end(), std::back_inserter(tmp));
    convex_hull = tmp;

    min_rectangle.clear();
    CGAL::min_rectangle_2(convex_hull.vertices_begin(), convex_hull.vertices_end(), std::back_inserter(min_rectangle));

    min_parallelogram.clear();
    CGAL::min_parallelogram_2(convex_hull.vertices_begin(), convex_hull.vertices_end(), std::back_inserter(min_parallelogram));

    std::vector<Point_2> center;
    double radius;
    if (points.size()>=P){
      CGAL::rectangular_p_center_2 (points.begin(), points.end(), std::back_inserter(center), radius, static_cast<int>(P));
      Vector_2 rvec(radius, radius);

      for(std::size_t i=0; i < center.size(); i++){
        p_center_iso_rectangle[i] = Iso_rectangle_2(center[i]-rvec, center[i]+rvec);
      }
    }
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
  } else {
    scene.removeEventFilter(pi);
  }
}


void
MainWindow::on_actionShowMinCircle_toggled(bool checked)
{
  cgi->setVisible(checked);
  Q_EMIT( changed());
}

void
MainWindow::on_actionShowMinEllipse_toggled(bool checked)
{
  egi->setVisible(checked);
  Q_EMIT( changed());
}


void
MainWindow::on_actionShowMinRectangle_toggled(bool checked)
{
  min_rectangle_gi->setVisible(checked);
  Q_EMIT( changed());
}

void
MainWindow::on_actionShowMinParallelogram_toggled(bool checked)
{
  min_parallelogram_gi->setVisible(checked);
  Q_EMIT( changed());
}

void
MainWindow::on_actionShowConvexHull_toggled(bool checked)
{
  convex_hull_gi->setVisible(checked);
  Q_EMIT( changed());
}

void
MainWindow::on_actionShowPCenter_toggled(bool checked)
{
  for(std::size_t i =0; i < P; i++){
    p_center[i]->setVisible(checked);
  }
  Q_EMIT( changed());
}

void
MainWindow::on_actionClear_triggered()
{
  mc.clear();
  me.clear();
  points.clear();
  convex_hull.clear();
  min_rectangle.clear();
  min_parallelogram.clear();
  for(std::size_t i=0; i < P;i++){
    p_center[i]->hide();
  }
  Q_EMIT( changed());
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
  for(int i = 0; i < number_of_points; ++i){
    Point_2 p = *pg++;
    mc.insert(p);
    me.insert(p);
    points.push_back(p);
  }

  update_from_points();

  // default cursor
  QApplication::restoreOverrideCursor();
  Q_EMIT( changed());
}


void
MainWindow::on_actionLoadPoints_triggered()
{
  QString fileName = QFileDialog::getOpenFileName(this,
                                                  tr("Open Points file"),
                                                  ".",
                                                  tr("CGAL files (*.pts.cgal);;"
                                                   #if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
                                                     "WKT files (*.WKT *.wkt);;"
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
  if(fileName.endsWith(".wkt", Qt::CaseInsensitive))
  {
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
    CGAL::read_multi_point_WKT(ifs, points);
    for(K::Point_2 p : points)
    {
      mc.insert(p);
      me.insert(p);
    }
#endif
  }
  else
  {
    K::Point_2 p;
    while(ifs >> p) {
      mc.insert(p);
      me.insert(p);
      points.push_back(p);
    }
  }
  update_from_points();

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
                                                     "WKT files (*.WKT *.wkt);;"
                                                   #endif
                                                     "All files (*)"));
  if(! fileName.isEmpty()){
    std::ofstream ofs(qPrintable(fileName));
    if(fileName.endsWith(".wkt", Qt::CaseInsensitive))
    {
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
      std::vector<K::Point_2> out_pts;
      out_pts.reserve(std::distance(mc.points_begin(),
                                    mc.points_end()));
      for(Min_circle::Point_iterator pit = mc.points_begin();
          pit != mc.points_end(); ++pit)
        out_pts.push_back(*pit);
      CGAL::write_multi_point_WKT(ofs, out_pts);
#endif
    }
    else
    {
      for(Min_circle::Point_iterator
          vit = mc.points_begin(),
          end = mc.points_end();
          vit!= end; ++vit)
      {
        ofs << *vit << std::endl;
      }
    }
  }
}


void
MainWindow::on_actionRecenter_triggered()
{
  this->graphicsView->setSceneRect(cgi->boundingRect());
  this->graphicsView->fitInView(cgi->boundingRect(), Qt::KeepAspectRatio);
}


#include "Bounding_volumes.moc"
#include <CGAL/Qt/resources.h>

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Bounding_volumes demo");

  // Import resources from libCGAL (Qt5).
  // See https://doc.qt.io/qt-5/qdir.html#Q_INIT_RESOURCE
  CGAL_QT_INIT_RESOURCES;

  MainWindow mainWindow;
  mainWindow.show();
  return app.exec();
}
