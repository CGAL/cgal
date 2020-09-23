#include <fstream>

// CGAL headers

#define USE_CORE_EXPR_KERNEL

#ifndef USE_CORE_EXPR_KERNEL
  #include <CGAL/Exact_circular_kernel_2.h>
  #include <CGAL/Hyperbolic_Delaunay_triangulation_CK_traits_2.h>
  #include <internal/Qt/HyperbolicPainterOstreamCK.h>
#else
  #include <CGAL/Cartesian.h>
  #include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
  #include <internal/Qt/HyperbolicPainterOstream.h>
#endif

#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/point_generators_2.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <QGraphicsEllipseItem>

// GraphicsView items and event filters (input classes)
#include <internal/Qt/TriangulationCircumcircle.h>
#include <internal/Qt/TriangulationConflictZone.h>
#include <internal/Qt/TriangulationPointInputAndConflictZone.h>
#include <internal/Qt/TriangulationGraphicsItem.h>
#include <internal/Qt/HyperbolicVoronoiGraphicsItem.h>

// for viewportsBbox
#include <CGAL/Qt/utility.h>

// the two base classes
#include <CGAL/Qt/DemosMainWindow.h>

#include "ui_HDT2.h"


#ifndef USE_CORE_EXPR_KERNEL
  typedef CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2<> K;
#else
  typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<> K;
#endif


typedef K::Point_2 Point_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;

typedef CGAL::Hyperbolic_Delaunay_triangulation_2<K> Delaunay;

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

  //CGAL::Qt::TriangulationMovingPoint<Delaunay> * mp;
  CGAL::Qt::TriangulationConflictZone<Delaunay> * cz;
  //CGAL::Qt::TriangulationRemoveVertex<Delaunay> * trv;
  CGAL::Qt::TriangulationPointInputAndConflictZone<Delaunay> * pi;
  CGAL::Qt::TriangulationCircumcircle<Delaunay> *tcc;
public:
  MainWindow();

public slots:

  void processInput(CGAL::Object o);

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

  disk = new QGraphicsEllipseItem(left_top_corner_x, left_top_corner_y, width, height);

  QPen pen;  // creates a default pen

  pen.setWidthF(0.025);
  pen.setBrush(Qt::black);
  disk->setPen(pen);

  scene.addItem(disk);

  // Add a GraphicItem for the Delaunay triangulation
  dgi = new CGAL::Qt::TriangulationGraphicsItem<Delaunay>(&dt);

  QObject::connect(this, SIGNAL(changed()),
                   dgi, SLOT(modelChanged()));

  QPen vpen;
  vpen.setStyle(::Qt::SolidLine);
  vpen.setWidth(15);
  vpen.setBrush(::Qt::red);
  vpen.setCapStyle(::Qt::RoundCap);
  vpen.setJoinStyle(::Qt::RoundJoin);
  dgi->setVerticesPen(vpen);

  QPen epen;
  epen.setWidthF(0.01);
  epen.setBrush(::Qt::black);
  dgi->setEdgesPen(epen);

  scene.addItem(dgi);

  // Add a GraphicItem for the Voronoi diagram
  vgi = new CGAL::Qt::VoronoiGraphicsItem<Delaunay>(&dt);

  QObject::connect(this, SIGNAL(changed()),
                   vgi, SLOT(modelChanged()));

  vgi->setEdgesPen(QPen(Qt::blue, 0.01, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(vgi);
  vgi->hide();

  // Setup input handlers. They get events before the scene gets them
  // and the input they generate is passed to the triangulation with
  // the signal/slot mechanism
  pi = new CGAL::Qt::TriangulationPointInputAndConflictZone<Delaunay>(&scene, &dt, this);

  QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
                   this, SLOT(processInput(CGAL::Object)));

  tcc = new CGAL::Qt::TriangulationCircumcircle<Delaunay>(&scene, &dt, this);
  tcc->setPen(QPen(Qt::red, 0.005, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

  cz = new CGAL::Qt::TriangulationConflictZone<Delaunay>(&scene, &dt, this);

  //
  // Manual handling of actions
  //

  QObject::connect(this->actionQuit, SIGNAL(triggered()),
                   this, SLOT(close()));

  // We put mutually exclusive actions in an QActionGroup
  QActionGroup* ag = new QActionGroup(this);
  ag->addAction(this->actionInsertPoint);
  //ag->addAction(this->actionMovingPoint);
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
  this->graphicsView->shear(230, 230);

  // Turn the vertical axis upside down
  this->graphicsView->transform().scale(1, -1);

  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_Hyperbolic_Delaunay_triangulation_2.html");
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
      dt.insert(p);
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
  } else {
    scene.removeEventFilter(pi);
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
  CGAL::Random_points_in_disc_2<Point_2> pg(1);
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
  QString fileName = QFileDialog::getSaveFileName(this,
                                                  tr("Save points"),
                                                  ".");
  if(! fileName.isEmpty()){
    std::ofstream ofs(qPrintable(fileName));
    for(Delaunay::All_vertices_iterator
          vit = dt.all_vertices_begin(),
          end = dt.all_vertices_end();
        vit!= end; ++vit)
    {
      ofs << dt.point(vit) << std::endl;
    }
  }
}


void
MainWindow::on_actionRecenter_triggered()
{
  this->graphicsView->setSceneRect(dgi->boundingRect());
  this->graphicsView->fitInView(dgi->boundingRect(), Qt::KeepAspectRatio);
}


#include "HDT2.moc"

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Delaunay_triangulation_2 demo");

  // Import resources from libCGALQt5
  // See http://doc.qt.io/qt-5/qdir.html#Q_INIT_RESOURCE
  CGAL_Qt_init_resources();// that function is in a DLL

  MainWindow mainWindow;
  mainWindow.show();

  QStringList args = app.arguments();
  args.removeAt(0);
  Q_FOREACH(QString filename, args) {
    mainWindow.open(filename);
  }

  return app.exec();
}
