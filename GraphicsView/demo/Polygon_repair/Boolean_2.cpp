#include <fstream>
// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Polygon_repair/repair.h>
#include <boost/config.hpp>
#include <boost/version.hpp>
#include <CGAL/IO/WKT.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QFileDialog>
#include <QInputDialog>
#include <QGraphicsLineItem>

// GraphicsView items and event filters (input classes)
#include <CGAL/Qt/GraphicsViewPolylineInput.h>
#include <CGAL/Qt/GraphicsViewPolygonWithHolesInput.h>
#include <CGAL/Qt/PolygonGraphicsItem.h>
#include <CGAL/Qt/PolygonWithHolesGraphicsItem.h>
#include <CGAL/Qt/MultipolygonWithHolesGraphicsItem.h>
#include <CGAL/Qt/LineGraphicsItem.h>

// the two base classes
#include "ui_Boolean_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef K::Line_2 Line_2;

typedef CGAL::Polygon_2<K> Polygon2;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;
typedef CGAL::Multipolygon_with_holes_2<K> Multipolygon_with_holes_2;

typedef std::shared_ptr<Polygon2> PolygonPtr ;

typedef std::vector<PolygonPtr> PolygonPtr_vector ;

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Polygon_2
{
  Q_OBJECT

private:

  bool add_2_A = true;
  CGAL::Qt::Converter<K> convert;
  Polygon2 poly;
  Polygon_with_holes_2 pwhA, pwhB;
  Multipolygon_with_holes_2 mpwhA, mpwhB, mpwhC;
  QGraphicsScene scene;

  CGAL::Qt::MultipolygonWithHolesGraphicsItem<Multipolygon_with_holes_2> *mpwhAgi, *mpwhBgi, *mpwhCgi;

  CGAL::Qt::GraphicsViewPolygonWithHolesInput<K> * pi;


public:
  MainWindow();

public Q_SLOTS:

  void processInput(CGAL::Object o);

  void on_actionClear_triggered();

  void on_actionLoadPolygon_triggered();
  void on_actionSavePolygon_triggered();

  void on_actionRecenter_triggered();

  void on_actionAdd_to_A_triggered();
  void on_actionAdd_to_B_triggered();

  void on_actionCreateInputPolygon_toggled(bool);

  void clear();

  virtual void open(QString);
Q_SIGNALS:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow()
{
  setupUi(this);

  this->graphicsView->setAcceptDrops(false);

  mpwhAgi = new CGAL::Qt::MultipolygonWithHolesGraphicsItem<Multipolygon_with_holes_2>(&mpwhA);
  mpwhAgi->setBrush(QBrush(::Qt::green));
  mpwhAgi->setZValue(3);
  QObject::connect(this, SIGNAL(changed()),
                   mpwhAgi, SLOT(modelChanged()));

  scene.addItem(mpwhAgi);

  mpwhBgi = new CGAL::Qt::MultipolygonWithHolesGraphicsItem<Multipolygon_with_holes_2>(&mpwhB);
  mpwhBgi->setZValue(4);
  mpwhBgi->setBrush(QBrush(QColor(255, 0, 0, 100)));
  QObject::connect(this, SIGNAL(changed()),
                   mpwhBgi, SLOT(modelChanged()));

  scene.addItem(mpwhBgi);

  mpwhCgi = new CGAL::Qt::MultipolygonWithHolesGraphicsItem<Multipolygon_with_holes_2>(&mpwhC);
  mpwhCgi->setZValue(5);
  mpwhCgi->setBrush(QBrush(QColor(211, 211, 211, 150)));
  QObject::connect(this, SIGNAL(changed()),
                   mpwhCgi, SLOT(modelChanged()));

  scene.addItem(mpwhCgi);


  // Setup input handlers. They get events before the scene gets them
  // pi = new CGAL::Qt::GraphicsViewPolylineInput<K>(this, &scene, 0, true);
  pi = new CGAL::Qt::GraphicsViewPolygonWithHolesInput<K>(this, &scene);
  pi->setZValue(10);

  this->actionCreateInputPolygon->setChecked(true);
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
  this->graphicsView->scale(1, -1);

  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_Polygon_2.html");
  this->addAboutCGAL();
  // this->setupExportSVG(action_Export_SVG, graphicsView);

  this->addRecentFiles(this->menuFile, this->actionQuit);
  connect(this, SIGNAL(openRecentFile(QString)),
          this, SLOT(open(QString)));
}


void
MainWindow::processInput(CGAL::Object o)
{
  if(add_2_A){
    if(assign(pwhA, o)){
      mpwhA.add_polygon_with_holes(pwhA);
    }
  }else{
    if(assign(pwhB, o)){
      mpwhB.add_polygon_with_holes(pwhB);
    }
  }
  if((! mpwhA.is_empty()) && (! mpwhB.is_empty())){
    mpwhC = CGAL::Polygon_repair::join(mpwhA, mpwhB);
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
MainWindow::on_actionAdd_to_A_triggered()
{
  this->actionAdd_to_A->setEnabled(false);
  this->actionAdd_to_B->setEnabled(true);
  add_2_A = true;
}
void
MainWindow::on_actionAdd_to_B_triggered()
{
  this->actionAdd_to_B->setEnabled(false);
  this->actionAdd_to_A->setEnabled(true);
  add_2_A = false;
}

void
MainWindow::on_actionClear_triggered()
{
  pwhA.clear();
  mpwhA.clear();
  pwhB.clear();
  mpwhB.clear();
  mpwhC.clear();
  clear();
  this->actionCreateInputPolygon->setChecked(true);
  Q_EMIT( changed());
}


void
MainWindow::on_actionLoadPolygon_triggered()
{
  QString fileName = QFileDialog::getOpenFileName(this,
                                                  tr("Open Polygon File"),
                                                  ".",
                                                  tr( "WKT files (*.wkt *.WKT);;"
                                                      "All file (*)"));
  if(! fileName.isEmpty()){
    open(fileName);
  }
}

void
MainWindow::open(QString fileName)
{
  this->actionCreateInputPolygon->setChecked(false);
  std::ifstream ifs(qPrintable(fileName));
  pwhA.clear();
  if(fileName.endsWith(".wkt", Qt::CaseInsensitive))
  {
    CGAL::IO::read_polygon_WKT(ifs, pwhA);
  }
  else
  {
    std::cout << "not supported" << std::endl;
  }
  clear();

  this->addToRecentFiles(fileName);
  Q_EMIT( changed());
}


void
MainWindow::on_actionSavePolygon_triggered()
{
  QString fileName = QFileDialog::getSaveFileName(this,
                                                  tr("Save Polygon"),
                                                  ".",
                                                  tr( "WKT files (*.wkt *.WKT);;"
                                                      "All file (*)"));
  if(! fileName.isEmpty()){
    std::ofstream ofs(qPrintable(fileName));
    if(fileName.endsWith(".wkt", Qt::CaseInsensitive))
    {
      CGAL::IO::write_polygon_WKT(ofs, pwhA);
    }
    else{
      std::cout << "not supported" << std::endl;
    }
  }
}


void
MainWindow::on_actionCreateInputPolygon_toggled(bool checked)
{
  //  poly.clear();
  clear();
  if(checked){
    scene.installEventFilter(pi);
  } else {
    scene.removeEventFilter(pi);
  }
  Q_EMIT( changed());
}

void
MainWindow::on_actionRecenter_triggered()
{
  this->graphicsView->setSceneRect(mpwhAgi->boundingRect());
  this->graphicsView->fitInView(mpwhAgi->boundingRect(), Qt::KeepAspectRatio);
}

void
MainWindow::clear()
{}


#include "Boolean_2.moc"
#include <CGAL/Qt/resources.h>

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("2D Boolean Operations");

  // Import resources from libCGAL (Qt6).
  // See https://doc.qt.io/qt-5/qdir.html#Q_INIT_RESOURCE
  CGAL_QT_INIT_RESOURCES;
  Q_INIT_RESOURCE(Boolean_2);

  MainWindow mainWindow;
  mainWindow.show();
  return app.exec();
}
