#include <fstream>

#include <boost/config.hpp>
#include <boost/version.hpp>
// CGAL headers
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Snap_rounding_traits_2.h>
#include <CGAL/Snap_rounding_2.h>
#include <CGAL/Snap_rounding_traits_2.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>

// GraphicsView items and event filters (input classes)
#include <CGAL/Qt/RegularGridGraphicsItem.h>
#include <CGAL/Qt/SegmentsGraphicsItem.h>
#include <CGAL/Qt/PolylinesGraphicsItem.h>
#include <CGAL/Qt/GraphicsViewPolylineInput.h>
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
#include <CGAL/IO/WKT.h>
#endif

// for viewportsBbox
#include <CGAL/Qt/utility.h>

// the two base classes
#include "ui_Snap_rounding_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Snap_rounding_traits_2<K>     Traits;

typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;



class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Snap_rounding_2
{
  Q_OBJECT

private:

  QGraphicsScene scene;

  CGAL::Qt::RegularGridGraphicsItem<K> * rgi;

  CGAL::Qt::GraphicsViewPolylineInput<K> * pi;

  std::list<Segment_2> input;
  std::list<std::list<Point_2> > output;

  typedef CGAL::Qt::SegmentsGraphicsItem<std::list<Segment_2> > InputSegmentsGraphicsItem;
  typedef CGAL::Qt::PolylinesGraphicsItem<std::list<std::list<Point_2> > > OutputPolylinesGraphicsItem;
  InputSegmentsGraphicsItem * isgi;
  OutputPolylinesGraphicsItem *plgi;
  double delta;

public:
  MainWindow();

  void resize(){
  this->graphicsView->setSceneRect(QRectF(0,0,20, 20));
  this->graphicsView->fitInView(0,0, 20, 20, Qt::KeepAspectRatio);
  }

public Q_SLOTS:

  void processInput(CGAL::Object o);

  void on_actionLoadSegments_triggered();

  void on_actionClear_triggered();

  void on_actionSaveSegments_triggered();

  void on_actionRecenter_triggered();

  void on_actionShowGrid_toggled(bool checked);
  void on_actionShowInput_toggled(bool checked);
  void on_actionShowSnappedSegments_toggled(bool checked);

  void deltaChanged(double);

  virtual void open(QString fileName);

Q_SIGNALS:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow(), delta(1.0)
{
  setupUi(this);

  this->graphicsView->setAcceptDrops(false);

  isgi = new InputSegmentsGraphicsItem(&input);
  scene.addItem(isgi);

  plgi = new OutputPolylinesGraphicsItem(&output);
  scene.addItem(plgi);

 // inputs polylines with 2 points
  pi = new CGAL::Qt::GraphicsViewPolylineInput<K>(this, &scene, 2, false);
  QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
                   this, SLOT(processInput(CGAL::Object)));

  scene.installEventFilter(pi);

  // Manual handling of actions
  //


  QObject::connect(this->doubleSpinBox, SIGNAL(valueChanged(double)),
                   this, SLOT(deltaChanged(double)));

  QObject::connect(this->actionQuit, SIGNAL(triggered()),
                   this, SLOT(close()));

  //
  // Setup the scene and the view
  //
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  this->graphicsView->setScene(&scene);
  // Turn the vertical axis upside down
  this->graphicsView->transform().scale(1, -1);
  this->graphicsView->setMouseTracking(true);

  rgi = new CGAL::Qt::RegularGridGraphicsItem<K>(delta, delta);

    QObject::connect(this, SIGNAL(changed()),
                     rgi, SLOT(modelChanged()));

    QObject::connect(this, SIGNAL(changed()),
                     isgi, SLOT(modelChanged()));

    QObject::connect(this, SIGNAL(changed()),
                     plgi, SLOT(modelChanged()));


  rgi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  rgi->setEdgesPen(QPen(Qt::gray, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(rgi);

  plgi->setEdgesPen(QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_Snap_rounding_2.html");
  this->addAboutCGAL();

  this->addRecentFiles(this->menuFile, this->actionQuit);
  connect(this, SIGNAL(openRecentFile(QString)),
          this, SLOT(open(QString)));
}



void
MainWindow::deltaChanged(double d)
{
  if(delta == d){
    return;
  }
  delta = d;
  output.clear();
  CGAL::snap_rounding_2<Traits,std::list<Segment_2>::const_iterator,std::list<std::list<Point_2> > >(input.begin(), input.end(), output, delta, true, false);
  rgi->setDelta(delta, delta);
  Q_EMIT( changed());
}


void
MainWindow::processInput(CGAL::Object o)
{

  std::list<Point_2> points;
  if(CGAL::assign(points, o)){
    if(points.size() == 2) {
      input.push_back(Segment_2(points.front(), points.back()));
      output.clear();
      CGAL::snap_rounding_2<Traits,std::list<Segment_2>::const_iterator,std::list<std::list<Point_2> > >(input.begin(), input.end(), output, delta, true, false);
    }
    else {
      std::cerr << points.size() << std::endl;
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
MainWindow::on_actionClear_triggered()
{
  input.clear();
  output.clear();
  Q_EMIT( changed());
}


void
MainWindow::on_actionShowGrid_toggled(bool checked)
{
  rgi->setVisible(checked);
  Q_EMIT( changed());
}

void
MainWindow::on_actionShowInput_toggled(bool checked)
{
  isgi->setVisible(checked);
  Q_EMIT( changed());
}



void
MainWindow::on_actionShowSnappedSegments_toggled(bool checked)
{
  plgi->setVisible(checked);
  Q_EMIT( changed());
}




void
MainWindow::on_actionLoadSegments_triggered()
{
  QString fileName = QFileDialog::getOpenFileName(this,
                                                  tr("Open segment file"),
                                                  ".",
                                                  tr("Edge files (*.edg);;"
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
  if(fileName.endsWith(".wkt", Qt::CaseInsensitive))
  {
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
    std::vector<std::vector<Point_2> > mls;
    CGAL::read_multi_linestring_WKT(ifs, mls);
    for(const std::vector<Point_2>& ls : mls)
    {
      if(ls.size() > 2)
        continue;
      Segment_2 seg(ls[0], ls[1]);
      input.push_back(seg);
    }
#endif
  }
  else {
    std::copy(std::istream_iterator<Segment_2>(ifs),
              std::istream_iterator<Segment_2>(),
              std::back_inserter(input));
  }
  output.clear();
  CGAL::snap_rounding_2<Traits,std::list<Segment_2>::const_iterator,std::list<std::list<Point_2> > >(input.begin(), input.end(), output, delta, true, false);
  ifs.close();
  // default cursor
  QApplication::restoreOverrideCursor();
  this->addToRecentFiles(fileName);
  on_actionRecenter_triggered();
  Q_EMIT( changed());
}

void
MainWindow::on_actionSaveSegments_triggered()
{
  QString fileName = QFileDialog::getSaveFileName(this,
                                                  tr("Save points"),
                                                  ".",
                                                  tr("Edge files (*.edg);;"
                                                   #if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
                                                     "WKT files (*.wkt *.WKT);;"
                                                   #endif
                                                     "All files (*)"));
  if(! fileName.isEmpty()){
    std::ofstream ofs(qPrintable(fileName));
    ofs.precision(12);
    if(fileName.endsWith(".wkt", Qt::CaseInsensitive))
    {
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
      std::vector<std::vector<Point_2> >mls;
      for(const Segment_2& seg : input)
      {
        std::vector<Point_2> ls(2);
        ls[0] = seg.source();
        ls[1] = seg.target();
        mls.push_back(ls);
      }
      CGAL::write_multi_linestring_WKT(ofs, mls);
#endif
    }
    else
      std::copy(input.begin(), input.end(),  std::ostream_iterator<Segment_2>(ofs, "\n"));
  }

}


void
MainWindow::on_actionRecenter_triggered()
{
  this->graphicsView->setSceneRect(isgi->boundingRect());
  this->graphicsView->fitInView(isgi->boundingRect(), Qt::KeepAspectRatio);
}


#include "Snap_rounding_2.moc"
#include <CGAL/Qt/resources.h>

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Snap_rounding_2 demo");

  // Import resources from libCGAL (Qt5).
  // See https://doc.qt.io/qt-5/qdir.html#Q_INIT_RESOURCE
  CGAL_QT_INIT_RESOURCES;
  Q_INIT_RESOURCE(Snap_rounding_2);

  MainWindow mainWindow;
  mainWindow.show();
  mainWindow.resize();
  return app.exec();
}
