//#define CGAL_USE_BOOST_BIMAP

#include <fstream>
#include <vector>

// CGAL headers

#include "svd-typedefs.h"
#include <boost/config.hpp>
#include <boost/version.hpp>
#include <CGAL/Timer.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QMessageBox>
#include <QGraphicsLineItem>

// GraphicsView items and event filters (input classes)
#include <CGAL/Qt/GraphicsViewPolylineInput.h>
#include <CGAL/Qt/SegmentDelaunayGraphGraphicsItem.h>
#include <CGAL/Constraints_loader.h>
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
#include <CGAL/IO/WKT.h>
#endif
//#include <CGAL/Qt/Converter.h>

// the two base classes
#include "ui_Segment_voronoi_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

// for viewportsBbox(QGraphicsScene*)
#include <CGAL/Qt/utility.h>

typedef Rep K;
typedef SDG_2 SVD;

typedef SVD::Vertex_handle Vertex_handle;



class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Segment_voronoi_2
{
  Q_OBJECT
  
private:  
  SVD svd; 
  QGraphicsScene scene;
  std::list<Point_2> seeds;

  CGAL::Qt::SegmentDelaunayGraphGraphicsItem<SVD> * sdggi;

  CGAL::Qt::GraphicsViewPolylineInput<K> * pi;

public:
  MainWindow();

private:
  template <typename Iterator> 
  void insert_polyline(Iterator b, Iterator e)
  {
    Point_2 p, q;
    SVD::Vertex_handle vh, wh;
    Iterator it = b;
    vh = svd.insert(*it);
    p = *it;
    ++it;
    for(; it != e; ++it){
      q = *it;
      if(p != q){
        wh = svd.insert(*it);
        svd.insert(vh,wh);
        vh = wh;
        p = q;
      } else {
        std::cout << "duplicate point: " << p << std::endl; 
      }
    }
    Q_EMIT( changed());
  }

protected Q_SLOTS:
 virtual void open(QString);

public Q_SLOTS:

  void processInput(CGAL::Object o);

  void on_actionInsertPolyline_toggled(bool checked);
  
  void on_actionClear_triggered();

  void on_actionRecenter_triggered();

  void on_actionLoadSegments_triggered();

  void loadPolygonConstraints(QString);

  void loadEdgConstraints(QString);
  void loadWKTConstraints(QString fileName);

Q_SIGNALS:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow()
{
  setupUi(this);

  this->graphicsView->setAcceptDrops(false);

  // Add a GraphicItem for the SVD triangulation
  sdggi = new CGAL::Qt::SegmentDelaunayGraphGraphicsItem<SVD>(&svd);
  QColor segmentColor(::Qt::blue);
  QColor voronoiColor(::Qt::black);
  segmentColor.setAlpha(150);
  sdggi->setSegmentPen(QPen(segmentColor,0));
  sdggi->setVoronoiPen(QPen(voronoiColor,0));
    
  QObject::connect(this, SIGNAL(changed()),
		   sdggi, SLOT(modelChanged()));

  sdggi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

  sdggi->setZValue(-1);
  scene.addItem(sdggi);

  // Setup input handlers. They get events before the scene gets them
  // and the input they generate is passed to the triangulation with 
  // the signal/slot mechanism    
  pi = new CGAL::Qt::GraphicsViewPolylineInput<K>(this, &scene, 0, true); // inputs polylines which are not closed
  QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
		   this, SLOT(processInput(CGAL::Object)));
    


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
  this->actionShowVoronoi->setChecked(true);
  this->actionShowConstraints->setChecked(true);

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
  this->addAboutDemo(":/cgal/help/about_Segment_voronoi_2.html");
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
      svd.insert(points.front());
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
MainWindow::on_actionInsertPolyline_toggled(bool checked)
{
  if(checked){
    scene.installEventFilter(pi);
  } else {
    scene.removeEventFilter(pi);
  }
}



void
MainWindow::on_actionClear_triggered()
{
  svd.clear();
  Q_EMIT( changed());
}


void 
MainWindow::open(QString fileName)
{
  if(! fileName.isEmpty()){
    if(fileName.endsWith(".polygons.cgal")){
      loadPolygonConstraints(fileName);
      this->addToRecentFiles(fileName);
    } else if(fileName.endsWith(".edg")){
      loadEdgConstraints(fileName);
      this->addToRecentFiles(fileName);
    } else if(fileName.endsWith(".wkt", Qt::CaseInsensitive)){
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
      loadWKTConstraints(fileName);
      this->addToRecentFiles(fileName);
#endif
    }
  }
}

void
MainWindow::on_actionLoadSegments_triggered()
{
  QString fileName = QFileDialog::getOpenFileName(this,
						  tr("Open Constraint File"),
						  ".",
						  tr("Edge files (*.edg);;"
                                                     "Polyline files (*.polygons.cgal);;"
                                                   #if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
                                                     "WKT files (*.wkt *.WKT)"
                                                   #endif
                                                     ));
  open(fileName);
}

void
MainWindow::loadPolygonConstraints(QString fileName)
{
  K::Point_2 p,q, first;
  SVD::Vertex_handle vp, vq, vfirst;
  std::ifstream ifs(qPrintable(fileName));
  int n;
  while(ifs >> n){
    ifs >> first;
    p = first;
    vfirst = vp = svd.insert(p);
    n--;
    while(n--){
      ifs >> q;
      vq = svd.insert(q, vp);
      svd.insert(vp,vq);
      p = q;
      vp = vq;
    }
    svd.insert(vp, vfirst);
  }
  
  
  Q_EMIT( changed());
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
  
  K::Point_2 p,q, qold(0,0); // Initialize qold, as otherwise some g++ issue a unjustified warning

  SVD::Vertex_handle vp, vq, vqold;
  while(ifs >> p) {
    ifs >> q;
    if(p == q){
      continue;
    }
    if((!first) && (p == qold)){
      vp = vqold;
    } else {
      vp = svd.insert(p);
    }
    vq = svd.insert(q, vp);
    svd.insert(vp,vq);
    qold = q;
    vqold = vq;
    first = false;
  }


  tim.stop();
  statusBar()->showMessage(QString("Insertion took %1 seconds").arg(tim.time()), 2000);
  // default cursor
  QApplication::restoreOverrideCursor();
  Q_EMIT( changed());
  actionRecenter->trigger();
}

void 
MainWindow::loadWKTConstraints(QString 
                               #if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
                               fileName
                               #endif
                               )
{
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
  typedef CGAL::Polygon_with_holes_2<K> Polygon;
  typedef std::vector<K::Point_2> LineString;
  
  //Multipolygon 
  K::Point_2 p,q, first;
  SVD::Vertex_handle vp, vq, vfirst;
  std::ifstream ifs(qPrintable(fileName));
  do{
    std::vector<Polygon> polygons;
    CGAL::read_multi_polygon_WKT(ifs, polygons);
    for(const Polygon& poly : polygons)
    {
      if(poly.outer_boundary().is_empty())
        continue;
      Polygon::General_polygon_2::const_iterator it
          =poly.outer_boundary().begin();
      first = *(it++);
      p = first;
      vfirst = vp = svd.insert(p);
      for(; it !=
          poly.outer_boundary().end();
          ++it){
        q = *it;
        vq = svd.insert(q, vp);
        svd.insert(vp,vq);
        p = q;
        vp = vq;
      }
      if(vp != vfirst)
        svd.insert(vp, vfirst);
    }
  }while(ifs.good() && !ifs.eof());
  //MultiLineString
  ifs.clear();
  ifs.seekg(ifs.beg);
  K::Point_2 qold(0,0); // Initialize qold, as otherwise some g++ issue a unjustified warning
  SVD::Vertex_handle vqold;
  do{
    std::vector<LineString > linestrings;
    CGAL::read_multi_linestring_WKT(ifs, linestrings);
    for(const LineString& ls : linestrings)
    {
      bool first_pass=true;      
      LineString::const_iterator it = ls.begin();
      for(; it!=ls.end(); ++it){
        p = *it++;
        q = *it;
        if(p == q){
          continue;
        }
        if((!first_pass) && (p == qold)){
          vp = vqold;
        } else {
          vp = svd.insert(p);
        }
        vq = svd.insert(q, vp);
        if(vp != vq)
          svd.insert(vp,vq);
        qold = q;
        vqold = vq;
        first_pass = false;
      }
    }
  }while(ifs.good() && !ifs.eof());
  Q_EMIT( changed());
  actionRecenter->trigger();
#endif
}

void
MainWindow::on_actionRecenter_triggered()
{
  this->graphicsView->setSceneRect(sdggi->boundingRect());
  this->graphicsView->fitInView(sdggi->boundingRect(), Qt::KeepAspectRatio);  
}

#include "Segment_voronoi_2.moc"
#include <CGAL/Qt/resources.h>

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Segment Voronoi 2 demo");

  // Import resources from libCGAL (Qt5).
  CGAL_QT_INIT_RESOURCES;

  MainWindow mainWindow;
  mainWindow.show();
  return app.exec();
}
