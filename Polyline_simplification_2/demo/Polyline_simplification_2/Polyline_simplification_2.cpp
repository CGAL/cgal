//#define CGAL_USE_BOOST_BIMAP

#include <fstream>
#include <vector>
#include <deque>
#include <boost/config.hpp>
#include <boost/version.hpp>

// CGAL headers
#include <CGAL/Bbox_2.h>
#include <CGAL/assertions_behaviour.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h>
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
#include <CGAL/IO/WKT.h>
#endif

//#define CGAL_TESTING_POLYLINE_SIMPLIFICATION
//#define CGAL_POLYLINE_SIMPLIFICATION_TRACE_LEVEL 15

bool lAppToLog = false ;
void Polyline_simplification_2_external_trace( std::string m )
{
  std::ofstream out("polysim_log.txt", ( lAppToLog ? std::ios::app | std::ios::ate : std::ios::trunc | std::ios::ate ) );
  out << std::setprecision(19) << m << std::endl << std::flush ;
  lAppToLog = true ;
}

void error_handler ( char const* what, char const* expr, char const* file, int line, char const* msg )
{
  std::cerr << "CGAL error: " << what << " violation!" << std::endl
       << "Expr: " << expr << std::endl
       << "File: " << file << std::endl
       << "Line: " << line << std::endl;
  if ( msg != 0)
      std::cerr << "Explanation:" << msg << std::endl;

  throw std::runtime_error("CGAL Error");
}


#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>
#include <CGAL/Polyline_simplification_2/Scaled_squared_distance_cost.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QSlider>
#include <QProgressBar>

// GraphicsView items and event filters (input classes)
#include <CGAL/Qt/GraphicsViewPolylineInput.h>
#include <CGAL/Qt/Polyline_simplification_2_graphics_item.h>
#include <CGAL/Qt/Converter.h>

// the two base classes
#include "ui_Polyline_simplification_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

// for viewportsBbox(QGraphicsScene*)
#include <CGAL/Qt/utility.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_2         Point_2;
typedef K::Segment_2       Segment_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;

typedef CGAL::Polyline_simplification_2::Vertex_base_2<K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>    Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>       TDS;
typedef CGAL::Exact_predicates_tag                        Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K,TDS, Itag> CDT;

typedef CGAL::Constrained_triangulation_plus_2<CDT>       PCT;

namespace PS2 = CGAL::Polyline_simplification_2 ;

enum Mode { ABS_P, ABS_E, REL_P, REL_E, MIXED_P, MIXED_E } ;

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Polyline_simplification_2
{
  Q_OBJECT

private:

  PCT                                                m_pct;
  QGraphicsScene                                    mScene;
  CGAL::Qt::PolylineSimplificationGraphicsItem<PCT> * mGI;
  CGAL::Qt::GraphicsViewPolylineInput<K> *          mPI;

private:

public:

  MainWindow();

protected:
  void dragEnterEvent(QDragEnterEvent *event);
  void dropEvent(QDropEvent *event);
private:

  void loadPoly(QString);
  void loadOSM (QString);
  void loadWKT (QString);

protected Q_SLOTS:
void open(QString);

public Q_SLOTS:

  void processInput(CGAL::Object o);

  void on_actionShowTriangulation_toggled(bool checked);

  void on_actionInsertPolyline_toggled(bool checked);

  void on_actionClear_triggered();

  void on_actionRecenter_triggered();

  void on_actionLoadConstraints_triggered();

  void on_actionSimplify_triggered();

  Mode getSimplificationMode() ;

  double getThreshold() ;

Q_SIGNALS:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow()
{
  CGAL::set_error_handler  (error_handler);
  CGAL::set_warning_handler(error_handler);

  setupUi(this);

  setAcceptDrops(true);

  // Add a GraphicItem for the PS triangulation
  mGI = new CGAL::Qt::PolylineSimplificationGraphicsItem<PCT>(&m_pct);

  QObject::connect(this, SIGNAL(changed()), mGI, SLOT(modelChanged()));

  mGI->setVerticesPen(QPen(Qt::black, 2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  mGI->setUnremovableVerticesPen(QPen(Qt::blue, 2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  mGI->setZValue(-1);
  mGI->setVisibleEdges(false);
  mGI->setVisibleConstraints(true);

  mScene.addItem(mGI);

  // Setup input handlers. They get events before the mScene gets them
  // and the input they generate is passed to the triangulation with
  // the signal/slot mechanism
  mPI = new CGAL::Qt::GraphicsViewPolylineInput<K>(this, &mScene, 0, true); // inputs polylines which are not closed
  this->on_actionInsertPolyline_toggled(true);
  QObject::connect(mPI, SIGNAL(generate(CGAL::Object)), this, SLOT(processInput(CGAL::Object)));


  //
  // Manual handling of actions
  //
  QObject::connect(this->actionQuit, SIGNAL(triggered()), this, SLOT(close()));

  // We put mutually exclusive actions in an QActionGroup
  //  QActionGroup* ag = new QActionGroup(this);
  //  ag->addAction(this->actionInsertPolyline);

  this->actionShowTriangulation->setChecked(false);

  //
  // Setup the mScene and the view
  //
  mScene.setItemIndexMethod(QGraphicsScene::NoIndex);
  mScene.setSceneRect(-100, -100, 100, 100);
  this->graphicsView->setScene(&mScene);
  this->graphicsView->setMouseTracking(true);

  // Turn the vertical axis upside down
  this->graphicsView->scale(1, -1);

  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_Polyline_simplification_2.html");
  this->addAboutCGAL();
  this->setupExportSVG(action_Export_SVG, graphicsView);

  this->addRecentFiles(this->menuFile, this->actionQuit);
  connect(this, SIGNAL(openRecentFile(QString)), this, SLOT(open(QString)));
}


void
MainWindow::dragEnterEvent(QDragEnterEvent *event)
{
  if (event->mimeData()->hasFormat("text/uri-list"))
    event->acceptProposedAction();
}

void
MainWindow::dropEvent(QDropEvent *event)
{
  QString filename = event->mimeData()->urls().at(0).path();
  open(filename);
  event->acceptProposedAction();
}

void
MainWindow::processInput(CGAL::Object o)
{
  std::list<Point_2> points;
  if(CGAL::assign(points, o))
  {
    if(points.size() >= 2)
    {
      m_pct.insert_constraint(points.begin(), points.end());
#if 0
      std::ofstream out("polygon.txt");
      out.precision(12);
      out << points.size() << std::endl;
      for(std::list<Point_2>::iterator it = points.begin(); it!= points.end(); ++it){
        out << *it << std::endl;
      }
#endif
      Q_EMIT( changed());
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
    mScene.installEventFilter(mPI);
  } else {
    mScene.removeEventFilter(mPI);
  }
}

void
MainWindow::on_actionShowTriangulation_toggled(bool checked)
{
  mGI->setVisibleEdges(checked);
  update();
}

void
MainWindow::on_actionClear_triggered()
{
  m_pct.clear();
  mGI->modelChanged();
  Q_EMIT( changed());
}

Mode MainWindow::getSimplificationMode()
{
  bool lP     = rbuttonPercentage->isChecked();
  int  lError = rbuttonAbsError  ->isChecked() ? 0 : ( rbuttonRelError  ->isChecked() ? 1 : 2 ) ;

  return lP ? lError == 0 ? ABS_P : ( lError == 1 ? REL_P : MIXED_P )
            : lError == 0 ? ABS_E : ( lError == 1 ? REL_E : MIXED_E ) ;
}

double MainWindow::getThreshold()
{
  double lValue = this->textValue->text().toDouble() ;
  return rbuttonPercentage->isChecked() ? lValue / 100.0
           : rbuttonAbsError->isChecked() ? lValue * lValue
              : lValue  ;
}

void MainWindow::on_actionSimplify_triggered()
{
  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);


  try
  {
    switch( getSimplificationMode() )
    {
    case ABS_P   : simplify(m_pct, PS2::Squared_distance_cost(),  PS2::Stop_below_count_ratio_threshold(getThreshold())   ) ;
 break ;
    case ABS_E   : simplify(m_pct, PS2::Squared_distance_cost(), PS2::Stop_above_cost_threshold(getThreshold()) ) ; break ;
    case REL_P   : simplify(m_pct, PS2::Scaled_squared_distance_cost(),  PS2::Stop_below_count_ratio_threshold(getThreshold()) ) ; break ;
    case REL_E   : simplify(m_pct, PS2::Scaled_squared_distance_cost(), PS2::Stop_above_cost_threshold(getThreshold()) ) ; break ;
    case MIXED_P : simplify(m_pct, PS2::Hybrid_squared_distance_cost<double>(1.0), PS2::Stop_below_count_ratio_threshold(getThreshold())) ; break ;
 case MIXED_E : simplify(m_pct, PS2::Hybrid_squared_distance_cost<double>(1.0), PS2::Stop_above_cost_threshold(getThreshold())) ; break ;

      break ;
    }

    statusBar()->showMessage(QString("Simplification done"));
  }
  catch(...)
  {
    statusBar()->showMessage(QString("Exception ocurred"));
  }

   // default cursor
  QApplication::restoreOverrideCursor();
  mGI->modelChanged();
  Q_EMIT( changed());
}


void
MainWindow::open(QString fileName)
{
  if(! fileName.isEmpty()){
    if(fileName.endsWith(".osm")){
      loadOSM(fileName);
      this->addToRecentFiles(fileName);
    }
    else if(fileName.endsWith(".wkt")){
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
      loadWKT(fileName);
      this->addToRecentFiles(fileName);
#endif
    }
  }
}

void
MainWindow::on_actionLoadConstraints_triggered()
{
  QString fileName = QFileDialog::getOpenFileName(this,
                                                  tr("Open Constraint File"),
                                                  "../data"
                                                #if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
                                                  ,tr("Polylines (*.osm *.wkt);;")
                                                #endif
                                                  );
  open(fileName);
}

std::string trim_right ( std::string str )
{
  if ( str.length() > 0 )
  {
    std::size_t pos = str.find_last_not_of( " " ) ;
    if ( pos != std::string::npos )
      return str.substr(0,pos+1);
  }

  return std::string("") ;
}

void MainWindow::loadWKT(QString
                         #if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
                         fileName
                         #endif
                         )
{
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
    typedef std::vector<Point_2> MultiPoint;

  typedef std::vector<Point_2> LineString;
  typedef std::deque<LineString> MultiLineString;

  typedef CGAL::Polygon_2<K> Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;
  typedef std::deque<Polygon_with_holes_2> MultiPolygon;

  std::ifstream ifs(qPrintable(fileName));
  MultiPoint points;
  MultiLineString polylines;
  MultiPolygon polygons;
  CGAL::read_WKT(ifs, points,polylines,polygons);

  m_pct.clear();
  mGI->modelChanged();

  if(! points.empty()){
    std::cout << "Ignore " << points.size() << " isolated points" << std::endl;
  }
  for(LineString poly : polylines){
    if ( poly.size() > 2 ){
      m_pct.insert_constraint(poly.begin(), poly.end());
    }
  }
 for(Polygon_with_holes_2 poly : polygons){
   m_pct.insert_constraint(poly.outer_boundary().vertices_begin(), poly.outer_boundary().vertices_end());
   for(Polygon_with_holes_2::Hole_const_iterator it = poly.holes_begin(); it != poly.holes_end(); ++it){
     const Polygon_2& hole = *it;
      m_pct.insert_constraint(hole);
    }
 }

  Q_EMIT( changed());

  actionRecenter->trigger();
#endif
}


void MainWindow::loadOSM(QString fileName)
{
  m_pct.clear();
  mGI->modelChanged();

  try
  {
    std::ifstream ifs(qPrintable(fileName));

    std::string line ;

    std::vector<Point_2> poly ;

    while ( std::getline(ifs,line) )
    {
      line = trim_right(line);

      if ( line.size() > 0 )
      {
        if ( line.find(':') != std::string::npos )
        {
          if ( poly.size() > 0 )
          {
            if ( poly.front() == poly.back() && poly.size() >= 4 )
            {
              if ( is_simple_2(poly.begin(), poly.end() - 1) )
                m_pct.insert_constraint(poly.begin(), poly.end() - 1 ) ;
            }
            else if ( poly.size() >= 2 )
            {
              if ( is_simple_2(poly.begin(), poly.end()) )
                m_pct.insert_constraint(poly.begin(), poly.end() ) ;
            }
          }
          poly.clear();
        }
        else
        {
          double x,y ;

          std::string::size_type pos = line.find(',');
          if ( pos != std::string::npos )
            line[pos]= ' ' ;

          std::istringstream ss(line);

          ss >> x >> y ;

          poly.push_back( Point_2(x,y) );
        }
      }
    }

    if ( poly.size() > 0 )
    {
      if ( poly.front() == poly.back() )
      {
        m_pct.insert_constraint(poly.begin(), poly.end() - 1 ) ;
      }
      else
      {
        m_pct.insert_constraint(poly.begin(), poly.end() ) ;
      }
    }
  }
  catch(...)
  {
    statusBar()->showMessage(QString("Exception ocurred"));
  }

  Q_EMIT( changed());

  actionRecenter->trigger();
}

void
MainWindow::on_actionRecenter_triggered()
{
  this->graphicsView->setSceneRect(mGI->boundingRect());
  this->graphicsView->fitInView(mGI->boundingRect(), Qt::KeepAspectRatio);
}

#include "Polyline_simplification_2.moc"

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Polyline_simplification_2 demo");

  // Import resources from libCGALQt5.
  CGAL_QT_INIT_RESOURCES;

  MainWindow mainWindow;
  mainWindow.show();
  return app.exec();
}
