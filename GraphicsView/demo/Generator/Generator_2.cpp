#include <fstream>
// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/function_objects.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/algorithm.h>
#ifndef Q_MOC_RUN
#include <CGAL/random_convex_hull_in_disc_2.h>
#endif


// Qt headers
#include <QtGui>
#include <QString>
#include <QInputDialog>
#include <QFileDialog>
#include <QGraphicsLineItem>

// GraphicsView items and event filters (input classes)
#include <CGAL/Qt/PointsGraphicsItem.h>
#include <CGAL/Qt/utility.h>
#include <CGAL/Qt/SegmentsGraphicsItem.h>
  
// the two base classes
#include "ui_Generator_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Vector_2 Vector_2;
typedef K::Segment_2 Segment_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Generator_2
{
  Q_OBJECT
  
private:



  CGAL::Qt::Converter<K> convert;
  std::vector<Point_2> points; 
  std::vector<Segment_2> segments; 
  QGraphicsScene scene;  

  CGAL::Qt::PointsGraphicsItem<std::vector<Point_2> > * pgi;
  CGAL::Qt::SegmentsGraphicsItem<std::vector<Segment_2> > * sgi;


  template <typename G>
  void
  on_actionGenerate_triggered()
  {
    QRectF rect = CGAL::Qt::viewportsBbox(&scene);
    CGAL::Qt::Converter<K> convert;  
    Iso_rectangle_2 isor = convert(rect);
    Point_2 center = CGAL::midpoint(isor[0], isor[2]);
    Vector_2 offset = center - CGAL::ORIGIN;
    double w = isor.xmax() - isor.xmin();
    double h = isor.ymax() - isor.ymin();
    double radius = (w<h) ? w/2 : h/2;

    G pg(radius);
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

    points.reserve(points.size() + number_of_points);
    for(int i = 0; i < number_of_points; ++i){
      points.push_back(*pg + offset);
      ++pg;
    }
    // default cursor
    QApplication::restoreOverrideCursor();
    Q_EMIT( changed());
  }

public:
  MainWindow();

public Q_SLOTS:

  void on_actionClear_triggered();


  void on_actionRecenter_triggered();
  void on_actionGeneratePointsOnCircle_triggered();
  void on_actionGeneratePointsInSquare_triggered();
  void on_actionGeneratePointsInDisc_triggered();
  void on_actionGenerateSegments_triggered();
  void on_actionGenerateSegmentFans_triggered();
  void on_actionGeneratePolytopeInDisc_triggered();
  void clear();

Q_SIGNALS:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow()
{
  setupUi(this);

  // Add a GraphicItem for the point set
  pgi = new CGAL::Qt::PointsGraphicsItem<std::vector<Point_2> >(&points);
  sgi = new CGAL::Qt::SegmentsGraphicsItem<std::vector<Segment_2> >(&segments);

  QObject::connect(this, SIGNAL(changed()),
		   pgi, SLOT(modelChanged()));


    QObject::connect(this, SIGNAL(changed()),
  		   sgi, SLOT(modelChanged()));

  pgi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  sgi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(pgi);
  scene.addItem(sgi);


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

  // Uncomment the following line to get antialiasing by default.
//   actionUse_Antialiasing->setChecked(true);

  // Turn the vertical axis upside down
  this->graphicsView->scale(1, -1);
                                                      
  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_Generator_2.html");
  this->addAboutCGAL();

}


/* 
 *  Qt Automatic Connections
 *  http://doc.qt.io/qt-5/designer-using-a-ui-file.html#automatic-connections
 * 
 *  setupUi(this) generates connections to the slots named
 *  "on_<action_name>_<signal_name>"
 */

void
MainWindow::on_actionClear_triggered()
{
  clear();
  Q_EMIT( changed());
}

void
MainWindow::on_actionRecenter_triggered()
{
  this->graphicsView->setSceneRect(pgi->boundingRect());
  this->graphicsView->fitInView(pgi->boundingRect(), Qt::KeepAspectRatio);  
}

void
MainWindow::on_actionGeneratePointsOnCircle_triggered()
{
  typedef CGAL::Random_points_on_circle_2<Point_2> Generator;
  on_actionGenerate_triggered<Generator>();
}


void
MainWindow::on_actionGeneratePointsInSquare_triggered()
{
  typedef CGAL::Random_points_in_square_2<Point_2> Generator;
  on_actionGenerate_triggered<Generator>();
}


void
MainWindow::on_actionGeneratePointsInDisc_triggered()
{
  typedef CGAL::Random_points_in_disc_2<Point_2> Generator;
  on_actionGenerate_triggered<Generator>();
}


void
MainWindow::on_actionGenerateSegments_triggered()
{
  segments.reserve(segments.size() + 200);
  
  // Prepare point generator for the horizontal segment, length 200.
  typedef  CGAL::Random_points_on_segment_2<Point_2>  Rpos_generator;
  Rpos_generator rpos( Point_2(-100,0), Point_2(100,0));
  
  // Prepare point generator for random points on circle, radius 250.
  typedef  CGAL::Random_points_on_circle_2<Point_2>  Rpoc_generator;
  Rpoc_generator rpoc( 250);
  
  // Create 200 segments.
  typedef CGAL::Creator_uniform_2< Point_2, Segment_2> Seg_creator;
  typedef CGAL::Join_input_iterator_2< Rpos_generator, Rpoc_generator, Seg_creator> Seg_iterator;
  Seg_iterator g( rpos, rpoc);
  CGAL::cpp11::copy_n( g, 200, std::back_inserter(segments));
  
  Q_EMIT( changed());
}


void
MainWindow::on_actionGenerateSegmentFans_triggered()
{
  typedef CGAL::Points_on_segment_2<Point_2>                PG;
  typedef CGAL::Creator_uniform_2< Point_2, Segment_2>        Creator;
  typedef CGAL::Join_input_iterator_2< PG, PG, Creator>   Segm_iterator;
  typedef CGAL::Counting_iterator<Segm_iterator,Segment_2>  Count_iterator;

  segments.reserve(segments.size() + 100);

  // A horizontal like fan.
  PG p1( Point_2(-250, -50), Point_2(-250, 50),50);   // Point generator.
  PG p2( Point_2( 250,-250), Point_2( 250,250),50);
  Segm_iterator  t1( p1, p2);                     // Segment generator.
  Count_iterator t1_begin( t1);                   // Finite range.
  Count_iterator t1_end(t1, 50);
  std::copy( t1_begin, t1_end, std::back_inserter(segments));
  
  // A vertical like fan.
  PG p3( Point_2( -50,-250), Point_2(  50,-250),50);
  PG p4( Point_2(-250, 250), Point_2( 250, 250),50);
  Segm_iterator  t2( p3, p4);
  Count_iterator t2_begin( t2);
  Count_iterator t2_end(t2, 50);
  std::copy( t2_begin, t2_end, std::back_inserter(segments));

  Q_EMIT( changed());
}
void
MainWindow::on_actionGeneratePolytopeInDisc_triggered()
{
    boost::mt19937 gen;
    gen.seed(time(0));
    std::vector<Point_2> points;
    QRectF rect = CGAL::Qt::viewportsBbox(&scene);
    CGAL::Qt::Converter<K> convert;
    Iso_rectangle_2 isor = convert(rect);
    Point_2 center = CGAL::midpoint(isor[0], isor[2]);
    Vector_2 offset = center - CGAL::ORIGIN;
    double w = isor.xmax() - isor.xmin();
    double h = isor.ymax() - isor.ymin();
    double radius = (w<h) ? w/2 : h/2;
    
    //G pg(radius);
    bool ok = false;
    const int number_of_points =
    QInputDialog::getInt(this,
                             tr("Number of random points in the disc"),
                             tr("Enter number of random points.\nThe polytope will be the convex hull of these points."),
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
    
    segments.reserve(segments.size() + 100);

    CGAL::random_convex_hull_in_disc_2(number_of_points,radius,gen,std::back_inserter(points),K());
    std::vector<Point_2>::iterator it2=points.begin();
    for(std::vector<Point_2>::iterator it=points.begin();it!=points.end();it++){
        it2++;
        if (it2==points.end()) it2=points.begin();
        Segment_2 p(*it+offset,*it2+offset);
        segments.push_back(p);
    }


    // default cursor
    QApplication::restoreOverrideCursor();
    
    Q_EMIT( changed());
}

void
MainWindow::clear()
{
  points.clear();
  segments.clear();
}


#include "Generator_2.moc"
#include <CGAL/Qt/resources.h>

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Generator_2 demo");

  // Import resources from libCGAL (Qt5).
  // See http://doc.qt.io/qt-5/qdir.html#Q_INIT_RESOURCE
  CGAL_QT_INIT_RESOURCES;
  Q_INIT_RESOURCE(Generator_2);

  MainWindow mainWindow;
  mainWindow.show();
  return app.exec();
}
