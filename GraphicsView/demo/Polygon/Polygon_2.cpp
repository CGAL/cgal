#include <fstream>
#include<boost/shared_ptr.hpp>
// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/Partition_traits_2.h>
#include<CGAL/create_straight_skeleton_2.h>
#include<CGAL/create_offset_polygons_2.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/extremal_polygon_2.h>
#include <CGAL/minkowski_sum_2.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QFileDialog>
#include <QInputDialog>
#include <QGraphicsLineItem>

// GraphicsView items and event filters (input classes)
#include <CGAL/Qt/GraphicsViewPolylineInput.h>
#include <CGAL/Qt/PolygonGraphicsItem.h>
#include <CGAL/Qt/PolygonWithHolesGraphicsItem.h>
#include <CGAL/Qt/LineGraphicsItem.h>
  
// the two base classes
#include "ui_Polygon_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef K::Line_2 Line_2;

typedef CGAL::Polygon_2<K,std::list< Point_2 > > Polygon2; // it must be a list for the partition
typedef CGAL::Polygon_with_holes_2<K,std::list< Point_2 > > Polygon_with_holes_2;

typedef CGAL::Straight_skeleton_2<K> Ss ;

typedef boost::shared_ptr<Ss> SsPtr ;

typedef boost::shared_ptr<Polygon2> PolygonPtr ;

typedef std::vector<PolygonPtr> PolygonPtr_vector ;

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Polygon_2
{
  Q_OBJECT
  
private:

  enum PartitionAlgorithm {YMonotone, ApproximateConvex, OptimalConvex} ;


  CGAL::Qt::Converter<K> convert;
  Polygon2 poly, kgon; 
  Polygon_with_holes_2 selfmink;
  QGraphicsScene scene;  

  CGAL::Qt::PolygonGraphicsItem<Polygon2> * pgi;

  CGAL::Qt::GraphicsViewPolylineInput<K> * pi;

  CGAL::Qt::PolygonWithHolesGraphicsItem<Polygon_with_holes_2> * minkgi;

  std::list<Polygon2> partitionPolygons;
  std::list<CGAL::Qt::PolygonGraphicsItem<Polygon2>* >  partitionGraphicsItems;
  std::list<QGraphicsLineItem* >  skeletonGraphicsItems;
  std::list<QGraphicsLineItem* >  offsetGraphicsItems;
  CGAL::Qt::LineGraphicsItem<K>* lgi;

  CGAL::Qt::PolygonGraphicsItem<Polygon2> * kgongi;

public:
  MainWindow();

public Q_SLOTS:

  void processInput(CGAL::Object o);

  void on_actionClear_triggered();

  void on_actionLoadPolygon_triggered();
  void on_actionSavePolygon_triggered();

  void on_actionSelfMinkowskiSum_triggered();
  void on_actionRecenter_triggered();
  void on_actionMaximumAreaKGon_triggered();
  void on_actionInnerSkeleton_triggered();
  void on_actionOuterOffset_triggered();
  void on_actionLinearLeastSquaresFitting_triggered();
  void on_actionLinearLeastSquaresFittingOfSegments_triggered();
  void on_actionCreateInputPolygon_toggled(bool);

  void on_actionYMonotonePartition_triggered();
  void on_actionApproximateConvexPartition_triggered();
  void on_actionOptimalConvexPartition_triggered();
  void partition(PartitionAlgorithm);

  void clearPartition();
  void clearMinkowski();
  void clearSkeleton();
  void clearOffset();
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

  minkgi = 0;
  // Add a GraphicItem for the Polygon2
  pgi = new CGAL::Qt::PolygonGraphicsItem<Polygon2>(&poly);

  QObject::connect(this, SIGNAL(changed()),
		   pgi, SLOT(modelChanged()));

  scene.addItem(pgi);

  kgongi =  new CGAL::Qt::PolygonGraphicsItem<Polygon2>(&kgon);
  kgongi->setEdgesPen(QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  kgongi->hide();
  scene.addItem(kgongi);
  
  
  lgi = new CGAL::Qt::LineGraphicsItem<K>();
  lgi->setPen(QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  lgi->hide();
  scene.addItem(lgi);
  assert(lgi->scene() == &scene);
  // Setup input handlers. They get events before the scene gets them
  pi = new CGAL::Qt::GraphicsViewPolylineInput<K>(this, &scene, 0, true);

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

  // Uncomment the following line to get antialiasing by default.
//   actionUse_Antialiasing->setChecked(true);

  // Turn the vertical axis upside down
  this->graphicsView->scale(1, -1);
                                                      
  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_Polygon_2.html");
  this->addAboutCGAL();
  this->setupExportSVG(action_Export_SVG, graphicsView);

  this->addRecentFiles(this->menuFile, this->actionQuit);
  connect(this, SIGNAL(openRecentFile(QString)),
	  this, SLOT(open(QString)));
}


void
MainWindow::processInput(CGAL::Object o)
{
  this->actionCreateInputPolygon->setChecked(false);
  std::list<Point_2> points;
  if(CGAL::assign(points, o)){
    if((points.size() == 1)&& poly.size()>0){
    
    } else {
      poly.clear();
      if(points.front() == points.back()){
	points.pop_back();
      }
      poly.insert(poly.vertices_begin(), points.begin(), points.end());
    }
    Q_EMIT( changed());
  }
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
  poly.clear();
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
						  tr( "Any file (*.*)"));
  if(! fileName.isEmpty()){
    open(fileName);
  }
}

void
MainWindow::open(QString fileName)
{
  this->actionCreateInputPolygon->setChecked(false);
  std::ifstream ifs(qPrintable(fileName));
  poly.clear();
  ifs >> poly;
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
						  tr("Any files (*.*)"));
  if(! fileName.isEmpty()){
    std::ofstream ofs(qPrintable(fileName));
    ofs << poly;
  }
}


void
MainWindow::on_actionCreateInputPolygon_toggled(bool checked)
{
  poly.clear();
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
  this->graphicsView->setSceneRect(pgi->boundingRect());
  this->graphicsView->fitInView(pgi->boundingRect(), Qt::KeepAspectRatio);  
}

void
MainWindow::on_actionSelfMinkowskiSum_triggered()
{
  if(poly.size()>0){
    if(! poly.is_simple()){
      return;
    }

    selfmink = minkowski_sum_2 (poly, poly);
    
    minkgi = new CGAL::Qt::PolygonWithHolesGraphicsItem<Polygon_with_holes_2>(&selfmink);
    
    minkgi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    scene.addItem(minkgi);
    minkgi->setEdgesPen(QPen(Qt::darkRed, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  }
}




void
MainWindow::on_actionInnerSkeleton_triggered()
{
  if(poly.size()>0){
    if(! poly.is_simple()){
      return;
    }
    clear();
    if(! poly.is_counterclockwise_oriented()){
      poly.reverse_orientation();
    }
    SsPtr iss = CGAL::create_interior_straight_skeleton_2(poly.vertices_begin(), poly.vertices_end());

    CGAL::Straight_skeleton_2<K> const& ss = *iss;

    typedef Ss::Vertex_const_handle     Vertex_const_handle ;
    typedef Ss::Halfedge_const_handle   Halfedge_const_handle ;
    typedef Ss::Halfedge_const_iterator Halfedge_const_iterator ;
  
    Halfedge_const_handle null_halfedge ;
    Vertex_const_handle   null_vertex ;

    for ( Halfedge_const_iterator i = ss.halfedges_begin(); i != ss.halfedges_end(); ++i )
      {
	if ( i->is_bisector() ){
	  Segment_2 s(i->opposite()->vertex()->point(), i->vertex()->point());
	  skeletonGraphicsItems.push_back(new QGraphicsLineItem(convert(s)));
	  scene.addItem(skeletonGraphicsItems.back());
	  skeletonGraphicsItems.back()->setPen(QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
	}
      }
      
  }

}

void
MainWindow::on_actionOuterOffset_triggered()
{
  if(poly.size()>0 && poly.is_simple())
  {
    clear();
    
    if(! poly.is_counterclockwise_oriented())
      poly.reverse_orientation();
    
    double w = poly.bbox().xmax() - poly.bbox().xmin() ;
    double h = poly.bbox().ymax() - poly.bbox().ymin() ;
    double s = (std::max)(w,h);
    double def_offset = s * 0.05 ;
    
    bool ok;
    const double off = QInputDialog::getDouble( this
                                              , tr("Offset distance")
                                              , tr("Offset distance:")
                                              , def_offset
                                              , 0.0
                                              , s
                                              , 2
                                              , &ok
                                              );  
    if ( ok )
    {
      PolygonPtr_vector lContours = create_exterior_skeleton_and_offset_polygons_2(off,poly) ;
      
      PolygonPtr_vector::const_iterator frame ;
      
      double max_area = 0.0 ;
      
      for( PolygonPtr_vector::const_iterator cit = lContours.begin(); cit != lContours.end(); ++ cit )
      {
        double area = (*cit)->area();
        if ( area > max_area )
        {
          max_area = area ;
          frame = cit ;
        }
      }
      
      for( PolygonPtr_vector::const_iterator cit = lContours.begin(); cit != lContours.end(); ++ cit )
      {
        if ( cit != frame )
        {
          Polygon2 const& lContour = **cit ;
          lContour.area();
          for ( Polygon2::Edge_const_iterator eit = lContour.edges_begin(); eit != lContour.edges_end(); ++ eit )
          {
            Segment_2 s(eit->source(), eit->target());
            offsetGraphicsItems.push_back(new QGraphicsLineItem(convert(s)));
            scene.addItem(offsetGraphicsItems.back());
            offsetGraphicsItems.back()->setPen(QPen(Qt::red, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
          } 
        }
      }
    }
  }  
}

void
MainWindow::on_actionMaximumAreaKGon_triggered()
{
  if( (poly.size()>2) && poly.is_convex()){
    clear();
    
    kgon.clear();
    std::vector<Point_2> points(poly.vertices_begin(),
                                poly.vertices_end());
    CGAL::maximum_area_inscribed_k_gon_2(points.begin(),
                                         points.end(),
					 3,
                                         std::back_inserter(kgon));
  
    kgongi->modelChanged();
    kgongi->show();
  } else {
    std::cout << "The polygon must be convex" << std::endl;
  }
}

void
MainWindow::on_actionLinearLeastSquaresFitting_triggered()
{
  if(poly.size()>2){
    clear();
    
    Line_2 line;
    CGAL::linear_least_squares_fitting_2(poly.vertices_begin(),
					 poly.vertices_end(),
					 line,
					 CGAL::Dimension_tag<0>());
  
    lgi->setLine(line);
    lgi->show();
  }
}

void
MainWindow::on_actionLinearLeastSquaresFittingOfSegments_triggered()
{
  if(poly.size()>2){
    clear();
    
    Line_2 line;
    CGAL::linear_least_squares_fitting_2(poly.edges_begin(),
					 poly.edges_end(),
					 line,
					 CGAL::Dimension_tag<1>());


    lgi->setLine(line);
    lgi->show();
  }
}

void
MainWindow::on_actionYMonotonePartition_triggered()
{
  partition(YMonotone);
}


void
MainWindow::on_actionOptimalConvexPartition_triggered()
{
  partition(OptimalConvex);
}


void
MainWindow::on_actionApproximateConvexPartition_triggered()
{
  partition(ApproximateConvex);
}


void 
MainWindow::partition(PartitionAlgorithm pa)
{
  if(poly.size()>0){
    clear();
    if(! poly.is_counterclockwise_oriented()){
      poly.reverse_orientation();
    }
    switch (pa) {
    case YMonotone :
      CGAL::y_monotone_partition_2(poly.vertices_begin(), poly.vertices_end(), std::back_inserter(partitionPolygons));
      break;
    case ApproximateConvex:
      CGAL::approx_convex_partition_2(poly.vertices_begin(), poly.vertices_end(), std::back_inserter(partitionPolygons));
      break;
    default:
      CGAL::optimal_convex_partition_2(poly.vertices_begin(), poly.vertices_end(), std::back_inserter(partitionPolygons));
      break;
    }
    for(std::list<Polygon2>::iterator it = partitionPolygons.begin();
	it != partitionPolygons.end();
	++it){
      partitionGraphicsItems.push_back(new CGAL::Qt::PolygonGraphicsItem<Polygon2>(&(*it)));
      scene.addItem(partitionGraphicsItems.back());
      partitionGraphicsItems.back()->setEdgesPen(QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    }
  }
}

void
MainWindow::clearPartition()
{
  partitionPolygons.clear();
  for(std::list<CGAL::Qt::PolygonGraphicsItem<Polygon2>* >::iterator it = partitionGraphicsItems.begin();
      it != partitionGraphicsItems.end();
      ++it){
    scene.removeItem(*it);
  }
  partitionGraphicsItems.clear();
}

void
MainWindow::clearMinkowski()
{
  if(minkgi != 0){
    scene.removeItem(minkgi);
    delete minkgi;
    minkgi = 0;
  }
}

void
MainWindow::clearSkeleton()
{ for(std::list<QGraphicsLineItem* >::iterator it = skeletonGraphicsItems.begin();
      it != skeletonGraphicsItems.end();
      ++it){
    scene.removeItem(*it);
  }
  skeletonGraphicsItems.clear();
}

void
MainWindow::clearOffset()
{ for(std::list<QGraphicsLineItem* >::iterator it = offsetGraphicsItems.begin();
      it != offsetGraphicsItems.end();
      ++it){
    scene.removeItem(*it);
  }
  offsetGraphicsItems.clear();
}

void
MainWindow::clear()
{
  clearPartition();
  clearMinkowski();
  clearSkeleton();
  clearOffset();
  lgi->hide();
  kgongi->hide();
}


#include "Polygon_2.moc"
#include <CGAL/Qt/resources.h>

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Polygon_2 demo");

  // Import resources from libCGAL (Qt5).
  // See http://doc.qt.io/qt-5/qdir.html#Q_INIT_RESOURCE
  CGAL_QT_INIT_RESOURCES;
  Q_INIT_RESOURCE(Polygon_2);

  MainWindow mainWindow;
  mainWindow.show();
  return app.exec();
}
