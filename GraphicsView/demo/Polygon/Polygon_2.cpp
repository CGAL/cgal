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

// Qt headers
#include <QtGui>
#include <QString>
#include <QFileDialog>
#include <QGraphicsLineItem>

// GraphicsView items and event filters (input classes)
#include <CGAL/Qt/GraphicsViewPolylineInput.h>
#include <CGAL/Qt/PolygonGraphicsItem.h>
#include <CGAL/Qt/LineGraphicsItem.h>
  
// the two base classes
#include "ui_Polygon_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef K::Line_2 Line_2;

typedef CGAL::Polygon_2<K,std::list<Point_2> > Polygon; // it must be a list for the partition

typedef CGAL::Straight_skeleton_2<K> Ss ;

typedef boost::shared_ptr<Ss> SsPtr ;

typedef boost::shared_ptr<Polygon> PolygonPtr ;

typedef std::vector<PolygonPtr> PolygonPtr_vector ;

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Polygon_2
{
  Q_OBJECT
  
private:

  enum PartitionAlgorithm {YMonotone, ApproximateConvex, OptimalConvex} ;


  CGAL::Qt::Converter<K> convert;
  Polygon poly; 
  QGraphicsScene scene;  

  CGAL::Qt::PolygonGraphicsItem<Polygon> * pgi;

  CGAL::Qt::GraphicsViewPolylineInput<K> * pi;

  std::list<Polygon> partitionPolygons;
  std::list<CGAL::Qt::PolygonGraphicsItem<Polygon>* >  partitionGraphicsItems;
  std::list<QGraphicsLineItem* >  skeletonGraphicsItems;
  std::list<QGraphicsLineItem* >  offsetGraphicsItems;
  CGAL::Qt::LineGraphicsItem<K>* lgi;

public:
  MainWindow();

public slots:

  void processInput(CGAL::Object o);

  void on_actionClear_triggered();

  void on_actionLoadPolygon_triggered();
  void on_actionSavePolygon_triggered();

  void on_actionRecenter_triggered();
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
  void clearSkeleton();
  void clearOffset();
  void clear();

  void open(const QString&);
signals:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow()
{
  setupUi(this);

  // Add a GraphicItem for the Polygon_2
  pgi = new CGAL::Qt::PolygonGraphicsItem<Polygon>(&poly);

  QObject::connect(this, SIGNAL(changed()),
		   pgi, SLOT(modelChanged()));

  pgi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(pgi);


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
    emit(changed());
  }
}

/* 
 *  Qt Automatic Connections
 *  http://doc.trolltech.com/4.4/designer-using-a-component.html#automatic-connections
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
  emit(changed());
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
MainWindow::open(const QString& fileName)
{
  this->actionCreateInputPolygon->setChecked(false);
  std::ifstream ifs(qPrintable(fileName));
  poly.clear();
  ifs >> poly;
  clear();

  this->addToRecentFiles(fileName);
  emit (changed());
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
  emit(changed());
}

void
MainWindow::on_actionRecenter_triggered()
{
  this->graphicsView->setSceneRect(pgi->boundingRect());
  this->graphicsView->fitInView(pgi->boundingRect(), Qt::KeepAspectRatio);  
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
                                              , def_offset
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
          Polygon const& lContour = **cit ;
          lContour.area();
          for ( Polygon::Edge_const_iterator eit = lContour.edges_begin(); eit != lContour.edges_end(); ++ eit )
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
    for(std::list<Polygon>::iterator it = partitionPolygons.begin();
	it != partitionPolygons.end();
	++it){
      partitionGraphicsItems.push_back(new CGAL::Qt::PolygonGraphicsItem<Polygon>(&(*it)));
      scene.addItem(partitionGraphicsItems.back());
      partitionGraphicsItems.back()->setEdgesPen(QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    }
  }
}

void
MainWindow::clearPartition()
{
  partitionPolygons.clear();
  for(std::list<CGAL::Qt::PolygonGraphicsItem<Polygon>* >::iterator it = partitionGraphicsItems.begin();
      it != partitionGraphicsItems.end();
      ++it){
    scene.removeItem(*it);
  }
  partitionGraphicsItems.clear();
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
  clearSkeleton();
  clearOffset();
  lgi->hide();
}


#include "Polygon_2.moc"

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Polygon_2 demo");

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  Q_INIT_RESOURCE(File);
  Q_INIT_RESOURCE(Polygon_2);
  Q_INIT_RESOURCE(Input);
  Q_INIT_RESOURCE(CGAL);

  MainWindow mainWindow;
  mainWindow.show();
  return app.exec();
}
