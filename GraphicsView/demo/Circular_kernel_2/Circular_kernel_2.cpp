
#include <fstream>

// CGAL headers
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Object.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>

// GraphicsView items and event filters (input classes)
#include <CGAL/Qt/GraphicsViewCircularArcInput.h>
#include "ArcsGraphicsItem.h"
  
// the two base classes
#include "ui_Circular_kernel_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

typedef CGAL::Quotient<CGAL::MP_Float>                       NT;
typedef CGAL::Cartesian<NT>                                 Linear_k;

typedef CGAL::Algebraic_kernel_for_circles_2_2<NT>          Algebraic_k;
typedef CGAL::Circular_kernel_2<Linear_k,Algebraic_k>       CircularKernel;

typedef CircularKernel::Point_2                                 Point_2;
typedef CircularKernel::Segment_2                               Segment_2;
typedef CircularKernel::Line_arc_2                              Line_arc_2;
typedef CircularKernel::Circular_arc_2                          Circular_arc_2;
typedef CircularKernel::Circular_arc_point_2                    Circular_arc_point_2;


typedef CGAL::Qt::ArcsGraphicsItem<CircularKernel>                 ArcsGraphicsItem;


typedef std::vector<CGAL::Object>                           ArcContainer;

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Circular_kernel_2
{
  Q_OBJECT
  
private:  
  ArcContainer arcs;
  ArcContainer intersections;
  QGraphicsScene scene;  

  ArcsGraphicsItem * agi;


  CGAL::Qt::GraphicsViewCircularArcInput<CircularKernel> * cai;

public:
  MainWindow();

public slots:

  virtual void open(QString);

  void processInput(CGAL::Object o);


  void on_actionInsertCircularArc_toggled(bool checked);
  
  void on_actionClear_triggered();

  void on_actionLoadLineAndCircularArcs_triggered();

  void on_actionRecenter_triggered();


signals:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow()
{
  setupUi(this);

  // Add a GraphicItem for the Circular triangulation
  agi = new CGAL::Qt::ArcsGraphicsItem<CircularKernel>(arcs, intersections);

  agi->setIntersectionsPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

  QObject::connect(this, SIGNAL(changed()),
		   agi, SLOT(modelChanged()));

  agi->setInputPen(QPen(Qt::black, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(agi);
  agi->hide();

  // Setup input handlers. They get events before the scene gets them
  // and the input they generate is passed to the triangulation with 
  // the signal/slot mechanism    
  cai = new CGAL::Qt::GraphicsViewCircularArcInput<CircularKernel>(this, &scene);

  QObject::connect(cai, SIGNAL(generate(CGAL::Object)),
		   this, SLOT(processInput(CGAL::Object)));

  // Manual handling of actions
  //
  QObject::connect(this->actionQuit, SIGNAL(triggered()), 
		   qApp, SLOT(quit()));

  // Check two actions 
  this->actionInsertCircularArc->setChecked(true);

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
  this->addAboutDemo(":/cgal/help/about_Circular_kernel_2.html");
  this->addAboutCGAL();
  this->addRecentFiles(this->menuFile, this->actionQuit);
  connect(this, SIGNAL(openRecentFile(QString)),
	  this, SLOT(open(QString)));
}


void
MainWindow::processInput(CGAL::Object o)
{
  Circular_arc_2 ca;
  Line_arc_2 la;
  bool is_circular = false;
  
  if(assign(ca, o)){
    is_circular = true;
  } else if(! assign(la, o)){
    std::cerr << "unknown object" << std::endl;
    return;
  }

  for(std::vector<CGAL::Object>::iterator it = arcs.begin(); it != arcs.end(); ++it){
    Circular_arc_2 vca;
    Line_arc_2 vla;
    if(assign(vca, *it)){
      if(is_circular){
	CGAL::intersection(ca, vca, std::back_inserter(intersections));
      } else {
	CGAL::intersection(la, vca, std::back_inserter(intersections));
      }
    } else if(assign(vla, *it)){
      if(is_circular){
	CGAL::intersection(ca, vla, std::back_inserter(intersections));
      } else {
	CGAL::intersection(la, vla, std::back_inserter(intersections));
      }
    }
  }
  arcs.push_back(o);
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
MainWindow::on_actionInsertCircularArc_toggled(bool checked)
{
  if(checked){
    scene.installEventFilter(cai);
  } else {
    scene.removeEventFilter(cai);
  }
}

void
MainWindow::on_actionClear_triggered()
{
  arcs.clear();
  intersections.clear();
  emit(changed());
}

void
MainWindow::on_actionLoadLineAndCircularArcs_triggered()
{
  QString fileName = QFileDialog::getOpenFileName(this,
						  tr("Open Line and Circular Arc File"),
						  ".",
						  tr("Edge files (*.arc)\n"));
  if(! fileName.isEmpty()){
    open(fileName);
    this->addToRecentFiles(fileName);
  }
}


void
MainWindow::open(QString fileName)
{
    std::ifstream ifs(qPrintable(fileName));
    char c;
    double x,y;
    Segment_2 s;
    
    while(ifs >> c){
      if(c == 's'){
	ifs >> x >> y;
	Point_2 p(x,y);
	ifs >> x >> y;
	Point_2 q(x,y);
	
	Line_arc_2 la(Segment_2(p,q));
	for(std::vector<CGAL::Object>::iterator it = arcs.begin(); it != arcs.end(); ++it){
	  Circular_arc_2 vca;
	  Line_arc_2 vla;
	  if(assign(vca, *it)){
	    CGAL::intersection(la, vca, std::back_inserter(intersections));
	  } else if(assign(vla, *it)){
	    CGAL::intersection(la, vla, std::back_inserter(intersections));
	  }
	}
	arcs.push_back(make_object(la));
      } else if(c == 'c'){
	ifs >> x >> y;
	Point_2 p(x,y);
	ifs >> x >> y;
	Point_2 q(x,y);
	ifs >> x >> y;
	Point_2 r(x,y);
	Circular_arc_2 ca(p,q,r);
	for(std::vector<CGAL::Object>::iterator it = arcs.begin(); it != arcs.end(); ++it){
	  Circular_arc_2 vca;
	  Line_arc_2 vla;
	  if(assign(vca, *it)){
	    CGAL::intersection(ca, vca, std::back_inserter(intersections));
	  } else if(assign(vla, *it)){
	    CGAL::intersection(ca, vla, std::back_inserter(intersections));
	  }
	}
	arcs.push_back(make_object(ca));
      }
    }
    emit (changed());
}

void
MainWindow::on_actionRecenter_triggered()
{
  this->graphicsView->setSceneRect(agi->boundingRect());
  this->graphicsView->fitInView(agi->boundingRect(), Qt::KeepAspectRatio);  
}


#include "Circular_kernel_2.moc"
#include <CGAL/Qt/resources.h>

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Circular_kernel_2 demo");

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  CGAL_QT4_INIT_RESOURCES;

  MainWindow mainWindow;
  mainWindow.show();
  return app.exec();
}
