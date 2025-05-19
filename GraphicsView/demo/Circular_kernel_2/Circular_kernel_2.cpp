#include <boost/config.hpp>
#include <boost/version.hpp>
#include <fstream>

// CGAL headers
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/IO/WKT.h>

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
typedef std::variant<Circular_arc_2, Line_arc_2 >               Arc_variant;
typedef std::variant<std::pair<Circular_arc_point_2, unsigned>,
                     Circular_arc_2, Line_arc_2 >               Inter_variant;

typedef CGAL::Qt::ArcsGraphicsItem<CircularKernel>                 ArcsGraphicsItem;


typedef std::vector<Arc_variant>                                ArcContainer;
typedef std::vector<Inter_variant>                              ArcIntersection;

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Circular_kernel_2
{
  Q_OBJECT

private:
  ArcContainer arcs;
  ArcIntersection intersections;
  QGraphicsScene scene;

  ArcsGraphicsItem * agi;


  CGAL::Qt::GraphicsViewCircularArcInput<CircularKernel> * cai;

public:
  MainWindow();

public Q_SLOTS:

  virtual void open(QString);

  void processInput(Arc_variant o);


  void on_actionInsertCircularArc_toggled(bool checked);

  void on_actionClear_triggered();

  void on_actionLoadLineAndCircularArcs_triggered();

  void on_actionRecenter_triggered();


Q_SIGNALS:
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

  QObject::connect(cai, SIGNAL(generate(Arc_variant)),
                   this, SLOT(processInput(Arc_variant)));

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
MainWindow::processInput(Arc_variant o)
{
  const Circular_arc_2* ca = nullptr;
  const Line_arc_2* la = nullptr;
  bool is_circular = false;

  if( (ca = std::get_if<Circular_arc_2>(&o)) ){
    is_circular = true;
  } else if(! (la = std::get_if<Line_arc_2>(&o))){
    std::cerr << "unknown object" << std::endl;
    return;
  }

  for(std::vector<Arc_variant>::iterator it = arcs.begin(); it != arcs.end(); ++it){
    Circular_arc_2 vca;
    Line_arc_2 vla;
    if(auto vca = std::get_if<Circular_arc_2>(&(*it))){
      if(is_circular){
        CGAL::intersection(*ca, *vca, std::back_inserter(intersections));
      } else {
        CGAL::intersection(*la, *vca, std::back_inserter(intersections));
      }
    } else if(auto vla = std::get_if<Line_arc_2>(&(*it))){
      if(is_circular){
        CGAL::intersection(*ca, *vla, std::back_inserter(intersections));
      } else {
        CGAL::intersection(*la, *vla, std::back_inserter(intersections));
      }
    }
  }
  arcs.push_back(o);
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
  Q_EMIT( changed());
}

void
MainWindow::on_actionLoadLineAndCircularArcs_triggered()
{
  QString fileName = QFileDialog::getOpenFileName(this,
                                                  tr("Open Line and Circular Arc File"),
                                                  ".",
                                                  tr("Edge files (*.arc)\n"
                                                     "WKT files (*.wkt *.WKT)\n"
                                                     ));
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
    if(fileName.endsWith(".wkt", Qt::CaseInsensitive))
    {
      //read pairs as Line_arc_2 and triplets as circular_arc_2
      do
      {
        std::vector<Point_2> multi_points;
        CGAL::IO::read_multi_point_WKT(ifs, multi_points);
        if(multi_points.size() == 2)
        {
          Line_arc_2 la(Segment_2(multi_points[0],
                        multi_points[1]));
          for(std::vector<Arc_variant>::iterator it = arcs.begin(); it != arcs.end(); ++it){
            if(auto vca = std::get_if<Circular_arc_2>(&(*it))){
              CGAL::intersection(la, *vca, std::back_inserter(intersections));
            } else if(auto vla = std::get_if<Line_arc_2>(&(*it))){
              CGAL::intersection(la, *vla, std::back_inserter(intersections));
            }
          }
          arcs.push_back(la);
        }
        else if(multi_points.size() == 3)
        {
          Circular_arc_2 ca(multi_points[0],
                            multi_points[1],
                            multi_points[2]);
          for(std::vector<Arc_variant>::iterator it = arcs.begin(); it != arcs.end(); ++it){
            Circular_arc_2 vca;
            Line_arc_2 vla;
            if(auto vca = std::get_if<Circular_arc_2>(&(*it))){
              CGAL::intersection(ca, *vca, std::back_inserter(intersections));
            } else if(auto vla = std::get_if<Line_arc_2>(&(*it))){
              CGAL::intersection(ca, *vla, std::back_inserter(intersections));
            }
          }
          arcs.push_back(ca);
        }
        else if(multi_points.size()>0)
        {
          std::cerr<<"unreadable object."<<std::endl;
        }
      }while(ifs.good() && !ifs.eof());
      ifs.close();
    }
    else
    {
      while(ifs >> c){
        if(c == 's'){
          ifs >> x >> y;
          Point_2 p(x,y);
          ifs >> x >> y;
          Point_2 q(x,y);

          Line_arc_2 la(Segment_2(p,q));
          for(std::vector<Arc_variant>::iterator it = arcs.begin(); it != arcs.end(); ++it){
            if(auto vca = std::get_if<Circular_arc_2>(&(*it))){
              CGAL::intersection(la, *vca, std::back_inserter(intersections));
            } else if(auto vla = std::get_if<Line_arc_2>(&(*it))){
              CGAL::intersection(la, *vla, std::back_inserter(intersections));
            }
          }
          arcs.push_back(la);
        } else if(c == 'c'){
          ifs >> x >> y;
          Point_2 p(x,y);
          ifs >> x >> y;
          Point_2 q(x,y);
          ifs >> x >> y;
          Point_2 r(x,y);
          Circular_arc_2 ca(p,q,r);
          for(std::vector<Arc_variant>::iterator it = arcs.begin(); it != arcs.end(); ++it){
            Circular_arc_2 vca;
            Line_arc_2 vla;
            if(auto vca = std::get_if<Circular_arc_2>(&(*it))){
              CGAL::intersection(ca, *vca, std::back_inserter(intersections));
            } else if(auto vla = std::get_if<Line_arc_2>(&(*it))){
              CGAL::intersection(ca, *vla, std::back_inserter(intersections));
            }
          }
          arcs.push_back(ca);
        }
      }
    }
    Q_EMIT( changed());
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

  // Import resources from libCGAL (Qt6).
  CGAL_QT_INIT_RESOURCES;

  MainWindow mainWindow;
  mainWindow.show();
  return app.exec();
}
