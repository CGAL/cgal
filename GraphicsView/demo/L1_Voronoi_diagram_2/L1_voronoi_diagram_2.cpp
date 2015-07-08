// Copyright (c) 2008  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:  $
// $Id:  $
// 
//
// Author(s)     : Ophir Setter <ophirset@post.tau.ac.il>
//                 

// CGAL headers
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/envelope_3.h>
#include <CGAL/L1_voronoi_traits_2.h>
#include <CGAL/point_generators_2.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>

#include <fstream>

// GraphicsView items and event filters (input classes)
#include <CGAL/Qt/ArrangementGraphicsItem.h>
#include <CGAL/Qt/ArrangementPointInput.h>
#include <CGAL/Qt/SetGraphicsItem.h>

// for viewportsBbox
#include <CGAL/Qt/utility.h>
  
// the two base classes
#include "ui_L1_voronoi_diagram_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::FT                                       Number_type;
typedef Kernel::Iso_rectangle_2                          Iso_rectangle_2;
typedef Kernel::Point_2                                  Point_2;
typedef std::vector<Point_2>                             Points;

typedef CGAL::L1_voronoi_traits_2<Kernel>                Traits_3;
typedef Traits_3::Surface_3                              Surface_3;
typedef CGAL::Envelope_diagram_2<Traits_3>               Envelope_diagram_2;

namespace CGAL {
  template <typename Kernel, typename T>
  Qt::PainterOstream<T>&
  operator<< (Qt::PainterOstream<T>& os,
              const Arr_linear_object_2<Kernel> &obj) {
    if (obj.is_segment())
      os << obj.segment();
    else if (obj.is_ray())
      os << obj.ray();
    else
      os << obj.line();
    
    return os;
  }
}

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::L1_voronoi_diagram_2
{
  Q_OBJECT
  
private:
  Points m_sites;
  Envelope_diagram_2 *m_envelope_diagram;
  QGraphicsScene m_scene;  

  CGAL::Qt::ArrangementGraphicsItem<Envelope_diagram_2> *m_graphics_item;
  CGAL::Qt::SetGraphicsItem<Points> * m_sites_graphics_item;
  CGAL::Qt::ArrangementPointInput<Envelope_diagram_2> * m_pi;


  void calculate_envelope();
  QRectF bounding_rect();

public:
  MainWindow();

public Q_SLOTS:

  void processInput(CGAL::Object o);

  void on_actionInsertPoint_toggled(bool checked);

  void on_actionInsertRandomPoints_triggered();

  void on_actionLoadPoints_triggered();

  void on_actionSavePoints_triggered();

  void on_actionClear_triggered();

  void on_actionRecenter_triggered();

  void on_actionShowVoronoi_toggled(bool checked);

  void open(QString fileName);

Q_SIGNALS:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow()
{
  setupUi(this);

  // Add a GraphicItem for the diagram
  m_envelope_diagram = new Envelope_diagram_2();
  m_graphics_item = new CGAL::Qt::ArrangementGraphicsItem<Envelope_diagram_2>
    (m_envelope_diagram);
  
  QObject::connect(this, SIGNAL(changed()),
		   m_graphics_item, SLOT(modelChanged()));

  m_graphics_item->
    setVerticesPen(QPen(Qt::red, 
                        3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  m_scene.addItem(m_graphics_item);
  
  // Add a GraphicItem for the sites
  m_sites_graphics_item = new CGAL::Qt::SetGraphicsItem<Points>
    (&m_sites);
  
  QObject::connect(this, SIGNAL(changed()),
		   m_sites_graphics_item, SLOT(modelChanged()));

  m_sites_graphics_item->
    setPen(QPen(Qt::blue, 
                5, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  m_scene.addItem(m_sites_graphics_item);

  
  // Setup input handlers. They get events before the scene gets them
  // and the input they generate is passed to the triangulation with 
  // the signal/slot mechanism
  // ophir
  m_pi = new CGAL::Qt::ArrangementPointInput<Envelope_diagram_2>(this);
  
  QObject::connect(m_pi, SIGNAL(generate(CGAL::Object)),
		   this, SLOT(processInput(CGAL::Object)));
  
  // 
  // Manual handling of actions
  //

  QObject::connect(this->actionQuit, SIGNAL(triggered()), 
		   this, SLOT(close()));

  // We put mutually exclusive actions in an QActionGroup
  QActionGroup* ag = new QActionGroup(this);
  ag->addAction(this->actionInsertPoint);

  // Check two actions 
  this->actionInsertPoint->setChecked(true);

  //
  // Setup the scene and the view
  //
  m_scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  m_scene.setSceneRect(-100, -100, 100, 100);
  this->graphicsView->setScene(&m_scene);
  this->graphicsView->setMouseTracking(true);

  // Turn the vertical axis upside down
  this->graphicsView->matrix().scale(1, -1);
                                                      
  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_L1_voronoi_diagram_2.html");
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
    m_sites.push_back(p);
    calculate_envelope();
  }
  Q_EMIT( changed());
}


/* 
 *  Qt Automatic Connections
 *  http://doc.qt.io/qt-5/designer-using-a-ui-file.html#automatic-connections
 * 
 *  setupUi(this) generates connections to the slots named
 *  "on_<action_name>_<signal_name>"
 */
void
MainWindow::on_actionInsertPoint_toggled(bool checked)
{
  if(checked){
    m_scene.installEventFilter(m_pi);
  } else {
    m_scene.removeEventFilter(m_pi);
  }
}

void
MainWindow::on_actionClear_triggered()
{
  m_sites.clear();
  calculate_envelope();
  Q_EMIT( changed());
}


void
MainWindow::on_actionInsertRandomPoints_triggered()
{
  QRectF rect = CGAL::Qt::viewportsBbox(&m_scene);
  CGAL::Qt::Converter<Kernel> convert;  
  Iso_rectangle_2 isor = convert(rect);
  CGAL::Random_points_in_iso_rectangle_2<Point_2> pg((isor.min)(), (isor.max)());
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
  for(int i = 0; i < number_of_points; ++i){
    m_sites.push_back(*pg++);
  }
  calculate_envelope();
  // default cursor
  QApplication::restoreOverrideCursor();
  Q_EMIT( changed());
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
  m_sites.clear();
  
  std::ifstream ifs(qPrintable(fileName));
  
  Kernel::Point_2 p;
  while(ifs >> p) {
    m_sites.push_back(p);
  }
  calculate_envelope();

  // default cursor
  QApplication::restoreOverrideCursor();
  this->addToRecentFiles(fileName);
  actionRecenter->trigger();
  Q_EMIT( changed());
}

void
MainWindow::on_actionSavePoints_triggered()
{
  QString fileName = QFileDialog::getSaveFileName(this,
						  tr("Save points"),
						  ".");
  if(! fileName.isEmpty()) {
    std::ofstream ofs(qPrintable(fileName));
    for(Points::iterator it = m_sites.begin();
        it != m_sites.end(); ++it)
      {
        ofs << *it << std::endl;
      }
  }
}


void
MainWindow::on_actionRecenter_triggered()
{
  this->graphicsView->setSceneRect(bounding_rect());
  this->graphicsView->fitInView(bounding_rect(), Qt::KeepAspectRatio);  
}

void
MainWindow::on_actionShowVoronoi_toggled(bool checked)
{
  this->m_graphics_item->setVisible(checked);
}


void
MainWindow::calculate_envelope() {
  if (m_envelope_diagram != NULL) {
    m_graphics_item->setArrangement(NULL);
    delete m_envelope_diagram;
    m_envelope_diagram = NULL;
  }
  
  m_envelope_diagram = new Envelope_diagram_2();
  m_graphics_item->setArrangement(m_envelope_diagram);
    
  CGAL::lower_envelope_3 (m_sites.begin(), m_sites.end(), *m_envelope_diagram);
}

QRectF
MainWindow::bounding_rect() {
  CGAL::Bbox_2 bbox(0, 0, 0, 0);

  if (m_envelope_diagram != NULL) {
    for (Envelope_diagram_2::Vertex_iterator it = 
           m_envelope_diagram->vertices_begin();
         it != m_envelope_diagram->vertices_end(); ++it) {
      double x = CGAL::to_double(it->point().x());
      double y = CGAL::to_double(it->point().y());
      CGAL::Bbox_2 temp(x, y, x, y);
      bbox = bbox + temp;
    }
  }
  
  QRectF rect(bbox.xmin(),
              bbox.ymin(), 
              bbox.xmax() - bbox.xmin(),
              bbox.ymax() - bbox.ymin());
  
  return rect;
}


#include "L1_voronoi_diagram_2.moc"
#include <CGAL/Qt/resources.h>

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("L1 Voronoi diagram_2 demo");

  // Import resources from libCGAL (Qt5).
  CGAL_QT_INIT_RESOURCES;

  MainWindow mainWindow;
  mainWindow.show();
  return app.exec();
}

