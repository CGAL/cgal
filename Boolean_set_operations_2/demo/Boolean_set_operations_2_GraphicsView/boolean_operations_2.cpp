// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
// 
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
// 
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Boolean_set_operations_2/demo/Boolean_set_operations_2/boolean_operations_2.cpp $
// $Id: boolean_operations_2.cpp 45454 2008-09-09 21:42:42Z lrineau $
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>



#include <CGAL/basic.h>
 

#include <fstream>
#include <string>

#include <CGAL/Bbox_2.h>
#include <CGAL/iterator.h>
 
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
//#include <CGAL/Qt/GraphicsViewPolylineInput.h>
//#include <CGAL/Qt/Polyline_simplification_2_graphics_item.h>
#include <CGAL/Qt/Converter.h>

// the two base classes
#include "ui_boolean_operations_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

// for viewportsBbox(QGraphicsScene*)
#include <CGAL/Qt/utility.h>

#include <CGAL/Bbox_2.h>
#include <CGAL/assertions_behaviour.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h> 

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

#include "typedefs.h"

#include <CGAL/IO/Dxf_bsop_reader.h>
//#include <qdeepcopy.h>

//global variable to aid naming windows 
int winsOpened=2;

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Boolean_operations_2
{
  Q_OBJECT
  
private:  

  QGraphicsScene                                    mScene;
  
private:  

  
public:

  MainWindow();

protected:
  void dragEnterEvent(QDragEnterEvent *event);
  void dropEvent(QDropEvent *event);
private:

public slots:
  
signals:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow()
{
  CGAL::set_error_handler  (error_handler);
  CGAL::set_warning_handler(error_handler);
  
  setupUi(this);

  setAcceptDrops(true);

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
  this->addAboutDemo(":/cgal/help/index.html");
  this->addAboutCGAL();

  this->addRecentFiles(this->menuFile, this->actionQuit);
  //connect(this, SIGNAL(openRecentFile(QString)), this, SLOT(open(QString)));
	  
}

void MainWindow::dragEnterEvent(QDragEnterEvent *event)
{
  if (event->mimeData()->hasFormat("text/uri-list"))
    event->acceptProposedAction();
}

void MainWindow::dropEvent(QDropEvent *event)
{
  QString filename = event->mimeData()->urls().at(0).path();
//  open(filename);
  event->acceptProposedAction();
}

#include "boolean_operations_2.moc"

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Boolean_operations_2 demo");

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  Q_INIT_RESOURCE(File);
  Q_INIT_RESOURCE(Input);
  Q_INIT_RESOURCE(CGAL);
  
  MainWindow mainWindow;
  mainWindow.show();
  return app.exec();
}


