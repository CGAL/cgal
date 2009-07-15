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


#include <fstream>
#include <string>
#include <list>

#include <boost/shared_ptr.hpp>

#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QSlider>
#include <QProgressBar>

#include <CGAL/basic.h>
#include <CGAL/Timer.h> 
#include <CGAL/Bbox_2.h>
#include <CGAL/iterator.h>
#include <CGAL/assertions_behaviour.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/approximated_offset_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Boolean_set_operations_2.h>

#ifdef CGAL_USE_GMP
  #include <CGAL/Gmpq.h>
#else
  #include <CGAL/MP_Float.h>
  #include <CGAL/Quotient.h>
#endif

#include <CGAL/Qt/BezierPolygonWithHolesGraphicsItem.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/DemosMainWindow.h>
#include <CGAL/Qt/utility.h>
#include <CGAL/IO/Dxf_bsop_reader.h>

// the two base classes
#include "ui_boolean_operations_2.h"

#include "typedefs.h"

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



//global variable to aid naming windows 
int winsOpened=2;

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Boolean_operations_2
{
  Q_OBJECT
  
private:  

  QGraphicsScene                  mScene;
  
  bool                            mBezier_red_active ;
  Bezier_GI*                      mBezier_red_GI;
  Bezier_GI*                      mBezier_blue_GI;
  Bezier_GI*                      mBezier_result_GI;
  Bezier_polygon_with_holes_list  mBezier_red_set ;
  Bezier_polygon_with_holes_list  mBezier_blue_set ;
  Bezier_polygon_with_holes_list  mBezier_result_set ;
  
private:  

  
public:

  MainWindow();

private:
  
  void dragEnterEvent(QDragEnterEvent *event);
  void dropEvent(QDropEvent *event);

public slots:
  
  void on_actionNew_triggered() {}
  void on_actionNewWindow_triggered() {}
  void on_actionOpenLinear_triggered() {}
  void on_actionOpenDXF_triggered() {}
  void on_actionOpenBezier_triggered() ;
  void on_actionPrint_triggered() {}
  void on_actionClose_triggered() {}
  void on_actionQuit_triggered() {}
  void on_actionHowTo_triggered() {}
  void on_actionAbout_triggered() {}
  void on_actionAboutCGAL_triggered() {}
  void on_actionUndo_triggered() {}
  void on_actionRedo_triggered() {}
  void on_actionMakeBlueActive_toggled(bool checked) {}
  void on_actionMakeRedActive_toggled(bool checked) {}
//   void on_actionInsertPolygon_triggered() {}
  void on_actionInsertCircle_triggered() {}
  void on_actionInsertBezier_triggered() {}
  void on_actionIntersection_triggered() {}
  void on_actionUnion_triggered() {}
  void on_actionBlueMinusRed_triggered() {}
  void on_actionRedMinusBlue_triggered() {}
  void on_actionSymmDiff_triggered() {}
  void on_actionMinkowskiSum_triggered() {}
  void on_actionBlueComplement_triggered() {}
  void on_actionRedComplement_triggered() {}
  void on_actionAllBlue_triggered() {}
  void on_actionAllRed_triggered() {}
  void on_actionDeleteBlue_triggered() {}
  void on_actionDeleteRed_triggered() {}
  void on_actionRefresh_triggered() {}
  
signals:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow()
  , mBezier_red_active(false)
{
  CGAL::set_error_handler  (error_handler);
  CGAL::set_warning_handler(error_handler);
  
  setupUi(this);

  setAcceptDrops(true);

  mBezier_red_GI    = new Bezier_GI(&mBezier_red_set);
  mBezier_blue_GI   = new Bezier_GI(&mBezier_blue_set);
  mBezier_result_GI = new Bezier_GI(&mBezier_result_set);
    
  QObject::connect(this, SIGNAL(changed()), mBezier_red_GI   , SLOT(modelChanged()));
  QObject::connect(this, SIGNAL(changed()), mBezier_blue_GI  , SLOT(modelChanged()));
  QObject::connect(this, SIGNAL(changed()), mBezier_result_GI, SLOT(modelChanged()));

  mBezier_red_GI   ->setPen(QPen(Qt::red   , 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  mBezier_blue_GI  ->setPen(QPen(Qt::blue  , 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  mBezier_result_GI->setPen(QPen(Qt::yellow, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  
  mScene.addItem(mBezier_red_GI);
  mScene.addItem(mBezier_blue_GI);
  mScene.addItem(mBezier_result_GI);
  
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

std::size_t read_Bezier_polygon_with_holes ( QString aFileName, Bezier_polygon_with_holes_list& rList )
{
  if( ! aFileName.isEmpty() && aFileName.endsWith(".dat") )
  {
    std::ifstream in_file (qPrintable(aFileName));
  
    if ( in_file )
    {
      // Red the number of bezier polygon with holes
      unsigned int n_regions ;
      in_file >> n_regions;
      
      for ( unsigned int r = 0 ; r < n_regions ; ++ r )
      {
        // Read the number of bezier curves.
        unsigned int n_curves;
        in_file >> n_curves;
      
        // Read the curves one by one, and construct the general polygon these
        // curve form (the outer boundary and the holes inside it).
        Bezier_traits_2                    traits;
        Bezier_traits_2::Make_x_monotone_2 make_x_monotone = traits.make_x_monotone_2_object();
        
        bool                               first = true;
        Bezier_rat_point                   p_0;
        std::list<Bezier_X_monotone_curve> xcvs;
        Bezier_rat_kernel                  ker;
        Bezier_rat_kernel::Equal_2         equal = ker.equal_2_object();
        Bezier_polygon_list                polygons ;
      
        for ( unsigned int k = 0; k < n_curves; ++ k ) 
        {
          // Read the current curve and subdivide it into x-monotone subcurves.
          
          Bezier_curve                            B;
          std::list<CGAL::Object>                 x_objs;
          std::list<CGAL::Object>::const_iterator xoit;
          Bezier_X_monotone_curve                 xcv;
      
          in_file >> B;
          make_x_monotone (B, std::back_inserter (x_objs));
          
          for (xoit = x_objs.begin(); xoit != x_objs.end(); ++xoit) 
          {
            if (CGAL::assign (xcv, *xoit))
              xcvs.push_back (xcv);
          }
          
          // Check if the current curve closes a polygon, namely whether it target
          // point (the last control point) equals the source of the first curve in
          // the current chain.
          if (! first) 
          {
            if (equal (p_0, B.control_point(B.number_of_control_points() - 1))) 
            {
              // Push a new polygon into the polygon list. Make sure that the polygon
              // is counterclockwise oriented if it represents the outer boundary
              // and clockwise oriented if it represents a hole.
              Bezier_polygon  pgn (xcvs.begin(), xcvs.end());
              CGAL::Orientation  orient = pgn.orientation();
              
              if ((polygons.empty() && orient == CGAL::CLOCKWISE) || (! polygons.empty() && orient == CGAL::COUNTERCLOCKWISE))
                pgn.reverse_orientation();
              
              polygons.push_back (pgn);
              xcvs.clear();
              first = true;
            }
          }
          else 
          {
            // This is the first curve in the chain - store its source point.
            p_0 = B.control_point(0);
            first = false;
          }
        }
      
        // Construct the polygon with holes.
        std::list<Bezier_polygon>::iterator  pit = polygons.begin();  ++pit;
        rList.push_back( Bezier_polygon_with_holes(polygons.front(), pit, polygons.end()) ) ;
      }
      
    }
  }
  
  return rList.size();
}

void MainWindow::on_actionOpenBezier_triggered()
{
  QString fileName = QFileDialog::getOpenFileName(this, tr("Open Bezier Curves file"), "../data", tr("Bezier Curve files (*.dat)") );
  
  if( read_Bezier_polygon_with_holes(fileName, mBezier_red_active ? mBezier_red_set : mBezier_blue_set ) > 0 )
  {
    emit(changed());
    this->addToRecentFiles(fileName);
  }
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


