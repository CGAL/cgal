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

#include <CGAL/Qt/PolygonWithHolesGraphicsItem.h>
#include <CGAL/Qt/BezierPolygonWithHolesGraphicsItem.h>
#include <CGAL/Qt/GraphicsViewPolygonWithHolesInput.h>
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

enum { RED, BLUE, RESULT } ;
enum { FIRST_KIND, BEZIER = FIRST_KIND , LAST_KIND } ;

QColor sPenColors  [] = { QColor(255,0,0)   , QColor(0,0,255)   , QColor(0,255,0)     } ;
QColor sBrushColors[] = { QColor(255,0,0,32), QColor(0,0,255,32), QColor(0,255,0,220) } ;

class Curve_set_base
{
public:
  
  virtual ~Curve_set_base() {}
  
  virtual CGAL::Qt::GraphicsItem* gi() const = 0 ;
  virtual CGAL::Qt::GraphicsItem* gi()       = 0 ;
  
  virtual QRectF bounding_rect() const { return gi()->boundingRect() ; }
  
  virtual bool is_empty() const = 0 ;
  
  virtual void clear               ()                               = 0 ;
  virtual void complement          ()                               = 0 ;
  virtual void assign              ( Curve_set_base const& aOther ) = 0 ;
  virtual void intersect           ( Curve_set_base const& aOther ) = 0 ;
  virtual void join                ( Curve_set_base const& aOther ) = 0 ;
  virtual void difference          ( Curve_set_base const& aOther ) = 0 ;
  virtual void symmetric_difference( Curve_set_base const& aOther ) = 0 ;
   
protected:
  
} ;

typedef boost::shared_ptr<Curve_set_base> Curve_set_base_ptr ;

typedef std::vector<Curve_set_base_ptr> Curve_set_vector ;

typedef Curve_set_vector::const_iterator Curve_set_const_iterator ;
typedef Curve_set_vector::iterator       Curve_set_iterator ;

template<class GI_, class Set_>
class Curve_set : public Curve_set_base
{
public:

  typedef GI_  GI  ;
  typedef Set_ Set ;
  
  typedef Curve_set<GI,Set> Self ;
  
  Curve_set ( int aGroup ) 
  {
    mGI = new GI(&mSet) ; 
    mGI->setPen  (QPen  ( sPenColors  [aGroup], 1, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    mGI->setBrush(QBrush( sBrushColors[aGroup] ));
  }
  
  Set const& set() const { return mSet ; }
  Set      & set()       { return mSet ; }
  
  virtual CGAL::Qt::GraphicsItem* gi() const { return mGI; }
  virtual CGAL::Qt::GraphicsItem* gi()       { return mGI; }
  
  virtual bool is_empty() const { return mSet.is_empty() ; }
  
  virtual void clear               ()                               { mSet.clear() ; }
  virtual void complement          ()                               { mSet.complement(); }
  virtual void assign              ( Curve_set_base const& aOther ) { mSet = cast(aOther).mSet; }
  virtual void intersect           ( Curve_set_base const& aOther ) { mSet.intersection        ( cast(aOther).mSet); }
  virtual void join                ( Curve_set_base const& aOther ) { mSet.join                ( cast(aOther).mSet); }
  virtual void difference          ( Curve_set_base const& aOther ) { mSet.difference          ( cast(aOther).mSet); }
  virtual void symmetric_difference( Curve_set_base const& aOther ) { mSet.symmetric_difference( cast(aOther).mSet); }
  
  static Self const& cast( Curve_set_base const& aOther ) { return dynamic_cast<Self const&>(aOther); }
  static Self      & cast( Curve_set_base      & aOther ) { return dynamic_cast<Self      &>(aOther); }
  
private:

  GI* mGI;
  Set mSet ;
} ;

class Bezier_curve_set : public Curve_set<Bezier_GI, Bezier_polygon_set>
{
  typedef Curve_set<Bezier_GI, Bezier_polygon_set> Base ;
  
public:
  
  Bezier_curve_set ( int aGroup ) : Base(aGroup) {} 
} ;

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Boolean_operations_2
{
  Q_OBJECT
  
private:  

  QGraphicsScene   mScene;
  bool             mBlue_active ;
  Curve_set_vector mCurve_sets ;
  
private:  

public:

  MainWindow();

private:
  
  void dragEnterEvent(QDragEnterEvent *event);
  void dropEvent(QDropEvent *event);
  void zoomToFit();
  
protected slots:
  
  void open( QString filename ) ;

public slots:
  
  void on_actionNew_triggered() ;
  void on_actionOpenLinear_triggered() {}
  void on_actionOpenDXF_triggered() {}
  void on_actionOpenBezier_triggered() ;
  void on_actionInsertPolygon_triggered() {}
  void on_actionInsertCircle_triggered() {}
  void on_actionInsertBezier_triggered() {}
  void on_actionIntersection_triggered() ;
  void on_actionUnion_triggered() ;
  void on_actionBlueMinusRed_triggered() ;
  void on_actionRedMinusBlue_triggered() ;
  void on_actionSymmDiff_triggered() ;
  void on_actionMinkowskiSum_triggered();
  void on_actionBlueComplement_triggered();
  void on_actionRedComplement_triggered();
  void on_actionAllBlue_triggered();
  void on_actionAllRed_triggered(); 
  void on_actionDeleteBlue_triggered();
  void on_actionDeleteRed_triggered();
  
  void on_checkboxShowBlue_toggled      (bool aChecked) { ToogleView(BLUE  ,aChecked); }
  void on_checkboxShowRed_toggled       (bool aChecked) { ToogleView(RED   ,aChecked); }
  void on_checkboxShowResult_toggled    (bool aChecked) { ToogleView(RESULT,aChecked); }
  
  void on_radioMakeBlueActive_toggled(bool aChecked) { mBlue_active =  aChecked ; }
  void on_radioMakeRedActive_toggled (bool aChecked) { mBlue_active = !aChecked ; }
  
signals:
  void changed();
  
private:
  
  Curve_set_base& set( int aKind, int aGroup ) { return *mCurve_sets[ (aKind*3) + aGroup ] ; }
  
  Curve_set_base& active_set( int aKind )    { return set(aKind, mBlue_active ? BLUE : RED) ; }
  
  Curve_set_base& result_set( int aKind )    { return set(aKind, RESULT) ; }

  void SetViewBlue  ( bool aChecked ) { checkboxShowBlue  ->setChecked(aChecked); }  
  void SetViewRed   ( bool aChecked ) { checkboxShowRed   ->setChecked(aChecked); }  
  void SetViewResult( bool aChecked ) { checkboxShowResult->setChecked(aChecked); }  

  void ToogleView( int aGROUP, bool aChecked );
};


MainWindow::MainWindow()
  : DemosMainWindow()
  , mBlue_active(true)
{
  CGAL::set_error_handler  (error_handler);
  CGAL::set_warning_handler(error_handler);
  
  setupUi(this);

  setAcceptDrops(true);

  mCurve_sets.push_back( Curve_set_base_ptr(new Bezier_curve_set(RED)    ) )  ;
  mCurve_sets.push_back( Curve_set_base_ptr(new Bezier_curve_set(BLUE)   ) ) ;
  mCurve_sets.push_back( Curve_set_base_ptr(new Bezier_curve_set(RESULT) ) ) ;
  
  for( Curve_set_iterator si = mCurve_sets.begin(); si != mCurve_sets.end() ; ++ si )
  {
    CGAL::Qt::GraphicsItem* lGI = (*si)->gi() ;
    QObject::connect(this, SIGNAL(changed()), lGI, SLOT(modelChanged()));
    mScene.addItem( lGI );
  }
  
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
  
  QObject::connect(this->actionQuit, SIGNAL(triggered()), this, SLOT(close()));
  QObject::connect(this, SIGNAL(openRecentFile(QString)), this, SLOT(open(QString)));
  
  QObject::connect(radioMakeBlueActive, SIGNAL(toggled(bool)), this, SLOT(on_radioMakeBlueActive_toggled (bool)));
  QObject::connect(radioMakeRedActive , SIGNAL(toggled(bool)), this, SLOT(on_radioMakeRedActive_toggled(bool)));
  
  QObject::connect(checkboxShowBlue   , SIGNAL(toggled(bool)), this, SLOT(on_checkboxShowBlue_toggled   (bool)));
  QObject::connect(checkboxShowRed    , SIGNAL(toggled(bool)), this, SLOT(on_checkboxShowRed_toggled    (bool)));
  QObject::connect(checkboxShowResult , SIGNAL(toggled(bool)), this, SLOT(on_checkboxShowResult_toggled (bool)));
  
	  
}

void MainWindow::on_actionNew_triggered() 
{
  for( Curve_set_iterator si = mCurve_sets.begin(); si != mCurve_sets.end() ; ++ si )
    (*si)->clear();
    
  SetViewBlue  (true);
  SetViewRed   (true);
  SetViewResult(true);
  
  radioMakeBlueActive->setChecked(true);
  
  emit(changed());
  
}

void MainWindow::dragEnterEvent(QDragEnterEvent *event)
{
  if (event->mimeData()->hasFormat("text/uri-list"))
    event->acceptProposedAction();
}

void MainWindow::dropEvent(QDropEvent *event)
{
  QString filename = event->mimeData()->urls().at(0).path();
  open(filename);
  event->acceptProposedAction();
}

bool read_Bezier_polygon_with_holes ( QString aFileName, Bezier_polygon_set& rSet )
{
  bool rOK = false ;
  
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
        Bezier_polygon_vector              polygons ;
      
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
        Bezier_polygon_vector::iterator  pit = polygons.begin();  ++pit;
        rSet.insert( Bezier_polygon_with_holes(polygons.front(), pit, polygons.end()) ) ;
        
        rOK = true ;
      }
      
    }
  }
  
  return rOK ;
}

void MainWindow::on_actionOpenBezier_triggered()
{
  open(QFileDialog::getOpenFileName(this, tr("Open Bezier Curves file"), "../data", tr("Bezier Curve files (*.dat)") ));
}

void MainWindow::on_actionIntersection_triggered() 
{
  bool lDone = false ;
  
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);
  
  for ( int k = FIRST_KIND; k != LAST_KIND ; ++ k )
  {
    if ( !set(k,BLUE).is_empty() && !set(k,RED).is_empty() )
    {
      set(k,RESULT).assign( set(k,RED) ) ;
      set(k,RESULT).intersect(set(k,BLUE));
      lDone = true ;
    }
  }
  
  this->setCursor(old);
  
  if ( lDone )
  {
    SetViewBlue  (false);
    SetViewRed   (false);
    SetViewResult(true);
    
    emit(changed());
  }
}

void MainWindow::on_actionUnion_triggered() 
{
  bool lDone = false ;
  
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);
  
  for ( int k = FIRST_KIND; k != LAST_KIND ; ++ k )
  {
    if ( !set(k,BLUE).is_empty() && !set(k,RED).is_empty() )
    {
      set(k,RESULT).assign( set(k,RED) ) ;
      set(k,RESULT).join(set(k,BLUE));
      lDone = true ;
    }
  }
  
  this->setCursor(old);
  
  if ( lDone )
  {
    SetViewBlue  (false);
    SetViewRed   (false);
    SetViewResult(true);
    
    emit(changed());
  }
}

void MainWindow::on_actionBlueMinusRed_triggered() 
{
  bool lDone = false ;
  
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);
  
  for ( int k = FIRST_KIND; k != LAST_KIND ; ++ k )
  {
    if ( !set(k,BLUE).is_empty() && !set(k,RED).is_empty() )
    {
      set(k,RESULT).assign( set(k,BLUE) ) ;
      set(k,RESULT).difference(set(k,RED));
      lDone = true ;
    }
  }
  
  this->setCursor(old);
    
  if ( lDone )
  {
    SetViewBlue  (false);
    SetViewRed   (false);
    SetViewResult(true);
    
    emit(changed());
  }
}

void MainWindow::on_actionRedMinusBlue_triggered() 
{
  bool lDone = false ;
  
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);
  
  for ( int k = FIRST_KIND; k != LAST_KIND ; ++ k )
  {
    if ( !set(k,BLUE).is_empty() && !set(k,RED).is_empty() )
    {
      set(k,RESULT).assign( set(k,RED) ) ;
      set(k,RESULT).difference(set(k,BLUE));
      lDone = true ;
    }
  }
  
  this->setCursor(old);
    
  if ( lDone )
  {
    SetViewBlue  (false);
    SetViewRed   (false);
    SetViewResult(true);
    
    emit(changed());
  }
}

void MainWindow::on_actionSymmDiff_triggered() 
{
  bool lDone = false ;
  
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);
  
  for ( int k = FIRST_KIND; k != LAST_KIND ; ++ k )
  {
    if ( !set(k,BLUE).is_empty() && !set(k,RED).is_empty() )
    {
      set(k,RESULT).assign( set(k,RED) ) ;
      set(k,RESULT).symmetric_difference(set(k,BLUE));
      lDone = true ;
    }
  }
  
  this->setCursor(old);
    
  
  if ( lDone )
  {
    SetViewBlue  (false);
    SetViewRed   (false);
    SetViewResult(true);
    
    emit(changed());
  }
}

void MainWindow::on_actionMinkowskiSum_triggered()
{
}

void MainWindow::on_actionBlueComplement_triggered()
{
  bool lDone = false ;
  
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);
  
  for ( int k = FIRST_KIND; k != LAST_KIND ; ++ k )
  {
    if ( !set(k,BLUE).is_empty() )
    {
      set(k,RESULT).assign( set(k,BLUE) ) ;
      set(k,RESULT).complement();
      lDone = true ;
    }
  }
  
  this->setCursor(old);
    
  if ( lDone )
  {
    SetViewBlue  (false);
    SetViewRed   (false);
    SetViewResult(true);
    
    emit(changed());
  }
}

void MainWindow::on_actionRedComplement_triggered()
{
  bool lDone = false ;
  
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);
  
  for ( int k = FIRST_KIND; k != LAST_KIND ; ++ k )
  {
    if ( !set(k,RED).is_empty() )
    {
      set(k,RESULT).assign( set(k,RED) ) ;
      set(k,RESULT).complement();
      lDone = true ;
    }
  }
  
  this->setCursor(old);
    
  if ( lDone )
  {
    SetViewBlue  (false);
    SetViewRed   (false);
    SetViewResult(true);
    
    emit(changed());
  }
}

void MainWindow::on_actionAllBlue_triggered()
{
  bool lDone = false ;
  
  for ( int k = FIRST_KIND; k != LAST_KIND ; ++ k )
  {
    bool lProceed = true ;
    
    if ( set(k,RESULT).is_empty() )
    {
      int answer = 0;
      answer = QMessageBox::warning(this, "Store result",
                                    QString( "Result is empty, all polygons will be deleted\n continue anyway?\n" ),
                                    "&Yes", "&No", QString::null, 1, 1 );
      lProceed = answer == 0 ;
    }
    
    if ( lProceed ) 
    {
      set(k,BLUE).assign( set(k,RESULT) ) ;
      set(k,RESULT).clear();
      radioMakeRedActive->setChecked(true);
      lDone = true ;
    }
  }
    
  if ( lDone )
  {
    SetViewBlue  (true);
    SetViewRed   (false);
    SetViewResult(true);
    
    emit(changed());
  }
}

void MainWindow::on_actionAllRed_triggered()
{
  bool lDone = false ;
  
  for ( int k = FIRST_KIND; k != LAST_KIND ; ++ k )
  {
    bool lProceed = true ;
    
    if ( set(k,RESULT).is_empty() )
    {
      int answer = 0;
      answer = QMessageBox::warning(this, "Store result",
                                    QString( "Result is empty, all polygons will be deleted\n continue anyway?\n" ),
                                    "&Yes", "&No", QString::null, 1, 1 );
      lProceed = answer == 0 ;
    }
    
    if ( lProceed ) 
    {
      set(k,RED).assign( set(k,RESULT) ) ;
      set(k,RESULT).clear();
      radioMakeBlueActive->setChecked(true);
      lDone = true ;
    }
  }
    
  if ( lDone )
  {
    SetViewBlue  (false);
    SetViewRed   (true);
    SetViewResult(true);
    
    emit(changed());
  }
}
void MainWindow::on_actionDeleteBlue_triggered()
{
  for ( int k = FIRST_KIND; k != LAST_KIND ; ++ k )
    set(k,BLUE).clear();
    
  SetViewBlue  (true);
  SetViewRed   (true);
  SetViewResult(true);
  
  emit(changed());
}

void MainWindow::on_actionDeleteRed_triggered()
{
  for ( int k = FIRST_KIND; k != LAST_KIND ; ++ k )
    set(k,RED).clear();
    
  SetViewBlue  (true);
  SetViewRed   (true);
  SetViewResult(true);
  
  emit(changed());
}

void MainWindow::open( QString fileName )
{
  if(! fileName.isEmpty())
  {
    if(fileName.endsWith(".dat"))
    {
      Bezier_curve_set& lSet = dynamic_cast<Bezier_curve_set&>(active_set(BEZIER));
      if( read_Bezier_polygon_with_holes(fileName,lSet.set()) )
      {
        emit(changed());
        zoomToFit();
        this->addToRecentFiles(fileName);
        
        if ( mBlue_active )
             radioMakeRedActive ->setChecked(true);
        else radioMakeBlueActive->setChecked(true);
        
      }
    }
  }  
}

void MainWindow::ToogleView( int aGROUP, bool aChecked )
{
  for ( int k = FIRST_KIND; k != LAST_KIND ; ++ k )
  {
    if ( aChecked )
         set(k,aGROUP).gi()->show();
    else set(k,aGROUP).gi()->hide();
  }
}


void MainWindow::zoomToFit()
{
  boost::optional<QRectF> lTotalRect ;
  
  for ( Curve_set_const_iterator si = mCurve_sets.begin() ; si != mCurve_sets.end() ; ++ si )
  {
    Curve_set_base const& lSet = **si ;
    
    if ( !lSet.is_empty() ) 
    {
      QRectF lRect = lSet.bounding_rect();
      if ( lTotalRect )
           lTotalRect = *lTotalRect | lRect ;
      else lTotalRect = lRect ;  
    }
  }
                   
  if ( lTotalRect )
  {
    this->graphicsView->setSceneRect(*lTotalRect);
    this->graphicsView->fitInView(*lTotalRect, Qt::KeepAspectRatio);  
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


