// Copyright (c) 2012  Tel-Aviv University (Israel).
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
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Apurva Bhatt <response2apurva@gmail.com>
//The demo contains no error handling

#include <QApplication>
#include <qmessagebox.h>
#include <fstream>

#include <QMainWindow>
#include <QGraphicsScene>
#include <QActionGroup>
#include <QtGui>
#include <QString>
#include <QFileDialog>
#include <QInputDialog>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QSlider>
#include <QProgressBar>
#include <QMessageBox>
#include <fstream>   
#include <string>
#include <sstream>
#include <iomanip>
#include <list>

#include <boost/shared_ptr.hpp>   

#include <CGAL/basic.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Timer.h> 
#include <CGAL/Bbox_2.h>
#include <CGAL/iterator.h>
#include <CGAL/assertions_behaviour.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/approximated_offset_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>

#ifdef CGAL_USE_GMP
  #include <CGAL/Gmpq.h>
#else
  #include <CGAL/MP_Float.h>
  #include <CGAL/Quotient.h>
#endif

#include <QT5/Circular_polygons.h>
#include <QT5/Linear_polygons.h>
#include <QT5/Graphics_view_circular_polygon_input.h>
#include <CGAL/Gps_segment_traits_2.h>
//#include <QT5/Graphics_view_linear_polygon_input.h>

//#include <QT5/Gps_segement_traits_2_apurva.h>
/*
#include <QT5/PiecewiseGraphicsItemBase.h>
#include <QT5/PiecewiseBoundaryGraphicsItem.h>
#include <QT5/PiecewiseRegionGraphicsItem.h>
#include <QT5/PiecewiseSetGraphicsItem.h>
*/
#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/DemosMainWindow.h>
#include <CGAL/Qt/utility.h>
#include <CGAL/IO/Dxf_bsop_reader.h>

#include <CGAL/Qt/GraphicsViewNavigation.h>

#include <QFileDialog>

#include "ui_Boolean_set_operations_2.h"

#include "Typedefs.h"

#define CGAL_POLYGON_EDGE_DEFAULT_COLOR "blue"
#define CGAL_POLYGON_FILL_DEFAULT_COLOR "blue"

using namespace std;

typedef CGAL::Qt::Circular_set_graphics_item<Circular_polygon_set> Circular_GI;
typedef CGAL::Qt::Linear_set_graphics_item<Linear_polygon_set>     Linear_GI;

void show_warning(std::string aS)
{
  QMessageBox::warning(NULL, "Warning", QString(aS.c_str()));
}

void show_error(std::string aS)
{
  QMessageBox::critical(NULL, "Critical Error", QString(aS.c_str()));
}

void error(std::string aS)
{
  show_error(aS);

  throw std::runtime_error(aS);
}

enum { BLUE_GROUP, RED_GROUP, RESULT_GROUP } ;

enum { CIRCULAR_TYPE, LINEAR_TYPE } ;

QPen   sPens   [] = { QPen(QColor(0,0,255),0,Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin)
                    , QPen(QColor(255,0,0),0,Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin)
                    , QPen(QColor(0,255,0),0,Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin) 
                    } ;
                    
QBrush sBrushes[] = { QBrush(QColor(0,0,255,32 ))
                    , QBrush(QColor(255,0,0,32 ))
                    , QBrush(QColor(0,255,0,220))
                    } ;
struct Rep_base
{
  virtual ~Rep_base() {}
  
  virtual int type () const = 0 ;
   
  virtual CGAL::Qt::GraphicsItem* gi() const = 0 ;
  virtual CGAL::Qt::GraphicsItem* gi()       = 0 ;
  
  virtual void set_pen  ( QPen   const& aPen   ) = 0 ;
  virtual void set_brush( QBrush const& aBrush ) = 0 ;
  
  virtual QRectF bounding_rect() const { return gi()->boundingRect() ; }
  
  virtual bool is_empty() const = 0 ;
  
  virtual void clear               ()                         = 0 ;
  virtual void complement          ()                         = 0 ;
  virtual void assign              ( Rep_base const& aOther ) = 0 ;
  virtual void intersect           ( Rep_base const& aOther ) = 0 ;
  virtual void join                ( Rep_base const& aOther ) = 0 ;
  virtual void difference          ( Rep_base const& aOther ) = 0 ;
  virtual void symmetric_difference( Rep_base const& aOther ) = 0 ;
  
} ;



template<class GI_, class Set_>
class Rep : public Rep_base
{
public:

  typedef GI_  GI  ;
  typedef Set_ Set ;
  
  typedef Rep<GI,Set> Self ;
  
  Rep() { mGI = new GI(&mSet) ; }
  
  Set const& set() const { return mSet ; }
  Set      & set()       { return mSet ; }
  
  virtual CGAL::Qt::GraphicsItem* gi() const { return mGI; }
  virtual CGAL::Qt::GraphicsItem* gi()       { return mGI; }
  
  virtual void set_pen  ( QPen   const& aPen   ) { mGI->setPen  (aPen);   } 
  virtual void set_brush( QBrush const& aBrush ) { mGI->setBrush(aBrush); }
      
  virtual bool is_empty() const { return mSet.is_empty() ; }
  
  virtual void clear()                         
  { 
    try
    {
      mSet.clear() ; 
    }
    catch(...)
    {
      show_error("Exception thrown during boolean operation");
    }
  }
  
  virtual void complement()                         
  { 
    try
    {
      mSet.complement(); 
    } 
    catch(...)
    {
      show_error("Exception thrown during boolean operation");
    } 
  }
  
  virtual void assign( Rep_base const& aOther ) 
  { 
    try
    {
      mSet = cast(aOther).mSet; 
    } 
    catch(...)
    {
      show_error("Exception thrown during boolean operation");
    } 
  }
  
  virtual void intersect( Rep_base const& aOther ) 
  { 
    try
    {
      mSet.intersection( cast(aOther).mSet); 
    } 
    catch(...)
    {
      show_error("Exception thrown during boolean operation");
    } 
  }
  
  virtual void join( Rep_base const& aOther ) 
  { 
    try
    {
      mSet.join( cast(aOther).mSet); 
    } 
    catch(...)
    {
      show_error("Exception thrown during boolean operation");
    } 
  }
  
  virtual void difference( Rep_base const& aOther ) 
  { 
    try
    {
      mSet.difference( cast(aOther).mSet); 
    } 
    catch(...)
    {
      show_error("Exception thrown during boolean operation");
    } 
  }
  
  virtual void symmetric_difference( Rep_base const& aOther ) 
  { 
    try
    {
      mSet.symmetric_difference( cast(aOther).mSet); 
    } 
    catch(...)
    {
      show_error("Exception thrown during boolean operation");
    } 
  }
  
  static Self const& cast( Rep_base const& aOther ) { return dynamic_cast<Self const&>(aOther); }
  static Self      & cast( Rep_base      & aOther ) { return dynamic_cast<Self      &>(aOther); }
  
private:

  GI* mGI;
  Set mSet ;
} ;

class Circular_rep : public Rep<Circular_GI, Circular_polygon_set>
{
  typedef Rep<Circular_GI, Circular_polygon_set> Base ;
  
public:
  
  Circular_rep () : Base() {} 
  
  virtual int type() const { return CIRCULAR_TYPE ; }
} ;

class Linear_rep : public Rep<Linear_GI, Linear_polygon_set>
{
  typedef Rep<Linear_GI, Linear_polygon_set> Base ;
public:
  
  Linear_rep () : Base() {} //error is here
  /*
  virtual int type() const { return LINEAR_TYPE ; }*/
} ;

class Curve_set
{
  typedef boost::shared_ptr<Rep_base> Rep_ptr ;
  
public:
  
  Curve_set( int aType, QPen aPen, QBrush aBrush ) 
  :
    mPen  (aPen)
  , mBrush(aBrush)
  {
    reset_type(aType);
  }
  
  void reset_type( int aType ) 
  {
    
    mRep = aType == CIRCULAR_TYPE ? Rep_ptr(new Circular_rep())
                                  : Rep_ptr(NULL);//new Linear_rep  ()) ;
    mRep->set_pen  (mPen);
    mRep->set_brush(mBrush);
  }
  
  CGAL::Qt::GraphicsItem const* gi() const { return mRep->gi() ; }
  CGAL::Qt::GraphicsItem*       gi()       { return mRep->gi() ; }
  
  QRectF bounding_rect() const { return mRep->bounding_rect() ; }
  
  bool is_empty() const { return !mRep || mRep->is_empty(); }
  
  void clear      () { mRep->clear() ; }
  void complement () { mRep->complement() ; }
  
  void assign ( Curve_set const& aOther ) 
  {
    if ( is_circular() && aOther.is_circular() )
    {
      get_circular_rep()->assign( *aOther.get_circular_rep() ) ;
    }/*
    else
    {
      get_linear_rep()->assign( *aOther.get_linear_rep() ) ;
    }*/
  }
  
  void intersect( Curve_set const& aOther ) 
  {
    if ( is_circular() && aOther.is_circular() )
    {
      get_circular_rep()->intersect( *aOther.get_circular_rep() ) ;
    }/*
    else
    {
      get_linear_rep()->intersect( *aOther.get_linear_rep() ) ;
    } */
  }
  
  void join ( Curve_set const& aOther ) 
  {
    if ( is_circular() && aOther.is_circular() )
    {
      get_circular_rep()->join( *aOther.get_circular_rep() ) ;
    }/*
    else
    {
      get_linear_rep()->join( *aOther.get_linear_rep() ) ;
    }  */
  }
  
  void difference( Curve_set const& aOther ) 
  {
    if ( is_circular() && aOther.is_circular() )
    {
      get_circular_rep()->difference( *aOther.get_circular_rep() ) ;
    }/*
    else
    {
      get_linear_rep()->difference( *aOther.get_linear_rep() ) ;
    } */
  }
  
  void symmetric_difference( Curve_set const& aOther ) 
  {
    if ( is_circular() && aOther.is_circular() )
    {
      get_circular_rep()->symmetric_difference( *aOther.get_circular_rep() ) ;
    }/*
    else
    {
      get_linear_rep()->symmetric_difference( *aOther.get_linear_rep() ) ;
    }  */
  }
   
  Rep_base const& rep() const { return *mRep ; }
  Rep_base&       rep()       { return *mRep ; }
  
  bool is_circular() const { return mRep->type() == CIRCULAR_TYPE ; }  
  //bool is_linear  () const { return mRep->type() == LINEAR_TYPE ; } // no need keep it for now 
  
  Circular_rep const* get_circular_rep() const { return dynamic_cast<Circular_rep const*>( boost::get_pointer(mRep) ); }
  Circular_rep      * get_circular_rep()       { return dynamic_cast<Circular_rep*      >( boost::get_pointer(mRep) ); }
  //Linear_rep   const* get_linear_rep  () const { return dynamic_cast<Linear_rep   const*>( boost::get_pointer(mRep) ); }
  //Linear_rep        * get_linear_rep  ()       { return dynamic_cast<Linear_rep  *      >( boost::get_pointer(mRep) ); }
  
  Circular_polygon_set const& circular() const { return get_circular_rep()->set(); }
  Circular_polygon_set      & circular()       { return get_circular_rep()->set(); }
  //Linear_polygon_set   const& linear  () const { return get_linear_rep  ()->set(); }
  //Linear_polygon_set        & linear  ()       { return get_linear_rep  ()->set(); }
  
private:

  QPen                        mPen ;
  QBrush                      mBrush ;
  boost::shared_ptr<Rep_base> mRep ;
  
} ;

typedef std::vector<Curve_set> Curve_set_container ;

typedef Curve_set_container::const_iterator Curve_set_const_iterator ;
typedef Curve_set_container::iterator       Curve_set_iterator ;


class MainWindow : public CGAL::Qt::DemosMainWindow ,public Ui::Boolean_set_operations_2
{
  Q_OBJECT// removing it gives error ui not declared
  
private:  

  QGraphicsScene                                                   mScene;
  bool                                                             mCircular_active ;
  bool                                                             mBlue_active ;
  Curve_set_container                                              mCurve_sets ;
  Circular_region_source_container                                 mBlue_circular_sources ;
  Circular_region_source_container                                 mRed_circular_sources ;
  //Linear_region_source_container                                   mBlue_linear_sources ; 
  //Linear_region_source_container                                   mRed_linear_sources ; 
  //CGAL::Qt::Graphics_view_linear_polygon_input<Kernel>*     mLinearInput ;
  CGAL::Qt::Graphics_view_circular_polygon_input<Kernel>* mCircularInput ;
   
public:

  MainWindow();

private:
  
  void dragEnterEvent(QDragEnterEvent *event);
  void dropEvent(QDropEvent *event);
  void zoomToFit();
  
protected slots:
  
  void open( QString filename ) ;

public slots:
  
  void processInput(CGAL::Object o);
  void on_actionNew_triggered() ;
  void on_actionRecenter_triggered();

  //void on_actionInsertLinear_toggled  (bool aChecked);
  void on_actionInsertCircular_triggered();
    
signals:

  void modelChanged();
  
private:/*
  void modelChanged()
  {
    emit(changed());
  }
  */
  bool ask_user_yesno( const char* aTitle, const char* aQuestion )
  {
    return QMessageBox::warning(this
                               ,aTitle
                               ,QString(aQuestion)
                               ,"&Yes"
                               ,"&No"
                               ,QString::null
                               , 1
                               , 1 
                               ) == 0 ;
  }
   
  Curve_set& set( int aGroup ) { return mCurve_sets[aGroup] ; }
  
  Curve_set& blue_set  () { return set(BLUE_GROUP)  ; }
  Curve_set& red_set   () { return set(RED_GROUP)   ; }
  Curve_set& result_set() { return set(RESULT_GROUP); }

  int active_group() const { return mBlue_active ? BLUE_GROUP : RED_GROUP ; }
  
  Curve_set& active_set()   { return set(active_group()) ; }

  Circular_region_source_container const& blue_circular_sources() const { return mBlue_circular_sources ; }
  Circular_region_source_container      & blue_circular_sources()       { return mBlue_circular_sources ; }

  Circular_region_source_container const& red_circular_sources () const { return mRed_circular_sources ; }
  Circular_region_source_container      & red_circular_sources ()       { return mRed_circular_sources ; }
  /*
  Linear_region_source_container const& blue_linear_sources() const { return mBlue_linear_sources ; }
  Linear_region_source_container      & blue_linear_sources()       { return mBlue_linear_sources ; }

  Linear_region_source_container const& red_linear_sources () const { return mRed_linear_sources ; }
  Linear_region_source_container      & red_linear_sources ()       { return mRed_linear_sources ; }

  Linear_region_source_container const& active_linear_sources() const { return mBlue_active ? mBlue_linear_sources : mRed_linear_sources ; }
  Linear_region_source_container      & active_linear_sources()       { return mBlue_active ? mBlue_linear_sources : mRed_linear_sources ; }
  */
  Circular_region_source_container const& active_circular_sources() const { return mBlue_active ? mBlue_circular_sources : mRed_circular_sources ; }
  Circular_region_source_container      & active_circular_sources()       { return mBlue_active ? mBlue_circular_sources : mRed_circular_sources ; }

  void ToogleView( int aGROUP, bool aChecked );
  
  void link_GI ( CGAL::Qt::GraphicsItem* aGI )
  {
    QObject::connect(this, SIGNAL(changed()), aGI, SLOT(modelChanged()));
    mScene.addItem( aGI );
  }
  
  void unlink_GI ( CGAL::Qt::GraphicsItem* aGI )
  {
    mScene.removeItem( aGI );
    QObject::disconnect(this, SIGNAL(changed()), aGI, SLOT(modelChanged()));
  }
  
  void switch_set_type( Curve_set& aSet, int aType );
  
  void switch_sets_type( int aType );
  
  bool ensure_circular_mode();
  
  //bool ensure_linear_mode();//see if it is need
};


MainWindow::MainWindow()
  : DemosMainWindow()
  , mCircular_active(true)
  , mBlue_active(true)
{
  setupUi(this);

  setAcceptDrops(true);
  cout<<"elementry setups"<<endl;
  mCurve_sets.push_back( Curve_set(CIRCULAR_TYPE, sPens[BLUE_GROUP]  , sBrushes[BLUE_GROUP]  ) ) ;
  mCurve_sets.push_back( Curve_set(CIRCULAR_TYPE, sPens[RED_GROUP]   , sBrushes[RED_GROUP]   ) ) ;
  mCurve_sets.push_back( Curve_set(CIRCULAR_TYPE, sPens[RESULT_GROUP], sBrushes[RESULT_GROUP]) ) ;
  cout<<"curve setups"<<endl;
  for( Curve_set_iterator si = mCurve_sets.begin(); si != mCurve_sets.end() ; ++ si )
    link_GI(si->gi()) ;
  
  //
  // Setup the mScene and the view
  //
  mScene.setItemIndexMethod(QGraphicsScene::NoIndex);
  mScene.setSceneRect(-100, -100, 100, 100);
  this->graphicsView->setScene(&mScene);
  this->graphicsView->setMouseTracking(true);

  // Turn the vertical axis upside down
  this->graphicsView->scale(1, -1);
    cout<<"UI setup"<<endl;                                                  
  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/index.html");
  this->addAboutCGAL();

  this->addRecentFiles(this->menuFile, this->actionQuit);
  cout<<"extra setup"<<endl;
  //mLinearInput  =new CGAL::Qt::Graphics_view_linear_polygon_input<Kernel>(this, &mScene);
  mCircularInput=new CGAL::Qt::Graphics_view_circular_polygon_input<Kernel>(this, &mScene);
  
  //QObject::connect(mLinearInput  , SIGNAL(generate(CGAL::Object)), this, SLOT(processInput(CGAL::Object)));
  QObject::connect(mCircularInput, SIGNAL(generate(CGAL::Object)), this, SLOT(processInput(CGAL::Object)));
  
  QObject::connect(this->actionQuit, SIGNAL(triggered()), this, SLOT(close()));
  QObject::connect(this->actionInsertCircular, SIGNAL(triggered()), this, SLOT(on_actionInsertCircular_triggered()));
  //QObject::connect(this, SIGNAL(openRecentFile(QString)), this, SLOT(open(QString)));//for file handling
  cout<<"connecting stuff"<<endl;
}

void MainWindow::on_actionNew_triggered() 
{
  for( Curve_set_iterator si = mCurve_sets.begin(); si != mCurve_sets.end() ; ++ si )
    si->clear();
    
  blue_circular_sources().clear();
    
  ToogleView(BLUE_GROUP  ,true);
  mCircular_active = true ;
  mBlue_active =  true ;
  modelChanged();
  
}

void MainWindow::on_actionRecenter_triggered()
{
  zoomToFit();
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

Circular_polygon linear_2_circ( Circular_Linear_polygon const& pgn )
{
  CGAL::Cartesian_converter<Kernel,Kernel> convert ;
  
  Circular_polygon rCP;
  
  for( Circular_Linear_polygon::Edge_const_iterator ei = pgn.edges_begin(); ei != pgn.edges_end(); ++ei )
  {
    if  ( ei->source() != ei->target() )
      rCP.push_back( Circular_X_monotone_curve( convert(ei->source()), convert(ei->target())) );
  }  

  return rCP;
}

Circular_polygon_with_holes linear_2_circ( Circular_Linear_polygon_with_holes const& pwh )
{
  Circular_polygon_with_holes rCP( linear_2_circ(pwh.outer_boundary()) ) ;
  
  for( Circular_Linear_polygon_with_holes::Hole_const_iterator hi = pwh.holes_begin(); hi != pwh.holes_end(); ++ hi )
    rCP.add_hole( linear_2_circ(*hi)  );

  return rCP;
}
//check out
void MainWindow::switch_set_type( Curve_set& aSet, int aType )
{
  unlink_GI( aSet.gi() ) ;
  
  aSet.reset_type(aType);
  
  link_GI( aSet.gi() ) ;
  
  modelChanged();
}

void MainWindow::switch_sets_type( int aType )
{
  switch_set_type( blue_set  (), aType ) ; 
  switch_set_type( red_set   (), aType ) ; 
  switch_set_type( result_set(), aType ) ; 
  
}

bool MainWindow::ensure_circular_mode()
{
  if ( ! mCircular_active )
  {
    bool lProceed = blue_set().is_empty() && red_set().is_empty() ;
    
    if ( ! lProceed )
      lProceed = ask_user_yesno("Linear/Circular mode switch"
                               ,"You are about to load a linear or circular poygon, but there are bezier curves already loaded.\n" \
                                "Both types are not interoperable. In order to proceed, the bezier curves must be removed first.\n" \
                                "OK to remove and proceed?\n"
                               ) ;
      
    if ( lProceed )
    {
      switch_sets_type(CIRCULAR_TYPE);
      mCircular_active = true ;
    }
  }
  return mCircular_active ;
}
//check out

void MainWindow::open( QString fileName )
{
  cout<<"To be done"<<endl;
    if(! fileName.isEmpty())
  {
    bool lRead = false ;
     
    if ( lRead )
    {
      modelChanged();
      zoomToFit();
      this->addToRecentFiles(fileName);
      
    }
  }  
}

void MainWindow::on_actionInsertCircular_triggered()
{
  cout<<"signal triggered"<<endl;
    bool aChecked=1;//temporality;
    if(aChecked)
       mScene.installEventFilter(mCircularInput);
  else mScene.removeEventFilter (mCircularInput);
}

void MainWindow::processInput(CGAL::Object o )
{
  Circular_polygon lCI ;
    mBlue_active =  true ;

  if ( CGAL::assign(lCI, o) )
  {
    if ( ensure_circular_mode() )
    {
      CGAL::Orientation o = lCI.orientation();
      if ( o == CGAL::CLOCKWISE )
        lCI.reverse_orientation();
        
      Circular_polygon_with_holes lCPWH(lCI);
      active_set().circular().join(lCPWH) ;  
      
      active_circular_sources().push_back(lCPWH);
    }
  }
  modelChanged();  
}

void MainWindow::ToogleView( int aGROUP, bool aChecked )
{
  if ( aChecked )
       set(aGROUP).gi()->show();
  else set(aGROUP).gi()->hide();
}


void MainWindow::zoomToFit()
{
  boost::optional<QRectF> lTotalRect ;
  
  for ( Curve_set_const_iterator si = mCurve_sets.begin() ; si != mCurve_sets.end() ; ++ si )
  {
    if ( !si->is_empty() ) 
    {
      QRectF lRect = si->bounding_rect();
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

#include "Boolean_set_operations_2.moc"
#include <CGAL/Qt/resources.h>
int main(int argc, char *argv[])
{
  //QApplication a(argc, argv);
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Boolean_operations_2 demo");
  CGAL_QT_INIT_RESOURCES;
  try
  {
//std::cout<<"hello";    
    MainWindow w;
    w.show();

    return app.exec();
  }
  catch (const std::exception e)
  {
    std::string s = e.what();
    show_error("Exception throne during run of the program:\n" + s);
  }
}
