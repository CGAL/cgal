// Copyright (c) 2009  GeometryFactory Sarl (France).
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
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>


//#define ENABLE_TRACE

#ifdef ENABLE_TRACE
#  define TRACE(m) { std::ostringstream ss ; ss << m << std::endl ; trace(ss.str()); }
#else
#  define TRACE(m) 
#endif


#include <fstream>
#include <string>
#include <sstream>
#include <list>

void trace( std::string s )
{
  static std::ofstream out("log.txt");
  out << s ;
}

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

#include <CGAL/Qt/BezierCurves.h>
#include <CGAL/Qt/CircularPolygons.h>
#include <CGAL/Qt/GraphicsViewBezierPolygonInput.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/DemosMainWindow.h>
#include <CGAL/Qt/utility.h>
#include <CGAL/IO/Dxf_bsop_reader.h>

// the two base classes
#include "ui_boolean_operations_2.h"

#include "typedefs.h"


void show_warning( std::string aS )
{
  QMessageBox::warning(NULL,"Warning",QString(aS.c_str()) ) ;
}

void show_error( std::string aS )
{
  QMessageBox::critical(NULL,"Critical Error",QString(aS.c_str()) ) ;
}

void error( std::string aS )
{
  show_error(aS);
  
  throw std::runtime_error(aS);  
}

void error_handler ( char const* what, char const* expr, char const* file, int line, char const* msg )
{
  std::ostringstream ss ;
  
  ss << "CGAL error: " << what << " violation!" << std::endl
     << "Expr: " << expr << std::endl
     << "File: " << file << std::endl
     << "Line: " << line << std::endl;
  if ( msg != 0)
    ss << "Explanation:" << msg << std::endl;
    
  error(ss.str());
    
}


enum { BLUE_GROUP, RED_GROUP, RESULT_GROUP } ;

enum { CIRCULAR_TYPE, BEZIER_TYPE } ;


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

class Bezier_rep : public Rep<Bezier_GI, Bezier_polygon_set>
{
  typedef Rep<Bezier_GI, Bezier_polygon_set> Base ;
  
public:
  
  
  Bezier_rep () : Base() {} 
  
  virtual int type() const { return BEZIER_TYPE ; }
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
                                  : Rep_ptr(new Bezier_rep  ()) ;
         
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
    }
    else if ( is_bezier() && aOther.is_bezier() )
    {
      get_bezier_rep()->assign( *aOther.get_bezier_rep() ) ;
    }  
  }
  
  void intersect( Curve_set const& aOther ) 
  {
    if ( is_circular() && aOther.is_circular() )
    {
      get_circular_rep()->intersect( *aOther.get_circular_rep() ) ;
    }
    else if ( is_bezier() && aOther.is_bezier() )
    {
      get_bezier_rep()->intersect( *aOther.get_bezier_rep() ) ;
    }  
  }
  
  void join ( Curve_set const& aOther ) 
  {
    if ( is_circular() && aOther.is_circular() )
    {
      get_circular_rep()->join( *aOther.get_circular_rep() ) ;
    }
    else if ( is_bezier() && aOther.is_bezier() )
    {
      get_bezier_rep()->join( *aOther.get_bezier_rep() ) ;
    }  
  }
  
  void difference( Curve_set const& aOther ) 
  {
    if ( is_circular() && aOther.is_circular() )
    {
      get_circular_rep()->difference( *aOther.get_circular_rep() ) ;
    }
    else if ( is_bezier() && aOther.is_bezier() )
    {
      get_bezier_rep()->difference( *aOther.get_bezier_rep() ) ;
    }  
  }
  
  void symmetric_difference( Curve_set const& aOther ) 
  {
    if ( is_circular() && aOther.is_circular() )
    {
      get_circular_rep()->symmetric_difference( *aOther.get_circular_rep() ) ;
    }
    else if ( is_bezier() && aOther.is_bezier() )
    {
      get_bezier_rep()->symmetric_difference( *aOther.get_bezier_rep() ) ;
    }  
  }
   
  Rep_base const& rep() const { return *mRep ; }
  Rep_base&       rep()       { return *mRep ; }
  
  bool is_circular() const { return mRep->type() == CIRCULAR_TYPE ; }  
  bool is_bezier  () const { return mRep->type() == BEZIER_TYPE ; }  
  
  Circular_rep const* get_circular_rep() const { return dynamic_cast<Circular_rep const*>( boost::get_pointer(mRep) ); }
  Circular_rep      * get_circular_rep()       { return dynamic_cast<Circular_rep*      >( boost::get_pointer(mRep) ); }
  Bezier_rep   const* get_bezier_rep  () const { return dynamic_cast<Bezier_rep   const*>( boost::get_pointer(mRep) ); }
  Bezier_rep        * get_bezier_rep  ()       { return dynamic_cast<Bezier_rep  *      >( boost::get_pointer(mRep) ); }
  
  Circular_polygon_set const& circular() const { return get_circular_rep()->set(); }
  Circular_polygon_set      & circular()       { return get_circular_rep()->set(); }
  Bezier_polygon_set   const& bezier  () const { return get_bezier_rep  ()->set(); }
  Bezier_polygon_set        & bezier  ()       { return get_bezier_rep  ()->set(); }
  
private:

  QPen                        mPen ;
  QBrush                      mBrush ;
  boost::shared_ptr<Rep_base> mRep ;
  
} ;

typedef std::vector<Curve_set> Curve_set_vector ;

typedef Curve_set_vector::const_iterator Curve_set_const_iterator ;
typedef Curve_set_vector::iterator       Curve_set_iterator ;

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Boolean_operations_2
{
  Q_OBJECT
  
private:  

  QGraphicsScene                                           mScene;
  bool                                                     mCircular_active ;
  bool                                                     mBlue_active ;
  Curve_set_vector                                         mCurve_sets ;
  CGAL::Qt::GraphicsViewBezierPolygonInput<Bezier_traits>* mBezierInput ;
    
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
  void on_actionOpenLinear_triggered() ;
  void on_actionOpenDXF_triggered() ;
  void on_actionOpenBezier_triggered() ;
  void on_actionSave_triggered() ;
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
  
  void on_actionInsertPWH_toggled(bool aChecked);
  
  void on_checkboxShowBlue_toggled      (bool aChecked) { ToogleView(BLUE_GROUP  ,aChecked); }
  void on_checkboxShowRed_toggled       (bool aChecked) { ToogleView(RED_GROUP   ,aChecked); }
  void on_checkboxShowResult_toggled    (bool aChecked) { ToogleView(RESULT_GROUP,aChecked); }
  
  void on_radioMakeBlueActive_toggled(bool aChecked) { mBlue_active =  aChecked ; }
  void on_radioMakeRedActive_toggled (bool aChecked) { mBlue_active = !aChecked ; }
  
signals:

  void changed();
  
private:
  
  void modelChanged()
  {
    emit(changed());
    zoomToFit();
  }
  
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
  
  void SetViewBlue  ( bool aChecked ) { checkboxShowBlue  ->setChecked(aChecked); }  
  void SetViewRed   ( bool aChecked ) { checkboxShowRed   ->setChecked(aChecked); }  
  void SetViewResult( bool aChecked ) { checkboxShowResult->setChecked(aChecked); }  

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
  
  bool ensure_bezier_mode();
};


MainWindow::MainWindow()
  : DemosMainWindow()
  , mCircular_active(true)
  , mBlue_active(true)
{
  CGAL::set_error_handler  (error_handler);
  CGAL::set_warning_handler(error_handler);
  
  setupUi(this);

  setAcceptDrops(true);

  mCurve_sets.push_back( Curve_set(CIRCULAR_TYPE, sPens[BLUE_GROUP]  , sBrushes[BLUE_GROUP]  ) ) ;
  mCurve_sets.push_back( Curve_set(CIRCULAR_TYPE, sPens[RED_GROUP]   , sBrushes[RED_GROUP]   ) ) ;
  mCurve_sets.push_back( Curve_set(CIRCULAR_TYPE, sPens[RESULT_GROUP], sBrushes[RESULT_GROUP]) ) ;
  
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
                                                      
  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/index.html");
  this->addAboutCGAL();

  this->addRecentFiles(this->menuFile, this->actionQuit);
  
  mBezierInput = new CGAL::Qt::GraphicsViewBezierPolygonInput<Bezier_traits>(this, &mScene);
  
  QObject::connect(mBezierInput, SIGNAL(generate(CGAL::Object)), this, SLOT(processInput(CGAL::Object)));

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
    si->clear();
    
  SetViewBlue  (true);
  SetViewRed   (true);
  SetViewResult(true);
  
  mCircular_active = true ;
  
  radioMakeBlueActive->setChecked(true);
  
  modelChanged();
  
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

Circular_polygon linear_2_circ( Linear_polygon const& pgn )
{
  CGAL::Cartesian_converter<Linear_kernel,Circular_kernel> convert ;
  
  Circular_polygon rCP;
  
  for( Linear_polygon::Edge_const_iterator ei = pgn.edges_begin(); ei != pgn.edges_end(); ++ei )
  {
    if  ( ei->source() != ei->target() )
      rCP.push_back( Circular_X_monotone_curve( convert(ei->source()), convert(ei->target())) );
  }  

  return rCP;
}

Circular_polygon_with_holes linear_2_circ( Linear_polygon_with_holes const& pwh )
{
  Circular_polygon_with_holes rCP( linear_2_circ(pwh.outer_boundary()) ) ;
  
  for( Linear_polygon_with_holes::Hole_const_iterator hi = pwh.holes_begin(); hi != pwh.holes_end(); ++ hi )
    rCP.add_hole( linear_2_circ(*hi)  );

  return rCP;
}

bool read_linear ( QString aFileName, Circular_polygon_set& rSet )
{
  bool rOK = false ;
  
  std::ifstream in_file (qPrintable(aFileName));

  if ( in_file )
  {
    unsigned int n_regions ;
    in_file >> n_regions;
    
    for ( unsigned int r = 0 ; r < n_regions ; ++ r )
    {
      unsigned int n_boundaries;
      in_file >> n_boundaries;
      
      Circular_polygon outer ;
      std::vector<Circular_polygon> holes ;
      
      for ( unsigned int r = 0 ; r < n_boundaries ; ++ r )
      {
        Linear_polygon p ;
        in_file >> p ;
        
        if ( r == 0 )
             outer = linear_2_circ(p); 
        else holes.push_back( linear_2_circ(p) );
      }
      
      rSet.join( Circular_polygon_with_holes(outer,holes.begin(),holes.end()) ) ;    
      rOK = true ;
    }
    
  }
  
  return rOK ;
}

bool read_dxf ( QString aFileName, Circular_polygon_set& rSet )
{
  bool rOK = false ;
  
  std::ifstream in_file (qPrintable(aFileName));

  if ( in_file )
  {
    CGAL::Dxf_bsop_reader<Circular_kernel>   reader;
    std::vector<Circular_polygon>            circ_polygons;
    std::vector<Circular_polygon_with_holes> circ_polygons_with_holes;
    
    reader(in_file
          ,std::back_inserter(circ_polygons)
          ,std::back_inserter(circ_polygons_with_holes)
          ,false
          );
          
    rSet.join( circ_polygons.begin()           , circ_polygons.end() 
             , circ_polygons_with_holes.begin(), circ_polygons_with_holes.end()
             ) ;
                      
    rOK = true ;
  }
  
  return rOK ;
}

Bezier_curve read_bezier_curve ( std::istream& is, bool aDoubleFormat )
{
  // Read the number of control points.
  unsigned int  n;

  is >> n;
  
  // Read the control points.
  std::vector<Bezier_rat_point> ctrl_pts;
  
  for ( unsigned int k = 0; k < n; k++)
  {
    Bezier_rat_point p ;
    if ( aDoubleFormat )
    {
      double x,y ;
      is >> x >> y ;
      p = Bezier_rat_point(x,y);
    }
    else
    {
      is >> p ;
    }
    
    if ( k == 0 || ctrl_pts[k-1] != p ) 
    {
      ctrl_pts.push_back(p) ;
    }
  }

  return Bezier_curve(ctrl_pts.begin(),ctrl_pts.end());
}

bool read_bezier ( QString aFileName, Bezier_polygon_set& rSet )
{
  
  bool rOK = false ;
  
  std::ifstream in_file (qPrintable(aFileName));
  
  if ( in_file )
  {
    try
    {
      
      std::string format ;
      std::getline(in_file,format);
      
      bool lDoubleFormat = ( format.length() >= 6 && format.substr(0,6) == "DOUBLE") ;
                              
      // Red the number of bezier polygon with holes
      unsigned int n_regions ;
      in_file >> n_regions;
      
      for ( unsigned int r = 0 ; r < n_regions ; ++ r )
      {
        Bezier_polygon_vector  polygons ;
        
        // Read the number of bezier curves.
        unsigned int n_boundaries;
        in_file >> n_boundaries;
      
        for ( unsigned int b = 0 ; b < n_boundaries ; ++ b )
        {
          // Read the number of bezier curves.
          unsigned int n_curves;
          in_file >> n_curves;
          
          // Read the curves one by one, and construct the general polygon these
          // curve form (the outer boundary and the holes inside it).
          
          std::list<Bezier_X_monotone_curve> xcvs;
        
          for ( unsigned int k = 0; k < n_curves; ++ k ) 
          {
            // Read the current curve and subdivide it into x-monotone subcurves.
            
            std::list<CGAL::Object>                 x_objs;
            std::list<CGAL::Object>::const_iterator xoit;
            Bezier_X_monotone_curve                 xcv;
            Bezier_traits                           traits;
            Bezier_traits::Make_x_monotone_2        make_x_monotone = traits.make_x_monotone_2_object();
        
            Bezier_curve B = read_bezier_curve(in_file, lDoubleFormat);
            if ( B.number_of_control_points() >= 2 )
            {
              //TRACE( "region " << r << " boundary " << b << " curve " << k );
                
              make_x_monotone (B, std::back_inserter (x_objs));
              
              for (xoit = x_objs.begin(); xoit != x_objs.end(); ++xoit) 
              {
                if (CGAL::assign (xcv, *xoit))
                {
                  //TRACE( " X montonote: " << xcv.source() << " -> " << xcv.target() << ( xcv.is_directed_right() ? " RIGHT":" LEFT") << ( xcv.is_vertical() ? " VERTICAL" : "")) ;
                  xcvs.push_back (xcv);
                }  
              }
            }
          }  
            
          Bezier_polygon  pgn (xcvs.begin(), xcvs.end());
          
          CGAL::Orientation  orient = pgn.orientation();
          //TRACE( "  Orientation: " << orient ) ;
            
          if (( b == 0 && orient == CGAL::CLOCKWISE) || ( b > 0 && orient == CGAL::COUNTERCLOCKWISE))
          {
            //TRACE( "Reversing orientation: " ) ;
            pgn.reverse_orientation();
          }
            
          polygons.push_back (pgn);
        }
      
        if ( polygons.size() > 0 )
        {
          Bezier_polygon_with_holes pwh(polygons.front());
          
          if ( polygons.size() > 1 )
          {
            for ( Bezier_polygon_vector::const_iterator it = CGAL::successor(polygons.begin())
                ; it != polygons.end()
                ; ++ it 
                )
              pwh.add_hole(*it);    
          }
          
          //if ( is_valid_polygon_with_holes(pwh, rSet.traits() ) )
          if ( true )
          {
            rSet.join(pwh) ;      
          }
          else
          {
            show_warning( "Bezier polygon is not valid" );
          }
        }
        
        rOK = true ;
      }
      
    }
    catch(...)
    {
      show_error("Exception ocurred during reading of bezier polygon set.");
    } 
  }
  
  return rOK ;
}

bool save_circular ( QString aFileName, Circular_polygon_set& rSet )
{
  bool rOK = false ;
  
  return rOK ;
}

void save_bezier_polygon( std::ostream& out_file, Bezier_polygon const& aBP )
{
  out_file << aBP.size() << std::endl ;
  
  for ( Bezier_polygon::Curve_const_iterator cit = aBP.curves_begin() ; cit != aBP.curves_end() ; ++ cit )
  {
    typedef std::vector<Linear_point> Linear_point_vector ;
    
    Linear_point_vector lQ ;

    CGAL::Qt::Bezier_helper::approximated_clip(*cit,std::back_inserter(lQ));  
    
    out_file << lQ.size() << std::endl ;
    
    for ( Linear_point_vector::const_iterator pit = lQ.begin() ; pit != lQ.end() ; ++ pit )
    {
      out_file << pit->x() << " " << pit->y() << std::endl ;
    }
  }
}
bool save_bezier ( QString aFileName, Bezier_polygon_set const& aSet )
{
  bool rOK = false ;
  
  std::ofstream out_file( qPrintable(aFileName) ) ;
  if ( out_file )
  {
    out_file << "DOUBLE" << std::endl ;
    
    std::vector<Bezier_polygon_with_holes> bpwh_container;
  
    aSet.polygons_with_holes( std::back_inserter(bpwh_container) ) ;
  
    out_file << bpwh_container.size() << std::endl ;
    
    for( std::vector<Bezier_polygon_with_holes>::const_iterator rit = bpwh_container.begin(); rit != bpwh_container.end() ; ++ rit )
    {
      Bezier_polygon_with_holes bpwh = *rit ;
      
      out_file << ( 1 + bpwh.number_of_holes() ) << std::endl ;
  
      save_bezier_polygon( out_file, bpwh.outer_boundary() ) ;
      
      for ( Bezier_polygon_with_holes::Hole_const_iterator hit = bpwh.holes_begin() ; hit != bpwh.holes_end() ; ++ hit )
        save_bezier_polygon(out_file, *hit);
      
      rOK = true ;
    }
  }
  
  return rOK ;
  
}

void MainWindow::on_actionOpenLinear_triggered()
{
  open(QFileDialog::getOpenFileName(this, tr("Open Linear Polygon"), "../data", tr("Linear Curve files (*.lps)") ));
}

void MainWindow::on_actionOpenDXF_triggered()
{
  open(QFileDialog::getOpenFileName(this, tr("Open DXF"), "../data", tr("DXF files (*.dxf)") ));
}

void MainWindow::on_actionOpenBezier_triggered()
{
  open(QFileDialog::getOpenFileName(this, tr("Open Bezier Polygon"), "../data", tr("Bezier Curve files (*.bps)") ));
}

void MainWindow::on_actionSave_triggered()
{
  if ( mCircular_active )
  {
    if ( !save_circular(QFileDialog::getSaveFileName(this, tr("Save Acive Circular Polygon Set"), "../data", tr("Linear Curve files (*.lps)") ) 
                       ,active_set().circular()
                       )
       )
    {
      show_error("Caanoit save circular polygon set.");
    }
       
  }
  else
  {
    if ( !save_bezier(QFileDialog::getSaveFileName(this, tr("Save Acive Bezier Polygon Set"), "../data", tr("Bezier Curve files (*.bps)") )
                     ,active_set().bezier() 
                     )
       )
    {
      show_error("Caanoit save bezier polygon set.");
    }
  }

}

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

bool MainWindow::ensure_bezier_mode()
{
  if ( mCircular_active )
  {
    bool lProceed = blue_set().is_empty() && red_set().is_empty() ;
    
    if ( ! lProceed )
      lProceed = ask_user_yesno("Bezier mode switch"
                               ,"You are about to load a Bezier curve, but there are linear and/or circular polygons already loaded.\n" \
                                "Both types are not interoperable. In order to proceed, the polygons must be removed first.\n" \
                                "OK to remove and proceed?\n"
                               ) ;
      
    if ( lProceed )
    {
      switch_sets_type(BEZIER_TYPE);
      mCircular_active = false ;
    }
  }
  return !mCircular_active ;
}

void MainWindow::open( QString fileName )
{
  if(! fileName.isEmpty())
  {
    bool lRead = false ;
    
    if(fileName.endsWith(".lps"))
    {
      if ( ensure_circular_mode() )
        lRead = read_linear(fileName,active_set().circular()) ;
    }
    else if (fileName.endsWith(".dxf"))
    {
      if ( ensure_circular_mode() )
        lRead = read_dxf(fileName,active_set().circular()) ;
    }
    else if (fileName.endsWith(".bps"))
    {
      if ( ensure_bezier_mode() )
        lRead = read_bezier(fileName,active_set().bezier()) ;
    }
     
    if ( lRead )
    {
      modelChanged();
      this->addToRecentFiles(fileName);
    }
  }  
}

void MainWindow::on_actionInsertPWH_toggled(bool aChecked)
{
  if(aChecked)
       mScene.installEventFilter(mBezierInput);
  else mScene.removeEventFilter (mBezierInput);
}

void MainWindow::processInput(CGAL::Object o )
{
  Bezier_polygon lBP ;
  if(CGAL::assign(lBP, o))
  {
    if ( ensure_bezier_mode() )
    {
      CGAL::Orientation o = lBP.orientation();
      if ( o == CGAL::CLOCKWISE )
        lBP.reverse_orientation();
      active_set().bezier().join( Bezier_polygon_with_holes(lBP) ) ;  
      
    }
  }
  modelChanged();  
}

void MainWindow::on_actionIntersection_triggered() 
{
  bool lDone = false ;
  
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);
  
  if ( !blue_set().is_empty() && !red_set().is_empty() )
  {
    result_set().assign( red_set() ) ;
    result_set().intersect(blue_set());
    lDone = true ;
  }
  
  this->setCursor(old);
  
  if ( lDone )
  {
    //SetViewBlue(false); SetViewRed(false); SetViewResult(true);
    
    modelChanged();
  }
}

void MainWindow::on_actionUnion_triggered() 
{
  bool lDone = false ;
  
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);
  
  if ( !blue_set().is_empty() && !red_set().is_empty() )
  {
    result_set().assign( red_set() ) ;
    result_set().join(blue_set());
    lDone = true ;
  }
  
  this->setCursor(old);
  
  if ( lDone )
  {
    //SetViewBlue(false);  SetViewRed(false);  SetViewResult(true);
    
    modelChanged();
  }
}

void MainWindow::on_actionBlueMinusRed_triggered() 
{
  bool lDone = false ;
  
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);
  
  if ( !blue_set().is_empty() && !red_set().is_empty() )
  {
    result_set().assign( blue_set() ) ;
    result_set().difference(red_set());
    lDone = true ;
  }
  
  this->setCursor(old);
    
  if ( lDone )
  {
    //SetViewBlue(false);  SetViewRed(false); SetViewResult(true);
    
    modelChanged();
  }
}

void MainWindow::on_actionRedMinusBlue_triggered() 
{
  bool lDone = false ;
  
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);
  
  if ( !blue_set().is_empty() && !red_set().is_empty() )
  {
    result_set().assign( red_set() ) ;
    result_set().difference(blue_set());
    lDone = true ;
  }
  
  this->setCursor(old);
    
  if ( lDone )
  {
    //SetViewBlue(false);  SetViewRed(false); SetViewResult(true);
    
    modelChanged();
  }
}

void MainWindow::on_actionSymmDiff_triggered() 
{
  bool lDone = false ;
  
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);
  
  if ( !blue_set().is_empty() && !red_set().is_empty() )
  {
    result_set().assign( red_set() ) ;
    result_set().symmetric_difference(blue_set());
    lDone = true ;
  }

  this->setCursor(old);
    
  
  if ( lDone )
  {
    //SetViewBlue(false); SetViewRed(false); SetViewResult(true);
    
    modelChanged();
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
  
  if ( !blue_set().is_empty() )
  {
    result_set().assign( blue_set() ) ;
    result_set().complement();
    lDone = true ;
  }
  
  this->setCursor(old);
    
  if ( lDone )
  {
    //SetViewBlue(false); SetViewRed(false); SetViewResult(true);
    
    modelChanged();
  }
}

void MainWindow::on_actionRedComplement_triggered()
{
  bool lDone = false ;
  
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);
  
  if ( !red_set().is_empty() )
  {
    result_set().assign( red_set() ) ;
    result_set().complement();
    lDone = true ;
  }
  
  this->setCursor(old);
    
  if ( lDone )
  {
    //SetViewBlue(false);  SetViewRed(false); SetViewResult(true);
    
    modelChanged();
  }
}

void MainWindow::on_actionAllBlue_triggered()
{
  bool lDone = false ;
  
  bool lProceed = result_set().is_empty() ? ask_user_yesno("Store result", "Result is empty, all polygons will be deleted\n continue anyway?\n")
                                          : true ;
                                          
  if ( lProceed ) 
  {
    blue_set().assign( result_set() ) ;
    result_set().clear();
    radioMakeRedActive->setChecked(true);
    lDone = true ;
  }
    
  if ( lDone )
  {
    //SetViewBlue(true);  SetViewRed(false); SetViewResult(true);
    
    modelChanged();
  }
}

void MainWindow::on_actionAllRed_triggered()
{
  bool lDone = false ;
  
  bool lProceed = result_set().is_empty() ? ask_user_yesno("Store result", "Result is empty, all polygons will be deleted\n continue anyway?\n")
                                          : true ;
  
  if ( lProceed ) 
  {
    red_set().assign( result_set() ) ;
    result_set().clear();
    radioMakeBlueActive->setChecked(true);
    lDone = true ;
  }
    
  if ( lDone )
  {
    //SetViewBlue(false); SetViewRed(true);  SetViewResult(true);
    
    modelChanged();
  }
}
void MainWindow::on_actionDeleteBlue_triggered()
{
  blue_set().clear();
    
  //SetViewBlue(true);SetViewRed(true); SetViewResult(true);
  
  modelChanged();
}

void MainWindow::on_actionDeleteRed_triggered()
{
  red_set().clear();
    
  //SetViewBlue(true); SetViewRed(true); SetViewResult(true);
  
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


