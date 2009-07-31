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
#include <CGAL/Circular_polygon_with_holes_sampler_2.h>
#include <CGAL/Bezier_polygon_with_holes_sampler_2.h>

#ifdef CGAL_USE_GMP
  #include <CGAL/Gmpq.h>
#else
  #include <CGAL/MP_Float.h>
  #include <CGAL/Quotient.h>
#endif

#include <CGAL/Qt/GeneralPolygonSetGraphicsItem.h>
//#include <CGAL/Qt/GraphicsViewPolygonWithHolesInput.h>
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
  
  virtual void clear               ()                         { mSet.clear() ; }
  virtual void complement          ()                         { mSet.complement(); }
  virtual void assign              ( Rep_base const& aOther ) { mSet = cast(aOther).mSet; }
  virtual void intersect           ( Rep_base const& aOther ) { mSet.intersection        ( cast(aOther).mSet); }
  virtual void join                ( Rep_base const& aOther ) { mSet.join                ( cast(aOther).mSet); }
  virtual void difference          ( Rep_base const& aOther ) { mSet.difference          ( cast(aOther).mSet); }
  virtual void symmetric_difference( Rep_base const& aOther ) { mSet.symmetric_difference( cast(aOther).mSet); }
  
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

  QGraphicsScene   mScene;
  bool             mCircular_active ;
  bool             mBlue_active ;
  Curve_set_vector mCurve_sets ;
  
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
  void on_actionOpenLinear_triggered() ;
  void on_actionOpenDXF_triggered() ;
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
  
  void on_checkboxShowBlue_toggled      (bool aChecked) { ToogleView(BLUE_GROUP  ,aChecked); }
  void on_checkboxShowRed_toggled       (bool aChecked) { ToogleView(RED_GROUP   ,aChecked); }
  void on_checkboxShowResult_toggled    (bool aChecked) { ToogleView(RESULT_GROUP,aChecked); }
  
  void on_radioMakeBlueActive_toggled(bool aChecked) { mBlue_active =  aChecked ; }
  void on_radioMakeRedActive_toggled (bool aChecked) { mBlue_active = !aChecked ; }
  
signals:

  void changed();
  
private:
  
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

Circular_polygon linear_2_circ( Linear_polygon const& pgn)
{
  Circular_polygon rCP;
  
  for( Linear_polygon::Edge_const_iterator ei = pgn.edges_begin(); ei != pgn.edges_end(); ++ei )
    rCP.push_back( Circular_X_monotone_curve(ei->source(), ei->target()) );

  return rCP;
}

bool read_linear ( QString aFileName, Circular_polygon_set& rSet )
{
  bool rOK = false ;
  
  std::ifstream in_file (qPrintable(aFileName));

  if ( in_file )
  {
    // Red the number of bezier polygon with holes
    unsigned int n_regions ;
    in_file >> n_regions;
    Circular_polygon outer ;
    std::vector<Circular_polygon> holes ;
    
    for ( unsigned int r = 0 ; r < n_regions ; ++ r )
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
  
  return rOK ;
}

bool read_dxf ( QString aFileName, Circular_polygon_set& rSet )
{
  bool rOK = false ;
  
  std::ifstream in_file (qPrintable(aFileName));

  if ( in_file )
  {
    CGAL::Dxf_bsop_reader<Kernel>            reader;
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

bool read_bezier ( QString aFileName, Bezier_polygon_set& rSet )
{
  bool rOK = false ;
  
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
      Bezier_traits                    traits;
      Bezier_traits::Make_x_monotone_2 make_x_monotone = traits.make_x_monotone_2_object();
      
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
      rSet.join( Bezier_polygon_with_holes(polygons.front(), pit, polygons.end()) ) ;
      
      rOK = true ;
    }
    
  }
  
  return rOK ;
}

void MainWindow::on_actionOpenLinear_triggered()
{
  open(QFileDialog::getOpenFileName(this, tr("Open Linear Polygon"), "../data", tr("Linear Curve files (*.poly)") ));
}

void MainWindow::on_actionOpenDXF_triggered()
{
  open(QFileDialog::getOpenFileName(this, tr("Open DXF"), "../data", tr("DXF files (*.dxf)") ));
}

void MainWindow::on_actionOpenBezier_triggered()
{
  open(QFileDialog::getOpenFileName(this, tr("Open Bezier Polygon"), "../data", tr("Bezier Curve files (*.bps)") ));
}

void MainWindow::switch_set_type( Curve_set& aSet, int aType )
{
  unlink_GI( aSet.gi() ) ;
  
  aSet.reset_type(aType);
  
  link_GI( aSet.gi() ) ;
  
  emit(changed());
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
    
    if(fileName.endsWith(".poly"))
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
      emit(changed());
      zoomToFit();
      this->addToRecentFiles(fileName);
    }
  }  
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
  
  if ( !blue_set().is_empty() && !red_set().is_empty() )
  {
    result_set().assign( red_set() ) ;
    result_set().join(blue_set());
    lDone = true ;
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
  
  if ( !blue_set().is_empty() && !red_set().is_empty() )
  {
    result_set().assign( blue_set() ) ;
    result_set().difference(red_set());
    lDone = true ;
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
  
  if ( !blue_set().is_empty() && !red_set().is_empty() )
  {
    result_set().assign( red_set() ) ;
    result_set().difference(blue_set());
    lDone = true ;
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
  
  if ( !blue_set().is_empty() && !red_set().is_empty() )
  {
    result_set().assign( red_set() ) ;
    result_set().symmetric_difference(blue_set());
    lDone = true ;
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
  
  if ( !blue_set().is_empty() )
  {
    result_set().assign( blue_set() ) ;
    result_set().complement();
    lDone = true ;
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
  
  if ( !red_set().is_empty() )
  {
    result_set().assign( red_set() ) ;
    result_set().complement();
    lDone = true ;
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
    SetViewBlue  (true);
    SetViewRed   (false);
    SetViewResult(true);
    
    emit(changed());
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
    SetViewBlue  (false);
    SetViewRed   (true);
    SetViewResult(true);
    
    emit(changed());
  }
}
void MainWindow::on_actionDeleteBlue_triggered()
{
  blue_set().clear();
    
  SetViewBlue  (true);
  SetViewRed   (true);
  SetViewResult(true);
  
  emit(changed());
}

void MainWindow::on_actionDeleteRed_triggered()
{
  red_set().clear();
    
  SetViewBlue  (true);
  SetViewRed   (true);
  SetViewResult(true);
  
  emit(changed());
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


