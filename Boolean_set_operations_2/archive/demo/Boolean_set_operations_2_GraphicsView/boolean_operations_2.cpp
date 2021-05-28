// Copyright (c) 2009  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
#include <iomanip>
#include <list>

void trace( std::string s )
{
  static std::ofstream out("log.txt");
  out << s ;
}

#include <memory>

#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QSlider>
#include <QProgressBar>
#include <QMessageBox>

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
#include <CGAL/Qt/GraphicsViewCircularPolygonInput.h>
//#include <CGAL/Qt/GraphicsViewGpsCircleInput.h>
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
  typedef std::shared_ptr<Rep_base> Rep_ptr ;

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
  std::shared_ptr<Rep_base> mRep ;

} ;

typedef std::vector<Curve_set> Curve_set_container ;

typedef Curve_set_container::const_iterator Curve_set_const_iterator ;
typedef Curve_set_container::iterator       Curve_set_iterator ;


class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Boolean_operations_2
{
  Q_OBJECT

private:

  QGraphicsScene                                                   mScene;
  bool                                                             mCircular_active ;
  bool                                                             mBlue_active ;
  Curve_set_container                                              mCurve_sets ;
  Circular_region_source_container                                 mBlue_circular_sources ;
  Circular_region_source_container                                 mRed_circular_sources ;
  Bezier_region_source_container                                   mBlue_bezier_sources ;
  Bezier_region_source_container                                   mRed_bezier_sources ;
  CGAL::Qt::GraphicsViewBezierPolygonInput<Bezier_traits>*         mBezierInput ;
  CGAL::Qt::GraphicsViewCircularPolygonInput<Gps_circular_kernel>* mCircularInput ;
  //CGAL::Qt::GraphicsViewGpsCircleSegmentInput<Circular_curve>* mCircularInput ;
  //CGAL::Qt::GraphicsViewGpsCircleInput<Circular_traits>*      mCircleInput ;

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
  void on_actionSaveBlue_triggered() ;
  void on_actionSaveRed_triggered() ;
  void on_actionSaveResult_triggered() ;
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
  void on_actionRecenter_triggered();

  void on_actionInsertBezier_toggled  (bool aChecked);
  void on_actionInsertCircular_toggled(bool aChecked);
  void on_actionInsertCircle_toggled  (bool aChecked);

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

  Circular_region_source_container const& blue_circular_sources() const { return mBlue_circular_sources ; }
  Circular_region_source_container      & blue_circular_sources()       { return mBlue_circular_sources ; }

  Circular_region_source_container const& red_circular_sources () const { return mRed_circular_sources ; }
  Circular_region_source_container      & red_circular_sources ()       { return mRed_circular_sources ; }

  Bezier_region_source_container const& blue_bezier_sources() const { return mBlue_bezier_sources ; }
  Bezier_region_source_container      & blue_bezier_sources()       { return mBlue_bezier_sources ; }

  Bezier_region_source_container const& red_bezier_sources () const { return mRed_bezier_sources ; }
  Bezier_region_source_container      & red_bezier_sources ()       { return mRed_bezier_sources ; }

  Bezier_region_source_container const& active_bezier_sources() const { return mBlue_active ? mBlue_bezier_sources : mRed_bezier_sources ; }
  Bezier_region_source_container      & active_bezier_sources()       { return mBlue_active ? mBlue_bezier_sources : mRed_bezier_sources ; }

  Circular_region_source_container const& active_circular_sources() const { return mBlue_active ? mBlue_circular_sources : mRed_circular_sources ; }
  Circular_region_source_container      & active_circular_sources()       { return mBlue_active ? mBlue_circular_sources : mRed_circular_sources ; }

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

  mBezierInput   = new CGAL::Qt::GraphicsViewBezierPolygonInput  <Bezier_traits>      (this, &mScene);
  mCircularInput = new CGAL::Qt::GraphicsViewCircularPolygonInput<Gps_circular_kernel>(this, &mScene);
  //mCircleInput   = new CGAL::Qt::GraphicsViewCircleInput       <Circular_traits>(this, &mScene);

  QObject::connect(mBezierInput  , SIGNAL(generate(CGAL::Object)), this, SLOT(processInput(CGAL::Object)));
  QObject::connect(mCircularInput, SIGNAL(generate(CGAL::Object)), this, SLOT(processInput(CGAL::Object)));
  //QObject::connect(mCircleInput  , SIGNAL(generate(CGAL::Object)), this, SLOT(processInput(CGAL::Object)));

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

  blue_circular_sources().clear();
  blue_bezier_sources  ().clear();
  red_circular_sources ().clear();
  red_bezier_sources   ().clear();

  SetViewBlue  (true);
  SetViewRed   (true);
  SetViewResult(true);

  mCircular_active = true ;

  radioMakeBlueActive->setChecked(true);

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

Circular_polygon linear_2_circ( Linear_polygon const& pgn )
{
  CGAL::Cartesian_converter<Linear_kernel,Gps_circular_kernel> convert ;

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

bool read_linear ( QString aFileName, Circular_polygon_set& rSet, Circular_region_source_container& rSources )
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

      Circular_polygon_with_holes pwh(outer,holes.begin(),holes.end());
      rSources.push_back(pwh);
      rSet.join(pwh) ;
      rOK = true ;
    }

  }

  return rOK ;
}

bool read_dxf ( QString aFileName, Circular_polygon_set& rSet, Circular_region_source_container& rSources )
{
  bool rOK = false ;

  std::ifstream in_file (qPrintable(aFileName));

  if ( in_file )
  {
    CGAL::Dxf_bsop_reader<Gps_circular_kernel>   reader;
    std::vector<Circular_polygon>            circ_polygons;
    std::vector<Circular_polygon_with_holes> circ_polygons_with_holes;

    reader(in_file
          ,std::back_inserter(circ_polygons)
          ,std::back_inserter(circ_polygons_with_holes)
          ,false
          );

    for ( std::vector<Circular_polygon>::iterator pit = circ_polygons.begin() ; pit != circ_polygons.end() ; ++ pit )
      circ_polygons_with_holes.push_back( Circular_polygon_with_holes(*pit) ) ;

    rSet.join( circ_polygons_with_holes.begin(), circ_polygons_with_holes.end() ) ;

    std::copy(circ_polygons_with_holes.begin(), circ_polygons_with_holes.end(), std::back_inserter(rSources) );

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
      Bezier_rational rx(static_cast<int> (1000 * x + 0.5), 1000);
      Bezier_rational ry(static_cast<int> (1000 * y + 0.5), 1000);
      p = Bezier_rat_point(rx,ry);
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

  std::vector<Bezier_rat_point> ctrl_pts2;

  typedef std::vector<Bezier_rat_point>::const_iterator cp_const_iterator ;

  cp_const_iterator beg  = ctrl_pts.begin();
  cp_const_iterator end  = ctrl_pts.end  ();
  cp_const_iterator last = end - 1 ;

  ctrl_pts2.push_back(*beg);

  if ( ctrl_pts.size() > 2 )
  {
    cp_const_iterator curr = beg ;
    cp_const_iterator next1 = curr  + 1 ;
    cp_const_iterator next2 = next1 + 1 ;

    do
    {
      CGAL::Orientation lOrient = orientation(*curr,*next1,*next2);

      if ( lOrient != CGAL::COLLINEAR )
        ctrl_pts2.push_back(*next1);

      ++ curr  ;
      ++ next1 ;
      ++ next2 ;

    }
    while ( next2 != end ) ;
  }

  ctrl_pts2.push_back(*last);

  return Bezier_curve(ctrl_pts2.begin(),ctrl_pts2.end());
}

bool read_bezier ( QString aFileName, Bezier_polygon_set& rSet, Bezier_region_source_container& rSources  )
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
        Bezier_polygon_vector bezier_polygons ;
        Bezier_region_source  br_source ;

        // Read the number of bezier curves.
        unsigned int n_boundaries;
        in_file >> n_boundaries;

        for ( unsigned int b = 0 ; b < n_boundaries ; ++ b )
        {
          Bezier_boundary_source bb_source ;

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

            Bezier_curve b = read_bezier_curve(in_file, lDoubleFormat);

            if ( b.number_of_control_points() >= 2 )
            {
              bb_source.push_back(b);
              //TRACE( "region " << r << " boundary " << b << " curve " << k );

              make_x_monotone (b, std::back_inserter (x_objs));

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

          br_source.push_back(bb_source);
          bezier_polygons.push_back (pgn);
        }

        if ( bezier_polygons.size() > 0 )
        {
          Bezier_polygon_with_holes pwh(bezier_polygons.front());

          if ( bezier_polygons.size() > 1 )
          {
            for ( Bezier_polygon_vector::const_iterator it = std::next(bezier_polygons.begin())
                ; it != bezier_polygons.end()
                ; ++ it
                )
              pwh.add_hole(*it);
          }

          if ( is_valid_polygon_with_holes(pwh, rSet.traits() ) )
          {
            rSet.join(pwh) ;
            rSources.push_back(br_source);
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
  typedef std::vector<Bezier_rat_point> Bezier_rat_point_vector ;

  int cc = aBP.size() ;
  int lc = cc - 1 ;

  out_file << "  " <<  cc << std::endl ;

  Bezier_rat_point lFirstP, lPrevP ;

  int i = 0 ;

  for ( Bezier_polygon::Curve_const_iterator cit = aBP.curves_begin() ; cit != aBP.curves_end() ; ++ cit, ++ i  )
  {
    Bezier_rat_point_vector lQ ;

    CGAL::Qt::Bezier_helper::clip(*cit,lQ);

    out_file << "   " << lQ.size() << std::endl ;

    if ( i == 0 )
      lFirstP = lQ.front();

    if ( i == lc )
      lQ.back() = lFirstP ;

    for ( Bezier_rat_point_vector::const_iterator pit = lQ.begin() ; pit != lQ.end() ; ++ pit )
    {
      Bezier_rat_point lP = pit == lQ.begin() && i > 0 ? lPrevP : *pit ;

      out_file << "    " << CGAL::to_double(lP.x()) << " " << CGAL::to_double(lP.y()) << std::endl ;

      lPrevP = lP ;
    }
  }
}

bool save_bezier_result ( QString aFileName, Bezier_polygon_set const& aSet )
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

      out_file << " " << ( 1 + bpwh.number_of_holes() ) << std::endl ;

      save_bezier_polygon( out_file, bpwh.outer_boundary() ) ;

      for ( Bezier_polygon_with_holes::Hole_const_iterator hit = bpwh.holes_begin() ; hit != bpwh.holes_end() ; ++ hit )
        save_bezier_polygon(out_file, *hit);

      rOK = true ;
    }
  }

  return rOK ;

}

bool save_bezier_sources ( QString aFileName, Bezier_region_source_container const& aSources )
{
  bool rOK = false ;

  std::ofstream out_file( qPrintable(aFileName) ) ;
  if ( out_file )
  {
    out_file << std::setprecision(19);

    out_file << "DOUBLE" << std::endl ;

    out_file << aSources.size() << std::endl ;

    for( Bezier_region_source_container::const_iterator rit = aSources.begin(); rit != aSources.end() ; ++ rit )
    {
      Bezier_region_source const& br = *rit ;

      out_file << "  " << br.size() << std::endl ;

      for( Bezier_region_source::const_iterator bit = br.begin(); bit != br.end() ; ++ bit )
      {
        Bezier_boundary_source const& bb = *bit ;

        out_file << "   " << bb.size() << std::endl ;

        for ( Bezier_boundary_source::const_iterator cit = bb.begin() ; cit != bb.end() ; ++ cit )
        {
          Bezier_curve const& bc = *cit ;

          out_file << "    " << bc.number_of_control_points() << std::endl ;

          for ( Bezier_curve::Control_point_iterator pit = bc.control_points_begin() ; pit != bc.control_points_end() ; ++ pit )
          {
            out_file << "     " << CGAL::to_double(pit->x()) << " " << CGAL::to_double(pit->y()) << std::endl ;
          }
        }
      }
    }

    rOK = true ;
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
void MainWindow::on_actionSaveBlue_triggered()
{
  if ( mCircular_active )
  {
    if ( !save_circular(QFileDialog::getSaveFileName(this, tr("Save 'Q' Circular Polygon Set"), "../data", tr("Linear Curve files (*.lps)") )
                       ,active_set().circular()
                       )
       )
    {
      show_error("Cannot save circular polygon set.");
    }

  }
  else
  {
    if ( !save_bezier_sources(QFileDialog::getSaveFileName(this, tr("Save 'Q' Bezier Polygon Set"), "../data", tr("Bezier Curve files (*.bps)") )
                             ,blue_bezier_sources()
                             )
       )
    {
      show_error("Cannot save bezier polygon set.");
    }
  }

}

void MainWindow::on_actionSaveRed_triggered()
{
  if ( mCircular_active )
  {
    if ( !save_circular(QFileDialog::getSaveFileName(this, tr("Save 'P' Circular Polygon Set"), "../data", tr("Linear Curve files (*.lps)") )
                       ,red_set().circular()
                       )
       )
    {
      show_error("Cannot save circular polygon set.");
    }

  }
  else
  {
    if ( !save_bezier_sources(QFileDialog::getSaveFileName(this, tr("Save 'P' Bezier Polygon Set"), "../data", tr("Bezier Curve files (*.bps)") )
                             ,red_bezier_sources()
                             )
       )
    {
      show_error("Cannot save bezier polygon set.");
    }
  }

}


void MainWindow::on_actionSaveResult_triggered()
{
  if ( mCircular_active )
  {
    if ( !save_circular(QFileDialog::getSaveFileName(this, tr("Save Result Circular Polygon Set"), "../data", tr("Linear Curve files (*.lps)") )
                       ,result_set().circular()
                       )
       )
    {
      show_error("Cannot save circular polygon set.");
    }

  }
  else
  {
    if ( !save_bezier_result(QFileDialog::getSaveFileName(this, tr("Save Result Bezier Polygon Set"), "../data", tr("Bezier Curve files (*.bps)") )
                            ,result_set().bezier()
                            )
       )
    {
      show_error("Cannot save bezier polygon set.");
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
        lRead = read_linear(fileName,active_set().circular(), active_circular_sources() ) ;
    }
    else if (fileName.endsWith(".dxf"))
    {
      if ( ensure_circular_mode() )
        lRead = read_dxf(fileName,active_set().circular(), active_circular_sources() ) ;
    }
    else if (fileName.endsWith(".bps"))
    {
      if ( ensure_bezier_mode() )
        lRead = read_bezier(fileName,active_set().bezier(), active_bezier_sources() ) ;
    }

    if ( lRead )
    {
      modelChanged();
      zoomToFit();
      this->addToRecentFiles(fileName);

    }
  }
}

void MainWindow::on_actionInsertBezier_toggled(bool aChecked)
{
  if(aChecked)
       mScene.installEventFilter(mBezierInput);
  else mScene.removeEventFilter (mBezierInput);
}

void MainWindow::on_actionInsertCircular_toggled(bool aChecked)
{
  if(aChecked)
       mScene.installEventFilter(mCircularInput);
  else mScene.removeEventFilter (mCircularInput);
}

void MainWindow::on_actionInsertCircle_toggled(bool aChecked)
{
//  if(aChecked)
//       mScene.installEventFilter(mCircleInput);
//  else mScene.removeEventFilter (mCircleInput);
}

void MainWindow::processInput(CGAL::Object o )
{
  std::pair<Bezier_polygon,Bezier_boundary_source>     lBI ;
  Circular_polygon lCI ;

  if(CGAL::assign(lBI, o))
  {
    if ( ensure_bezier_mode() )
    {
      CGAL::Orientation o = lBI.first.orientation();
      if ( o == CGAL::CLOCKWISE )
        lBI.first.reverse_orientation();

      active_set().bezier().join( Bezier_polygon_with_holes(lBI.first) ) ;

      Bezier_region_source br ; br.push_back (lBI.second);

      active_bezier_sources().push_back(br);

    }
  }
  else if ( CGAL::assign(lCI, o) )
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
  blue_set             ().clear();
  blue_circular_sources().clear();
  blue_bezier_sources  ().clear();

  //SetViewBlue(true);SetViewRed(true); SetViewResult(true);

  modelChanged();
}

void MainWindow::on_actionDeleteRed_triggered()
{
  red_set             ().clear();
  red_circular_sources().clear();
  red_bezier_sources  ().clear();

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
#include <CGAL/Qt/resources.h>

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Boolean_operations_2 demo");

  // Import resources from libCGALQt5.
  CGAL_QT_INIT_RESOURCES;

  MainWindow mainWindow;
  mainWindow.show();
  return app.exec();
}


