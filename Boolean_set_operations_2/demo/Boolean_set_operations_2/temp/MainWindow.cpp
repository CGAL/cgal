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
#include <CGAL/Arr_Bezier_curve_traits_2.h>
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

//#include <QT5/BezierCurves.h>
#include <QT5/CircularPolygons.h>
#include <QT5/LinearPolygons.h>
//#include <QT5/GraphicsViewBezierPolygonInput.h>
#include <QT5/GraphicsViewCircularPolygonInput.h>
#include <QT5/GraphicsViewLinearPolygonInput.h>
//#include <CGAL/Qt/GraphicsViewGpsCircleInput.h>

//Boundary_pieces_graphics_item
#include <QT5/PiecewiseGraphicsItemBase.h>
#include <QT5/PiecewiseBoundaryGraphicsItem.h>
#include <QT5/PiecewiseRegionGraphicsItem.h>

#include <QT5/PiecewiseSetGraphicsItem.h>

#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/DemosMainWindow.h>
#include <CGAL/Qt/utility.h>
#include <CGAL/IO/Dxf_bsop_reader.h>
#include "Typedefs.h"

#include <CGAL/Qt/GraphicsViewNavigation.h>

#include <QFileDialog>

//#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "MainWindow.moc"

#include "Typedefs.h"

#define CGAL_POLYGON_EDGE_DEFAULT_COLOR "blue"
#define CGAL_POLYGON_FILL_DEFAULT_COLOR "blue"

using namespace std;

typedef CGAL::Qt::Circular_set_graphics_item<Circular_polygon_set> Circular_GI;
typedef CGAL::Qt::Linear_set_graphics_item<Linear_polygon_set>     Linear_GI;

void show_error( std::string aS )
{
  QMessageBox::critical(NULL,"Critical Error",QString(aS.c_str()) ) ;
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
  
  Linear_rep () : Base() {} 
  
  virtual int type() const { return LINEAR_TYPE ; }
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
                                  : Rep_ptr(new Linear_rep  ()) ;
         
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
    else
    {
      get_linear_rep()->assign( *aOther.get_linear_rep() ) ;
    }  
  }
  
  void intersect( Curve_set const& aOther ) 
  {
    if ( is_circular() && aOther.is_circular() )
    {
      get_circular_rep()->intersect( *aOther.get_circular_rep() ) ;
    }
    else
    {
      get_linear_rep()->intersect( *aOther.get_linear_rep() ) ;
    }  
  }
  
  void join ( Curve_set const& aOther ) 
  {
    if ( is_circular() && aOther.is_circular() )
    {
      get_circular_rep()->join( *aOther.get_circular_rep() ) ;
    }
    else
    {
      get_linear_rep()->join( *aOther.get_linear_rep() ) ;
    }  
  }
  
  void difference( Curve_set const& aOther ) 
  {
    if ( is_circular() && aOther.is_circular() )
    {
      get_circular_rep()->difference( *aOther.get_circular_rep() ) ;
    }
    else
    {
      get_linear_rep()->difference( *aOther.get_linear_rep() ) ;
    }  
  }
  
  void symmetric_difference( Curve_set const& aOther ) 
  {
    if ( is_circular() && aOther.is_circular() )
    {
      get_circular_rep()->symmetric_difference( *aOther.get_circular_rep() ) ;
    }
    else
    {
      get_linear_rep()->symmetric_difference( *aOther.get_linear_rep() ) ;
    }  
  }
   
  Rep_base const& rep() const { return *mRep ; }
  Rep_base&       rep()       { return *mRep ; }
  
  bool is_circular() const { return mRep->type() == CIRCULAR_TYPE ; }  
  bool is_linear  () const { return mRep->type() == LINEAR_TYPE ; }  
  
  Circular_rep const* get_circular_rep() const { return dynamic_cast<Circular_rep const*>( boost::get_pointer(mRep) ); }
  Circular_rep      * get_circular_rep()       { return dynamic_cast<Circular_rep*      >( boost::get_pointer(mRep) ); }
  Linear_rep   const* get_linear_rep  () const { return dynamic_cast<Linear_rep   const*>( boost::get_pointer(mRep) ); }
  Linear_rep        * get_linear_rep  ()       { return dynamic_cast<Linear_rep  *      >( boost::get_pointer(mRep) ); }
  
  Circular_polygon_set const& circular() const { return get_circular_rep()->set(); }
  Circular_polygon_set      & circular()       { return get_circular_rep()->set(); }
  Linear_polygon_set   const& linear  () const { return get_linear_rep  ()->set(); }
  Linear_polygon_set        & linear  ()       { return get_linear_rep  ()->set(); }
  
private:

  QPen                        mPen ;
  QBrush                      mBrush ;
  boost::shared_ptr<Rep_base> mRep ;
  
} ;

typedef std::vector<Curve_set> Curve_set_container ;

typedef Curve_set_container::const_iterator Curve_set_const_iterator ;
typedef Curve_set_container::iterator       Curve_set_iterator ;


class MainWindow :  public CGAL::Qt::DemosMainWindow,  public Ui::Boolean_operations_2
{
  Q_OBJECT// removing it gives error
  
private:  

  QGraphicsScene                                                   mScene;
  bool                                                             mCircular_active ;
  bool                                                             mBlue_active ;
  Curve_set_container                                              mCurve_sets ;
  Circular_region_source_container                                 mBlue_circular_sources ;
  Circular_region_source_container                                 mRed_circular_sources ;
  Linear_region_source_container                                   mBlue_linear_sources ; 
  Linear_region_source_container                                   mRed_linear_sources ; 
  CGAL::Qt::GraphicsViewLinearPolygonInput<Linear_traits>*         mLinearInput ;
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
  void on_actionOpenLinear_triggered() ;
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

  void on_actionInsertLinear_toggled  (bool aChecked);
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

  Linear_region_source_container const& blue_linear_sources() const { return mBlue_linear_sources ; }
  Linear_region_source_container      & blue_linear_sources()       { return mBlue_linear_sources ; }

  Linear_region_source_container const& red_linear_sources () const { return mRed_linear_sources ; }
  Linear_region_source_container      & red_linear_sources ()       { return mRed_linear_sources ; }

  Linear_region_source_container const& active_linear_sources() const { return mBlue_active ? mBlue_linear_sources : mRed_linear_sources ; }
  Linear_region_source_container      & active_linear_sources()       { return mBlue_active ? mBlue_linear_sources : mRed_linear_sources ; }

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
  
  bool ensure_linear_mode();//see its need
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
  
  mLinearInput   = new CGAL::Qt::GraphicsViewLinearPolygonInput  <Linear_traits>      (this, &mScene);
  mCircularInput = new CGAL::Qt::GraphicsViewCircularPolygonInput<Gps_circular_kernel>(this, &mScene);
  //mCircleInput   = new CGAL::Qt::GraphicsViewCircleInput       <Circular_traits>(this, &mScene);
  
  QObject::connect(mLinearInput  , SIGNAL(generate(CGAL::Object)), this, SLOT(processInput(CGAL::Object)));
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
  blue_linear_sources  ().clear();
  red_circular_sources ().clear();
  red_linear_sources   ().clear();
    
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

Linear_curve read_linear_curve ( std::istream& is, bool aDoubleFormat )
{
  // Read the number of control points.
  unsigned int  n;

  is >> n;
  
  // Read the control points.
  std::vector<Linear_rat_point> ctrl_pts;
  
  for ( unsigned int k = 0; k < n; k++)
  {
    Linear_rat_point p ;
    if ( aDoubleFormat )
    {
      double x,y ;
      is >> x >> y ;
      Linear_rational rx(static_cast<int> (1000 * x + 0.5), 1000);
      Linear_rational ry(static_cast<int> (1000 * y + 0.5), 1000); 
      p = Linear_rat_point(rx,ry);
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

  std::vector<Linear_rat_point> ctrl_pts2;

  typedef std::vector<Linear_rat_point>::const_iterator cp_const_iterator ;

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

  return Linear_curve(ctrl_pts2.begin(),ctrl_pts2.end());
}

bool save_circular ( QString aFileName, Circular_polygon_set& rSet )
{
  bool rOK = false ;
  
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
//check out
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
       
  }/*
  else
  {
    if ( !save_bezier_sources(QFileDialog::getSaveFileName(this, tr("Save 'Q' Bezier Polygon Set"), "../data", tr("Bezier Curve files (*.bps)") )
                             ,blue_bezier_sources() 
                             )
       )
    {
      show_error("Cannot save bezier polygon set.");
    }
  }*/

}
///check out
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
       
  }/*
  else
  {
    if ( !save_bezier_sources(QFileDialog::getSaveFileName(this, tr("Save 'P' Bezier Polygon Set"), "../data", tr("Bezier Curve files (*.bps)") )
                             ,red_bezier_sources() 
                             )
       )
    {
      show_error("Cannot save bezier polygon set.");
    }
  }*/

}

//check out
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
       
  }/*
  else
  {
    if ( !save_bezier_result(QFileDialog::getSaveFileName(this, tr("Save Result Bezier Polygon Set"), "../data", tr("Bezier Curve files (*.bps)") )
                            ,result_set().bezier() 
                            )
       )
    {
      show_error("Cannot save bezier polygon set.");
    }
  }*/

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
//check out/*
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
*/
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
    /*else if (fileName.endsWith(".bps"))
    {
      if ( ensure_bezier_mode() )
        lRead = read_bezier(fileName,active_set().bezier(), active_bezier_sources() ) ;
    }*/
     
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
  /*
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
  }*/
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
