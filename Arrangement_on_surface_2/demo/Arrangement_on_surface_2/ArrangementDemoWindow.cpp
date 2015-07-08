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
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#include "ArrangementDemoWindow.h"
#include "NewTabDialog.h"
#include "OverlayDialog.h"
#include "ArrangementDemoPropertiesDialog.h"
#include "ArrangementDemoTab.h"
#include "Conic_reader.h"

#include "DeleteCurveMode.h"
#include "ArrangementGraphicsItem.h"

#include <boost/math/special_functions/fpclassify.hpp>

#include <QActionGroup>
#include <QFileDialog>
#include <QMessageBox>
#include <QColorDialog>

#include <CGAL/IO/Arr_with_history_iostream.h>
#include <CGAL/IO/Arr_text_formatter.h>
#include <CGAL/IO/Arr_with_history_text_formatter.h>

ArrangementDemoWindow::ArrangementDemoWindow(QWidget* parent) :
  CGAL::Qt::DemosMainWindow( parent ),
  lastTabIndex(static_cast<unsigned int>(-1)),
  ui( new Ui::ArrangementDemoWindow )
{
  this->setupUi( );

  // set up the demo window
  // ArrangementDemoTabBase* demoTab =
  this->makeTab( SEGMENT_TRAITS );
  this->setupStatusBar( );
  this->setupOptionsMenu( );
  this->addAboutDemo( ":/help/about.html" );
  this->addAboutCGAL( );

  // set up callbacks
  QObject::connect( this->modeGroup, SIGNAL( triggered( QAction* ) ),
                    this, SLOT( updateMode( QAction* ) ) );
  QObject::connect( this->envelopeGroup, SIGNAL( triggered( QAction* ) ),
                    this, SLOT( updateEnvelope( QAction* ) ) );
  QObject::connect( this->snapGroup, SIGNAL( triggered( QAction* ) ),
                    this, SLOT( updateSnapping( QAction* ) ) );

  QObject::connect( this->conicTypeGroup, SIGNAL( triggered( QAction* ) ),
                    this, SLOT( updateConicType( QAction* ) ) );
}

ArrangementDemoWindow::~ArrangementDemoWindow() {}

ArrangementDemoTabBase* ArrangementDemoWindow::makeTab( TraitsType tt )
{
  static int tabLabelCounter = 1;
  QString tabLabel;

  ArrangementDemoTabBase* demoTab;
  Seg_arr* seg_arr;
  Pol_arr* pol_arr;

#ifdef CGAL_USE_CORE
  Conic_arr* conic_arr;
#endif

  Lin_arr* lin_arr;
  Arc_arr* arc_arr;
  // Alg_seg_arr* alg_seg_arr;
  CGAL::Object arr;

  switch ( tt )
  {
   default:
   case SEGMENT_TRAITS:
    seg_arr = new Seg_arr;
    demoTab = new ArrangementDemoTab< Seg_arr >( seg_arr, 0 );
    arr = CGAL::make_object( seg_arr );
    tabLabel = QString( "%1 - Segment" ).arg( tabLabelCounter++ );
    break;
   case POLYLINE_TRAITS:
    pol_arr = new Pol_arr;
    demoTab = new ArrangementDemoTab< Pol_arr >( pol_arr, 0 );
    arr = CGAL::make_object( pol_arr );
    tabLabel = QString( "%1 - Polyline" ).arg( tabLabelCounter++ );
    break;

#ifdef CGAL_USE_CORE
   case CONIC_TRAITS:
    conic_arr = new Conic_arr;
    demoTab = new ArrangementDemoTab< Conic_arr >( conic_arr, 0 );
    arr = CGAL::make_object( conic_arr );
    tabLabel = QString( "%1 - Conic" ).arg( tabLabelCounter++ );
    break;
#endif

   case LINEAR_TRAITS:
    lin_arr = new Lin_arr;
    demoTab = new ArrangementDemoTab< Lin_arr >( lin_arr, 0 );
    arr = CGAL::make_object( lin_arr );
    tabLabel = QString( "%1 - Linear" ).arg( tabLabelCounter++ );
    break;
   case CIRCULAR_ARC_TRAITS:
    arc_arr = new Arc_arr;
    demoTab = new ArrangementDemoTab< Arc_arr >( arc_arr, 0 );
    arr = CGAL::make_object( arc_arr );
    tabLabel = QString( "%1 - Circular Arc" ).arg( tabLabelCounter++ );
    break;
   // case ALGEBRAIC_TRAITS:
   //  alg_seg_arr = new Alg_seg_arr;
   //  demoTab = new ArrangementDemoTab< Alg_seg_arr >( alg_seg_arr, 0 );
   //  arr = CGAL::make_object( alg_seg_arr );
   //  tabLabel = QString( "%1 - Algebraic" ).arg( tabLabelCounter++ );
   //  break;
  }

  this->arrangements.push_back( arr );
  this->tabs.push_back( demoTab );

  QGraphicsView* view = demoTab->getView( );
  this->addNavigation( view );
  this->ui->tabWidget->addTab( demoTab, tabLabel );
  this->lastTabIndex = this->ui->tabWidget->currentIndex( );
  this->ui->tabWidget->setCurrentWidget( demoTab );

  this->resetCallbackState( this->ui->tabWidget->currentIndex( ) );
  this->removeCallback( this->ui->tabWidget->currentIndex( ) );
  this->updateMode( this->modeGroup->checkedAction( ) );
  this->updateFillColorSwatch( );

  return demoTab;
}

ArrangementDemoTabBase* ArrangementDemoWindow::getTab( unsigned int tabIndex )
  const
{
  if (tabIndex > this->tabs.size()) return NULL;
  return this->tabs[tabIndex];
}

ArrangementDemoTabBase* ArrangementDemoWindow::getCurrentTab( ) const
{
  int currentIndex = this->ui->tabWidget->currentIndex( );
  if ( currentIndex == -1 )
    return NULL;

  ArrangementDemoTabBase* res = this->tabs[ currentIndex ];
  return res;
}

std::vector< QString > ArrangementDemoWindow::getTabLabels( ) const
{
  std::vector< QString > res;
  for ( int i = 0; i < this->ui->tabWidget->count( ); ++i )
  {
    res.push_back( this->ui->tabWidget->tabText( i ) );
  }
  return res;
}

std::vector< CGAL::Object > ArrangementDemoWindow::getArrangements( ) const
{
  std::vector< CGAL::Object > res;
  for ( unsigned int i = 0; i < this->arrangements.size( ); ++i )
  {
    res.push_back( this->arrangements[ i ] );
  }
  return res;
}

void ArrangementDemoWindow::setupUi( )
{
  this->ui->setupUi( this );

  this->modeGroup = new QActionGroup( this );
  this->modeGroup->addAction( this->ui->actionDrag );
  this->modeGroup->addAction( this->ui->actionInsert );
  this->modeGroup->addAction( this->ui->actionDelete );
  this->modeGroup->addAction( this->ui->actionPointLocation );
  this->modeGroup->addAction( this->ui->actionRayShootingUp );
  this->modeGroup->addAction( this->ui->actionRayShootingDown );
  this->modeGroup->addAction( this->ui->actionMerge );
  this->modeGroup->addAction( this->ui->actionSplit );
  this->modeGroup->addAction( this->ui->actionFill );
  this->activeModes.push_back( this->ui->actionInsert );

  this->envelopeGroup = new QActionGroup( this );
  this->envelopeGroup->addAction( this->ui->actionLowerEnvelope );
  this->envelopeGroup->addAction( this->ui->actionUpperEnvelope );
  this->envelopeGroup->setExclusive( false );

  this->snapGroup = new QActionGroup( this );
  this->snapGroup->addAction( this->ui->actionSnapMode );
  this->snapGroup->addAction( this->ui->actionGridSnapMode );
  this->snapGroup->setExclusive( false );
  this->ui->actionGridSnapMode->setEnabled( false );

  this->conicTypeGroup = new QActionGroup( this );
  this->conicTypeGroup->addAction( this->ui->actionConicSegment );
  this->conicTypeGroup->addAction( this->ui->actionConicCircle );
  this->conicTypeGroup->addAction( this->ui->actionConicEllipse );
  this->conicTypeGroup->addAction( this->ui->actionConicThreePoint );
  this->conicTypeGroup->addAction( this->ui->actionConicFivePoint );
  this->conicTypeGroup->addAction( this->ui->actionCurveRay );
  this->conicTypeGroup->addAction( this->ui->actionCurveLine );

  this->updateFillColorSwatch( );
}

void ArrangementDemoWindow::updateMode( QAction* newMode )
{
  // QWidget* widget = this->ui->tabWidget->currentWidget( );
  // ArrangementDemoTabBase* demoTab =
  //   static_cast< ArrangementDemoTabBase* >( widget );
  const unsigned int TabIndex = this->ui->tabWidget->currentIndex( );
  if (TabIndex == static_cast<unsigned int>(-1)) return;
  ArrangementDemoTabBase* activeTab = this->tabs[ TabIndex ];
  QGraphicsScene* activeScene = activeTab->getScene( );
  QGraphicsView* activeView = activeTab->getView( );

  this->resetCallbackState( TabIndex );
  this->removeCallback( TabIndex );

  // update the active mode
  this->activeModes.at( 0 ) = newMode;

  // hook up the new active mode
  if ( newMode == this->ui->actionInsert )
  {
    activeScene->installEventFilter( activeTab->getCurveInputCallback( ) );
  }
  else if ( newMode == this->ui->actionDrag )
  {
    activeView->setDragMode( QGraphicsView::ScrollHandDrag );
  }
  else if ( newMode == this->ui->actionDelete )
  {
    activeScene->installEventFilter( activeTab->getDeleteCurveCallback( ) );
  }
  else if ( newMode == this->ui->actionPointLocation )
  {
    activeScene->installEventFilter( activeTab->getPointLocationCallback( ) );
  }
  else if ( newMode == this->ui->actionRayShootingUp )
  {
    // -y is up for Qt, so we shoot down
    activeTab->getVerticalRayShootCallback( )->setShootingUp( true );
    activeScene->installEventFilter( activeTab->getVerticalRayShootCallback());
  }
  else if ( newMode == this->ui->actionRayShootingDown )
  {
    // the bottom of the viewport for Qt is +y, so we shoot up
    activeTab->getVerticalRayShootCallback( )->setShootingUp( false );
    activeScene->installEventFilter( activeTab->getVerticalRayShootCallback());
  }
  else if ( newMode == this->ui->actionMerge )
  {
    activeScene->installEventFilter( activeTab->getMergeEdgeCallback( ) );
  }
  else if ( newMode == this->ui->actionSplit )
  {
    activeScene->installEventFilter( activeTab->getSplitEdgeCallback( ) );
  }
  else if ( newMode == this->ui->actionFill )
  {
    activeScene->installEventFilter( activeTab->getFillFaceCallback( ) );
  }
  this->updateFillColorSwatch( );
}

void ArrangementDemoWindow::resetCallbackState( unsigned int tabIndex )
{
  if (tabIndex == static_cast<unsigned int>(-1)
      || tabIndex >= this->tabs.size( )) return;

  ArrangementDemoTabBase* activeTab = this->tabs[ tabIndex ];

  QAction* activeMode = this->activeModes.at( 0 );

  // unhook the old active mode
  if ( activeMode == this->ui->actionInsert )
  {  }
  else if ( activeMode == this->ui->actionDrag )
  {  }
  else if ( activeMode == this->ui->actionDelete )
  {
    activeTab->getDeleteCurveCallback( )->reset( );
  }
  else if ( activeMode == this->ui->actionPointLocation )
  {
    activeTab->getPointLocationCallback( )->reset( );
  }
  else if ( activeMode == this->ui->actionRayShootingUp )
  {
    activeTab->getVerticalRayShootCallback( )->reset( );
  }
  else if ( activeMode == this->ui->actionRayShootingDown )
  {
    activeTab->getVerticalRayShootCallback( )->reset( );
  }
  else if ( activeMode == this->ui->actionMerge )
  {
    activeTab->getMergeEdgeCallback( )->reset( );
  }
  else if ( activeMode == this->ui->actionSplit )
  {
    activeTab->getSplitEdgeCallback( )->reset( );
  }
  else if ( activeMode == this->ui->actionFill )
  {
    activeTab->getFillFaceCallback( )->reset( );
  }
}

void ArrangementDemoWindow::removeCallback( unsigned int tabIndex )
{
  if (tabIndex == static_cast<unsigned int>(-1)) return;

  ArrangementDemoTabBase* activeTab = this->tabs[ tabIndex ];
  QGraphicsScene* activeScene = activeTab->getScene( );
  QGraphicsView* activeView = activeTab->getView( );
#if 0
  QAction* activeMode = this->activeModes[ tabIndex ];
#endif

  activeScene->removeEventFilter( activeTab->getCurveInputCallback( ) );
  activeView->setDragMode( QGraphicsView::NoDrag );
  activeScene->removeEventFilter( activeTab->getDeleteCurveCallback( ) );
  activeScene->removeEventFilter( activeTab->getPointLocationCallback( ) );
  activeScene->removeEventFilter( activeTab->getVerticalRayShootCallback( ) );
  activeScene->removeEventFilter( activeTab->getMergeEdgeCallback( ) );
  activeScene->removeEventFilter( activeTab->getSplitEdgeCallback( ) );
  activeScene->removeEventFilter( activeTab->getVerticalRayShootCallback( ) );
  activeScene->removeEventFilter( activeTab->getFillFaceCallback( ) );

  // unhook the old active mode
#if 0
  if ( activeMode == this->ui->actionInsert )
  {
    activeScene->removeEventFilter( activeTab->getCurveInputCallback( ) );
  }
  else if ( activeMode == this->ui->actionDrag )
  {
    activeView->setDragMode( QGraphicsView::NoDrag );
  }
  else if ( activeMode == this->ui->actionDelete )
  {
    activeScene->removeEventFilter( activeTab->getDeleteCurveCallback( ) );
  }
  else if ( activeMode == this->ui->actionPointLocation )
  {
    activeScene->removeEventFilter( activeTab->getPointLocationCallback( ) );
  }
  else if ( activeMode == this->ui->actionRayShootingUp )
  {
    activeScene->removeEventFilter( activeTab->getVerticalRayShootCallback( ) );
  }
  else if ( activeMode == this->ui->actionRayShootingDown )
  {
    activeScene->removeEventFilter( activeTab->getVerticalRayShootCallback( ) );
  }
  else if ( activeMode == this->ui->actionMerge )
  {
    activeScene->removeEventFilter( activeTab->getMergeEdgeCallback( ) );
  }
  else if ( activeMode == this->ui->actionSplit )
  {
    activeScene->removeEventFilter( activeTab->getSplitEdgeCallback( ) );
  }
#endif
}

void ArrangementDemoWindow::updateFillColorSwatch( )
{
  unsigned int currentTabIndex = this->ui->tabWidget->currentIndex( );
  if (currentTabIndex == static_cast<unsigned int>(-1)) return;
  ArrangementDemoTabBase* currentTab = this->tabs[ currentTabIndex ];
  FillFaceCallbackBase* fillFaceCallback = currentTab->getFillFaceCallback( );
  QColor fillColor = fillFaceCallback->getColor( );
  if ( !fillColor.isValid( ) )
  {
    fillColor = ::Qt::black;
  }

  QPixmap fillColorPixmap( 16, 16 );
  fillColorPixmap.fill( fillColor );
  QIcon fillColorIcon( fillColorPixmap );
  this->ui->actionFillColor->setIcon( fillColorIcon );
}

void ArrangementDemoWindow::openArrFile( QString filename )
{
  int index = this->ui->tabWidget->currentIndex( );
  if ( index == -1 )
  {
    QMessageBox::information( this, "Oops", "Create a new tab first" );
    return;
  }
  if ( filename.isNull( ) )
  {
    return;
  }

  std::ifstream ifs( filename.toStdString( ).c_str( ) );
  CGAL::Object arr = this->arrangements[ index ];
  Seg_arr* seg;
  Pol_arr* pol;

#ifdef CGAL_USE_CORE
  Conic_arr* conic;
#endif

  // Alg_seg_arr* alg;
  if ( CGAL::assign( seg, arr ) )
  {
    typedef CGAL::Arr_text_formatter<Seg_arr>           Seg_text_formatter;
    typedef CGAL::Arr_with_history_text_formatter<Seg_text_formatter>
      ArrFormatter;
    typedef ArrangementDemoTab<Seg_arr>                 TabType;

    ArrFormatter arrFormatter;
    CGAL::read( *seg, ifs, arrFormatter );
    this->arrangements[ index ] = CGAL::make_object( seg );
    TabType* tab = static_cast< TabType* >( this->tabs[ index ] );
    tab->setArrangement( seg );
  }
  else if ( CGAL::assign( pol, arr ) )
  {
    typedef CGAL::Arr_text_formatter< Pol_arr >         Pol_text_formatter;
    typedef CGAL::Arr_with_history_text_formatter<Pol_text_formatter>
      ArrFormatter;
    typedef ArrangementDemoTab< Pol_arr >               TabType;

    ArrFormatter arrFormatter;
    CGAL::read( *pol, ifs, arrFormatter );
    this->arrangements[ index ] = CGAL::make_object( pol );
    TabType* tab = static_cast< TabType* >( this->tabs[ index ] );
    tab->setArrangement( pol );
  }

#ifdef CGAL_USE_CORE
  else if (CGAL::assign(conic, arr)) {
#if 0
    typedef CGAL::Arr_text_formatter<Conic_arr>         Conic_text_formatter;
    typedef CGAL::Arr_with_history_text_formatter<Conic_text_formatter>
      ArrFormatter;
    ArrFormatter                                        arrFormatter;
    CGAL::read( *conic, ifs, arrFormatter );
    this->arrangements[ index ] = CGAL::make_object( conic );
    tab->setArrangement( conic );
#endif
    typedef ArrangementDemoTab< Conic_arr > TabType;
    Conic_reader< Conic_arr::Geometry_traits_2 > conicReader;
    std::vector< Conic_arr::Curve_2 > curve_list;
    CGAL::Bbox_2 bbox;
    conicReader.read_data( filename.toStdString( ).c_str( ),
                           std::back_inserter( curve_list ), bbox );
    this->arrangements[ index ] = CGAL::make_object( conic );
    TabType* tab = static_cast< TabType* >( this->tabs[ index ] );

    CGAL::insert( *conic, curve_list.begin(), curve_list.end() );
    tab->setArrangement( conic );
    //QMessageBox::information( this, "Oops",
    //  "Reading conic arrangement not supported" );
  }
#endif

  // else if ( CGAL::assign( alg, arr ) )
  // {
  //   typedef CGAL::Arr_text_formatter< Alg_seg_arr >     Seg_text_formatter;
  //   typedef CGAL::Arr_with_history_text_formatter<Seg_text_formatter>
  //     ArrFormatter;
  //   typedef ArrangementDemoTab< Alg_seg_arr >           TabType;

  //   ArrFormatter arrFormatter;
  //   CGAL::read( *alg, ifs, arrFormatter );
  //   this->arrangements[ index ] = CGAL::make_object( alg );
  //   TabType* tab = static_cast< TabType* >( this->tabs[ index ] );
  //   tab->setArrangement( alg );
  // }
  ifs.close( );
}

void ArrangementDemoWindow::openDatFile( QString filename )
{
  int index = this->ui->tabWidget->currentIndex( );
  if ( index == -1 )
  {
    QMessageBox::information( this, "Oops", "Create a new tab first" );
    return;
  }
  if ( filename.isNull( ) )
  {
    return;
  }

  std::ifstream inputFile( filename.toStdString( ).c_str( ) );
  CGAL::Object arr = this->arrangements[ index ];
  Seg_arr* seg;
  Pol_arr* pol;

#ifdef CGAL_USE_CORE
  Conic_arr* conic;
#endif

  // Alg_seg_arr* alg;

  // Creates an ofstream object named inputFile
  if (! inputFile.is_open() ) // Always test file open
  {
    std::cerr << "Error opening input file" << std::endl;
    return;
  }

  Pol_traits traits;
  Pol_traits::Construct_curve_2 poly_const =
    traits.construct_curve_2_object();


  if ( CGAL::assign( pol, arr ) )
  {
    pol->clear( );

    std::vector<Arr_pol_point_2> points;

    unsigned int num_polylines;
    inputFile >> num_polylines;
    std::list<Arr_pol_2> pol_list;

    unsigned int i;
    for (i = 0; i < num_polylines; i++)
    {
      unsigned int num_segments;
      inputFile >> num_segments;
      points.clear();
      unsigned int j;
      for (j = 0; j < num_segments; j++)
      {
        int ix, iy;
        inputFile >> ix >> iy;
        points.push_back (Arr_pol_point_2(NT(ix),NT(iy)));
      }

      Arr_pol_2 curve = poly_const(points.begin(), points.end());
      pol_list.push_back(curve);
    }
    CGAL::insert(*pol, pol_list.begin(), pol_list.end());

    typedef ArrangementDemoTab< Pol_arr > TabType;
    TabType* tab = static_cast< TabType* >( this->tabs[ index ] );
    tab->setArrangement( pol );
  }
  else if ( CGAL::assign( seg, arr ) )
  {
    seg->clear( );

    int count;
    inputFile >> count;
    int i;
    std::list<Arr_seg_2> seg_list;
    for (i = 0; i < count; i++)
    {
      NT x0, y0, x1, y1;
      inputFile >> x0 >> y0 >> x1 >> y1;

      Arr_seg_point_2 p1(x0, y0);
      Arr_seg_point_2 p2(x1, y1);

      Arr_seg_2 curve(p1, p2);

      seg_list.push_back(curve);
    }

    CGAL::insert(*(seg), seg_list.begin(), seg_list.end());

    typedef ArrangementDemoTab< Seg_arr > TabType;
    TabType* tab = static_cast< TabType* >( this->tabs[ index ] );
    tab->setArrangement( seg );
  }

#ifdef CGAL_USE_CORE
  else if (CGAL::assign( conic, arr)) {
    conic->clear( );
    Conic_reader< Conic_arr::Geometry_traits_2 > reader;
    std::list<Arr_conic_2> curve_list;
    CGAL::Bbox_2 bbox;
    reader.read_data( filename.toStdString().c_str(),
                      std::back_inserter(curve_list), bbox );
    CGAL::insert (*conic, curve_list.begin(), curve_list.end());

    typedef ArrangementDemoTab< Conic_arr > TabType;
    TabType* tab = static_cast< TabType* >( this->tabs[ index ] );
    tab->setArrangement( conic );
  }
#endif

  inputFile.close();
}

void ArrangementDemoWindow::updateEnvelope( QAction* newMode )
{
  if ( this->ui->tabWidget->currentIndex( ) == -1 ) return;
  ArrangementDemoTabBase* activeTab =
    this->tabs[ this->ui->tabWidget->currentIndex( ) ];
  // QGraphicsScene* activeScene = activeTab->getScene( );
  // QGraphicsView* activeView = activeTab->getView( );

  bool show = newMode->isChecked( );
  if ( newMode == this->ui->actionLowerEnvelope )
  {
    activeTab->getEnvelopeCallback( )->showLowerEnvelope( show );
  }
  else if ( newMode == this->ui->actionUpperEnvelope )
  {
    activeTab->getEnvelopeCallback( )->showUpperEnvelope( show );
  }
}

void ArrangementDemoWindow::updateSnapping( QAction* newMode )
{
  ArrangementDemoTabBase* activeTab =
    this->tabs[ this->ui->tabWidget->currentIndex( ) ];
  QGraphicsScene* activeScene = activeTab->getScene( );
  ArrangementDemoGraphicsView* activeView = activeTab->getView( );

  bool enabled = newMode->isChecked( );
  if ( newMode == this->ui->actionSnapMode )
  {
    activeTab->getCurveInputCallback( )->setSnappingEnabled( enabled );
    activeTab->getSplitEdgeCallback( )->setSnappingEnabled( enabled );
    if ( ! enabled )
    {
      this->ui->actionGridSnapMode->setChecked( false );
      this->ui->actionGridSnapMode->setEnabled( false );
      activeTab->getCurveInputCallback( )->setSnapToGridEnabled( false );
      activeTab->getSplitEdgeCallback( )->setSnapToGridEnabled( false );
    }
    else
    {
      this->ui->actionGridSnapMode->setEnabled( true );
    }
  }
  else if ( newMode == this->ui->actionGridSnapMode )
  {
    activeTab->getCurveInputCallback( )->setSnapToGridEnabled( enabled );
    activeTab->getSplitEdgeCallback( )->setSnapToGridEnabled( enabled );
    activeView->setShowGrid( enabled );
  }
  activeScene->update( );
}

void ArrangementDemoWindow::updateConicType( QAction* newType )
{
  ArrangementDemoTabBase* activeTab =
    this->tabs[ this->ui->tabWidget->currentIndex( ) ];
  // QGraphicsScene* activeScene = activeTab->getScene( );
  // ArrangementDemoGraphicsView* activeView = activeTab->getView( );

#ifdef CGAL_USE_CORE
  Conic_arr* conic_arr;
  bool isConicArr =
    CGAL::assign( conic_arr,
                  this->arrangements[ this->ui->tabWidget->currentIndex( ) ] );
#endif

  Lin_arr* lin_arr;
  bool isLinearArr =
    CGAL::assign( lin_arr,
                  this->arrangements[ this->ui->tabWidget->currentIndex( ) ] );

  if ( isLinearArr )
  {
    typedef Lin_arr::Geometry_traits_2         Line_geom_traits;
    typedef CGAL::Qt::GraphicsViewCurveInput<Line_geom_traits>
      LinearCurveInputCallback;
    LinearCurveInputCallback* curveInputCallback =
      ( LinearCurveInputCallback* ) activeTab->getCurveInputCallback( );
    if ( newType == this->ui->actionConicSegment )
    {
      curveInputCallback->setCurveType( LinearCurveInputCallback::SEGMENT );
    }
    else if ( newType == this->ui->actionCurveRay )
  {
      curveInputCallback->setCurveType( LinearCurveInputCallback::RAY );
    }
    else if ( newType == this->ui->actionCurveLine )
    {
      curveInputCallback->setCurveType( LinearCurveInputCallback::LINE );
    }
    //std::cout << "do nothing" << std::endl;
  }

#ifdef CGAL_USE_CORE
  else if (isConicArr) {
    // std::cout << "do something conic arr related" << std::endl;
    typedef Conic_arr::Geometry_traits_2       Conic_geom_traits;
    typedef CGAL::Qt::GraphicsViewCurveInput<Conic_geom_traits>
      ConicCurveInputCallback;
    ConicCurveInputCallback* curveInputCallback =
      ( ConicCurveInputCallback* ) activeTab->getCurveInputCallback( );
    if ( newType == this->ui->actionConicSegment )
    {
      curveInputCallback->setConicType(ConicCurveInputCallback::CONIC_SEGMENT);
    }
    else if ( newType == this->ui->actionConicCircle )
    {
      curveInputCallback->setConicType(ConicCurveInputCallback::CONIC_CIRCLE);
    }
    else if ( newType == this->ui->actionConicEllipse )
    {
      curveInputCallback->setConicType(ConicCurveInputCallback::CONIC_ELLIPSE);
    }
    else if ( newType == this->ui->actionConicThreePoint )
    {
      curveInputCallback->setConicType( ConicCurveInputCallback::
                                        CONIC_THREE_POINT );
    }
    else if ( newType == this->ui->actionConicFivePoint )
    {
      curveInputCallback->setConicType( ConicCurveInputCallback::
                                        CONIC_FIVE_POINT );
    }
  }
#endif
}

void ArrangementDemoWindow::on_actionSaveAs_triggered( )
{
  int index = this->ui->tabWidget->currentIndex( );
  if ( index == -1 )
    return;
  QString filename =
    QFileDialog::getSaveFileName( this, tr( "Save file" ),
                                  "", "Arrangement (*.arr)" );
  if ( filename.isNull( ) )
    return;

  std::ofstream ofs( filename.toStdString( ).c_str( ) );
  CGAL::Object arr = this->arrangements[ index ];
  Seg_arr* seg;
  Pol_arr* pol;

#ifdef CGAL_USE_CORE
  Conic_arr* conic;
#endif

  if ( CGAL::assign( seg, arr ) )
  {
    typedef CGAL::Arr_text_formatter<Seg_arr>           Seg_text_formatter;
    typedef CGAL::Arr_with_history_text_formatter<Seg_text_formatter>
      ArrFormatter;
    ArrFormatter                                        arrFormatter;
    CGAL::write( *seg, ofs, arrFormatter );
  }
  else if ( CGAL::assign( pol, arr ) )
  {
    typedef CGAL::Arr_text_formatter<Pol_arr>           Pol_text_formatter;
    typedef CGAL::Arr_with_history_text_formatter<Pol_text_formatter>
      ArrFormatter;
    ArrFormatter                                        arrFormatter;
    CGAL::write( *pol, ofs, arrFormatter );
  }

#ifdef CGAL_USE_CORE
  else if (CGAL::assign(conic, arr)) {
#if 0
    typedef CGAL::Arr_text_formatter<Conic_arr>         Conic_text_formatter;
    typedef CGAL::Arr_with_history_text_formatter<Conic_text_formatter>
      ArrFormatter;
    ArrFormatter                                        arrFormatter;
    CGAL::write( *conic, ofs, arrFormatter );
#endif
    ofs << conic->number_of_curves( ) << std::endl;
    for ( Conic_arr::Curve_iterator it = conic->curves_begin( );
          it != conic->curves_end( ); ++it )
    {
      if ( it->is_full_conic( ) )
      {
        ofs << "F ";
        ofs << it->r( ) << " ";
        ofs << it->s( ) << " ";
        ofs << it->t( ) << " ";
        ofs << it->u( ) << " ";
        ofs << it->v( ) << " ";
        ofs << it->w( ) << " ";
        ofs << std::endl;
      }
      else if ( it->orientation( ) == CGAL::COLLINEAR )
      {
        ofs << "S ";
        ofs << it->source( ) << " ";
        ofs << it->target( ) << " ";
        ofs << std::endl;
      }
      else
      {
        ofs << "A ";
        ofs << it->r( ) << " ";
        ofs << it->s( ) << " ";
        ofs << it->t( ) << " ";
        ofs << it->u( ) << " ";
        ofs << it->v( ) << " ";
        ofs << it->w( ) << " ";
        if ( it->orientation( ) == CGAL::COUNTERCLOCKWISE )
          ofs << "1 ";
        else if ( it->orientation( ) == CGAL::CLOCKWISE )
          ofs << "-1 ";
        else
          ofs << "0 ";
        ofs << it->source( ) << " ";
        ofs << it->target( ) << " ";
        ofs << std::endl;
      }
    }
  }
#endif

  ofs.close( );
}

void ArrangementDemoWindow::on_actionOpen_triggered( )
{
  int index = this->ui->tabWidget->currentIndex( );
  if ( index == -1 )
  {
    QMessageBox::information( this, "Oops", "Create a new tab first" );
    return;
  }
  QString filename =
    QFileDialog::getOpenFileName( this, tr( "Open file" ),
                                  "", "Arrangement files (*.arr *.dat);;All files (*.*)" );
  if ( filename.isNull( ) )
    return;

  if ( filename.endsWith( ".arr" ) )
  {
    this->openArrFile( filename );
  }
  else
  {
    this->openDatFile( filename );
  }

  ArrangementDemoTabBase* currentTab = this->tabs[ index ];
  CGAL::Qt::ArrangementGraphicsItemBase* agi =
    currentTab->getArrangementGraphicsItem( );
  QRectF bb = agi->boundingRect( );
  QGraphicsView* view = currentTab->getView( );
  // std::cout << bb.left( ) << " " << bb.bottom( ) << ", " << bb.right( )
  //           << " " << bb.top( ) << std::endl;
  if ( boost::math::isinf(bb.left( )) ||
       boost::math::isinf(bb.right( )) ||
       boost::math::isinf(bb.top( )) ||
       boost::math::isinf(bb.bottom( )) )
  {
    // std::cout << "unbounded; using default bb" << std::endl;
    bb = QRectF( -100, -100, 200, 200 );
    view->setSceneRect( bb );
  }
  else
  {
    view->fitInView( bb, ::Qt::KeepAspectRatio );
    view->setSceneRect( bb );
  }
#if 0
  view->centerOn( bb.center( ) );
#endif
}

void ArrangementDemoWindow::on_actionQuit_triggered( )
{
  qApp->exit( );
}

void ArrangementDemoWindow::on_actionNewTab_triggered( )
{
  NewTabDialog* newTabDialog = new NewTabDialog;
  if ( newTabDialog->exec( ) == QDialog::Accepted )
  {
    int id = newTabDialog->checkedId( );
    if ( id == SEGMENT_TRAITS )
    {
      this->makeTab( SEGMENT_TRAITS );
    }
    else if ( id == POLYLINE_TRAITS )
    {
      this->makeTab( POLYLINE_TRAITS );
    }
    else if ( id == CONIC_TRAITS )
    {
      this->makeTab( CONIC_TRAITS );
    }
    else if ( id == LINEAR_TRAITS )
    {
      this->makeTab( LINEAR_TRAITS );
    }
    else if ( id == CIRCULAR_ARC_TRAITS )
    {
      this->makeTab( CIRCULAR_ARC_TRAITS );
    }
    // else if ( id == ALGEBRAIC_TRAITS )
    // {
    //   this->makeTab( ALGEBRAIC_TRAITS );
    // }
    else
    {
      std::cout << "Sorry, this trait is not yet supported" << std::endl;
    }
  }
  delete newTabDialog;
}

void ArrangementDemoWindow::on_tabWidget_currentChanged( )
{
  // std::cout << "Tab changed" << std::endl;
  // disable the callback for the previously active tab
  this->resetCallbackState( this->lastTabIndex );
  this->removeCallback( this->lastTabIndex );
  this->lastTabIndex = this->ui->tabWidget->currentIndex( );

  this->updateMode( this->modeGroup->checkedAction( ) );

  CGAL::Object arr;
  if ( this->ui->tabWidget->currentIndex( ) != -1 )
    arr = this->arrangements[ this->ui->tabWidget->currentIndex( ) ];

  // Seg_arr* seg;
  // Pol_arr* pol;
  Lin_arr* lin;

#ifdef CGAL_USE_CORE
  Conic_arr* conic;
#endif

  if ( CGAL::assign( lin, arr ) )
  {
    this->ui->actionConicSegment->setChecked( true );

    this->ui->actionCurveRay->setVisible( true );
    this->ui->actionCurveLine->setVisible( true );

    this->ui->actionConicCircle->setVisible( false );
    this->ui->actionConicEllipse->setVisible( false );
    this->ui->actionConicThreePoint->setVisible( false );
    this->ui->actionConicFivePoint->setVisible( false );

    this->conicTypeGroup->setEnabled( true );
  }

#ifdef CGAL_USE_CORE
  else if (CGAL::assign( conic, arr)) {
    this->ui->actionConicSegment->setChecked( true );

    this->ui->actionCurveRay->setVisible( false );
    this->ui->actionCurveLine->setVisible( false );

    this->ui->actionConicCircle->setVisible( true );
    this->ui->actionConicEllipse->setVisible( true );
    this->ui->actionConicThreePoint->setVisible( true );
    this->ui->actionConicFivePoint->setVisible( true );

    this->conicTypeGroup->setEnabled( true );
  }
#endif

  else { // segment or polyline
    this->ui->actionConicSegment->setChecked( true );

    this->ui->actionCurveRay->setVisible( false );
    this->ui->actionCurveLine->setVisible( false );

    this->ui->actionConicCircle->setVisible( false );
    this->ui->actionConicEllipse->setVisible( false );
    this->ui->actionConicThreePoint->setVisible( false );
    this->ui->actionConicFivePoint->setVisible( false );

    this->conicTypeGroup->setEnabled( true );
  }
}

void ArrangementDemoWindow::on_actionOverlay_triggered( )
{
  OverlayDialog* overlayDialog = new OverlayDialog( this );
  if ( overlayDialog->exec( ) == QDialog::Accepted )
  {
    std::vector< CGAL::Object > arrs = overlayDialog->selectedArrangements( );
    if ( arrs.size( ) == 2 )
    {
      Seg_arr* seg_arr;
      Seg_arr* seg_arr2;
      Pol_arr* pol_arr;
      Pol_arr* pol_arr2;

#ifdef CGAL_USE_CORE
      Conic_arr* conic_arr;
      Conic_arr* conic_arr2;
#endif

      Lin_arr* lin_arr;
      Lin_arr* lin_arr2;
      Arc_arr* arc_arr;
      Arc_arr* arc_arr2;
      // Alg_seg_arr* alg_arr;
      // Alg_seg_arr* alg_arr2;
      if ( CGAL::assign( seg_arr, arrs[ 0 ] ) &&
           CGAL::assign( seg_arr2, arrs[ 1 ] ) )
      {
        this->makeOverlayTab( seg_arr, seg_arr2 );
      }
      if ( CGAL::assign( pol_arr, arrs[ 0 ] ) &&
           CGAL::assign( pol_arr2, arrs[ 1 ] ) )
      {
        this->makeOverlayTab( pol_arr, pol_arr2 );
      }

#ifdef CGAL_USE_CORE
      if (CGAL::assign(conic_arr, arrs[0]) && CGAL::assign(conic_arr2, arrs[1]))
      {
        this->makeOverlayTab( conic_arr, conic_arr2 );
      }
#endif

      if ( CGAL::assign( lin_arr, arrs[ 0 ] ) &&
           CGAL::assign( lin_arr2, arrs[ 1 ] ) )
      {
        this->makeOverlayTab( lin_arr, lin_arr2 );
      }
      if ( CGAL::assign( arc_arr, arrs[ 0 ] ) &&
           CGAL::assign( arc_arr2, arrs[ 1 ] ) )
      {
        this->makeOverlayTab( arc_arr, arc_arr2 );
      }
      // if ( CGAL::assign( alg_arr, arrs[ 0 ] ) &&
      //      CGAL::assign( alg_arr2, arrs[ 1 ] ) )
      // {
      //   this->makeOverlayTab( alg_arr, alg_arr2 );
      // }

#if 0
      if ( CGAL::assign( conic_arr, arrs[ 0 ] ) ||
           CGAL::assign( conic_arr, arrs[ 1 ] ) )
      {
        this->makeOverlayTab< Conic_arr >( arrs[ 0 ], arrs[ 1 ] );
      }
      else if ( CGAL::assign( pol_arr, arrs[ 0 ] ) ||
                CGAL::assign( pol_arr, arrs[ 1 ] ) )
      {
        this->makeOverlayTab< Pol_arr >( arrs[ 0 ], arrs[ 1 ] );
      }
      else
      {
        this->makeOverlayTab< Seg_arr >( arrs[ 0 ], arrs[ 1 ] );
      }
#endif
    }
  }
  delete overlayDialog;
}

void ArrangementDemoWindow::on_actionCloseTab_triggered( )
{
  unsigned int currentTabIndex = this->ui->tabWidget->currentIndex( );
  if (! this->ui->tabWidget->count() ||
      (currentTabIndex == static_cast<unsigned int>(-1)))
    return;

  // delete the tab
  this->ui->tabWidget->removeTab( currentTabIndex );
  this->tabs.erase( this->tabs.begin( ) + currentTabIndex );

  // delete the arrangement
  this->arrangements.erase( this->arrangements.begin( ) + currentTabIndex );
}

void ArrangementDemoWindow::on_actionPrintConicCurves_triggered( )
{
  // int currentTabIndex = this->ui->tabWidget->currentIndex( );
  // Conic_arr* arr;
  // if (currentTabIndex == static_cast<unsigned int>(-1)) return;
  // CGAL::Object o = this->arrangements[ currentTabIndex ];
  // if ( ! CGAL::assign( arr, o ) )
  //     return;
  // typedef typename Conic_arr::Curve_iterator Curve_iterator;
  // // std::cout << arr->number_of_curves( ) << std::endl;
  // Curve_iterator it;
  // for (it = arr->curves_begin(); it != arr->curves_end(); ++it) {
  //     std::cout << *it << std::endl;
  //     if ( it->is_full_conic( ) )
  //     {
  //         std::cout << "F ";
  //         std::cout << it->r( ) << " ";
  //         std::cout << it->s( ) << " ";
  //         std::cout << it->t( ) << " ";
  //         std::cout << it->u( ) << " ";
  //         std::cout << it->v( ) << " ";
  //         std::cout << it->w( ) << " ";
  //         std::cout << std::endl;
  //     }
  //     else if ( it->orientation( ) == CGAL::COLLINEAR )
  //     {
  //         std::cout << "S ";
  //         std::cout << it->source( ) << " ";
  //         std::cout << it->target( ) << " ";
  //         std::cout << std::endl;
  //     }
  //     else
  //     {
  //         std::cout << "A ";
  //         std::cout << it->r( ) << " ";
  //         std::cout << it->s( ) << " ";
  //         std::cout << it->t( ) << " ";
  //         std::cout << it->u( ) << " ";
  //         std::cout << it->v( ) << " ";
  //         std::cout << it->w( ) << " ";
  //         if ( it->orientation( ) == CGAL::COUNTERCLOCKWISE )
  //             std::cout << "1 ";
  //         else if ( it->orientation( ) == CGAL::CLOCKWISE )
  //             std::cout << "-1 ";
  //         else
  //             std::cout << "0 ";
  //         std::cout << it->source( ) << " ";
  //         std::cout << it->target( ) << " ";
  //         std::cout << std::endl;
  //     }
  // }
}

void ArrangementDemoWindow::on_actionZoomIn_triggered( )
{
  unsigned int currentTabIndex = this->ui->tabWidget->currentIndex( );
  if (currentTabIndex == static_cast<unsigned int>(-1)) return;
  ArrangementDemoTabBase* currentTab = this->tabs[ currentTabIndex ];
  QGraphicsView* view = currentTab->getView( );
  view->scale( 2.0, 2.0 );
}

void ArrangementDemoWindow::on_actionZoomOut_triggered( )
{
  unsigned int currentTabIndex = this->ui->tabWidget->currentIndex( );
  if (currentTabIndex == static_cast<unsigned int>(-1)) return;
  ArrangementDemoTabBase* currentTab = this->tabs[ currentTabIndex ];
  QGraphicsView* view = currentTab->getView( );
  view->scale( 0.5, 0.5 );
}

void ArrangementDemoWindow::on_actionPreferences_triggered( )
{
  unsigned int currentTabIndex = this->ui->tabWidget->currentIndex( );
  if (currentTabIndex == static_cast<unsigned int>(-1)) return;
  ArrangementDemoTabBase* currentTab = this->tabs[ currentTabIndex ];
  CGAL::Qt::ArrangementGraphicsItemBase* agi =
    currentTab->getArrangementGraphicsItem( );
  ArrangementDemoGraphicsView* view = currentTab->getView( );
  EnvelopeCallbackBase* envelopeCallback = currentTab->getEnvelopeCallback( );
  VerticalRayShootCallbackBase* verticalRayShootCallback =
    currentTab->getVerticalRayShootCallback( );
  SplitEdgeCallbackBase* splitEdgeCallback = currentTab->getSplitEdgeCallback( );

#if 0
  QPen vertexPen = agi->getVerticesPen( );
  QPen edgePen = agi->getEdgesPen( );
  QBrush vertexPenBrush = vertexPen.brush( );
  QColor vertexColor = vertexPenBrush.color( );
  QBrush edgePenBrush = edgePen.brush( );
  QColor edgeColor = edgePenBrush.color( );
#endif

  ArrangementDemoPropertiesDialog* dialog =
    new ArrangementDemoPropertiesDialog( this );
  if ( dialog->exec( ) == QDialog::Accepted )
  {
    typedef ArrangementDemoPropertiesDialog Dialog;

    QColor edgeColor =  dialog->property(Dialog::EDGE_COLOR_KEY).value<QColor>();

    unsigned int edgeWidth = dialog->property(Dialog::EDGE_WIDTH_KEY).value<unsigned int>();

    QColor vertexColor = dialog->property(Dialog::VERTEX_COLOR_KEY).value<QColor>();

    unsigned int vertexRadius = dialog->property(Dialog::VERTEX_RADIUS_KEY).value<unsigned int>();

    QColor envelopeEdgeColor = dialog->property(Dialog::ENVELOPE_EDGE_COLOR_KEY).value<QColor>();

    unsigned int envelopeEdgeWidth = dialog->property(Dialog::ENVELOPE_EDGE_WIDTH_KEY).value<unsigned int>();

    QColor envelopeVertexColor = dialog->property(Dialog::ENVELOPE_VERTEX_COLOR_KEY).value<QColor>();

    unsigned int envelopeVertexRadius = dialog->property(Dialog::ENVELOPE_VERTEX_RADIUS_KEY).value<unsigned int>();

    QColor verticalRayEdgeColor = dialog->property(Dialog::VERTICAL_RAY_EDGE_COLOR_KEY).value<QColor>();

    unsigned int verticalRayEdgeWidth = dialog->property(Dialog::VERTICAL_RAY_EDGE_WIDTH_KEY).value<unsigned int>();

    DeleteCurveMode mode = dialog->property(Dialog::DELETE_CURVE_MODE_KEY).value<DeleteCurveMode>();

    unsigned int gridSize = dialog->property(Dialog::GRID_SIZE_KEY).value<unsigned int>();

    QColor gridColor = dialog->property(Dialog::GRID_COLOR_KEY).value<QColor>();
    //end new for Qt5 version !


    QPen edgesPen(QBrush(edgeColor), edgeWidth);
    QPen verticesPen(QBrush(vertexColor), vertexRadius);
    agi->setEdgesPen( edgesPen );
    agi->setVerticesPen( verticesPen );
    agi->modelChanged( );
    view->setGridSize( gridSize );
    view->setGridColor( gridColor );
    envelopeCallback->setEnvelopeEdgeColor( envelopeEdgeColor );
    envelopeCallback->setEnvelopeEdgeWidth( envelopeEdgeWidth );
    envelopeCallback->setEnvelopeVertexColor( envelopeVertexColor );
    envelopeCallback->setEnvelopeVertexRadius( envelopeVertexRadius );
    verticalRayShootCallback->setEdgeColor( verticalRayEdgeColor );
    verticalRayShootCallback->setEdgeWidth( verticalRayEdgeWidth );
    splitEdgeCallback->setColor( edgeColor );

#if 0
    std::cout << edgeColor.name( ).toStdString( ) << std::endl;
    std::cout << edgeWidth << std::endl;
    std::cout << vertexColor.name( ).toStdString( ) << std::endl;
    std::cout << vertexRadius << std::endl;
    std::cout << DeleteCurveMode::ToString( mode ).toStdString( ) << std::endl;
#endif

  }
}

void ArrangementDemoWindow::on_actionFillColor_triggered( )
{
  unsigned int currentTabIndex = this->ui->tabWidget->currentIndex( );
  if (currentTabIndex == static_cast<unsigned int>(-1)) return;
  ArrangementDemoTabBase* currentTab = this->tabs[ currentTabIndex ];
  FillFaceCallbackBase* fillFaceCallback = currentTab->getFillFaceCallback( );
  QColor fillColor = fillFaceCallback->getColor( );

  QColor selectedColor = QColorDialog::getColor( fillColor );
  if ( selectedColor.isValid( ) )
  {
    fillFaceCallback->setColor( selectedColor );
    this->updateFillColorSwatch( );
  }
}
