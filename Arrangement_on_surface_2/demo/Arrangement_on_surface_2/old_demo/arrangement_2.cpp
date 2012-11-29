// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/features/gsoc2012-Arrangement_on_surface_2-demo-atsui/Arrangement_on_surface_2/demo/Arrangement_on_surface_2/arrangement_2.cpp $
// $Id: arrangement_2.cpp 67117 2012-01-13 18:14:48Z lrineau $
//
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#include <CGAL/basic.h>

#include "arrangement_2.h"
#include "forms.h"
#include "qt_layer.h"
#include "demo_tab.h"

#include <CGAL/IO/pixmaps/hand.xpm>
#include <CGAL/IO/pixmaps/movepoint.xpm>
#include <CGAL/IO/pixmaps/point.xpm>
#include <CGAL/IO/pixmaps/line.xpm>
#include <CGAL/IO/pixmaps/arrow.xpm>
#include <CGAL/IO/pixmaps/alpha_shape.xpm>
#include "icons/polyline.xpm"
#include "icons/conic.xpm"
#include "icons/ray_shooting2.xpm"
#include "icons/conic_types.xpm"

#include "icons/demo_insert.xpm"
#include "icons/demo_delete.xpm"
#include "icons/demo_snapgrid.xpm"
#include "icons/demo_snapvertex.xpm"
#include "icons/demo_merge.xpm"
#include "icons/demo_split.xpm"
#include "icons/demo_zoomout.xpm"
#include "icons/demo_zoomin.xpm"
#include "icons/demo_pointlocation.xpm"

#include "icons/demo_conic_3points.xpm"
#include "icons/demo_conic_5points.xpm"
#include "icons/demo_conic_circle.xpm"
#include "icons/demo_conic_ellipse.xpm"
#include "icons/demo_conic_segment.xpm"
#include "icons/demo_rayshoot_down.xpm"
#include "icons/demo_rayshoot_up.xpm"
#include "icons/lower_env_xpm.xpm"
#include "icons/upper_env_xpm.xpm"
#include "icons/demo_fill.xpm"
#include "icons/demo_colors.xpm"






//////////////////////////////////////////////////////////////////////////////
/*! MyWindow constructor
 * \param w - window width
 * \param h - window hight
 */
MyWindow::MyWindow(int w, int h) : num_of_colors(18)
{
  myBar = new QTabWidget(this);
  setCentralWidget(myBar);
  m_width = w;
  m_height = h;
  tab_number = 0;
  number_of_tabs = 0;
  testlayer = new Qt_layer( myBar );
  colors_flag = true;
  statusBar();

  m_scailing_factor = 2;

  // Traits Group
  QActionGroup *traitsGroup = new QActionGroup( this ); // Connected later
  traitsGroup->setExclusive( TRUE );

  setSegmentTraits = new QAction("Segment Traits",
                                 QPixmap( (const char**)line_xpm ),
                                 "&Segment Traits", 0 ,traitsGroup,
                                 "Segment Traits" );
  setSegmentTraits->setToggleAction( TRUE );

  setPolylineTraits = new QAction("Polyline Traits",
                                  QPixmap( (const char**)polyline_xpm ),
                                  "&Polyline Traits", 0 , traitsGroup,
                                  "Polyline Traits" );
  setPolylineTraits->setToggleAction( TRUE );

#ifdef CGAL_USE_CORE
  setConicTraits = new QAction("Conic Traits",
                               QPixmap( (const char**)conic_xpm ),
                               "&Conic Traits", 0 , traitsGroup,
                               "Conic Traits" );
  setConicTraits->setToggleAction( TRUE );
#endif

  // Snap Mode Group

  setSnapMode = new QAction("Snap Mode", QPixmap( (const char**)snapvertex_xpm ),
                            "&Snap Mode", 0 , this, "Snap Mode" );
  setSnapMode->setToggleAction( TRUE );

  setGridSnapMode = new QAction("Grid Snap Mode",
                                QPixmap( (const char**)snapgrid_xpm ),
                                "&Grid Snap Mode", 0 , this,
                                "Grid Snap Mode" );
  setGridSnapMode->setToggleAction( TRUE );

  // insert - delete - point_location Mode Group
  QActionGroup *modeGroup = new QActionGroup( this ); // Connected later
  modeGroup->setExclusive( TRUE );

  insertMode = new QAction("Insert", QPixmap( (const char**)insert_xpm ),
                           "&Insert", 0 , modeGroup, "Insert" );
  insertMode->setToggleAction( TRUE );

  deleteMode = new QAction("Delete", QPixmap( (const char**)delete_xpm ),
                           "&Delete", 0 , modeGroup, "Delete" );
  deleteMode->setToggleAction( TRUE );

  pointLocationMode = new QAction("PointLocation",
                                  QPixmap( (const char**)pointlocation_xpm ),
                                  "&Point Location", 0 , modeGroup,
                                  "Point Location" );
  pointLocationMode->setToggleAction( TRUE );

  rayShootingUpMode = new QAction("RayShootingUp",
                                QPixmap( (const char**)demo_rayshoot_up_xpm ),
                                "&Ray Shooting Up", 0 , modeGroup,
                                "Ray Shooting Up" );
  rayShootingUpMode->setToggleAction( TRUE );

  rayShootingDownMode = new QAction("RayShootingDown",
                                QPixmap( (const char**)demo_rayshoot_down_xpm ),
                                "&Ray Shooting Down", 0 , modeGroup,
                                "Ray Shooting Down" );
  rayShootingDownMode->setToggleAction( TRUE );

  dragMode = new QAction("Drag", QPixmap( (const char**)hand_xpm ),
                         "&Drag", 0 , modeGroup, "Drag" );
  dragMode->setToggleAction( TRUE );

  mergeMode = new QAction("Merge", QPixmap( (const char**)merge_xpm ),
                         "&Merge", 0 , modeGroup, "Merge" );
  mergeMode->setToggleAction( TRUE );

  splitMode = new QAction("Split", QPixmap( (const char**)split_xpm ),
                         "&Split", 0 , modeGroup, "Split" );
  splitMode->setToggleAction( TRUE );

   fillfaceMode = new QAction("Fill", QPixmap( (const char**)demo_fill_xpm ),
                         "&Fill", 0 , modeGroup, "Fill" );
  fillfaceMode->setToggleAction( TRUE );



  // zoom in
  zoominBt = new QAction("Zoom in", QPixmap( (const char**)zoomin_xpm ),
                         "&Zoom in", 0 , this, "Zoom in" );
  // zoom out
  zoomoutBt = new QAction("Zoom out", QPixmap( (const char**)zoomout_xpm ),
                          "&Zoom out", 0 , this, "Zoom out" );

   // color dialog
  color_dialog_bt = new QAction("Choose color", QPixmap( (const char**)demo_colors_xpm ),
                         "&choose color", 0 , this, "choose color" );

  lower_env_dialog_bt = new QAction("Lower envelope", QPixmap( (const char**)lower_env_xpm ),
                             "&lower envelope", 0, this, "Lower envelop" );
  lower_env_dialog_bt->setToggleAction( TRUE );

  upper_env_dialog_bt = new QAction("Upper envelope", QPixmap( (const char**)upper_env_xpm ),
                             "&upper envelope", 0, this, "Upper envelop" );
  upper_env_dialog_bt->setToggleAction( TRUE );

#ifdef CGAL_USE_CORE
  // Conic Type Group
  QActionGroup *conicTypeGroup = new QActionGroup( this ); // Connected later
  conicTypeGroup->setExclusive( TRUE );

  setCircle = new QAction("Circle",
                                 QPixmap( (const char**)demo_conic_circle_xpm ),
                                 "&Circle", 0 ,conicTypeGroup,
                                 "Circle" );
  setCircle->setToggleAction( TRUE );
  setSegment = new QAction("Segment",
                                 QPixmap( (const char**)demo_conic_segment_xpm ),
                                 "&Segment", 0 ,conicTypeGroup,
                                 "Segment" );
  setSegment->setToggleAction( TRUE );
  setEllipse = new QAction("Ellipse",
                                 QPixmap( (const char**)demo_conic_ellipse_xpm ),
                                 "&Ellipse", 0 ,conicTypeGroup,
                                 "Ellipse" );
  setEllipse->setToggleAction( TRUE );
  setParabola = new QAction("3 Points Arc",
                                 QPixmap( (const char**)demo_conic_3points_xpm ),
                                 "&3 Points Arc", 0 ,conicTypeGroup,
                                 "3 Points Arc" );
  setParabola->setToggleAction( TRUE );
  setHyperbola = new QAction("5 Points Arc",
                                 QPixmap( (const char**)demo_conic_5points_xpm ),
                                 "&5 Points Arc", 0 ,conicTypeGroup,
                                 "5 Points Arc" );
  setHyperbola->setToggleAction( TRUE );
#endif

  //create a timer for checking if somthing changed
  QTimer *timer = new QTimer( this );
  connect( timer, SIGNAL(timeout()),
           this, SLOT(timer_done()) );
  timer->start( 200, FALSE );

  // file menu
  QPopupMenu * file = new QPopupMenu( this );
  menuBar()->insertItem( "&File", file );
  file->insertItem("&Open Segment File...", this, SLOT(fileOpenSegment()));
  file->insertItem("&Open Polyline File...", this, SLOT(fileOpenPolyline()));\

#ifdef CGAL_USE_CORE
  file->insertItem("&Open Conic File...", this, SLOT(fileOpenConic()));
#endif

  file->insertItem("&Open Segment Arr File...", this, SLOT(fileOpenSegmentPm()));
  file->insertItem("&Open Polyline Arr File...", this, SLOT(fileOpenPolylinePm()));
  //file->insertItem("&Open Conic Pm File", this, SLOT(fileOpenConicPm()));
  file->insertItem("&Save...", this, SLOT(fileSave()));
  file->insertItem("&Save As...", this, SLOT(fileSaveAs()));
  //file->insertItem("&Save to ps...", this, SLOT(fileSave_ps()));
  file->insertSeparator();
  file->insertItem("&Print...", this , SLOT(print()));
  file->insertSeparator();
  file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
  file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), CTRL+Key_Q );
  menuBar()->insertSeparator();

  // tab menu
  QPopupMenu * tab = new QPopupMenu( this );
  menuBar()->insertItem( "&Tab", tab );
  tab->insertItem("Add &Segment Tab", this, SLOT(add_segment_tab()));
  tab->insertItem("Add &Polyline Tab", this, SLOT(add_polyline_tab()));

#ifdef CGAL_USE_CORE
  tab->insertItem("Add &Conic Tab", this, SLOT(add_conic_tab()));
  tab->insertSeparator();
#endif

  tab->insertItem("Remove &Tab", this, SLOT(remove_tab()));
  menuBar()->insertSeparator();

  // mode menu
  QPopupMenu * mode = new QPopupMenu( this );
  menuBar()->insertItem( "&Mode", mode );
  insertMode->addTo( mode );
  deleteMode->addTo( mode );
  pointLocationMode->addTo( mode );
  rayShootingUpMode->addTo( mode );
  rayShootingDownMode->addTo( mode );
  dragMode->addTo( mode );
  mergeMode->addTo( mode );
  splitMode->addTo( mode );
  fillfaceMode->addTo( mode );
  menuBar()->insertSeparator();

  // snap mode menu
  QPopupMenu * snap_mode = new QPopupMenu( this );
  menuBar()->insertItem( "&Snap mode", snap_mode );
  setSnapMode->addTo(snap_mode);
  setGridSnapMode->addTo(snap_mode);
  menuBar()->insertSeparator();

  // traits menu
  QPopupMenu * traits = new QPopupMenu( this );
  menuBar()->insertItem( "&Traits Type", traits );
  setSegmentTraits->addTo(traits);
  setPolylineTraits->addTo(traits);
#ifdef CGAL_USE_CORE
  setConicTraits->addTo(traits);
#endif

  // options menu
  QPopupMenu * options = new QPopupMenu( this );
  menuBar()->insertItem( "&Options", options );
  options->insertSeparator();
  options->insertItem("Overlay...", this, SLOT(overlay_pm()));
  options->insertSeparator();
  options->insertItem("Properties...", this, SLOT(properties()));
  options->insertSeparator();
  options->insertItem("Show Grid", this, SLOT(showGrid()));
  options->insertItem("Hide Grid", this, SLOT(hideGrid()));
  options->insertSeparator();
  //options->insertItem("Conic Type", this, SLOT(conicType()));
  //options->insertSeparator();
  options->insertItem("Unbounded Face Color...", this, SLOT(backGroundColor()));
  options->insertSeparator();
  options->insertItem("Edge Color...", this, SLOT(changeEdgeColor()));
  options->insertSeparator();
  options->insertItem("Vertex Color...", this, SLOT(changeVertexColor()));
  options->insertSeparator();
  options->insertItem("Point-Locaiton Strategy....", this ,
                                         SLOT(pointLocationStrategy()));


  // help menu
  QPopupMenu * help = new QPopupMenu( this );
  menuBar()->insertItem( "&Help", help );
  help->insertItem("How To...", this, SLOT(howto()), Key_F1);
  help->insertSeparator();
  help->insertItem("&About...", this, SLOT(about()), CTRL+Key_A );
  help->insertItem("About &Qt...", this, SLOT(aboutQt()) );

  QToolBar *modeTools = new QToolBar( this, "mode operations" );
  modeTools->setLabel( "Mode Operations" );
  insertMode->addTo( modeTools );
  deleteMode->addTo( modeTools );
  dragMode->addTo( modeTools );
  pointLocationMode->addTo( modeTools );
  rayShootingUpMode->addTo( modeTools );
  rayShootingDownMode->addTo( modeTools );
  mergeMode->addTo( modeTools );
  splitMode->addTo( modeTools );
  fillfaceMode->addTo( modeTools );
  modeTools->addSeparator();

  QToolBar *snapModeTools = new QToolBar( this, "snapMode operations" );
  snapModeTools->setLabel( "Snap Mode Operations" );
  snapModeTools->addSeparator();
  setSnapMode->addTo( snapModeTools );
  setGridSnapMode->addTo( snapModeTools );
  snapModeTools->addSeparator();

  QToolBar *traitsTool = new QToolBar( this, "traits type" );
  traitsTool->setLabel( "Traits Type" );
  traitsTool->addSeparator();
  setSegmentTraits->addTo( traitsTool );
  setPolylineTraits->addTo( traitsTool );
#ifdef CGAL_USE_CORE
  setConicTraits->addTo( traitsTool );
#endif
  traitsTool->addSeparator();

  QToolBar *zoomTool = new QToolBar( this, "zoom" );
  zoomTool->setLabel( "Zoom" );
  zoomTool->addSeparator();
  zoomoutBt->addTo( zoomTool );
  zoominBt->addTo( zoomTool );
  zoomTool->addSeparator();

  QToolBar *colorTool = new QToolBar( this, "color" );
  colorTool->addSeparator();
  colorTool->setLabel("Choose color");
  color_dialog_bt->addTo(colorTool);
  colorTool->addSeparator();

  QToolBar *envelopeTool = new QToolBar( this, "envelopes" );
  envelopeTool->addSeparator();
  envelopeTool->setLabel("Envelopes");
  lower_env_dialog_bt->addTo(envelopeTool);
  upper_env_dialog_bt->addTo(envelopeTool);
  envelopeTool->addSeparator();


#ifdef CGAL_USE_CORE
  conicTypeTool = new QToolBar( this, "conic type" );
  conicTypeTool->setLabel( "Conic Type" );
  conicTypeTool->addSeparator();
  setSegment->addTo( conicTypeTool );
  setCircle->addTo( conicTypeTool );
  setEllipse->addTo( conicTypeTool );
  setParabola->addTo( conicTypeTool );
  setHyperbola->addTo( conicTypeTool );
#endif

  connect( zoomoutBt, SIGNAL( activated () ) ,
       this, SLOT( zoomout() ) );

  connect( zoominBt, SIGNAL( activated () ) ,
       this, SLOT( zoomin() ) );

  connect (color_dialog_bt , SIGNAL( activated()) ,
          this , SLOT(openColorDialog() ) );

  connect (lower_env_dialog_bt, SIGNAL(toggled(bool)) ,
           this, SLOT(lowerEnvelope(bool) ));

  connect (upper_env_dialog_bt, SIGNAL(toggled(bool)) ,
           this, SLOT(upperEnvelope(bool) ));
  // connect mode group
  connect( modeGroup, SIGNAL( selected(QAction*) ),
           this, SLOT( updateMode(QAction*) ) );

  // connect Traits Group
  connect( traitsGroup, SIGNAL( selected(QAction*) ),
           this, SLOT( updateTraitsType(QAction*) ) );

#ifdef CGAL_USE_CORE
  // connect Conic Type Group
  connect( conicTypeGroup, SIGNAL( selected(QAction*) ),
           this, SLOT( updateConicType(QAction*) ) );
#endif

  // connect Snap Mode

  connect( setSnapMode, SIGNAL( toggled( bool ) ) ,
           this, SLOT( updateSnapMode( bool ) ) );

  connect( setGridSnapMode, SIGNAL( toggled( bool ) ) ,
       this, SLOT( updateGridSnapMode( bool ) ) );

  // connect the change of current tab
  connect( myBar, SIGNAL( currentChanged(QWidget * )  ),
           this, SLOT( update() ) );

  colors = new QColor[num_of_colors];
  colors[0]  =  Qt::blue;
  colors[1]  =  Qt::gray;
  colors[2]  =  Qt::green;
  colors[3]  =  Qt::cyan;
  colors[4]  =  Qt::magenta;
  colors[5]  =  Qt::darkRed;
  colors[6]  =  Qt::darkGreen;
  colors[7]  =  Qt::darkBlue;
  colors[8]  =  Qt::darkMagenta;
  colors[9]  =  Qt::darkCyan;
  colors[10] =  Qt::yellow;
  colors[11] =  Qt::white;
  colors[12] =  Qt::darkGray;
  colors[13] =  Qt::gray;
  colors[14] =  Qt::red;
  colors[15] =  Qt::cyan;
  colors[16] =  Qt::darkYellow;
  colors[17] =  Qt::lightGray;

  //state flag
  old_state = 0;
  add_segment_tab();
  resize(m_width,m_height);
}

/*! distructor */
MyWindow::~MyWindow()
{
  delete []colors;
}

#include "arrangement_2.moc"


/*! main */
int main(int argc, char **argv)
{
  const QString my_title_string("Arrangement Demo with CGAL Qt_widget");
  QApplication app( argc, argv );
  MyWindow widget(707,707); // physical window size
  app.setMainWidget(&widget);
  widget.setCaption(my_title_string);
  widget.setMouseTracking(TRUE);
  widget.show();
  return app.exec();
}

