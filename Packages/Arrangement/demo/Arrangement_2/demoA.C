#include <CGAL/basic.h>

#ifdef CGAL_USE_QT

#include "demo1.h"
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
#include "icons/draw.xpm"
#include "icons/conic_types.xpm"

#include "icons/demo_insert.xpm"
#include "icons/demo_delete.xpm"
#include "icons/demo_snapgrid.xpm"
#include "icons/demo_rayshoot.xpm"
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
#include "icons/demo_arrow_down.xpm"
#include "icons/demo_arrow_up.xpm"




//////////////////////////////////////////////////////////////////////////////
/*! MyWindow constructor
 * \param w - window width
 * \param h - window hight
 */
MyWindow::MyWindow(int w, int h) 
{
  myBar = new QTabWidget(this);
  setCentralWidget(myBar);
  m_width = w;
  m_height = h;
  tab_number = 1;
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
  
  setConicTraits = new QAction("Conic Traits",
                               QPixmap( (const char**)conic_xpm ),
                               "&Conic Traits", 0 , traitsGroup,
                               "Conic Traits" );
  setConicTraits->setToggleAction( TRUE );
  
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
  
  rayShootingMode = new QAction("RayShooting",
                                QPixmap( (const char**)demo_rayshoot_up_xpm ),
                                "&Ray Shooting", 0 , modeGroup,
                                "Ray Shooting" );
  rayShootingMode->setToggleAction( TRUE );

  dragMode = new QAction("Drag", QPixmap( (const char**)hand_xpm ),
                         "&Drag", 0 , modeGroup, "Drag" );
  dragMode->setToggleAction( TRUE );

  mergeMode = new QAction("Merge", QPixmap( (const char**)merge_xpm ),
                         "&Merge", 0 , modeGroup, "Merge" );
  mergeMode->setToggleAction( TRUE );

  splitMode = new QAction("Split", QPixmap( (const char**)split_xpm ),
                         "&Split", 0 , modeGroup, "Split" );
  splitMode->setToggleAction( TRUE );

  // zoom in
  zoominBt = new QAction("Zoom in", QPixmap( (const char**)zoomin_xpm ),
                         "&Zoom in", 0 , this, "Zoom in" );
  // zoom out
  zoomoutBt = new QAction("Zoom out", QPixmap( (const char**)zoomout_xpm ),
                          "&Zoom out", 0 , this, "Zoom out" );
  
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
  
  

  //create a timer for checking if somthing changed
  QTimer *timer = new QTimer( this );
  connect( timer, SIGNAL(timeout()),
           this, SLOT(timer_done()) );
  timer->start( 200, FALSE );
  
  // file menu
  QPopupMenu * file = new QPopupMenu( this );
  menuBar()->insertItem( "&File", file );
  file->insertItem("&Open Segment File...", this, SLOT(fileOpenSegment()));
  file->insertItem("&Open Polyline File...", this, SLOT(fileOpenPolyline()));
  file->insertItem("&Open Conic File...", this, SLOT(fileOpenConic()));
  file->insertItem("&Open Segment Pm File...", this, SLOT(fileOpenSegmentPm()));
  file->insertItem("&Open Polyline Pm File...", this, SLOT(fileOpenPolylinePm()));
  //file->insertItem("&Open Conic Pm File", this, SLOT(fileOpenConicPm()));
  file->insertItem("&Save...", this, SLOT(fileSave()));
  file->insertItem("&Save As...", this, SLOT(fileSaveAs()));
  file->insertItem("&Save to ps...", this, SLOT(fileSave_ps()));
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
  tab->insertItem("Add &Conic Tab", this, SLOT(add_conic_tab()));
  tab->insertSeparator();
  tab->insertItem("Remove &Tab", this, SLOT(remove_tab()));
  menuBar()->insertSeparator();
  
  // mode menu
  QPopupMenu * mode = new QPopupMenu( this );
  menuBar()->insertItem( "&Mode", mode );
  insertMode->addTo( mode );
  deleteMode->addTo( mode );
  pointLocationMode->addTo( mode );
  rayShootingMode->addTo( mode );
  dragMode->addTo( mode );
  mergeMode->addTo( mode );
  splitMode->addTo( mode );
  menuBar()->insertSeparator();
  
  // snap mode menu
  QPopupMenu * snap_mode = new QPopupMenu( this );
  menuBar()->insertItem( "&Snap mode", snap_mode );
  setSnapMode->addTo(snap_mode); 
  setGridSnapMode->addTo(snap_mode);
  menuBar()->insertSeparator();
  
  // options menu
  QPopupMenu * options = new QPopupMenu( this );
  menuBar()->insertItem( "&Options", options );
  setSegmentTraits->addTo(options); 
  setPolylineTraits->addTo(options); 
  setConicTraits->addTo(options); 
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
  options->insertItem("Background Color...", this, SLOT(backGroundColor()));
  options->insertSeparator(); 
  options->insertItem("Planar Map Color...", this, SLOT(changePmColor()));
  options->insertSeparator();
  options->insertItem("Ray-Shooting Direction...", this, 
	                                       SLOT(rayShootingDirection()));
  
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
  rayShootingMode->addTo( modeTools );
  mergeMode->addTo( modeTools );
  splitMode->addTo( modeTools );
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
  setConicTraits->addTo( traitsTool );
  traitsTool->addSeparator();
  
  QToolBar *zoomTool = new QToolBar( this, "zoom" );
  zoomTool->setLabel( "Zoom" );
  zoomTool->addSeparator();
  zoomoutBt->addTo( zoomTool );
  zoominBt->addTo( zoomTool );
  zoomTool->addSeparator();

  conicTypeTool = new QToolBar( this, "conic type" );
  conicTypeTool->setLabel( "Conic Type" );
  conicTypeTool->addSeparator();
  setSegment->addTo( conicTypeTool );
  setCircle->addTo( conicTypeTool );
  setEllipse->addTo( conicTypeTool );
  setParabola->addTo( conicTypeTool );
  setHyperbola->addTo( conicTypeTool );

  
  connect( zoomoutBt, SIGNAL( activated () ) , 
       this, SLOT( zoomout() ) );
  
  connect( zoominBt, SIGNAL( activated () ) , 
       this, SLOT( zoomin() ) );
  
  // connect mode group
  connect( modeGroup, SIGNAL( selected(QAction*) ), 
           this, SLOT( updateMode(QAction*) ) );
  
  // connect Traits Group
  connect( traitsGroup, SIGNAL( selected(QAction*) ),
           this, SLOT( updateTraitsType(QAction*) ) );
  
  // connect Conic Type Group
  connect( conicTypeGroup, SIGNAL( selected(QAction*) ),
           this, SLOT( updateConicType(QAction*) ) );

  // connect Snap Mode 
  
  connect( setSnapMode, SIGNAL( toggled( bool ) ) , 
           this, SLOT( updateSnapMode( bool ) ) );
  
  connect( setGridSnapMode, SIGNAL( toggled( bool ) ) , 
       this, SLOT( updateGridSnapMode( bool ) ) );
    
  // connect the change of current tab
  connect( myBar, SIGNAL( currentChanged(QWidget * )  ),
           this, SLOT( update() ) );
  
  //state flag 
  old_state = 0;
  add_segment_tab();
  resize(m_width,m_height);  
}

/*! distructor */
MyWindow::~MyWindow()
{}

/*! something_changed - change the current page's current_state
 *  and thats makes him redraw
 */
void MyWindow::something_changed()
{
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*, 
  // as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab     *w_demo_p = 
    dynamic_cast<Qt_widget_base_tab  *> (myBar->currentPage());
  
  w_demo_p->current_state++;
}

/*! get_new_object - get a point object from current page 
 * \param obj - a CGAL object
 */
void MyWindow::get_new_object(CGAL::Object obj)
{
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*, 
  // as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab    *w_demo_p = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  
  // point location
  Coord_point p;
  if(CGAL::assign(p,obj)) {
    
    w_demo_p->pl_point = p;
    //something_changed();
  }
  something_changed();
}

/*! about - message box about the demo */
void MyWindow::about()
{
  QMessageBox::about( this, "About",
                      "This is a demo for the Arrangement package\n"
                      "Copyright CGAL @2003");
}

/*! aboutQt - message box about Qt */
void MyWindow::aboutQt()
{
  QMessageBox::aboutQt( this, "About Qt" );
}

/*! howto - help menue */
void MyWindow::howto()
{
  QString home;
  home = "help/index.html";
  CGAL::Qt_help_window * help =
    new CGAL::Qt_help_window(home, ".", 0, "help viewer");
  help->resize(400, 400);
  help->setCaption("Demo HowTo");
  help->show();
}

/*! initialize the widget */
void MyWindow::init(Qt_widget_base_tab *widget)
{
  connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),
          this, SLOT(get_new_object(CGAL::Object)));
  widget->attach(testlayer);
  widget->setCursor(QCursor( QPixmap( (const char**)small_draw_xpm)));
  rayShootingMode->setIconSet(QPixmap((const char**)demo_rayshoot_up_xpm ));
  QColor c = Qt::white;
  c.setRgb(236,236,236);
  widget->setBackgroundColor(c);
  tab_number++;
  number_of_tabs++;
  // add the new widget to myBar
  myBar->insertTab( widget, 
                    QString("Arr " + QString::number( widget->index ) ),
                    widget->index );
  myBar->setCurrentPage(myBar->indexOf(widget));
 
  update();
  something_changed();
  
}

/*! add a tab widget with segment traits */
void MyWindow::add_segment_tab()
{
  Qt_widget_demo_tab<Segment_tab_traits> *widget = 
    new Qt_widget_demo_tab<Segment_tab_traits>
    (SEGMENT_TRAITS , this, tab_number);
  init(widget);
}

/*! add a tab widget with polyline traits */
void MyWindow::add_polyline_tab()
{
  Qt_widget_demo_tab<Polyline_tab_traits> *widget = 
    new Qt_widget_demo_tab<Polyline_tab_traits>
    (POLYLINE_TRAITS , this, tab_number);
  init(widget);
}

/*! add a tab widget with conic traits */
void MyWindow::add_conic_tab()
{
  Qt_widget_demo_tab<Conic_tab_traits> *widget = 
    new Qt_widget_demo_tab<Conic_tab_traits>
    (CONIC_TRAITS , this , tab_number);
  init(widget);
}

/*! add a tab widget with polyline traits */
void MyWindow::remove_tab()
{
  if (number_of_tabs > 1)
  {
    myBar->removePage(myBar->currentPage());
    number_of_tabs--;
  }
  else
  {
    QMessageBox::information( this, "Remove Tab",
                              "Can not remove last tab");
  }
}

/*! redraw the widget when timer ends */
void MyWindow::timer_done()
{
  // We peform downcasting from QWigdet* to Qt_widget_demo_tab*, 
  // as we know that only
  // Qt_widget_demo_tab objects are stored in the tab pages.
  Qt_widget_base_tab    *w_demo_p = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  
  if(old_state!=w_demo_p->current_state){
    w_demo_p->redraw();
    old_state = w_demo_p->current_state;
  }
}    

/*! update traits type
 * \param action the new traits type
 */
void MyWindow::updateTraitsType( QAction *action )
{
  Qt_widget_base_tab *old_widget = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  Qt_widget_base_tab * widget = 0;
  
  if (action == setSegmentTraits)
  {
    if (old_widget->traits_type == SEGMENT_TRAITS) return;
    widget = new Qt_widget_demo_tab<Segment_tab_traits>
      (SEGMENT_TRAITS , this);
  }
  else if (action == setPolylineTraits)
  {
    if (old_widget->traits_type == POLYLINE_TRAITS) return;
    widget = new Qt_widget_demo_tab<Polyline_tab_traits>
      (POLYLINE_TRAITS , this);
  }
  else if (action == setConicTraits)
  {
    if (old_widget->traits_type == CONIC_TRAITS) return;
    widget = new Qt_widget_demo_tab<Conic_tab_traits>
      (CONIC_TRAITS , this);
  }
  
  if( !old_widget->empty ) // pm is not empty
  {
    switch( QMessageBox::warning( this, "Update Traits Type",
        "This action will destroy the current planar map.\n"
        "Do you want to continue ?",
        "Yes",
        "No", 0, 0, 1 ) ) {
      case 0: 
          // continue
          break;
      case 1: // The user clicked the Quit or pressed Escape
		  update();
          something_changed();
          return;
          break;
    }
  }

  int old_index = old_widget->index;
  int index = myBar->currentPageIndex();
  QString label = myBar->label(index);
  myBar->removePage(myBar->currentPage());
  
  //initialize the new tab widget
  *widget << CGAL::LineWidth(2); // << CGAL::BackgroundColor (CGAL::WHITE);
  QColor c = Qt::white;
  c.setRgb(236,236,236);
  widget->setBackgroundColor(c);
  widget->set_window(-10, 10, -10, 10);
  widget->setMouseTracking(TRUE);
  connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),
          this, SLOT(get_new_object(CGAL::Object)));
  widget->attach(testlayer);
  widget->m_line_width = 2;
  widget->index = old_index;
  widget->pm_color = widget->colors[old_index];
  widget->setCursor(QCursor( QPixmap( (const char**)small_draw_xpm)));
  rayShootingMode->setIconSet(QPixmap((const char**)demo_rayshoot_up_xpm ));

  // add the new widget to myBar
  myBar->insertTab( widget, label , index );
  
  myBar->setCurrentPage(index);
    
  update();
  something_changed();
  
}

/*! change the buttons stste according to the traits type */
void MyWindow::setTraits( TraitsType t )
{
  switch ( t ) {
   case SEGMENT_TRAITS:
    setSegmentTraits->setOn( TRUE );
	conicTypeTool->hide();
    break;
   case POLYLINE_TRAITS:
    setPolylineTraits->setOn( TRUE );
	conicTypeTool->hide();
    break;
   case CONIC_TRAITS:
    setConicTraits->setOn( TRUE );
	conicTypeTool->show();
    break;
  }
}

/*! update widget conic type
 * \param action - the new conic type 
 */
void MyWindow::updateConicType( QAction *action )
{
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*,
  // as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab    *w_demo_p = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  if ( action == setCircle ) 
    w_demo_p->conic_type = CIRCLE;
  else if ( action == setSegment ) 
    w_demo_p->conic_type = SEGMENT;
  else if ( action == setEllipse ) 
    w_demo_p->conic_type = ELLIPSE;
  else if ( action == setParabola ) 
    w_demo_p->conic_type = PARABOLA;
  else if ( action == setHyperbola ) 
    w_demo_p->conic_type = HYPERBOLA;
}

/*! change the buttons stste according to the traits type */
void MyWindow::setConicType( ConicType t )
{
  switch ( t ) {
   case CIRCLE:
    setCircle->setOn( TRUE );
    break;
   case SEGMENT:
    setSegment->setOn( TRUE );
    break;
   case ELLIPSE:
    setEllipse->setOn( TRUE );
    break;
   case PARABOLA:
    setParabola->setOn( TRUE );
    break;
   case HYPERBOLA:
    setHyperbola->setOn( TRUE );
    break;
  }
}

/*! update snap mode - change the snap mode and the relevant buttons 
 *  state after the snap mode button was clicked
 * \param on - true if the user activate the button.
 */
void MyWindow::updateSnapMode( bool on )
{
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*, 
  // as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab    *w_demo_p = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  
  if (on)
  {
    setGridSnapMode->setEnabled( TRUE );
	if (w_demo_p->snap_mode == GRID)
	  setGridSnapMode->setOn( TRUE );
	else
	{
      setGridSnapMode->setOn( FALSE );
      w_demo_p->snap_mode = POINT;
	}
    w_demo_p->snap = true;
  }
  else
  {
    SnapMode old = w_demo_p->snap_mode;
    setGridSnapMode->setOn( FALSE );
	setSnapMode->setOn( FALSE );
    setGridSnapMode->setEnabled( FALSE );
    w_demo_p->snap = false;
    w_demo_p->snap_mode = NONE;
    if (old == GRID)
      something_changed();
  }
}

/*! update grid snap mode - we have 2 states of snap
 *  - to grid points 
 *  - to closest point in the planar map
 * \param on - true when we choos grid mode
 */
void MyWindow::updateGridSnapMode( bool on )
{
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*,
  // as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab    *w_demo_p = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  
  if (on && w_demo_p->snap)
    w_demo_p->snap_mode = GRID;
  else
    w_demo_p->snap_mode = POINT;
  
  something_changed();
}

/*! update widget mode
 * \param action - the new widget mode
 */
void MyWindow::updateMode( QAction *action )
{
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*,
  // as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab    *w_demo_p = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  
  if ( action == insertMode ) 
  {
    w_demo_p->mode = INSERT;
    w_demo_p->setCursor(QCursor( QPixmap( (const char**)small_draw_xpm)));
    something_changed();
  }
  else if ( action == deleteMode ) 
  {
    w_demo_p->mode = DELETE;
    w_demo_p->setCursor(QCursor( QPixmap( (const char**)delete_xpm)));
    something_changed();
  }
  
  else if ( action == pointLocationMode ) 
  {
    w_demo_p->mode = POINT_LOCATION;
    w_demo_p->setCursor(Qt::CrossCursor);
    //w_demo_p->setCursor(QCursor( QPixmap( (const char**)point_xpm)));
    current_label = point_location_label;
  }
  else if ( action == rayShootingMode ) 
  {
    w_demo_p->mode = RAY_SHOOTING;
	if (w_demo_p->ray_shooting_direction)
      w_demo_p->setCursor(
	  QCursor(QPixmap((const char**)demo_arrow_up_xpm)));
	else
      w_demo_p->setCursor(
	  QCursor(QPixmap((const char**)demo_arrow_down_xpm)));
  }
  else if ( action == dragMode ) 
  {
    w_demo_p->mode = DRAG;
    w_demo_p->setCursor(QCursor( QPixmap( (const char**)hand_xpm)));
    something_changed();
  }
  else if ( action == mergeMode ) 
  {
    w_demo_p->mode = MERGE;
    w_demo_p->setCursor(Qt::IbeamCursor );
    something_changed();
  }
  else if ( action == splitMode ) 
  {
    w_demo_p->mode = SPLIT;
    w_demo_p->setCursor(Qt::SplitHCursor  );
    something_changed();
  }
}

/*! a dialog form to set the rayShooting Diraction. */
void MyWindow::rayShootingDirection()
{
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*, 
  // as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab * w_demo_p = 
    dynamic_cast<Qt_widget_base_tab  *> (myBar->currentPage());
  RayShootingOptionsForm *form = new RayShootingOptionsForm();
  if ( form->exec() ) 
  {
    QString type = form->arrComboBox1->currentText();
    if (strcmp(type,"Up") == 0)
	{
      w_demo_p->ray_shooting_direction = true;
      rayShootingMode->setIconSet(QPixmap((const char**)demo_rayshoot_up_xpm ));
	  if (w_demo_p->mode == RAY_SHOOTING)
	    w_demo_p->setCursor(
	    QCursor(QPixmap((const char**)demo_arrow_up_xpm)));
	}
    else if (strcmp(type,"Down") == 0)
	{
      w_demo_p->ray_shooting_direction = false;
      rayShootingMode->setIconSet( 
		  QPixmap((const char**)demo_rayshoot_down_xpm ));
	  if (w_demo_p->mode == RAY_SHOOTING)
	    w_demo_p->setCursor(
	    QCursor(QPixmap((const char**)demo_arrow_down_xpm)));
	}
  }
}

#endif
