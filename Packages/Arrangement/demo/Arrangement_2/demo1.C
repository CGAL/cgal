#include "demo1.h"
#include "forms.h"
#include "qt_layer.h"
#include "demo_tab.h"

////////////////////////////////////////////////////////////////////////
#ifndef CGAL_USE_QT
#include <iostream>
int main(int, char*)
{
  std::cout << "Sorry, this demo needs QT...";
  std::cout << std::endl;
  
  return 0;
}

#else

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


const QString my_title_string("Arrangement Demo with CGAL Qt_widget");

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
  
  colors[1] = Qt::blue;
  colors[2] = Qt::gray;
  colors[3] = Qt::green;
  colors[4] = Qt::cyan;
  colors[5] = Qt::magenta;
  colors[6] = Qt::black;
  colors[7] = Qt::darkGreen;
  colors[8] = Qt::darkBlue;
  colors[9] = Qt::darkMagenta;
  colors[10] = Qt::darkCyan;
  colors[11] = Qt::yellow;

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
                                QPixmap( (const char**)rayshoot_xpm ),
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
                                 QPixmap( (const char**)mycircle_xpm ),
                                 "&Circle", 0 ,conicTypeGroup,
                                 "Circle" );
  setCircle->setToggleAction( TRUE );
  setSegment = new QAction("Segment",
                                 QPixmap( (const char**)segment_xpm ),
                                 "&Segment", 0 ,conicTypeGroup,
                                 "Segment" );
  setSegment->setToggleAction( TRUE );
  setEllipse = new QAction("Ellipse",
                                 QPixmap( (const char**)ellipse_xpm ),
                                 "&Ellipse", 0 ,conicTypeGroup,
                                 "Ellipse" );
  setEllipse->setToggleAction( TRUE );
  setParabola = new QAction("Parabola",
                                 QPixmap( (const char**)parabula_xpm ),
                                 "&Parabola", 0 ,conicTypeGroup,
                                 "Parabola" );
  setParabola->setToggleAction( TRUE );
  setHyperbola = new QAction("Hyperbola",
                                 QPixmap( (const char**)hyperbula_xpm ),
                                 "&Hyperbola", 0 ,conicTypeGroup,
                                 "Hyperbola" );
  setHyperbola->setToggleAction( TRUE );
  
  

  //create a timer for checking if somthing changed
  QTimer *timer = new QTimer( this );
  connect( timer, SIGNAL(timeout()),
           this, SLOT(timer_done()) );
  timer->start( 200, FALSE );
  
  // file menu
  QPopupMenu * file = new QPopupMenu( this );
  menuBar()->insertItem( "&File", file );
  file->insertItem("&Open Segment File", this, SLOT(fileOpenSegment()));
  file->insertItem("&Open Polyline File", this, SLOT(fileOpenPolyline()));
  file->insertItem("&Open Conic File", this, SLOT(fileOpenConic()));
  file->insertItem("&Open Segment Pm File", this, SLOT(fileOpenSegmentPm()));
  file->insertItem("&Open Polyline Pm File", this, SLOT(fileOpenPolylinePm()));
  //file->insertItem("&Open Conic Pm File", this, SLOT(fileOpenConicPm()));
  file->insertItem("&Save", this, SLOT(fileSave()));
  file->insertItem("&Save As", this, SLOT(fileSaveAs()));
  file->insertItem("&Save to ps", this, SLOT(fileSave_ps()));
  file->insertSeparator();
  file->insertItem("&Print", this , SLOT(print()));
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
  options->insertItem("Overlay", this, SLOT(overlay_pm()));
  options->insertSeparator();
  options->insertItem("Properties", this, SLOT(properties()));
  options->insertSeparator();
  options->insertItem("Show Grid", this, SLOT(showGrid()));
  options->insertItem("Hide Grig", this, SLOT(hideGrid()));
  options->insertSeparator();
  options->insertItem("Conic Type", this, SLOT(conicType()));
  options->insertSeparator();
  options->insertItem("Background Color", this, SLOT(backGroundColor()));
  options->insertSeparator();
  options->insertItem("RayShooting Diraction", this, 
	                                       SLOT(rayShootingDirection()));
  
  // help menu
  QPopupMenu * help = new QPopupMenu( this );
  menuBar()->insertItem( "&Help", help );
  help->insertItem("How To", this, SLOT(howto()), Key_F1);
  help->insertSeparator();
  help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
  help->insertItem("About &Qt", this, SLOT(aboutQt()) );
  
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
  QMessageBox::about( this, my_title_string,
                      "This is a demo for the Arrangement package\n"
                      "Copyright CGAL @2003");
}

/*! aboutQt - message box about Qt */
void MyWindow::aboutQt()
{
  QMessageBox::aboutQt( this, my_title_string );
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
  rayShootingMode->setIconSet(QPixmap((const char**)ray_shooting_xpm ));
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
  resize(m_width,m_height);
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
    QMessageBox::information( this, my_title_string,
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
  Qt_widget_base_tab *widget;
  
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
  widget->setCursor(QCursor( QPixmap( (const char**)small_draw_xpm)));
  rayShootingMode->setIconSet(QPixmap((const char**)ray_shooting_xpm ));

  // add the new widget to myBar
  myBar->insertTab( widget, label , index );
  
  myBar->setCurrentPage(index);
  
  resize(m_width,m_height);
    
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
    setGridSnapMode->setOn( FALSE );
    w_demo_p->snap_mode = POINT;
    w_demo_p->snap = true;
  }
  else
  {
    SnapMode old = w_demo_p->snap_mode;
    setGridSnapMode->setOn( FALSE );
    setGridSnapMode->setEnabled( FALSE );
    w_demo_p->snap = false;
    w_demo_p->snap_mode = NONE;
    if (old == GRID)
      something_changed();
  }
}

/*! update grid snap mode - we have to states of snap
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
	  QCursor(QPixmap((const char**)small_ray_shooting_xpm)));
	else
      w_demo_p->setCursor(
	  QCursor(QPixmap((const char**)small_ray_shooting_down_xpm)));
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

/*! set the buttons state according to the current mode */
void MyWindow::setMode( Mode m )
{
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*, 
  // as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab    *w_demo_p = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  
  w_demo_p->mode = m;
  switch ( m ) {
   case INSERT: insertMode->setOn( TRUE ); break;
   case DELETE: deleteMode->setOn( TRUE ); break;
   case POINT_LOCATION: pointLocationMode->setOn( TRUE ); break;
   case RAY_SHOOTING: rayShootingMode->setOn( TRUE ); break;
   case DRAG: dragMode->setOn( TRUE ); break;
   case MERGE: mergeMode->setOn( TRUE ); break;
   case SPLIT: splitMode->setOn( TRUE ); break;
  }
}

/*! update all modes */
void MyWindow::update()
{
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*, 
  // as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab    *w_demo_p = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  setMode( w_demo_p->mode );
  updateSnapMode( false );
  setTraits( w_demo_p->traits_type );
  setConicType( w_demo_p->conic_type ); 
}

/*! zoom in - enlarge the picture */
void MyWindow::zoomin()
{
  Qt_widget_base_tab    *w_demo_p = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  w_demo_p->zoom(m_scailing_factor);
}

/*! zoom out - lessen the picture */
void MyWindow::zoomout()
{
  Qt_widget_base_tab    *w_demo_p = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  w_demo_p->zoom(1/m_scailing_factor);
}

/*! properties form dialog */
void MyWindow::properties()
{
  Qt_widget_base_tab    *w_demo_p = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  
  PropertiesForm *optionsForm = 
    new PropertiesForm( myBar , this ,number_of_tabs );
  
  if ( optionsForm->exec() ) 
  {    
    m_width = optionsForm->box1->value();
    m_height = optionsForm->box2->value();
    w_demo_p->m_line_width = optionsForm->box3->value();
	w_demo_p->m_vertex_width = optionsForm->box8->value();
    double new_factor = static_cast<double> (optionsForm->box4->value());
    QString paint_mode = optionsForm->box5->currentText();
    w_demo_p->cube_size = optionsForm->box6->value();
    if (strcmp(paint_mode,"Different Colors At Overlay") == 0)
      colors_flag = true;
    else
      colors_flag = false;
	QString remove_mode = optionsForm->box7->currentText();
	if (strcmp(remove_mode,"Remove All Original Curve") == 0)
      w_demo_p->remove_org_curve = true;
    else
      w_demo_p->remove_org_curve = false;
    m_scailing_factor = (new_factor/10);
    resize(m_width,m_height);
    w_demo_p->redraw();
    something_changed();
  }
  delete optionsForm;
}
/*! open a segment file and add new tab */
void MyWindow::fileOpenSegment()
{
  Qt_widget_base_tab *w_demo_p = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  bool flag = (w_demo_p->traits_type == SEGMENT_TRAITS);
  if( w_demo_p->empty ) // pm is empty
  {
    updateTraitsType( setSegmentTraits );
    fileOpen(true);
  }
  else
  {
    FileOpenOptionsForm
         *form = new FileOpenOptionsForm(flag);
    if ( form->exec() ) 
    {
        int id = form->buttonGroup->id(form->buttonGroup->selected());
		switch ( id ) 
		{
          case 0: // open file in a new tab
            add_segment_tab();
            fileOpen();      
            break;
          case 1: // open file in current tab (delete current Pm)
            updateTraitsType( setSegmentTraits );
            fileOpen(true);
          break;
          case 2: // merge file into current tab
            fileOpen();
          break;
		}// switch
	}// if
  }
}// fileOpenSegment

/*! open a segment file and add new tab */
void MyWindow::fileOpenSegmentPm()
{
  Qt_widget_base_tab *w_demo_p = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  if( w_demo_p->empty ) // pm is empty
  {    
    updateTraitsType( setSegmentTraits );
    fileOpenPm();
  }
  else
  {
    FileOpenOptionsForm
         *form = new FileOpenOptionsForm(false);
    if ( form->exec() ) 
    {
        int id = form->buttonGroup->id(form->buttonGroup->selected());
		switch ( id ) 
		{
          case 0: // open file in a new tab
            add_segment_tab();
            fileOpenPm();      
            break;
          case 1: // open file in current tab (delete current Pm)
            updateTraitsType( setSegmentTraits );
            fileOpenPm();
          break;          
		}// switch
	}// if
  }
}// fileOpenSegment

/*! open a polyline file and add new tab */
void MyWindow::fileOpenPolyline()
{
  Qt_widget_base_tab *w_demo_p = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  bool flag = (w_demo_p->traits_type == POLYLINE_TRAITS);
  if( w_demo_p->empty ) // pm is empty
  {    
    updateTraitsType( setPolylineTraits );
    fileOpen(true);	
  }
  else
  {
    FileOpenOptionsForm
         *form = new FileOpenOptionsForm(flag);
    if ( form->exec() ) 
    {
        int id = form->buttonGroup->id(form->buttonGroup->selected());
		switch ( id ) 
		{
          case 0: // open file in a new tab
            add_polyline_tab();
            fileOpen();      
            break;
          case 1: // open file in current tab (delete current Pm)
            updateTraitsType( setPolylineTraits );
            fileOpen(true);
          break;
          case 2: // merge file into current tab
            fileOpen();
          break;
		}// switch
	}// if
  } 
}// fileOpenPolyline

/*! open a polyline file and add new tab */
void MyWindow::fileOpenPolylinePm()
{
  Qt_widget_base_tab *w_demo_p = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  if( w_demo_p->empty ) // pm is empty
  {
    updateTraitsType( setPolylineTraits );
    fileOpenPm();
  }
  else
  {
    FileOpenOptionsForm
         *form = new FileOpenOptionsForm(false);
    if ( form->exec() ) 
    {
        int id = form->buttonGroup->id(form->buttonGroup->selected());
		switch ( id ) 
		{
          case 0: // open file in a new tab
            add_polyline_tab();
            fileOpenPm();      
            break;
          case 1: // open file in current tab (delete current Pm)
            updateTraitsType( setPolylineTraits );
            fileOpenPm();
          break;          
		}// switch
	}// if
  }
}// fileOpenPolylinePm

/*! open a polyline file and add new tab */
void MyWindow::fileOpenConic()
{
  Qt_widget_base_tab *w_demo_p = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  bool flag = (w_demo_p->traits_type == CONIC_TRAITS);
  if( w_demo_p->empty ) // pm is empty
  {
    updateTraitsType( setConicTraits );
    fileOpen(true);
  }
  else
  {
    FileOpenOptionsForm
         *form = new FileOpenOptionsForm(flag);
    if ( form->exec() ) 
    {
        int id = form->buttonGroup->id(form->buttonGroup->selected());
		switch ( id ) 
		{
          case 0: // open file in a new tab
            add_conic_tab();
            fileOpen();      
            break;
          case 1: // open file in current tab (delete current Pm)
            updateTraitsType( setConicTraits );
            fileOpen(true);
          break;
          case 2: // merge file into current tab
            fileOpen();
          break;
		}// switch
	}// if
  }// else  
}// fileOpenConic

/*! open a file */
void MyWindow::fileOpen( bool clear_flag )
{
  
  QString filename =
    QFileDialog::getOpenFileName(QString::null, 0, this,
                                 "file open", "Demo -- File Open" );
  if ( !filename.isEmpty() )
    load( filename , clear_flag);
  else
    statusBar()->message( "File Open abandoned", 2000 );
}

/*! open a Pm file */
void MyWindow::fileOpenPm()
{  
  QString filename =
    QFileDialog::getOpenFileName(QString::null, 0, this,
                                 "file open", "Demo -- File Open" );
  if ( filename.isEmpty() )
  {
    statusBar()->message( "File Open abandoned", 2000 );
	return;
  }
  std::ifstream inputFile(filename);
  // Creates an ifstream object named inputFile
  if (! inputFile.is_open()) // Always test file open
  {
    std::cout << "Error opening input file" << std::endl;
    return;
  }
 
  Qt_widget_base_tab    *w_demo_p1 = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());

  QCursor old = w_demo_p1->cursor();
  w_demo_p1->setCursor(Qt::WaitCursor);

  switch ( w_demo_p1->traits_type ) {
   case SEGMENT_TRAITS:
    {
      Qt_widget_demo_tab<Segment_tab_traits> *w_demo_p = 
       static_cast<Qt_widget_demo_tab<Segment_tab_traits> *> 
       (myBar->currentPage());
      w_demo_p->m_curves_arr.read(inputFile);
	  if( w_demo_p->m_curves_arr.number_of_vertices() == 0 )
    	w_demo_p->empty = false;
      else 
        w_demo_p->empty = true;
     break;
    }
   case POLYLINE_TRAITS: // dosen't work !!
    {
      Qt_widget_demo_tab<Polyline_tab_traits> *w_demo_p = 
       static_cast<Qt_widget_demo_tab<Polyline_tab_traits> *> 
       (myBar->currentPage());
       w_demo_p->m_curves_arr.read(inputFile);   
	   if( w_demo_p->m_curves_arr.number_of_vertices() == 0 )
	     w_demo_p->empty = false;
       else 
         w_demo_p->empty = true;
	   break;
    }
   case CONIC_TRAITS: // dosen't work !!
    {
  //   Qt_widget_demo_tab<Conic_tab_traits> *w_demo_p = 
  //     static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> 
  //     (myBar->currentPage());
	 //w_demo_p->m_curves_arr.read(inputFile);
     break;
    }
  }  
  
  inputFile.close();
  
  w_demo_p1->set_window(w_demo_p1->bbox.xmin() , w_demo_p1->bbox.xmax() , 
                     w_demo_p1->bbox.ymin() , w_demo_p1->bbox.ymax());
  
  inputFile.close();
  w_demo_p1->setCursor(old);
  something_changed();

  setCaption( QString( "Planar Map -- %1" ).arg( m_filename ) );
  statusBar()->message( QString( "Opened \'%1\'" ).arg( m_filename ), 2000 );

}

/*! open a polyline or conic file 
 * \param filename - name of the file
 */
void MyWindow::load( const QString& filename , bool clear_flag )
{
  std::ifstream inputFile(filename);
  // Creates an ofstream object named inputFile
  if (! inputFile.is_open()) // Always test file open
  {
    std::cout << "Error opening input file" << std::endl;
    return;
  }
  
  Qt_widget_base_tab    *w_demo = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());

  QCursor old = w_demo->cursor();
  w_demo->setCursor(Qt::WaitCursor);
   
  if (w_demo->traits_type == CONIC_TRAITS)
  {
    Qt_widget_demo_tab<Conic_tab_traits>    *w_demo_p = 
      static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> 
      (myBar->currentPage());
    if (clear_flag)
        w_demo_p->m_curves_arr.clear();
    char dummy[256];
    Pm_base_conic_2* cv;
    int count;
    inputFile >> count;
    inputFile.getline(dummy, sizeof(dummy));
    for (int i = 0; i < count; i++) 
    {
      cv = new Pm_base_conic_2();
      ReadCurve(inputFile, *cv);
      
      Curve_conic_data cd;
      cd.m_type = Curve_conic_data::LEAF;
      cd.m_index = w_demo_p->index;
      cd.m_ptr.m_curve = cv;
      
      w_demo_p->m_curves_arr.insert(Pm_conic_2( *cv , cd));
      
      CGAL::Bbox_2 curve_bbox = cv->bbox();
      if (i == 0)
        w_demo->bbox = curve_bbox;
      else
        w_demo->bbox = w_demo->bbox + curve_bbox;
	}
	if( w_demo_p->m_curves_arr.number_of_vertices() == 0 )
	  w_demo_p->empty = false;
    else 
      w_demo_p->empty = true;    
  }
  
  else if (w_demo->traits_type == POLYLINE_TRAITS)
  {
    Qt_widget_demo_tab<Polyline_tab_traits>    *w_demo_p = 
      static_cast<Qt_widget_demo_tab<Polyline_tab_traits> *> 
      (myBar->currentPage());
    if (clear_flag)
        w_demo_p->m_curves_arr.clear();
    
    int num_polylines, num_segments;
    int ix, iy;
    std::vector<Pm_pol_point_2> points;
    int i, j;
    
    inputFile >> num_polylines;
    for (i = 0; i < num_polylines; i++) 
    {
      inputFile >> num_segments;
      points.clear();
      for (j = 0; j < num_segments; j++)
      {
        inputFile >> ix >> iy;
        points.push_back (Pm_pol_point_2(NT(ix),NT(iy)));
      }
      
      Pm_base_pol_2 *base_polyline = 
        new Pm_base_pol_2(points.begin(), points.end());
      
      CGAL::Bbox_2 curve_bbox = base_polyline->bbox();
      if (i == 0)
        w_demo->bbox = curve_bbox;
      else
        w_demo->bbox = w_demo->bbox + curve_bbox;
      
      Curve_pol_data cd;
      cd.m_type = Curve_pol_data::LEAF;
      cd.m_index = w_demo_p->index;
      cd.m_ptr.m_curve = base_polyline;
      
      w_demo_p->m_curves_arr.insert(Pm_pol_2( *base_polyline , cd));
    }
    if( w_demo_p->m_curves_arr.number_of_vertices() == 0 )
      w_demo_p->empty = false;
    else 
      w_demo_p->empty = true;
  }
  
  else if (w_demo->traits_type == SEGMENT_TRAITS)
  {
    Qt_widget_demo_tab<Segment_tab_traits>    *w_demo_p = 
      static_cast<Qt_widget_demo_tab<Segment_tab_traits> *> 
      (myBar->currentPage());
    if (clear_flag)
        w_demo_p->m_curves_arr.clear();
    
    int count;
    inputFile >> count;
    
    int i;
    for (i = 0; i < count; i++) {
      NT x0, y0, x1, y1;
      inputFile >> x0 >> y0 >> x1 >> y1;
    
      Pm_seg_point_2 p1(x0, y0);
      Pm_seg_point_2 p2(x1, y1);

	  Pm_base_seg_2 *base_seg =
		  new Pm_base_seg_2(p1, p2);
	
      CGAL::Bbox_2 curve_bbox = base_seg->bbox();
      if (i == 0)
        w_demo->bbox = curve_bbox;
      else
        w_demo->bbox = w_demo->bbox + curve_bbox;
      
      Curve_data cd;
      cd.m_type = Curve_data::LEAF;
      cd.m_index = w_demo_p->index;
      cd.m_ptr.m_curve = base_seg;
      
      w_demo_p->m_curves_arr.insert(Pm_seg_2( *base_seg , cd));

    }
	if( w_demo_p->m_curves_arr.number_of_vertices() == 0 )
	  w_demo_p->empty = false;
    else 
      w_demo_p->empty = true;
      
  }
  w_demo->set_window(w_demo->bbox.xmin() , w_demo->bbox.xmax() , 
                     w_demo->bbox.ymin() , w_demo->bbox.ymax());
  
  inputFile.close();
  w_demo->setCursor(old);
  
  something_changed();
}

/*! read a conic curve 
 * \param is - input file stream
 * \param cv - will hold the reading curve
 */
void MyWindow::ReadCurve(std::ifstream & is, Pm_base_conic_2 & cv)
{
  // Read a line from the input file.
  char one_line[128];
  
  skip_comments (is, one_line);
  std::string stringvalues(one_line);
  std::istringstream str_line (stringvalues, std::istringstream::in);
  
  // Get the arc type.
  char     type;
  bool     is_circle = false;       // Is this a circle.
  Pm_conic_circle_2 circle;
  CONIC_NT       r, s, t, u, v, w;  // The conic coefficients.
  
  str_line >> type;
  
  // An ellipse (full ellipse or a partial ellipse):
  if (type == 'f' || type == 'F' || type == 'e' || type == 'E')
  {  
    // Read the ellipse (using the format "a b x0 y0"):
    //
    //     x - x0   2      y - y0   2
    //  ( -------- )  + ( -------- )  = 1
    //       a               b
    //
    CONIC_NT     a, b, x0, y0;
    
    str_line >> a >> b >> x0 >> y0;
    
    CONIC_NT     a_sq = a*a;
    CONIC_NT     b_sq = b*b;
    
    if (a == b)
    {
      is_circle = true;
      circle = Pm_conic_circle_2 
        (Pm_conic_point_2 (x0, y0), a*b, CGAL::CLOCKWISE);
    }
    else
    {
      r = b_sq;
      s = a_sq;
      t = 0;
      u = -2*x0*b_sq;
      v = -2*y0*a_sq;
      w = x0*x0*b_sq + y0*y0*a_sq - a_sq*b_sq;
    }
    
    if (type == 'f' || type == 'F')
    {
      // Create a full ellipse (or circle).
      if (is_circle)
        cv = Pm_base_conic_2 (circle);
      else
        cv = Pm_base_conic_2 (r, s, t, u, v, w);
      
      return;
    }
  }
  else if (type == 'h' || type == 'H')
  {
    // Read the hyperbola (using the format "a b x0 y0"):
    //
    //     x - x0   2      y - y0   2
    //  ( -------- )  - ( -------- )  = 1
    //       a               b
    //
    CONIC_NT     a, b, x0, y0;
    
    str_line >> a >> b >> x0 >> y0;
    
    CONIC_NT     a_sq = a*a;
    CONIC_NT     b_sq = b*b;
    
    r = b_sq;
    s= -a_sq;
    t = 0;
    u = -2*x0*b_sq;
    v = 2*y0*a_sq;
    w = x0*x0*b_sq - y0*y0*a_sq - a_sq*b_sq;  
  }
  else if (type == 'p' || type == 'P')
  {
    // Read the parabola (using the format "c x0 y0"):
    //
    //                        2
    //  4c*(y - y0) = (x - x0)
    //
    CONIC_NT     c, x0, y0;
    
    str_line >> c >> x0 >> y0;
    
    r = 1;
    s = 0;
    t = 0;
    u = -2*x0;
    v = -4*c;
    w = x0*x0 + 4*c*y0;
  }
  else if (type == 'c' || type == 'C' || type == 'a' || type == 'A')
  {
    // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
    str_line >> r >> s >> t >> u >> v >> w;
    
    if (type == 'c' || type == 'C')
    {
      // Create a full conic (should work only for ellipses).
      cv = Pm_base_conic_2 (r, s, t, u, v, w);
      return;
    }
  }
  else if (type == 's' || type == 'S')
  {
    // Read a segment, given by its endpoints (x1,y1) and (x2,y2);
    CONIC_NT      x1, y1, x2, y2;
    
    str_line >> x1 >> y1 >> x2 >> y2;
    
    Pm_conic_point_2   source (x1, y1);
    Pm_conic_point_2   target (x2, y2);
    Pm_conic_segment_2 segment (source, target);
    
    // Create the segment.
    cv = Pm_base_conic_2(segment);
    return;
  }
  else if (type == 'i' || type == 'I')
  {
    // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
    str_line >> r >> s >> t >> u >> v >> w;
    
    // Read the approximated source, along with a general conic 
    // <r_1,s_1,t_1,u_1,v_1,w_1> whose intersection with <r,s,t,u,v,w>
    // defines the source.
    CONIC_NT     r1, s1, t1, u1, v1, w1;
    CONIC_NT     x1, y1;
    
    str_line >> x1 >> y1;
    str_line >> r1 >> s1 >> t1 >> u1 >> v1 >> w1;
    
    Pm_conic_point_2   app_source (x1, y1);
    
    // Read the approximated target, along with a general conic 
    // <r_2,s_2,t_2,u_2,v_2,w_2> whose intersection with <r,s,t,u,v,w>
    // defines the target.
    CONIC_NT     r2, s2, t2, u2, v2, w2;
    CONIC_NT     x2, y2;
    
    str_line >> x2 >> y2;
    str_line >> r2 >> s2 >> t2 >> u2 >> v2 >> w2;
    
    Pm_conic_point_2   app_target (x2, y2);
    
    // Create the conic arc.
    cv = Pm_base_conic_2 (r, s, t, u, v ,w,
                          app_source, r1, s1, t1, u1, v1, w1,
                          app_target, r2, s2, t2, u2, v2, w2);
    return;
  }
  else
  {
    std::cerr << "Illegal conic type specification: " << type << "."
              << std::endl;
    return;
  }
  
  // Read the end points of the arc and create it.
  CONIC_NT    x1, y1, x2, y2;
  
  str_line >> x1 >> y1 >> x2 >> y2;
  
  Pm_conic_point_2 source (x1, y1);
  Pm_conic_point_2 target (x2, y2);
  
  // Create the conic (or circular) arc.
  if (is_circle)
  {
    cv = Pm_base_conic_2 (circle,
                          source, target);
  }
  else
  {
    cv = Pm_base_conic_2 (r, s, t, u, v, w,
                          source, target);
  }
  
  return;
}

/*! calls from readCurve to skip comments in the input file */
void MyWindow::skip_comments( std::ifstream& is, char* one_line )
{
  while( !is.eof() )
  {
    is.getline( one_line, 128 );
    if( one_line[0] != '#' )
    {
      break;
    }
  }
}

/*! save file in a different name */
void MyWindow::fileSaveAs()
{
  QString filename =
    QFileDialog::getSaveFileName(QString::null, "Planar Map (*.pm)", this,
                                 "file save as", "Planar Map -- File Save As");
  if ( !filename.isEmpty() ) 
  {
    int answer = 0;
    if ( QFile::exists( filename ) )
      answer = QMessageBox::warning(
                                    this, "Overwrite File",
                                    QString( "Overwrite\n\'%1\'?" ).
                                    arg( filename ),
                                    "&Yes", "&No", QString::null, 1, 1 );
    if ( answer == 0 ) 
    {
      m_filename = filename;
      //updateRecentFiles( filename );
      fileSave();
      return;
    }
  }
  statusBar()->message( "Saving abandoned", 2000 );
}

/*! save file */
void MyWindow::fileSave()
{
  if ( m_filename.isEmpty() ) {
    fileSaveAs();
    return;
  }
  
  std::ofstream outFile(m_filename);
  // Creates an ofstream object named outFile
  if (! outFile.is_open()) // Always test file open
  {
    std::cout << "Error opening input file" << std::endl;
    return;
  }
  
  Qt_widget_base_tab    *w_demo_p1 = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  
  switch ( w_demo_p1->traits_type ) {
   case SEGMENT_TRAITS:
    {
     Qt_widget_demo_tab<Segment_tab_traits> *w_demo_p = 
       static_cast<Qt_widget_demo_tab<Segment_tab_traits> *> 
       (myBar->currentPage());
     outFile << w_demo_p->m_curves_arr;
     break;
    }
   case POLYLINE_TRAITS:
    {
     Qt_widget_demo_tab<Polyline_tab_traits> *w_demo_p = 
       static_cast<Qt_widget_demo_tab<Polyline_tab_traits> *> 
       (myBar->currentPage());
     outFile << w_demo_p->m_curves_arr;
     break;
    }
   case CONIC_TRAITS:
    {
     Qt_widget_demo_tab<Conic_tab_traits> *w_demo_p = 
       static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> 
       (myBar->currentPage());
     outFile << w_demo_p->m_curves_arr;
     break;
    }
  }  
  
  outFile.close();
  
  setCaption( QString( "Planar Map -- %1" ).arg( m_filename ) );
  statusBar()->message( QString( "Saved \'%1\'" ).arg( m_filename ), 2000 );
  //m_changed = FALSE;
}

/*! save planar map to post script */
void MyWindow::fileSave_ps()
{
#if 0
  
  Qt_widget_base_tab    *w_demo_p1 = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  
  switch ( w_demo_p1->traits_type ) {
   case SEGMENT_TRAITS:
    {
     Qt_widget_demo_tab<Segment_tab_traits> *w_demo_p = 
       static_cast<Qt_widget_demo_tab<Segment_tab_traits> *> 
       (myBar->currentPage());
     // Print to Postscript file:
     CGAL::Postscript_file_stream  LPF(m_width, m_height ,"pm.ps");
     LPF.init(-3,3,-3);
     LPF.set_line_width(1);
     LPF << w_demo_p->m_curves_arr;
     break;
    }
   case POLYLINE_TRAITS:
    {
     //Qt_widget_demo_tab<Polyline_tab_traits> *w_demo_p = 
     //static_cast<Qt_widget_demo_tab<Polyline_tab_traits> *> 
     //  (myBar->currentPage());
     //outFile << w_demo_p->m_curves_arr;
     break;
    }
   case CONIC_TRAITS:
    {
     //Qt_widget_demo_tab<Conic_tab_traits> *w_demo_p = 
     //static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> 
      // (myBar->currentPage());
     //outFile << w_demo_p->m_curves_arr;
     break;
    }
  }  
 #endif 
}

/*! print planar map */
void MyWindow::print()
{
  Qt_widget_base_tab    *w_demo_p1 = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  w_demo_p1->print_to_ps();
}

/*! open the overlay dialog form and read its output */
void MyWindow::overlay_pm()
{
  OverlayForm *form = new OverlayForm( myBar , this ,tab_number - 1);
  if ( form->exec() ) 
  {    
    unsigned int i = 2;
    if (form->listBox2->count() < i)
    {
      QMessageBox::information( this, my_title_string,
          "Please!!! you need more than one planar map to make an overlay...");
      return;
    }
    i = 12;
    if (form->listBox2->count() > i)
    {
      QMessageBox::information( this, my_title_string,
                              "Max number of Planar Maps to overlay is 12!!!");
      return;
    }
    std::list<int> indexes;
    TraitsType t;
    Qt_widget_base_tab *w_demo_p;
    int index,real_index;
    
    for (unsigned int i = 0; i < form->listBox2->count(); i++)
    {
      form->listBox2->setCurrentItem(i);
      char s[100];
      strcpy(s, form->listBox2->currentText()); 
      char * pch;
      pch = strtok(s," ");
      pch = strtok(NULL, " ");
      index = atoi(pch);
      real_index = realIndex(index);
      indexes.push_back(real_index);
    }
    
    w_demo_p = static_cast<Qt_widget_base_tab *> (myBar->page( real_index ));
    t = w_demo_p->traits_type;

    QCursor old = w_demo_p->cursor();
    w_demo_p->setCursor(Qt::WaitCursor);
    make_overlay( indexes , t);
    w_demo_p->setCursor(old);
    
  }
  delete form;
}

/*! make the overlay
 * \param indexes - list of the planar maps indexes in the overlay
 * \param t - the traits type of the overlay
 */
void MyWindow::make_overlay( std::list<int> indexes , TraitsType t)
{
  switch ( t ) 
  {
   case SEGMENT_TRAITS:
    {
     add_segment_tab();
     Qt_widget_demo_tab<Segment_tab_traits> *w_demo_p_new = 
       static_cast<Qt_widget_demo_tab<Segment_tab_traits> *> 
       (myBar->currentPage());
     Qt_widget_demo_tab<Segment_tab_traits> *w_demo_p;
     
     std::list<Pm_seg_2> seg_list;
     Pm_seg_const_iter itp;
     *w_demo_p_new << CGAL::LineWidth(2);
     int current;
     
     while (! indexes.empty())
     {
       current = indexes.front();
       indexes.pop_front();
       w_demo_p_new->setColor(colors[current-1]);
       w_demo_p = 
         static_cast<Qt_widget_demo_tab<Segment_tab_traits> *> 
         (myBar->page( current ));
       
       Seg_arr::Edge_iterator hei;
       for (hei = w_demo_p->m_curves_arr.edges_begin(); 
            hei != w_demo_p->m_curves_arr.edges_end(); ++hei) 
       {
         Pm_xseg_2 xcurve = hei->curve();
         Pm_seg_2 curve(xcurve, xcurve.get_data());
         seg_list.push_back(curve);
       }
       
     }
     
     w_demo_p_new->m_curves_arr.insert(seg_list.begin(),seg_list.end());
     
     if (!colors_flag)
     {
       // update new planner map indexes
       Seg_arr::Halfedge_iterator hei;
       for (hei = w_demo_p_new->m_curves_arr.halfedges_begin(); 
            hei != w_demo_p_new->m_curves_arr.halfedges_end(); ++hei) 
       {
         Curve_data cd = hei->curve().get_data();
         cd.m_type = Curve_data::INTERNAL;
         cd.m_index = w_demo_p_new->index;
         hei->curve().set_data( cd );
       }
       
     }
     
     break;
    }
   case POLYLINE_TRAITS:
    {
     add_polyline_tab();
     Qt_widget_demo_tab<Polyline_tab_traits> *w_demo_p_new = 
       static_cast<Qt_widget_demo_tab<Polyline_tab_traits> *> 
       (myBar->currentPage());
     Qt_widget_demo_tab<Polyline_tab_traits> *w_demo_p;
     
     std::list<Pm_pol_2> pol_list;
     Pm_pol_const_iter itp;
     int current;
     
     while (! indexes.empty())
     {
       current = indexes.front();
       indexes.pop_front();
       w_demo_p_new->setColor(colors[current-1]);
       w_demo_p = 
         static_cast<Qt_widget_demo_tab<Polyline_tab_traits> *> 
         (myBar->page( current ));
       
       Pol_arr::Halfedge_iterator hei;
       for (hei = w_demo_p->m_curves_arr.halfedges_begin(); 
            hei != w_demo_p->m_curves_arr.halfedges_end(); ++hei) 
       {
         Pm_xpol_2 xcurve = hei->curve();
         Pm_pol_2 curve(xcurve, xcurve.get_data());
         pol_list.push_back(curve);
       }
       
     }
     
     w_demo_p_new->m_curves_arr.insert(
                                       pol_list.begin(),pol_list.end());
     
     
     if (!colors_flag)
     {
       // update new planner map indexes
       Pol_arr::Halfedge_iterator hei;
       for (hei = w_demo_p_new->m_curves_arr.halfedges_begin(); 
            hei != w_demo_p_new->m_curves_arr.halfedges_end(); ++hei) 
       {
         Curve_pol_data cd = hei->curve().get_data();
         cd.m_type = Curve_pol_data::INTERNAL;
         cd.m_index = w_demo_p_new->index;
         hei->curve().set_data( cd );
       }
     }
     

     break;
    }
   case CONIC_TRAITS:
    {
     add_conic_tab();
     Qt_widget_demo_tab<Conic_tab_traits> *w_demo_p_new = 
       static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> 
       (myBar->currentPage());
     Qt_widget_demo_tab<Conic_tab_traits> *w_demo_p;
     
     Pm_xconic_const_iter itp;
     *w_demo_p_new << CGAL::LineWidth(3);
     int current;
     
     while (! indexes.empty())
     {
       current = indexes.front();
       indexes.pop_front();
       w_demo_p_new->setColor(colors[current-1]);
       w_demo_p = static_cast<Qt_widget_demo_tab<Conic_tab_traits> *>
         (myBar->page( current ));
       
       Conic_arr::Halfedge_iterator hei;
       for (hei = w_demo_p->m_curves_arr.halfedges_begin(); 
            hei != w_demo_p->m_curves_arr.halfedges_end(); ++hei) 
       {
         Pm_xconic_2 xcurve = hei->curve();
         Pm_conic_2 curve(xcurve, xcurve.get_data());
         //seg_list.push_back(curve);
         w_demo_p_new->m_curves_arr.insert(curve);
       }
     }
     
     if (!colors_flag)
     {
       // update new planner map indexes
       Conic_arr::Halfedge_iterator hei;
       for (hei = w_demo_p_new->m_curves_arr.halfedges_begin(); 
            hei != w_demo_p_new->m_curves_arr.halfedges_end(); ++hei) 
       {
         Curve_conic_data cd = hei->curve().get_data();
         cd.m_type = Curve_conic_data::INTERNAL;
         cd.m_index = w_demo_p_new->index;
         hei->curve().set_data( cd );
       }
     }
     
     break;
    }
    
  }
}

/*! real index - finds the tab bar index of a tab
 * \param index - the tab index (the same one it was 
 *                initialized with)
 *\ return the tab bar index of this tab
 */
int MyWindow::realIndex(int index)
{
  Qt_widget_base_tab * w_demo_p;
  for (int i = 0; i < tab_number-1; i++)
  {
    if ( myBar->isTabEnabled( myBar->page(i) ) )
    {
      // We peform downcasting from QWigdet* to 
      // Qt_widget_base_tab*, as we know that only
      // Qt_widget_base_tab objects are stored in the tab pages.
      w_demo_p = static_cast<Qt_widget_base_tab *> (myBar->page(i));
      if (w_demo_p->index == index)
        return i;
    }
  }
  return -1;
}

/*! show grid without been in a grid snap mode */
void MyWindow::showGrid()
{
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*, 
  // as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab     *w_demo_p = 
    dynamic_cast<Qt_widget_base_tab  *> (myBar->currentPage());
  w_demo_p->grid = true;
  something_changed();
}

/*! hide the grid (can be applied only when we 
 *  are not in a grid snap mode.
 */
void MyWindow::hideGrid()
{
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*, 
  // as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab     *w_demo_p = 
    dynamic_cast<Qt_widget_base_tab  *> (myBar->currentPage());
  w_demo_p->grid = false;
  something_changed();
}

/*! a dialog form to set the conic type for conics insertion. */
void MyWindow::conicType()
{
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*, 
  // as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab     *w_demo_p = 
    dynamic_cast<Qt_widget_base_tab  *> (myBar->currentPage());
  OptionsForm *form = new OptionsForm();
  if ( form->exec() ) 
  {
    QString type = form->arrComboBox1->currentText();
    //std::cout << type << std::endl;
    if (strcmp(type,"Circle") == 0)
      w_demo_p->conic_type = CIRCLE;
    else if (strcmp(type,"Segment") == 0)
      w_demo_p->conic_type = SEGMENT;
	else if (strcmp(type,"Ellipse") == 0)
      w_demo_p->conic_type = ELLIPSE;
	else if (strcmp(type,"Parabola") == 0)
      w_demo_p->conic_type = PARABOLA;
	else if (strcmp(type,"Hyperbola") == 0)
      w_demo_p->conic_type = HYPERBOLA;
  }
}

void MyWindow::backGroundColor()
{
  QColor c = QColorDialog::getColor();
  Qt_widget_base_tab     *w_demo_p = 
    dynamic_cast<Qt_widget_base_tab  *> (myBar->currentPage());
  
  w_demo_p->setBackgroundColor( c ); 
  
  something_changed();
}

/*! a dialog form to set the rayShooting Diraction. */
void MyWindow::rayShootingDirection()
{
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*, 
  // as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab     *w_demo_p = 
    dynamic_cast<Qt_widget_base_tab  *> (myBar->currentPage());
  RayShootingOptionsForm *form = new RayShootingOptionsForm();
  if ( form->exec() ) 
  {
    QString type = form->arrComboBox1->currentText();
    if (strcmp(type,"Up") == 0)
	{
      w_demo_p->ray_shooting_direction = true;
      rayShootingMode->setIconSet(QPixmap((const char**)ray_shooting_xpm ));
	  if (w_demo_p->mode == RAY_SHOOTING)
	    w_demo_p->setCursor(
	    QCursor(QPixmap((const char**)small_ray_shooting_xpm)));
	}
    else if (strcmp(type,"Down") == 0)
	{
      w_demo_p->ray_shooting_direction = false;
      rayShootingMode->setIconSet( 
		  QPixmap((const char**)ray_shooting_down_xpm ));
	  if (w_demo_p->mode == RAY_SHOOTING)
	    w_demo_p->setCursor(
	    QCursor(QPixmap((const char**)small_ray_shooting_down_xpm)));
	}
  }
}

/*! main */
int main(int argc, char **argv)
{
  QApplication app( argc, argv );
  MyWindow widget(700,700); // physical window size
  app.setMainWidget(&widget);
  widget.setCaption(my_title_string);
  widget.setMouseTracking(TRUE);
  widget.show();
  return app.exec();  
}

#endif // CGAL_USE_QT

//std::cout << xseg << std::endl;
//std::fflush(stdout);
