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

#include <CGAL/IO/pixmaps/zoom_out.xpm>
#include <CGAL/IO/pixmaps/zoom_in.xpm>
#include <CGAL/IO/pixmaps/hand.xpm>
#include <CGAL/IO/pixmaps/movepoint.xpm>
#include <CGAL/IO/pixmaps/point.xpm>
#include <CGAL/IO/pixmaps/line.xpm>
#include <CGAL/IO/pixmaps/arrow.xpm>
#include <CGAL/IO/pixmaps/alpha_shape.xpm>
#include "icons/polyline.xpm"
#include "icons/insert.xpm"
#include "icons/delete.xpm"
#include "icons/grid.xpm"
#include "icons/conic.xpm"
#include "icons/none.xpm"
#include "icons/po.xpm"
#include "icons/ray_shooting2.xpm"
#include "icons/draw.xpm"

const QString my_title_string("Arrangement Demo with CGAL Qt_widget");

//////////////////////////////////////////////////////////////////////////////

MyWindow::MyWindow(int w, int h){

  myBar = new QTabWidget(this);
    setCentralWidget(myBar);
  m_width = w;
  m_height = h;
  tab_number = 0;
  number_of_tabs = 0;
  overlay_flag = false; // flag for add_conic_tab
  testlayer = new Qt_layer( myBar );

  list_of_colors.push_back(Qt::green);
  list_of_colors.push_back(Qt::blue);
  list_of_colors.push_back(Qt::black);
  list_of_colors.push_back(Qt::yellow);
  list_of_colors.push_back(Qt::cyan);
  list_of_colors.push_back(Qt::magenta);
  list_of_colors.push_back(Qt::gray);
  list_of_colors.push_back(Qt::darkGreen);
  list_of_colors.push_back(Qt::darkBlue);
  list_of_colors.push_back(Qt::darkYellow);
  list_of_colors.push_back(Qt::darkCyan);
  list_of_colors.push_back(Qt::darkMagenta);

  // status bar labels:
  insert_label = new QLabel("Insert" , this);
  delete_label = new QLabel("Delete" , this);
  drag_label = new QLabel("Drag" , this);
  ray_shooting_label = new QLabel("Ray Shooting" , this);
  point_location_label = new QLabel("Point Location" , this);
  current_label = insert_label;
  m_scailing_factor = 2;

  // Traits Group
  QActionGroup *traitsGroup = new QActionGroup( this ); // Connected later
  traitsGroup->setExclusive( TRUE );

  setSegmentTraits = new QAction(
      "Segment Traits", QPixmap( (const char**)line_xpm ),
      "&Segment Traits", 0 ,traitsGroup, "Segment Traits" );
  setSegmentTraits->setToggleAction( TRUE );

  setPolylineTraits = new QAction(
      "Polyline Traits", QPixmap( (const char**)polyline_xpm ),
      "&Polyline Traits", 0 , traitsGroup, "Polyline Traits" );
  setPolylineTraits->setToggleAction( TRUE );

  setConicTraits = new QAction(
      "Conic Traits", QPixmap( (const char**)conic_xpm ),
      "&Conic Traits", 0 , traitsGroup, "Conic Traits" );
  setConicTraits->setToggleAction( TRUE );

  // Snap Mode Group
  QActionGroup *snapModeGroup = new QActionGroup( this ); // Connected later
  snapModeGroup->setExclusive( TRUE );

  setNoneSnapMode = new QAction(
      "None Snap Mode", QPixmap( (const char**)none_xpm ),
      "&None Snap Mode", 0 , snapModeGroup, "None Snap Mode" );
  setNoneSnapMode->setToggleAction( TRUE );

  setGridSnapMode = new QAction(
      "Grid Snap Mode", QPixmap( (const char**)grid_xpm ),
      "&Grid Snap Mode", 0 , snapModeGroup, "Grid Snap Mode" );
  setGridSnapMode->setToggleAction( TRUE );

  setPointSnapMode = new QAction(
      "Point Snap Mode", QPixmap( (const char**)po_xpm ),
      "&Point Snap Mode", 0 , snapModeGroup, "Point Snap Mode" );
  setPointSnapMode->setToggleAction( TRUE );

  // insert - delete - point_location Mode Group
  QActionGroup *modeGroup = new QActionGroup( this ); // Connected later
  modeGroup->setExclusive( TRUE );

  insertMode = new QAction(
      "Insert", QPixmap( (const char**)draw_xpm ),
      "&Insert", 0 , modeGroup, "Insert" );
  insertMode->setToggleAction( TRUE );

  deleteMode = new QAction(
      "Delete", QPixmap( (const char**)delete_xpm ),
      "&Delete", 0 , modeGroup, "Delete" );
  deleteMode->setToggleAction( TRUE );

  pointLocationMode = new QAction(
      "PointLocation", QPixmap( (const char**)point_xpm ),
      "&Point Location", 0 , modeGroup, "Point Location" );
  pointLocationMode->setToggleAction( TRUE );

  rayShootingMode = new QAction(
      "RayShooting", QPixmap( (const char**)ray_shooting_xpm ),
      "&Ray Shooting", 0 , modeGroup, "Ray Shooting" );
  rayShootingMode->setToggleAction( TRUE );
 
  dragMode = new QAction(
    "Drag", QPixmap( (const char**)hand_xpm ),
    "&Drag", 0 , modeGroup, "Drag" );
  dragMode->setToggleAction( TRUE );
  
  //create a timer for checking if somthing changed
  QTimer *timer = new QTimer( this );
  connect( timer, SIGNAL(timeout()),
           this, SLOT(timer_done()) );
  timer->start( 200, FALSE );

    // file menu
  QPopupMenu * file = new QPopupMenu( this );
  menuBar()->insertItem( "&File", file );
  file->insertItem("&Open Conic File", this, SLOT(add_conic_tab()));
  file->insertItem("&Open Polyline File", this, SLOT(fileOpenPolyline()));
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
  menuBar()->insertSeparator();

  // snap mode menu
  QPopupMenu * snap_mode = new QPopupMenu( this );
  menuBar()->insertItem( "&Snap mode", snap_mode );
  setNoneSnapMode->addTo(snap_mode); 
  setGridSnapMode->addTo(snap_mode);
  setPointSnapMode->addTo(snap_mode);
  menuBar()->insertSeparator();

  // options menu
  QPopupMenu * options = new QPopupMenu( this );
  menuBar()->insertItem( "&Options", options );
  setSegmentTraits->addTo(options); 
  setPolylineTraits->addTo(options); 
  setConicTraits->addTo(options); 
  options->insertSeparator();
  options->insertItem("Overlay", this, SLOT(overlay_pm()));
  menuBar()->insertSeparator();
  options->insertItem("Properties", this, SLOT(properties()));

  // help menu
  QPopupMenu * help = new QPopupMenu( this );
  menuBar()->insertItem( "&Help", help );
  help->insertItem("How To", this, SLOT(howto()), Key_F1);
  help->insertSeparator();
  help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
  help->insertItem("About &Qt", this, SLOT(aboutQt()) );
   
  // options toolbar
  QToolBar *optionsTools = new QToolBar( this, "options operations" );
  optionsTools->setLabel( "Options Operations" );
  insertMode->addTo( optionsTools );
  deleteMode->addTo( optionsTools );
  dragMode->addTo( optionsTools );
  pointLocationMode->addTo( optionsTools );
  rayShootingMode->addTo( optionsTools );
  optionsTools->addSeparator();
  setNoneSnapMode->addTo( optionsTools );
  setGridSnapMode->addTo( optionsTools );
  setPointSnapMode->addTo( optionsTools );
  optionsTools->addSeparator();
  setSegmentTraits->addTo( optionsTools );
  setPolylineTraits->addTo( optionsTools );
  setConicTraits->addTo( optionsTools );
  optionsTools->addSeparator();
  optionsTools->addSeparator();
  
  // zoom button
  QToolButton* zoominBt =
    new QToolButton(QPixmap(zoomin_xpm ),
                    "Zoom in", 
                    0, 
                    this, 
                    SLOT(zoomin()), 
                    optionsTools, 
                    "Zoom in");
  zoominBt->setTextLabel("Scaling factor ");

  QToolButton* zoomoutBt = 
    new QToolButton(QPixmap(zoomout_xpm ),
                    "Zoom out", 
                    0, 
                    this, 
                    SLOT(zoomout()), 
                    optionsTools, 
                    "Zoom out");
  zoomoutBt->setTextLabel("Scaling factor ");

  // connect mode group
  connect( modeGroup, SIGNAL( selected(QAction*) ), 
           this, SLOT( updateMode(QAction*) ) );

  // connect Traits Group
  connect( traitsGroup, SIGNAL( selected(QAction*) ),
           this, SLOT( updateTraitsType(QAction*) ) );

  // connect Snap Mode Group
  connect( snapModeGroup, SIGNAL( selected(QAction*) ),
           this, SLOT( updateSnapMode(QAction*) ) );
  
  //state flag 
  old_state = 0;
   
  add_segment_tab();

  // connect the change of current tab
  connect( myBar, SIGNAL( currentChanged(QWidget * )  ),
           this, SLOT( update() ) );
  
  statusBar();
 }

MyWindow::~MyWindow()
{}

void MyWindow::something_changed()
{
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*, as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab   *w_demo_p = dynamic_cast<Qt_widget_base_tab  *> (myBar->currentPage());

  w_demo_p->current_state++;
}
 
void MyWindow::get_new_object(CGAL::Object obj)
  {
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*, as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab  *w_demo_p = static_cast<Qt_widget_base_tab *> (myBar->currentPage());

  // point location
    Coord_point p;
    if(CGAL::assign(p,obj)) {

      w_demo_p->pl_point = p;
      //something_changed();
    }
  something_changed();
  }

  void MyWindow::about()
  {
    QMessageBox::about( this, my_title_string,
    "This is a demo for the Arrangement package\n"
      "Copyright CGAL @2003");
  }

  void MyWindow::aboutQt()
  {
    QMessageBox::aboutQt( this, my_title_string );
  }

  void MyWindow::howto(){
    QString home;
    home = "help/index.html";
    CGAL::Qt_help_window * help =
      new CGAL::Qt_help_window(home, ".", 0, "help viewer");
    help->resize(400, 400);
    help->setCaption("Demo HowTo");
    help->show();
  }

  void MyWindow::add_segment_tab()
  {
  Qt_widget_demo_tab<Segment_tab_traits> *widget = new Qt_widget_demo_tab<Segment_tab_traits>(SEGMENT_TRAITS , this, tab_number);
  connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),this, SLOT(get_new_object(CGAL::Object)));
  widget->attach(testlayer);
  widget->setCursor(QCursor( QPixmap( (const char**)small_draw_xpm)));
  tab_number++;
  number_of_tabs++;
  // add the new widget to myBar
  myBar->insertTab( widget, QString("Arr " + QString::number( widget->index ) ) , widget->index );
  myBar->setCurrentPage(myBar->indexOf(widget));
  resize(m_width,m_height);
  something_changed();
  }

  void MyWindow::add_polyline_tab()
  {
  Qt_widget_demo_tab<Polyline_tab_traits> *widget = new Qt_widget_demo_tab<Polyline_tab_traits>(POLYLINE_TRAITS , this, tab_number);
  connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),this, SLOT(get_new_object(CGAL::Object)));
  widget->attach(testlayer);
  widget->setCursor(QCursor( QPixmap( (const char**)small_draw_xpm)));
  tab_number++;
  number_of_tabs++;
  // add the new widget to myBar
  myBar->insertTab( widget, QString("Arr " + QString::number( widget->index ) ) , widget->index );
  myBar->setCurrentPage(myBar->indexOf(widget));
  resize(m_width,m_height);
  something_changed();
  }
  
  void MyWindow::add_conic_tab()
  {
  Qt_widget_demo_tab<Conic_tab_traits> *widget = new Qt_widget_demo_tab<Conic_tab_traits>(CONIC_TRAITS , this , tab_number);
  connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),this, SLOT(get_new_object(CGAL::Object)));
  widget->attach(testlayer);
  widget->setCursor(QCursor( QPixmap( (const char**)small_draw_xpm)));
  tab_number++;
  number_of_tabs++;
  // add the new widget to myBar
  myBar->insertTab( widget, QString("Arr " + QString::number( widget->index ) ) , widget->index );
  myBar->setCurrentPage(myBar->indexOf(widget));
  if (!overlay_flag)
    fileOpen();
  resize(m_width,m_height);
  something_changed();
  }

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


  void MyWindow::timer_done()
  {
  // We peform downcasting from QWigdet* to Qt_widget_demo_tab*, as we know that only
  // Qt_widget_demo_tab objects are stored in the tab pages.
  Qt_widget_base_tab  *w_demo_p = static_cast<Qt_widget_base_tab *> (myBar->currentPage());

    if(old_state!=w_demo_p->current_state){
      w_demo_p->redraw();
      old_state = w_demo_p->current_state;
    }
  }  

void MyWindow::updateTraitsType( QAction *action )
{
  //if (tab_number == 0) return;
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*, as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab *old_widget = static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  Qt_widget_base_tab *widget;

  if (action == setSegmentTraits)
  {
    if (old_widget->traits_type == SEGMENT_TRAITS) return;
    widget = new Qt_widget_demo_tab<Segment_tab_traits>(SEGMENT_TRAITS , this);
  }
  else if (action == setPolylineTraits)
  {
    if (old_widget->traits_type == POLYLINE_TRAITS) return;
    widget = new Qt_widget_demo_tab<Polyline_tab_traits>(POLYLINE_TRAITS , this);
  }
  else if (action == setConicTraits)
  {
    if (old_widget->traits_type == CONIC_TRAITS) return;
    widget = new Qt_widget_demo_tab<Conic_tab_traits>(CONIC_TRAITS , this);
  }

  int old_index = old_widget->index;
  int index = myBar->currentPageIndex();
  QString label = myBar->label(index);
  myBar->removePage(myBar->currentPage());

  //initialize the new tab widget
  *widget << CGAL::LineWidth(2) << CGAL::BackgroundColor (CGAL::WHITE);
  widget->set_window(-10, 10, -10, 10);
  widget->setMouseTracking(TRUE);
  connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),this, SLOT(get_new_object(CGAL::Object)));
  widget->attach(testlayer);
  widget->m_line_width = 2;
  widget->index = old_index;
  widget->setCursor(QCursor( QPixmap( (const char**)small_draw_xpm)));

  // add the new widget to myBar
  myBar->insertTab( widget, label , index );

  myBar->setCurrentPage(index);
    
  resize(m_width,m_height);
  
  if (action == setConicTraits)
    fileOpen();

  something_changed();

}

void MyWindow::setTraits( TraitsType t )
{
  switch ( t ) {
  case SEGMENT_TRAITS:
    setSegmentTraits->setOn( TRUE );
    insertMode->setEnabled( TRUE );
    deleteMode->setEnabled( TRUE );
    break;
  case POLYLINE_TRAITS:
    setPolylineTraits->setOn( TRUE );
    insertMode->setEnabled( TRUE );
    deleteMode->setEnabled( TRUE );
    break;
  case CONIC_TRAITS:
    setConicTraits->setOn( TRUE );
    insertMode->setEnabled( FALSE );
    deleteMode->setEnabled( FALSE );
    break;
  }
}

void MyWindow::updateSnapMode( QAction *action )
{
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*, as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab  *w_demo_p = static_cast<Qt_widget_base_tab *> (myBar->currentPage());

  if ( action == setNoneSnapMode ) {
  w_demo_p->snap_mode = NONE;
  }
  else if ( action == setGridSnapMode ) {
  w_demo_p->snap_mode = GRID;
  }
  else if ( action == setPointSnapMode ) {
  w_demo_p->snap_mode = POINT;
  }

  something_changed();
}

void MyWindow::setSnapMode( SnapMode m )
{
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*, as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab  *w_demo_p = static_cast<Qt_widget_base_tab *> (myBar->currentPage());

  w_demo_p->snap_mode = m;
  switch ( m ) {
  case NONE:
    setNoneSnapMode->setOn( TRUE );
    something_changed();
    break;
  case GRID:
    setGridSnapMode->setOn( TRUE );
    something_changed();
    break;
  case POINT:
    setPointSnapMode->setOn( TRUE );
    something_changed();
    break;

  }
}

void MyWindow::updateMode( QAction *action )
{
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*, as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab  *w_demo_p = static_cast<Qt_widget_base_tab *> (myBar->currentPage());

  if ( action == insertMode ) 
  {
        w_demo_p->mode = INSERT;
    //w_demo_p->setCursor(QCursor(Qt::ArrowCursor));
    w_demo_p->setCursor(QCursor( QPixmap( (const char**)small_draw_xpm)));
    statusBar()->removeWidget(current_label);
    statusBar()->addWidget(insert_label);
    current_label = insert_label;
    something_changed();
  }
  else if ( action == deleteMode ) 
  {
    w_demo_p->mode = DELETE;
    w_demo_p->setCursor(QCursor( QPixmap( (const char**)delete_xpm)));
    statusBar()->removeWidget(current_label);
    statusBar()->addWidget(delete_label);
    current_label = delete_label;
    something_changed();
  }
  else if ( action == pointLocationMode ) 
  {
    w_demo_p->mode = POINT_LOCATION;
    w_demo_p->setCursor(QCursor( QPixmap( (const char**)point_xpm)));
    statusBar()->removeWidget(current_label);
    statusBar()->addWidget(point_location_label);
    current_label = point_location_label;
  }
  else if ( action == rayShootingMode ) 
  {
    w_demo_p->mode = RAY_SHOOTING;
    w_demo_p->setCursor(QCursor( QPixmap( (const char**)small_ray_shooting_xpm)));
    statusBar()->removeWidget(current_label);
    statusBar()->addWidget(ray_shooting_label);
    current_label = ray_shooting_label;
  }
  else if ( action == dragMode ) 
  {
    w_demo_p->mode = DRAG;
    w_demo_p->setCursor(QCursor( QPixmap( (const char**)hand_xpm)));
    statusBar()->removeWidget(current_label);
    statusBar()->addWidget(drag_label);
    current_label = drag_label;
    something_changed();
  }
}

void MyWindow::setMode( Mode m )
{
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*, as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab  *w_demo_p = static_cast<Qt_widget_base_tab *> (myBar->currentPage());

  w_demo_p->mode = m;
  switch ( m ) {
  case INSERT:
    insertMode->setOn( TRUE );
    break;
  case DELETE:
    deleteMode->setOn( TRUE );
    break;
  case POINT_LOCATION:
    pointLocationMode->setOn( TRUE );
    break;
  case RAY_SHOOTING:
    rayShootingMode->setOn( TRUE );
    break;
  case DRAG:
    dragMode->setOn( TRUE );
    break;
  }
}

void MyWindow::update()
{
  // We peform downcasting from QWigdet* to Qt_widget_base_tab*, as we know that only
  // Qt_widget_base_tab objects are stored in the tab pages.
  Qt_widget_base_tab  *w_demo_p = static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  setMode( w_demo_p->mode );
  setSnapMode( w_demo_p->snap_mode );
  setTraits( w_demo_p->traits_type );
      
}

void MyWindow::zoomin()
{
    Qt_widget_base_tab  *w_demo_p = static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  w_demo_p->zoom(m_scailing_factor);
}

void MyWindow::zoomout()
{
    Qt_widget_base_tab  *w_demo_p = static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  w_demo_p->zoom(1/m_scailing_factor);
}

void MyWindow::properties()
{
    PropertiesForm *optionsForm = new PropertiesForm( myBar , this ,number_of_tabs );
    
    if ( optionsForm->exec() ) 
  {  
    Qt_widget_base_tab  *w_demo_p1 = static_cast<Qt_widget_base_tab *> (myBar->currentPage());
    m_width = optionsForm->box1->value();
    m_height = optionsForm->box2->value();
    w_demo_p1->m_line_width = optionsForm->box3->value();
    double new_factor = static_cast<double> (optionsForm->box4->value());
    m_scailing_factor = (new_factor/10);
    std::cout << optionsForm->box4->value() << " " << m_scailing_factor << std::endl;
    resize(m_width,m_height);
    w_demo_p1->redraw();
    something_changed();
    }
    delete optionsForm;
}

void MyWindow::fileOpenPolyline()
{
  add_polyline_tab();
  fileOpen();
}

void MyWindow::fileOpen()
{
    
    QString filename = QFileDialog::getOpenFileName(
          QString::null, 0, this,
          "file open", "Demo -- File Open" );
    if ( !filename.isEmpty() )
  load( filename );
    else
  statusBar()->message( "File Open abandoned", 2000 );
}

void MyWindow::load( const QString& filename )
{
    std::ifstream inputFile(filename);
        // Creates an ofstream object named inputFile
    if (! inputFile.is_open()) // Always test file open
    {
      std::cout << "Error opening input file" << std::endl;
      return;
    }

    
        Qt_widget_base_tab  *w_demo = static_cast<Qt_widget_base_tab *> (myBar->currentPage());

    if (w_demo->traits_type == CONIC_TRAITS)
    {
      Qt_widget_demo_tab<Conic_tab_traits>  *w_demo_p = static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());
      char dummy[256];
      Pm_base_conic_2* cv;
      int count;
      //std::list<Pm_conic_2> curve_list;
    
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

        //curve_list.push_back(Pm_conic_2( *cv , cd));
        w_demo_p->m_curves_arr.insert(Pm_conic_2( *cv , cd));
        
        CGAL::Bbox_2 curve_bbox = cv->bbox();
        if (i == 0)
          w_demo->bbox = curve_bbox;
        else
          w_demo->bbox = w_demo->bbox + curve_bbox;
      }
      
      //w_demo_p->m_curves_arr.insert(curve_list.begin() , curve_list.end());
      // insert xcurve into xcurve list

      Conic_arr::Edge_iterator ei;

      for (ei = w_demo_p->m_curves_arr.edges_begin(); ei != w_demo_p->m_curves_arr.edges_end(); ++ei) 
      {
        Pm_xconic_2 *xseg = new Pm_xconic_2(ei->curve());
        w_demo_p->m_curves_list.push_back(xseg);
      }

    }

    else if (w_demo->traits_type == POLYLINE_TRAITS)
    {
      Qt_widget_demo_tab<Polyline_tab_traits>  *w_demo_p = static_cast<Qt_widget_demo_tab<Polyline_tab_traits> *> (myBar->currentPage());
      w_demo_p->m_curves_list.clear();
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

        Pm_base_pol_2 *base_polyline = new Pm_base_pol_2(points.begin(), points.end());
        
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
        w_demo_p->m_curves_list.push_back(new Pm_pol_2( *base_polyline , cd));
    
      }

    }

    w_demo->set_window(w_demo->bbox.xmin() , w_demo->bbox.xmax() , w_demo->bbox.ymin() , w_demo->bbox.ymax());

    inputFile.close();
    something_changed();
}


void MyWindow::ReadCurve(std::ifstream & is, Pm_base_conic_2 & cv)
{
      // Read a line from the input file.
      char one_line[128];
      
      skip_comments (is, one_line);
      std::string stringvalues(one_line);
      std::istringstream str_line (stringvalues, std::istringstream::in);
      
      // Get the arc type.
      char     type;
      bool     is_circle = false;              // Is this a circle.
      Pm_conic_circle_2 circle;
      CONIC_NT       r, s, t, u, v, w;               // The conic coefficients.
      
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
              circle = Pm_conic_circle_2 (Pm_conic_point_2 (x0, y0), a*b, CGAL::CLOCKWISE);
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

void MyWindow::fileSaveAs()
{
    QString filename = QFileDialog::getSaveFileName(
          QString::null, "Planar Map (*.pm)", this,
          "file save as", "Planar Map -- File Save As" );
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

    Qt_widget_base_tab  *w_demo_p1 = static_cast<Qt_widget_base_tab *> (myBar->currentPage());

  switch ( w_demo_p1->traits_type ) {
    case SEGMENT_TRAITS:
    {
      Qt_widget_demo_tab<Segment_tab_traits> *w_demo_p = static_cast<Qt_widget_demo_tab<Segment_tab_traits> *> (myBar->currentPage());
      outFile << w_demo_p->m_curves_arr;
      break;
    }
    case POLYLINE_TRAITS:
    {
      Qt_widget_demo_tab<Polyline_tab_traits> *w_demo_p = static_cast<Qt_widget_demo_tab<Polyline_tab_traits> *> (myBar->currentPage());
      outFile << w_demo_p->m_curves_arr;
      break;
    }
    case CONIC_TRAITS:
    {
      Qt_widget_demo_tab<Conic_tab_traits> *w_demo_p = static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());
      outFile << w_demo_p->m_curves_arr;
      break;
    }
  }  

  outFile.close();

    setCaption( QString( "Planar Map -- %1" ).arg( m_filename ) );
    statusBar()->message( QString( "Saved \'%1\'" ).arg( m_filename ), 2000 );
    //m_changed = FALSE;
}

void MyWindow::fileSave_ps()
{
#if 0
  CGAL::Postscript_file_stream ps_stream(m_width, m_height ,"pm.ps");

    Qt_widget_base_tab  *w_demo_p1 = static_cast<Qt_widget_base_tab *> (myBar->currentPage());

  switch ( w_demo_p1->traits_type ) {
    case SEGMENT_TRAITS:
    {
      Qt_widget_demo_tab<Segment_tab_traits> *w_demo_p = static_cast<Qt_widget_demo_tab<Segment_tab_traits> *> (myBar->currentPage());
      
      //ps_stream.set_line_width(w_demo_p->m_line_width);
      CGAL::Pm_drawer<Seg_arr, CGAL::Postscript_file_stream> drawer(ps_stream);
      ps_stream << CGAL::BLUE;
      drawer.draw_halfedges(w_demo_p->m_curves_arr.halfedges_begin(), w_demo_p->m_curves_arr.halfedges_end());
      ps_stream << CGAL::RED;
      drawer.draw_vertices(w_demo_p->m_curves_arr.vertices_begin(), w_demo_p->m_curves_arr.vertices_end());
      break;
    }
    case POLYLINE_TRAITS:
    {
      //Qt_widget_demo_tab<Polyline_tab_traits> *w_demo_p = static_cast<Qt_widget_demo_tab<Polyline_tab_traits> *> (myBar->currentPage());
      //outFile << w_demo_p->m_curves_arr;
      break;
    }
    case CONIC_TRAITS:
    {
      //Qt_widget_demo_tab<Conic_tab_traits> *w_demo_p = static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());
      //outFile << w_demo_p->m_curves_arr;
      break;
    }
  }  

#endif 
}

void MyWindow::print()
{
    Qt_widget_base_tab  *w_demo_p1 = static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  w_demo_p1->print_to_ps();
}


void MyWindow::overlay_pm()
{
    OverlayForm *form = new OverlayForm( myBar , this ,tab_number );
    
    if ( form->exec() ) 
  {  
    unsigned int i = 2;
    if (form->listBox2->count() < i)
    {
      QMessageBox::information( this, my_title_string,"Please!!! you need more than one planar map to make an overlay...");
      return;
    }
    std::list<int> indexes;
    TraitsType t;
    Qt_widget_base_tab  *w_demo_p;
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

    w_demo_p = static_cast<Qt_widget_base_tab *> (myBar->page( index ));
    t = w_demo_p->traits_type;

    make_overlay( indexes , t);
  
    }
    delete form;
}

void MyWindow::make_overlay( std::list<int> indexes , TraitsType t)
{
  switch ( t ) 
  {
    case SEGMENT_TRAITS:
    {
      add_segment_tab();
      Qt_widget_demo_tab<Segment_tab_traits> *w_demo_p_new = static_cast<Qt_widget_demo_tab<Segment_tab_traits> *> (myBar->currentPage());
            Qt_widget_demo_tab<Segment_tab_traits> *w_demo_p;

      std::list<Pm_seg_2> seg_list;
      Pm_seg_const_iter itp;
      *w_demo_p_new << CGAL::LineWidth(3);
      int current;
      std::list<QColor>::const_iterator it = list_of_colors.begin();

      while (! indexes.empty())
      {
        current = indexes.front();
        indexes.pop_front();
        w_demo_p_new->setColor(*it);
        it++;
        w_demo_p = static_cast<Qt_widget_demo_tab<Segment_tab_traits> *> (myBar->page( current ));

        for (itp = w_demo_p->m_curves_list.begin(); itp != w_demo_p->m_curves_list.end(); ++itp)
        {
          Pm_seg_2 * seg = new Pm_seg_2( **itp );
          w_demo_p_new->m_curves_list.push_back(seg);
          seg_list.push_back(*seg);
          *w_demo_p_new << *(*itp);
        }
      }
      
      w_demo_p_new->m_curves_arr.insert(seg_list.begin(),seg_list.end());

      *w_demo_p_new << CGAL::RED;
      *w_demo_p_new << CGAL::DISC;
      // Go over all vertices and for each vertex print the ID numbers of the
      // base curves that go through it.
      Seg_arr::Vertex_iterator   vit;
      for (vit = w_demo_p_new->m_curves_arr.vertices_begin(); vit != w_demo_p_new->m_curves_arr.vertices_end(); vit++)
      {
        Seg_arr::Halfedge_around_vertex_circulator 
        eit, first = (*vit).incident_halfedges();

        eit = first;

        int ind1;
        int ind2 = (*eit).curve().get_data().m_index;
        do 
        {
          ind1 = (*eit).curve().get_data().m_index;

          // Keep track of IDs we haven't seen before.
          if (ind1 != ind2)
          {
            const Pm_seg_point_2& p = (*vit).point();
            *w_demo_p_new << p;
            break;
          }

          eit++;

        } while (eit != first);
      }

      // allow the user to see the overlay
      w_demo_p_new->current_state = old_state;
  
      // update new planner map index
      Pm_seg_iter iter;
      for (iter = w_demo_p_new->m_curves_list.begin(); iter != w_demo_p_new->m_curves_list.end(); ++iter)
      {
        Curve_data cd = (**iter).get_data();
        cd.m_type = Curve_data::INTERNAL;
        cd.m_index = w_demo_p_new->index;
        (**iter).set_data( cd );
      }

      break;
    }
    case POLYLINE_TRAITS:
    {
      add_polyline_tab();
      Qt_widget_demo_tab<Polyline_tab_traits> *w_demo_p_new = static_cast<Qt_widget_demo_tab<Polyline_tab_traits> *> (myBar->currentPage());
            Qt_widget_demo_tab<Polyline_tab_traits> *w_demo_p;

      std::list<Pm_pol_2> pol_list;
      Pm_pol_const_iter itp;
      *w_demo_p_new << CGAL::LineWidth(3);
      int current;
      std::list<QColor>::const_iterator it = list_of_colors.begin();

      while (! indexes.empty())
      {
        current = indexes.front();
        indexes.pop_front();
        w_demo_p_new->setColor(*it);
        it++;
        w_demo_p = static_cast<Qt_widget_demo_tab<Polyline_tab_traits> *> (myBar->page( current ));

        for (itp = w_demo_p->m_curves_list.begin(); itp != w_demo_p->m_curves_list.end(); ++itp)
        {
          Pm_pol_2 * seg = new Pm_pol_2( **itp );
          w_demo_p_new->m_curves_list.push_back(seg);
          pol_list.push_back(*seg);
          *w_demo_p_new << *(*itp);
        }
      }
      
      w_demo_p_new->m_curves_arr.insert(pol_list.begin(),pol_list.end());

      *w_demo_p_new << CGAL::RED;
      *w_demo_p_new << CGAL::DISC;
      // Go over all vertices and for each vertex print the ID numbers of the
      // base curves that go through it.
      Pol_arr::Vertex_iterator   vit;
      for (vit = w_demo_p_new->m_curves_arr.vertices_begin(); vit != w_demo_p_new->m_curves_arr.vertices_end(); vit++)
      {
        Pol_arr::Halfedge_around_vertex_circulator 
        eit, first = (*vit).incident_halfedges();

        eit = first;

        int ind1;
        int ind2 = (*eit).curve().get_data().m_index;
        do 
        {
          ind1 = (*eit).curve().get_data().m_index;

          // Keep track of IDs we haven't seen before.
          if (ind1 != ind2)
          {
            const Pm_pol_point_2& p = (*vit).point();
            *w_demo_p_new << p;
            break;
          }

          eit++;

        } while (eit != first);
      }

      // allow the user to see the overlay
      w_demo_p_new->current_state = old_state;
  
      // update new planner map index
      Pm_pol_iter iter;
      for (iter = w_demo_p_new->m_curves_list.begin(); iter != w_demo_p_new->m_curves_list.end(); ++iter)
      {
        Curve_pol_data cd = (**iter).get_data();
        cd.m_type = Curve_pol_data::INTERNAL;
        cd.m_index = w_demo_p_new->index;
        (**iter).set_data( cd );
      }

      break;
    }
    case CONIC_TRAITS:
    {
      overlay_flag = true;
      add_conic_tab();
      overlay_flag = false;
      Qt_widget_demo_tab<Conic_tab_traits> *w_demo_p_new = static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());
            Qt_widget_demo_tab<Conic_tab_traits> *w_demo_p;

      //std::list<Pm_xconic_2> conic_list;
      Pm_xconic_const_iter itp;
      *w_demo_p_new << CGAL::LineWidth(3);
      int current;
      std::list<QColor>::const_iterator it = list_of_colors.begin();

      while (! indexes.empty())
      {
        current = indexes.front();
        indexes.pop_front();
        w_demo_p_new->setColor(*it);
        it++;
        w_demo_p = static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->page( current ));

        for (itp = w_demo_p->m_curves_list.begin(); itp != w_demo_p->m_curves_list.end(); ++itp)
        {
          Pm_xconic_2 * seg = new Pm_xconic_2( **itp );
          w_demo_p_new->m_curves_list.push_back(seg);
          //conic_list.push_back(*seg);
          Pm_base_conic_2* org_conic = w_demo_p_new->m_tab_traits.get_origin_curve(*seg);
          w_demo_p_new->m_curves_arr.insert(Pm_conic_2( *org_conic , seg->get_data() ));

          w_demo_p_new->draw_curve( **itp);
        }
      }
      
      //w_demo_p_new->m_curves_arr.insert(conic_list.begin(),conic_list.end());

      *w_demo_p_new << CGAL::RED;
      *w_demo_p_new << CGAL::DISC;
      // Go over all vertices and for each vertex print the ID numbers of the
      // base curves that go through it.
      Conic_arr::Vertex_iterator   vit;
      for (vit = w_demo_p_new->m_curves_arr.vertices_begin(); vit != w_demo_p_new->m_curves_arr.vertices_end(); vit++)
      {
        Conic_arr::Halfedge_around_vertex_circulator 
        eit, first = (*vit).incident_halfedges();

        eit = first;

        int ind1;
        int ind2 = (*eit).curve().get_data().m_index;
        do 
        {
          ind1 = (*eit).curve().get_data().m_index;

          // Keep track of IDs we haven't seen before.
          if (ind1 != ind2)
          {
            const Pm_conic_point_2& p = (*vit).point();
            *w_demo_p_new << p;
            break;
          }

          eit++;

        } while (eit != first);
      }

      // allow the user to see the overlay
      w_demo_p_new->current_state = old_state;
  
      // update new planner map index
      Pm_xconic_iter iter;
      for (iter = w_demo_p_new->m_curves_list.begin(); iter != w_demo_p_new->m_curves_list.end(); ++iter)
      {
        Curve_conic_data cd = (**iter).get_data();
        cd.m_type = Curve_conic_data::INTERNAL;
        cd.m_index = w_demo_p_new->index;
        (**iter).set_data( cd );
      }

      break;
    }

  }
}

int MyWindow::realIndex(int index)
{
  Qt_widget_base_tab * w_demo_p;
  for (int i = 0; i <= tab_number; i++)
  {
    if ( myBar->isTabEnabled( myBar->page(i) ) )
    {
      // We peform downcasting from QWigdet* to Qt_widget_base_tab*, as we know that only
      // Qt_widget_base_tab objects are stored in the tab pages.
      w_demo_p = static_cast<Qt_widget_base_tab *> (myBar->page(i));
      if (w_demo_p->index == index)
        return i;
    }
  }
  return -1;
}


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



