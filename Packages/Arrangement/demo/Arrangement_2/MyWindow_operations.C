#include <CGAL/basic.h>

#ifdef CGAL_USE_QT

#include "MyWindow.h"
#include "forms.h"
#include "qt_layer.h"
#include "demo_tab.h"

#include "icons/demo_delete.xpm"
#include "icons/draw.xpm"

#include "icons/demo_rayshoot_down.xpm"
#include "icons/demo_rayshoot_up.xpm"
#include "icons/demo_arrow_down.xpm"
#include "icons/demo_arrow_up.xpm"

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
  Qt_widget_base_tab *w_demo_p = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  
  // point location
  Coord_point p;
  if(CGAL::assign(p,obj)) 
    w_demo_p->pl_point = p;
  
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
    new PropertiesForm( myBar , this ,number_of_tabs ,w_demo_p , 
		m_scailing_factor , colors_flag);
  
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
	if (strcmp(remove_mode,"Remove entire original curve") == 0)
      w_demo_p->remove_org_curve = true;
    else
      w_demo_p->remove_org_curve = false;
	QString draw_vertex_mode = optionsForm->box9->currentText();
    if (strcmp(draw_vertex_mode,"Draw") == 0)
      w_demo_p->draw_vertex = true;
    else
      w_demo_p->draw_vertex = false;
    m_scailing_factor = (new_factor/10);
    resize(m_width,m_height);
    w_demo_p->redraw();
    something_changed();
  }
  delete optionsForm;
}


/*! print planar map */
void MyWindow::print()
{
  Qt_widget_base_tab    *w_demo_p1 = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  w_demo_p1->print_to_ps();
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


void MyWindow::backGroundColor()
{
  QColor c = QColorDialog::getColor();
  Qt_widget_base_tab     *w_demo_p = 
    dynamic_cast<Qt_widget_base_tab  *> (myBar->currentPage());
  
  w_demo_p->setBackgroundColor( c ); 
  
  something_changed();
}

void MyWindow::changePmColor()
{
  QColor c = QColorDialog::getColor();
  Qt_widget_base_tab     *w_demo_p = 
    dynamic_cast<Qt_widget_base_tab  *> (myBar->currentPage());
  
  w_demo_p->pm_color = c;
  w_demo_p->change_pm_color = true;
  something_changed();
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

/*! a dialog form to set the pointLocationStrategy */
void MyWindow::pointLocationStrategy()
{
  //TODO - implement
  /* Qt_widget_base_tab * w_demo_p = 
    dynamic_cast<Qt_widget_base_tab  *> (myBar->currentPage());*/
   PointLocationStrategyForm *form = new PointLocationStrategyForm();
   if ( form->exec() )
   {
     QString type = form->arrComboBox1->currentText();
     if(! strcmp(type,"Naive"))
       strategy = NAIVE ;
     else 
       if(!strcmp(type,"Simple"))
         strategy = SIMPLE;  
       else
         if(!strcmp(type,"Trapezoiedal"))
           strategy = TRAP;
         else
           if(!strcmp(type,"Walk"))
             strategy = WALK;
   }
}




/*! initialize the widget */
void MyWindow::init(Qt_widget_base_tab *widget)
{
  connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),
          this, SLOT(get_new_object(CGAL::Object)));
  widget->attach(testlayer);
  widget->setCursor(QCursor( QPixmap( (const char**)small_draw_xpm)));
  rayShootingMode->setIconSet(QPixmap((const char**)demo_rayshoot_up_xpm ));
  widget->setBackgroundColor(def_bg_color);
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
    (SEGMENT_TRAITS ,strategy , this, tab_number);
  init(widget);
  widget->draw();
}

/*! add a tab widget with polyline traits */
void MyWindow::add_polyline_tab()
{
  Qt_widget_demo_tab<Polyline_tab_traits> *widget = 
    new Qt_widget_demo_tab<Polyline_tab_traits>
    (POLYLINE_TRAITS ,strategy, this, tab_number);
  init(widget);
  widget->draw();
}

/*! add a tab widget with conic traits */
void MyWindow::add_conic_tab()
{
  Qt_widget_demo_tab<Conic_tab_traits> *widget = 
    new Qt_widget_demo_tab<Conic_tab_traits>
    (CONIC_TRAITS ,strategy, this , tab_number);
  init(widget);
  widget->draw();

  //widget->set_window(widget->x_pixel(widget->x_min()), widget->x_pixel(widget->x_max()), 
    //                 widget->x_pixel(widget->y_min()), widget->x_pixel(widget->y_max()));

}

/*! remove the current page (tab) from myBar */
void MyWindow::remove_tab()
{
  if (number_of_tabs > 1)
  {
    QWidget *w_demo_p = myBar->currentPage(); // w_demo_p is a pointer to Qt_widget_demo_tab object
    myBar->removePage(w_demo_p);
    delete w_demo_p;  //the destructor of Qt_widget_demo_tab will be called (virtual...) 
    number_of_tabs--;
  }
  else
  {
    QMessageBox::information( this, "Remove Tab",
                              "Can not remove last tab");
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
      (SEGMENT_TRAITS ,strategy, this);
  }
  else if (action == setPolylineTraits)
  {
    if (old_widget->traits_type == POLYLINE_TRAITS) return;
    widget = new Qt_widget_demo_tab<Polyline_tab_traits>
      (POLYLINE_TRAITS ,strategy , this);
  }
  else if (action == setConicTraits)
  {
    if (old_widget->traits_type == CONIC_TRAITS) return;
    widget = new Qt_widget_demo_tab<Conic_tab_traits>
      (CONIC_TRAITS ,strategy , this);
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
  *widget << CGAL::LineWidth(2); 
  widget->setBackgroundColor(def_bg_color);
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

/*! open color dialog for faces color */
void MyWindow::openColorDialog()
{
	Qt_widget_base_tab    *w_demo_p = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
    QColor c = QColorDialog::getColor();
	if( c.isValid())
	  w_demo_p->fill_face_color = c;
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
  else if ( action == fillfaceMode ) 
  {
    w_demo_p->mode = FILLFACE;
    w_demo_p->setCursor(Qt::CrossCursor  );	
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
   case FILLFACE: fillfaceMode->setOn( TRUE ); break;
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
  //updateSnapMode( w_demo_p->snap );
  setTraits( w_demo_p->traits_type );
  setConicType( w_demo_p->conic_type ); 

  if ( w_demo_p->snap )
  {
	setGridSnapMode->setEnabled( TRUE );
	setSnapMode->setOn( TRUE );
	if ( w_demo_p->snap_mode == GRID)
	  setGridSnapMode->setOn( TRUE );
	else
	  setGridSnapMode->setOn( FALSE );
  }
  else
  {
	setSnapMode->setOn( FALSE );
	setGridSnapMode->setOn( FALSE );
    setGridSnapMode->setEnabled( FALSE );
  }
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


#endif // CGAL_USE_QT
