#include "demo1.h"
#include "forms.h"
#include "qt_layer.h"
#include "demo_tab.h"


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
	if (strcmp(remove_mode,"Remove All Original Curve") == 0)
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
    FileOpenOptionsForm *form = new FileOpenOptionsForm(flag);
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
    FileOpenOptionsForm * form = new FileOpenOptionsForm(false);
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
    FileOpenOptionsForm * form = new FileOpenOptionsForm(flag);
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
  {
    load( filename , clear_flag);
	updateMode( dragMode );
  }
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
    	w_demo_p->empty = true;
      else 
        w_demo_p->empty = false;
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
  updateMode( dragMode );
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
	  w_demo_p->empty = true;
    else 
      w_demo_p->empty = false;    
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
    std::list<Pm_pol_2> pol_list;

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
      
	  Pm_pol_2 curve(*base_polyline, cd);
      pol_list.push_back(curve);
      //w_demo_p->m_curves_arr.insert(Pm_pol_2( *base_polyline , cd));
    }
	 w_demo_p->m_curves_arr.insert(pol_list.begin(),pol_list.end());

    if( w_demo_p->m_curves_arr.number_of_vertices() == 0 )
      w_demo_p->empty = true;
    else 
      w_demo_p->empty = false;
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
    std::list<Pm_seg_2> seg_list;
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
      
	  Pm_seg_2 curve(*base_seg, cd);
      seg_list.push_back(curve);
      //w_demo_p->m_curves_arr.insert(Pm_seg_2( *base_seg , cd));
    }
	w_demo_p->m_curves_arr.insert(seg_list.begin(),seg_list.end());

	if( w_demo_p->m_curves_arr.number_of_vertices() == 0 )
	  w_demo_p->empty = true;
    else 
      w_demo_p->empty = false;
      
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
