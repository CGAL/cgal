#ifndef MY_WINDOW_H
#define MY_WINDOW_H

//************************************
// my window
//************************************
class My_Window : public QMainWindow {
  Q_OBJECT

  friend class Layers_toolbar;
private:
  CGAL::Qt_widget *widget;
  Layers_toolbar *layers_toolbar;
  File_toolbar   *file_toolbar;
  CGAL::Qt_widget_standard_toolbar *stoolbar;
  CGAL::Qt_widget_get_segment<Rep> get_segment;
  CGAL::Qt_widget_get_point<Rep> get_point;
  CGAL::Qt_widget_get_simple_polygon<Polygon_2> get_polygon;
  Input_mode input_mode;
  bool is_remove_mode;
  bool is_snap_mode;

public:
  My_Window(int x, int y)
  {

    //******************
    num_selected = 0;

    is_remove_mode = false;
    input_mode = SVD_SEGMENT;
    is_snap_mode = false;

    widget = new CGAL::Qt_widget(this);
    setCentralWidget(widget);

    *widget << CGAL::BackgroundColor(CGAL::BLACK);
    resize(x,y);
    widget->set_window(0, x, 0, y);
    widget->show();

    //    setUsesBigPixmaps(TRUE);

    //How to attach the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar(widget, this,
						    this, FALSE, "");

    file_toolbar = new File_toolbar("File operations",
				    this, this, FALSE,
				    "File operations");

    layers_toolbar = new Layers_toolbar(widget, svd,
					"Geometric Operations",
					this, this, FALSE,
					"Geometric Operations");

    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)), this,
	    SLOT(get_object(CGAL::Object)));

    connect(layers_toolbar, SIGNAL(inputModeChanged(Input_mode)), this,
    	    SLOT(get_input_mode(Input_mode)));

    connect(layers_toolbar, SIGNAL(insertModeChanged(bool)), this,
    	    SLOT(get_remove_mode(bool)));

    connect(layers_toolbar, SIGNAL(snapModeChanged(bool)), this,
    	    SLOT(get_snap_mode(bool)));

    connect(file_toolbar, SIGNAL(fileToRead(const QString&)), this,
	    SLOT(read_from_file(const QString&)));

    connect(file_toolbar, SIGNAL(fileToWrite(const QString&)), this,
	    SLOT(write_to_file(const QString&)));

    connect(file_toolbar, SIGNAL(printScreen()), this,
	    SLOT(print_screen()));

    connect(file_toolbar, SIGNAL(clearAll()), this,
	    SLOT(remove_all()));

    setMouseTracking(true);
    widget->setMouseTracking(true);

    widget->attach(&get_point);
    widget->attach(&get_segment);
    widget->attach(&get_polygon);

    get_segment.activate();
    get_point.deactivate();
    get_polygon.deactivate();
  }

  ~My_Window(){}

  void set_window(double xmin, double xmax,
		  double ymin, double ymax)
  {
    widget->set_window(xmin, xmax, ymin, ymax);
  }

private:
  void get_object_remove_mode(CGAL::Object obj)
  {
    std::cout << "in remove mode" << std::endl;
#if 1
    if ( svd.number_of_vertices() == 0 ) { return; }

    Point_2 p;
    if ( CGAL::assign(p, obj) ) {
      SVD_2::Vertex_handle v = svd.nearest_neighbor(p);
#if 1
      std::cout << "degree: " << v->degree() << std::endl;
      if ( v->site().is_segment() &&
	   !v->site().is_input() ) {
	std::cout << "site: " << v->site() << std::endl;
	std::cout << "supporting segment: "
		  << v->site().supporting_site().segment() << std::endl;
	if ( !v->site().is_input(0) ) {
	  std::cout << "crossing segment for source: "
		    << v->site().crossing_site(0).segment() << std::endl;
	}
	if ( !v->site().is_input(1) ) {
	  std::cout << "crossing segment for target: "
		    << v->site().crossing_site(1).segment() << std::endl;
	}
	SVD_2::Vertex_circulator vc = svd.incident_vertices(v);
	SVD_2::Vertex_circulator vc_start = vc;
	do {
	  SVD_2::Vertex_handle vv(vc);
	  if ( !svd.is_infinite(vc) &&
	       vv->site().is_point() &&
	       (vv->site().point() == v->site().source() ||
		vv->site().point() == v->site().target()) ) {
	    std::cout << "degree of endpoint " << vv->site()
		      << " : " << vv->degree() << std::endl;
	  }
	  ++vc;
	} while ( vc_start != vc );
      }

      widget->redraw();

      SVD_2::Vertex_circulator vc = svd.incident_vertices(v);
      SVD_2::Vertex_circulator vc_start = vc;
      *widget << CGAL::GREEN;
      do {
	SVD_2::Vertex_handle vv(vc);
	if ( !svd.is_infinite(vc) ) {
	  SVD_2::Site_2 site = vv->site();
	  if ( site.is_segment() ) {
	    *widget << site.segment();
	  } else {
	    *widget << site.point();
	  }
	}
	++vc;
      } while ( vc_start != vc );
#else
      *widget << CGAL::BLACK;
      if ( v->is_point() ) {
	std::cout << v->point() << std::endl;
	*widget << v->point();
      } else {
	std::cout << v->segment() << std::endl;
	*widget << v->segment();
      }
#endif
    }
#else
    if ( svd.number_of_vertices() == 0 ) { return; }

    Point_2 q;

    if ( CGAL::assign(q, obj) ) {
      SVD_2::Vertex_handle v = svd.nearest_neighbor(q);
      svd.remove(v, is_snap_mode);
      widget->redraw();
    }
#endif
  }

private slots:
  void get_object(CGAL::Object obj)
  {
    if ( is_remove_mode ) {
      get_object_remove_mode(obj);
      return;
    }

    CGAL::Timer timer;
    char msg[100];

    if ( input_mode == SVD_POINT ) {
      if ( is_snap_mode ) {
	Point_2 p;
	if ( CGAL::assign(p, obj) ) {
	  SVD_2::Vertex_handle v;
	  v = svd.nearest_neighbor(p);
	}
	return;
      }

      Point_2 p;
      if ( CGAL::assign(p, obj) ) {
	timer.start();
	insert_point(svd, p);
	timer.stop();

	CGAL_CLIB_STD::sprintf(msg, "Insertion time: %f", timer.time());
	statusBar()->message(msg);
      }
    } else if ( input_mode == SVD_SEGMENT ) {
      Segment s;

      if ( is_snap_mode ) {	
	if( CGAL::assign(s, obj) ) {
	  SVD_2::Vertex_handle v1, v2;
	  v1 = svd.nearest_neighbor(s.source());
	  v2 = svd.nearest_neighbor(s.target());
	  if ( v1 != NULL && v1->is_point() &&
	       v2 != NULL && v2->is_point() ) {
	    timer.start();
	    insert_segment(svd, v1->site().point(), v2->site().point() );
	    timer.stop();

	    CGAL_CLIB_STD::sprintf(msg,	"Insertion time: %f", timer.time());
	    statusBar()->message(msg);
	  }
	}
      } else {
	if( CGAL::assign(s, obj) ) {
	  timer.start();
	  insert_segment(svd, s.source(), s.target());
	  timer.stop();

	  CGAL_CLIB_STD::sprintf(msg,	"Insertion time: %f", timer.time());
	  statusBar()->message(msg);
	}
      }
    } else if ( input_mode == SVD_POLYGON ) {
      Polygon_2 pgn;
      if ( CGAL::assign(pgn, obj) ) {
	timer.start();
	insert_polygon(svd, pgn);
	timer.stop();

	CGAL_CLIB_STD::sprintf(msg,	"Insertion time: %f", timer.time());
	statusBar()->message(msg);
      }
    }

    //    svd.is_valid(true,1);
    //    std::cout << std::endl;
    widget->redraw();
  }

  void get_remove_mode(bool b)
    {
      is_remove_mode = b;

      if ( is_remove_mode ) {	
	get_point.activate();
	get_segment.deactivate();
	get_polygon.deactivate();
      } else {
	if ( input_mode == SVD_SEGMENT ) {
	  get_point.deactivate();
	  get_segment.activate();
	} else if ( input_mode == SVD_POLYGON ) {
	  get_point.deactivate();
	  get_polygon.activate();
	}
      }

    }


  void get_input_mode(Input_mode im) {
    input_mode = im;

    if ( input_mode == SVD_POINT ) {
      get_point.activate();
      get_segment.deactivate();
      get_polygon.deactivate();
    } else if ( input_mode == SVD_SEGMENT ) {
      get_point.deactivate();
      get_segment.activate();
      get_polygon.deactivate();
    } else if ( input_mode == SVD_POLYGON ) {
      get_point.deactivate();
      get_segment.deactivate();
      get_polygon.activate();
    }
  }

  void get_snap_mode(bool b) {
    is_snap_mode = b;
  }

  void read_from_file(const QString& fileName)
  {
    typedef SVD_2::Vertex_handle Vertex_handle;
    CGAL::Timer timer;
    svd.clear();

    std::ifstream f(fileName);
    assert( f );

    int counter = 0;
    timer.start();

    char msg[100];
    bool bbox_empty = true;
    CGAL::Bbox_2 bbox;
    char type;
    while (f >> type) {
      CGAL::Bbox_2 tbox;
      if (type == 'p') {
        Point_2 p;
        f >> p;
        insert_point(svd, p);
	tbox = p.bbox();
	counter++;
      } else if (type == 's') {
        Point_2 p1, p2;
        f >> p1 >> p2;
        insert_segment(svd, p1, p2);
        tbox = Segment(p1,p2).bbox();
	counter++;
      } else if (type == 'l') {
	Vertex_handle vh;
        int nr_of_points;
        f >> nr_of_points;
        Point_2 p1, p2;
        f >> p1;
	tbox = p1.bbox();
	bool got_location = false;
        while(--nr_of_points!=0){
	  f >> p2;
	  if(!got_location){
	    vh = insert_segment(svd, p1, p2);
	    got_location = true;
	  } else
	    vh = insert_segment(svd, p1, p2, vh);
	  tbox = tbox + Segment(p1,p2).bbox();
	  counter++;
	  p1 = p2;
        }
      }
      
      if(bbox_empty) {
	bbox = tbox;
	bbox_empty = false;
      } else {
	bbox = bbox + tbox;
      }

      if ( counter % 500 == 0 ) {
	std::cout << "\r" << counter
		  << " sites haved been inserted..." << std::flush;

	sprintf(msg, "%d sites have been inserted...", counter);
	statusBar()->message(msg);
      }
    }//endwhile

    std::cout << "\r" << counter
	      << " sites haved been inserted... Done!" << std::endl;

    timer.stop();
    std::cout << "Insertion time: " << timer.time() << std::endl;

    svd.is_valid(true, 1);
    std::cerr << std::endl;

    CGAL_CLIB_STD::sprintf(msg,
			   "%d sites inserted. Insertion time: %f",
			   counter, timer.time());
    statusBar()->message(msg);

    double width =  bbox.xmax() - bbox.xmin();
    double height =  bbox.ymax() - bbox.ymin();
    double s = 0.1;
    set_window(bbox.xmin() - s * width,	bbox.xmax() + s * width,
	       bbox.ymin() - s * height, bbox.ymax() + s * height);
    widget->redraw();
    //    widget->clear_history();
  }

  void write_to_file(const QString& fileName)
  {
    std::ofstream f(fileName);
    assert( f );
    f.precision(18);
    SVD_2::Finite_vertices_iterator vit;
    for (vit = svd.finite_vertices_begin();
	 vit != svd.finite_vertices_end(); ++vit) {
      f << vit->site() << std::endl;
    }
    //    svd.write_sites(f);
  }

  void print_screen()
  {
    widget->print_to_ps();
  }

  void remove_all()
    {
      sitelist.clear();
      num_selected = 0;
      svd.clear();
      widget->redraw();
    }

};


#endif //  MY_WINDOW_H
