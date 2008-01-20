// Copyright (c) 2004,2005  INRIA Sophia-Antipolis (France) and
// Notre Dame University (U.S.A.).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>

#ifndef MY_WINDOW_H
#define MY_WINDOW_H

#include <cassert>

#include <qlayout.h>
#include <qlabel.h>

#include "Which_diagram.h"

//************************************
// Layout_widget
//************************************

class Layout_widget
  : public QWidget
{
 public:
  Layout_widget(QWidget *parent, const char *name=0)
    : QWidget(parent, name)
  {
    QBoxLayout* topLayout = new QVBoxLayout(this, QBoxLayout::TopToBottom);

      // create/initialize the label
    label = new QLabel(this, "label");
    label->setText("");

    // create/initialize Qt_widget
    widget = new CGAL::Qt_widget(this);

    // add widgets to layout
    topLayout->addWidget(widget, 1);
    topLayout->addWidget(label, 0);
  }

  ~Layout_widget() {}

  CGAL::Qt_widget* get_qt_widget() { return widget; }
  QLabel*          get_label() { return label; }

  // methods to access functionality of CGAL::Qt_widget;
  template<class T>
  void attach(const T& t) { widget->attach(t); }
  void redraw() { widget->redraw(); }
  void print_to_ps() { widget->print_to_ps(); }
  void set_window(double xmin, double xmax, double ymin, double ymax) {
    widget->set_window(xmin, xmax, ymin, ymax);
  }

 private:
  CGAL::Qt_widget *widget;
  QLabel          *label;
};

template<class T>
Layout_widget& operator<<(Layout_widget& lw, const T& t)
{
  *lw.get_qt_widget() << t;
  return lw;
}


//************************************
// my window
//************************************
class My_Window : public QMainWindow {
  Q_OBJECT

  friend class Layers_toolbar;
private:
  Layout_widget  *widget;
  Layers_toolbar *layers_toolbar;
  File_toolbar   *file_toolbar;
  CGAL::Qt_widget_standard_toolbar *stoolbar;
  CGAL::Qt_widget_get_segment<Rep> get_segment;
  CGAL::Qt_widget_get_point<Rep> get_point;
  CGAL::Qt_widget_get_simple_polygon<Polygon_2> get_polygon;
  Input_mode input_mode;
  bool is_remove_mode;
  bool is_snap_mode;
  QString title_;
  QString qmsg;
  char msg[300];
  bool is_pvd_;

public:
  My_Window(int x, int y, bool is_pvd = false)
  {
    is_pvd_ = is_pvd;

    //******************
    num_selected = 0;

    is_remove_mode = false;
    if ( is_pvd_ ) {
      input_mode = SDG_POLYGON;
    } else {
      input_mode = SDG_SEGMENT;
    }
    is_snap_mode = false;

    widget = new Layout_widget(this);
    setCentralWidget(widget);

    *widget << CGAL::BackgroundColor(CGAL::BLACK);
    resize(x,y);
    widget->set_window(0, x, 0, y);
    widget->show();

    //    setUsesBigPixmaps(TRUE);

    //How to attach the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar(widget->get_qt_widget(),
						    this, this, FALSE, "");

    file_toolbar = new File_toolbar("File operations",
				    this, this, FALSE,
				    "File operations");

    layers_toolbar = new Layers_toolbar(widget->get_qt_widget(), sdg,
					"Geometric Operations",
					this, this, FALSE,
					"Geometric Operations", 0, is_pvd_);

    this->addToolBar(stoolbar, Top, FALSE);
    this->addToolBar(file_toolbar, Top, FALSE);
    this->addToolBar(layers_toolbar, Top, FALSE);

    connect(widget->get_qt_widget(), SIGNAL(new_cgal_object(CGAL::Object)),
	    this, SLOT(get_object(CGAL::Object)));

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

    get_point.deactivate();
    if ( is_pvd_ ) {
      get_segment.deactivate();
      get_polygon.activate();
    } else {
      get_segment.activate();
      get_polygon.deactivate();
    }

    // Adding menus

    // file menu
    QPopupMenu* file = new QPopupMenu(this);
    menuBar()->insertItem("&File", file);
    file->insertItem("&Clear", this, SLOT(remove_all()), CTRL+Key_C);
    file->insertSeparator();
    if ( is_pvd_ ) {
      file->insertItem("&Load polygon Delaunay graph", this,
		       SLOT(open_from_file()), CTRL+Key_O);
    } else {
      file->insertItem("&Load segment Delaunay graph", this,
		       SLOT(open_from_file()), CTRL+Key_O);
    }
    if ( is_pvd_ ) {
      file->insertItem("&Save polygon Delaunay graph", this,
		       SLOT(save_to_file()), CTRL+Key_S);
    } else {
      file->insertItem("&Save segment Delaunay graph", this,
		       SLOT(save_to_file()), CTRL+Key_S);
    }
    file->insertSeparator();
    file->insertItem("&Read input data", this,
		     SLOT(read_input_from_file()), CTRL+Key_R);
    file->insertItem("&Save output data", this,
		     SLOT(write_output_to_file()), CTRL+Key_W);
    file->insertSeparator();
    file->insertItem("Print", this, SLOT(print_screen()), CTRL+Key_P);
    //    file->insertSeparator();
    //    file->insertItem("&Quit",this,SLOT(closeAll()),CTRL+Key_Q);

    // view menu
    QPopupMenu* view = new QPopupMenu(this);
    QPopupMenu* bcolor = new QPopupMenu(view);
    menuBar()->insertItem("&View", view);
    view->insertItem("Background Color", bcolor);
    bcolor->insertItem("White", this, SLOT(change_bcolor_to_white()));
    bcolor->insertItem("Black", this, SLOT(change_bcolor_to_black()));
    bcolor->insertItem("Yellow", this, SLOT(change_bcolor_to_yellow()));

    // about menu
    QPopupMenu* about = new QPopupMenu(this);
    menuBar()->insertItem("&About", about);
    about->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
    about->insertItem("About &Qt", this, SLOT(aboutQt()) );

    if ( is_pvd_ ) {
      title_ = tr("Polygon Delaunay graph 2");
    } else {
      title_ = tr("Segment Delaunay graph 2");
    }
  }

  ~My_Window(){}

  void set_window(double xmin, double xmax,
		  double ymin, double ymax)
  {
    widget->set_window(xmin, xmax, ymin, ymax);
  }

  const QString& get_title() const { return title_; }
  void set_title(const QString& title) { title_ = title; }

private:
  void set_msg(const QString& str) {
    widget->get_label()->setText(str);
  }

  void get_object_remove_mode(CGAL::Object obj)
  {
    if ( sdg.number_of_vertices() == 0 ) { return; }

    if ( !is_snap_mode ) {
      Point_2 q;

      if ( CGAL::assign(q, obj) ) {
	CGAL::Timer timer;
	SDG_2::Vertex_handle v = sdg.nearest_neighbor(q);
	//      sdg.remove(v, is_snap_mode);
	timer.start();
	bool success = sdg.remove(v);
	timer.stop();
	std::sprintf(msg, "Removal time: %f", timer.time());
	if ( success ) {
	  set_msg(QString("Removal was successful") + " - " + msg);
	} else {
	  set_msg("Removal was unsuccessful");
	}
	if ( success ) {
	  sdg.is_valid(true, 1);
	  std::cerr << std::endl;
	  widget->redraw();
	}
      }
      return;
    }

    if ( sdg.number_of_vertices() == 0 ) { return; }

    Point_2 p;
    if ( CGAL::assign(p, obj) ) {
      SDG_2::Vertex_handle v = sdg.nearest_neighbor(p);

      std::cerr << "degree: " << sdg.data_structure().degree(v) << std::endl;

      if ( sdg.dimension() == 2 ) {
	SDG_2::Face_circulator fc_start = sdg.incident_faces(v);
	SDG_2::Face_circulator fc = fc_start;
	do {
	  int id = fc->index(v);
	  SDG_2::Vertex_handle vopp = sdg.tds().mirror_vertex(fc, id);
	  std::cerr << "Testing incircle: " << std::endl;
	  if ( !sdg.is_infinite(fc->vertex(0)) &&
	       !sdg.is_infinite(fc->vertex(1)) &&
	       !sdg.is_infinite(fc->vertex(2)) &&
	       !sdg.is_infinite(vopp) ) {
	    SDG_2::Geom_traits::Vertex_conflict_2 incircle =
	      sdg.geom_traits().vertex_conflict_2_object();
	    std::cerr << "  vertices of face: " << std::endl;
	    for (int i = 0; i < 3; i++) {
	      std::cerr << "     " << fc->vertex(i)->site() << std::endl;
	    }
	    std::cerr << "  other vertex: " << vopp->site() << std::endl;
	    CGAL::Sign s = incircle(fc->vertex(0)->site(),
				    fc->vertex(1)->site(),
				    fc->vertex(2)->site(),
				    vopp->site());
	    std::cerr << "  incircle: " << int(s) << std::endl;
	  }
	} while ( ++fc != fc_start );
      }


      if ( v->site().is_segment() &&
	   !v->site().is_input() ) {
	std::cerr << "site: " << v->site() << std::endl;
	std::cerr << "supporting segment: "
		  << v->site().supporting_site().segment() << std::endl;
	if ( !v->site().is_input(0) ) {
	  std::cerr << "crossing segment for source: "
		    << v->site().crossing_site(0).segment() << std::endl;
	}
	if ( !v->site().is_input(1) ) {
	  std::cerr << "crossing segment for target: "
		    << v->site().crossing_site(1).segment() << std::endl;
	}
	SDG_2::Vertex_circulator vc = sdg.incident_vertices(v);
	SDG_2::Vertex_circulator vc_start = vc;
	do {
	  SDG_2::Vertex_handle vv(vc);
	  if ( !sdg.is_infinite(vc) &&
	       vv->site().is_point() &&
	       (vv->site().point() == v->site().source() ||
		vv->site().point() == v->site().target()) ) {
	    std::cerr << "degree of endpoint " << vv->site()
		      << " : " << sdg.data_structure().degree(vv)
		      << std::endl;
	  }
	  ++vc;
	} while ( vc_start != vc );
      }

      widget->redraw();

      SDG_2::Vertex_circulator vc = sdg.incident_vertices(v);
      SDG_2::Vertex_circulator vc_start = vc;
      *widget << CGAL::GREEN;
      do {
	SDG_2::Vertex_handle vv(vc);
	if ( !sdg.is_infinite(vc) ) {
	  SDG_2::Site_2 site = vv->site();
	  if ( site.is_segment() ) {
	    *widget << site.segment();
	  } else {
	    *widget << site.point();
	  }
	}
	++vc;
      } while ( vc_start != vc );
    }

  }

private slots:
  void get_object(CGAL::Object obj)
  {
    set_msg("");
    if ( is_remove_mode ) {
      get_object_remove_mode(obj);
      return;
    }

    CGAL::Timer timer;

    bool is_polygon = false;
    if ( input_mode == SDG_POINT ) {
      if ( is_snap_mode ) {
	Point_2 p;
	if ( CGAL::assign(p, obj) ) {
	  SDG_2::Vertex_handle v;
	  v = sdg.nearest_neighbor(p);
	}
	return;
      }

      Point_2 p;
      if ( CGAL::assign(p, obj) ) {
	timer.start();
	insert_point(sdg, p);
	timer.stop();

	std::sprintf(msg, "Insertion time: %f", timer.time());
	set_msg(msg);
      }
    } else if ( input_mode == SDG_SEGMENT ) {
      Segment s;

      if ( is_snap_mode ) {
	if( CGAL::assign(s, obj) ) {
	  SDG_2::Vertex_handle v1, v2;
	  v1 = sdg.nearest_neighbor(s.source());
	  v2 = sdg.nearest_neighbor(s.target());
	  if ( v1 != NULL && v1->is_point() &&
	       v2 != NULL && v2->is_point() ) {
	    timer.start();
	    insert_segment(sdg, v1->site().point(), v2->site().point() );
	    timer.stop();

	    std::sprintf(msg,	"Insertion time: %f", timer.time());
	    set_msg(msg);
	  }
	}
      } else {
	if( CGAL::assign(s, obj) ) {
	  timer.start();
	  insert_segment(sdg, s.source(), s.target());
	  timer.stop();

	  std::sprintf(msg, "Insertion time: %f", timer.time());
	  set_msg(msg);
	}
      }
    } else if ( input_mode == SDG_POLYGON ) {
      is_polygon = true;
      Polygon_2 pgn;
      if ( CGAL::assign(pgn, obj) ) {
	timer.start();
	insert_polygon(sdg, pgn);
	timer.stop();

	std::sprintf(msg, "Insertion time: %f", timer.time());
	set_msg(msg);
      }
    }

    sdg.is_valid(true,1);
    std::cerr << std::endl;
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
	if ( input_mode == SDG_SEGMENT ) {
	  get_point.deactivate();
	  get_segment.activate();
	} else if ( input_mode == SDG_POLYGON ) {
	  get_point.deactivate();
	  get_polygon.activate();
	}
      }

    }


  void get_input_mode(Input_mode im) {
    input_mode = im;

    if ( input_mode == SDG_POINT ) {
      get_point.activate();
      get_segment.deactivate();
      get_polygon.deactivate();
    } else if ( input_mode == SDG_SEGMENT ) {
      get_point.deactivate();
      get_segment.activate();
      get_polygon.deactivate();
    } else if ( input_mode == SDG_POLYGON ) {
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
    set_msg("");
    typedef SDG_2::Vertex_handle Vertex_handle;
    CGAL::Timer timer;
    sdg.clear();

    std::ifstream f(fileName);
    assert( f );

    int counter = 0;
    timer.start();

    bool bbox_empty = true;
    CGAL::Bbox_2 bbox;
    char type;
    while (f >> type) {
      CGAL::Bbox_2 tbox;
      if (type == 'p') {
        Point_2 p;
        f >> p;
        insert_point(sdg, p);
	tbox = p.bbox();
	counter++;
      } else if (type == 's') {
        Point_2 p1, p2;
        f >> p1 >> p2;
        insert_segment(sdg, p1, p2);
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
	    vh = insert_segment(sdg, p1, p2);
	    got_location = true;
	  } else
	    vh = insert_segment(sdg, p1, p2, vh);
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
	std::sprintf(msg, "%d sites have been inserted...",
			       counter);
	set_msg(msg);
      }
    } // end while

    timer.stop();

    std::sprintf(msg, "%d sites inserted. Insertion time: %f",
			   counter, timer.time());
    qmsg = QString(msg) + " - Validating diagram...";
    set_msg(qmsg);

    sdg.is_valid(true, 1);
    std::cerr << std::endl;

    qmsg = qmsg + " done!";
    set_msg(qmsg);

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
    set_msg("");
    std::ofstream f(fileName);
    assert( f );
    f.precision(18);

    QString qmsg = "Writing input sites to file...";
    set_msg(qmsg);
    SDG_2::Input_sites_iterator sit;
    for (sit = sdg.input_sites_begin();
	 sit != sdg.input_sites_end(); ++sit) {
      f << (*sit) << std::endl;
    }
    qmsg = qmsg + " done!";
    set_msg(qmsg);
  }

  void read_input_from_file()
  {
    set_msg("");
    QString fileName =
      QFileDialog::getOpenFileName(QString::null, QString::null,
				   this, "Open file...");

    if ( !fileName.isNull() ) {
      read_from_file(fileName);
    }
  }

  void write_output_to_file()
  {
    set_msg("");
    QString fileName =
      QFileDialog::getSaveFileName(tr("data.cin"), QString::null,
				   this, "Save as...");

    if ( !fileName.isNull() ) {
      write_to_file(fileName);
    }
  }

  void open_from_file()
  {
    set_msg("");
    QString fileName =
      QFileDialog::getOpenFileName(QString::null, QString::null,
				   this, "Open file...");

    if ( !fileName.isNull() ) {
      if ( is_pvd_ ) {
	qmsg = "Reading polygon Delaunay graph from file...";
      } else {
	qmsg = "Reading segment Delaunay graph from file...";
      }
      set_msg(qmsg);

      std::ifstream f(fileName);
      assert(f);

      CGAL::Timer timer;
      timer.start();
      f >> sdg;
      timer.stop();

      int n_sites = static_cast<int>(sdg.number_of_input_sites());
      std::sprintf(msg,
			     "%d sites inserted. Insertion time: %f",
			     n_sites, timer.time());
      qmsg = qmsg + " done! " + msg;
      set_msg(qmsg);
      widget->redraw();
    }
  }

  QString extension(CGAL::Tag_true)  const { return "hsdg"; }
  QString extension(CGAL::Tag_false) const { return "sdg"; }

  void save_to_file()
  {
    set_msg("");

    QString ext = extension(CGAL::Which_diagram<SDG_2>::Is_hierarchy());
    QString fileName =
      QFileDialog::getSaveFileName(QString("data.") + ext,
				   QString::null,
				   this, "Save as...");

    if ( !fileName.isNull() ) {
      if ( is_pvd_ ) {
	qmsg = "Saving polygon Delaunay graph to file...";
      } else {
	qmsg = "Saving segment Delaunay graph to file...";
      }
      set_msg(qmsg);
      std::ofstream f(fileName);
      assert(f);
      f << sdg;
      qmsg = qmsg + " done!";
      set_msg(qmsg);
    }
  }

  void print_screen()
  {
    set_msg("");
    widget->print_to_ps();
  }

  void remove_all()
    {
      set_msg("");
      sitelist.clear();
      num_selected = 0;
      sdg.clear();
      widget->redraw();
    }

  void about()
  {
    QString vd_type;
    if ( is_pvd_ ) {
      vd_type = "polygon";
    } else {
      vd_type = "segment";
    }
    QMessageBox::about( this, get_title(),
			QString("This is a demo for the 2D ") + vd_type
			+ " Delaunay graph\n\n" +
			+ "Author: Menelaos Karavelas <mkaravel@tem.uoc.gr>"
			+ "\n\n"
			+ "Copyright(c) INRIA and University of Notre Dame"
			+ " 2003,2004,2005");
  }

  void aboutQt()
  {
    QMessageBox::aboutQt( this, get_title() );
  }

  void change_bcolor_to_white()  { change_bcolor(CGAL::WHITE); }
  void change_bcolor_to_black()  { change_bcolor(CGAL::BLACK); }
  void change_bcolor_to_yellow() { change_bcolor(CGAL::YELLOW); }

  void change_bcolor(const CGAL::Color& c)
  {
    *widget << CGAL::BackgroundColor(c);
    widget->redraw();
  }

};


#endif //  MY_WINDOW_H
