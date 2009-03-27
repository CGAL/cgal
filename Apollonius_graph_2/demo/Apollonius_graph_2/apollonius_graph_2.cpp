// Copyright (c) 2003,2004,2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
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

#include <CGAL/basic.h>

#include <iostream>
#include <fstream>
#include <cassert>

#include <qapplication.h>
#include <qmainwindow.h>
#include <qpixmap.h>
#include <qlabel.h>
#include <qlayout.h>

//#include <CGAL/IO/Color.h>
#include <CGAL/Timer.h>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_get_circle.h>
#include <CGAL/IO/Qt_widget_get_point.h>
#include <CGAL/IO/Qt_widget_Apollonius_site_2.h>
#include <CGAL/IO/pixmaps/demoicon.xpm>


#include "typedefs.h"

#include "qt_file_toolbar.h"
#include "qt_layers_toolbar.h"
#include "qt_layers.h"
#include "edit_vertex_layer.h"


//************************************
// global variables
//************************************
AG_2 ag;

//************************************
// conversion functions
//************************************
inline Apollonius_site_2
to_site(const Circle_2 &c)
{
  double r = CGAL_NTS sqrt(CGAL_NTS to_double(c.squared_radius()));
  return  Apollonius_site_2(c.center(), Rep::RT(r));
}

inline Circle_2
to_circle(const Apollonius_site_2 &wp)
{
  return Circle_2(wp.point(), CGAL_NTS square(wp.weight()) );
}


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
  CGAL::Qt_widget_get_circle<Rep> get_circle;
  CGAL::Qt_widget_get_point<Rep> get_point;
  bool is_edit_mode;
  bool is_remove_mode;
  bool is_insert_point_mode;
  QString title_;
  QString qmsg;
  char msg[300];

public:
  My_Window(int x, int y)
  {
    is_edit_mode = false;
    is_remove_mode = false;
    is_insert_point_mode = false;

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

    file_toolbar = new File_toolbar("File Operations",
				    this, this, FALSE,
				    "File Operations");

    layers_toolbar = new Layers_toolbar(widget->get_qt_widget(), ag,
					"Geometric Operations",
					this, this, FALSE,
					"Geometric Operations");

    this->addToolBar(stoolbar, Top, FALSE);
    this->addToolBar(file_toolbar, Top, FALSE);
    this->addToolBar(layers_toolbar, Top, FALSE);

    connect(widget->get_qt_widget(), SIGNAL(new_cgal_object(CGAL::Object)),
	    this, SLOT(get_object(CGAL::Object)));

    connect(layers_toolbar, SIGNAL(inputModeChanged(bool)), this,
	    SLOT(get_input_mode(bool)));

    connect(layers_toolbar, SIGNAL(removeModeChanged(bool)), this,
	    SLOT(get_remove_mode(bool)));

    connect(layers_toolbar, SIGNAL(editModeChanged(bool)), this,
    	    SLOT(get_edit_mode(bool)));

    connect(file_toolbar, SIGNAL(fileToRead(const QString&)), this,
	    SLOT(read_from_file(const QString&)));

    connect(file_toolbar, SIGNAL(fileToWrite(const QString&)), this,
	    SLOT(write_to_file(const QString&)));

    connect(file_toolbar, SIGNAL(printScreen()), this,
	    SLOT(print_screen()));
    connect(file_toolbar, SIGNAL(clearAll()), this,
	    SLOT(remove_all()));

    widget->attach(&get_circle);

    setMouseTracking(true);
    widget->setMouseTracking(true);

    widget->attach(&get_point);

    get_circle.activate();
    get_point.deactivate();


    // Adding menus

    // file menu
    QPopupMenu* file = new QPopupMenu(this);
    menuBar()->insertItem("&File", file);
    file->insertItem("&Clear", this, SLOT(remove_all()), CTRL+Key_C);
    file->insertSeparator();
    file->insertItem("&Load Apollonius graph", this,
		     SLOT(open_from_file()), CTRL+Key_O);
    file->insertItem("&Save Apollonius graph", this,
		     SLOT(save_to_file()), CTRL+Key_S);
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

    title_ = tr("Apollonius graph 2");

    //    get_circle.setMouseTracking(true);
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

private slots:
  void get_object(CGAL::Object obj)
  {
    set_msg("");
    if ( is_edit_mode ) { return; }
    if ( is_remove_mode ) {
      if ( ag.number_of_vertices() == 0 ) { return; }
      Point_2 p;
      if ( CGAL::assign(p, obj) ) {
	Vertex_handle v = ag.nearest_neighbor(p);
	if ( v != Vertex_handle() ) {
	  qmsg = "Removing input site...";
	  set_msg(qmsg);
	  ag.remove(v);
	  qmsg = qmsg + " done! - Validating Apollonius diagram...";
	  set_msg(qmsg);
	  assert( ag.is_valid(false,1) );
	  qmsg = qmsg + " done!";
	  set_msg(qmsg);
	}
      }
      widget->redraw();
      return;
    }
    Circle_2 c;
    Point_2 p;

    CGAL::Timer timer;
    if ( CGAL::assign(c, obj) ) {
      Apollonius_site_2 wp = to_site(c);
      ag.insert(wp);
    } else if ( CGAL::assign(p, obj) ) {
      Apollonius_site_2 wp(p, Weight(0));
      ag.insert(wp);
    }
    std::sprintf(msg, "Insertion time: %f.", timer.time());
    qmsg = QString(msg) + " - Validating Apollonius graph...";
    set_msg(qmsg);
    assert( ag.is_valid(false, 1) );
    qmsg = qmsg + " done!";
    set_msg(qmsg);
    widget->redraw();
  }

  void get_input_mode(bool b)
  {
    is_insert_point_mode = b;

    if ( !is_remove_mode && !is_edit_mode ) {
      if ( is_insert_point_mode ) {
	get_point.activate();
	get_circle.deactivate();
      } else {
	get_point.deactivate();
	get_circle.activate();
      }
    }
  }

  void get_remove_mode(bool b)
  {
    is_remove_mode = b;

    if ( is_remove_mode ) {
      get_point.activate();
      get_circle.deactivate();
    } else {
      if ( !is_insert_point_mode ) {
	get_point.deactivate();
	get_circle.activate();
      }
    }
  }

  void get_edit_mode(bool b)
  {
    is_edit_mode = b;

    if ( is_edit_mode ) {
      get_point.activate();
      get_circle.deactivate();
    } else {
      if ( !is_insert_point_mode ) {
	get_point.deactivate();
	get_circle.activate();
      }
    }
  }

  void read_from_file(const QString& fileName)
  {
    set_msg("");
    std::ifstream f(fileName);
    assert( f );

    Apollonius_site_2 wp;

    int counter = 0;
    std::cout << std::endl;

    CGAL::Timer timer;
    timer.start();
    while ( f >> wp ) {
      ag.insert(wp);
      counter++;
      if ( counter % 500 == 0 ) {
	std::sprintf(msg, "%d have been inserted...", counter);
	set_msg(msg);
      }
    }
    timer.stop();

    std::sprintf(msg,
			   "%d sites inserted. Insertion time: %f "
			   "- Validating Apollonius graph...",
			   counter, timer.time());
    set_msg(msg);

    assert( ag.is_valid(false, 1) );

    qmsg = QString(msg) + " done!";
    set_msg(qmsg);
    widget->redraw();
  }

  void write_to_file(const QString& fileName)
  {
    set_msg("");
    std::ofstream f(fileName);
    assert( f );
    f.precision(18);

    QString qmsg = "Writing input sites to file...";
    set_msg(qmsg);
    for( AG_2::Sites_iterator it = ag.sites_begin();
	 it != ag.sites_end(); it++ ) {
      f << (*it) << std::endl;
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
      QFileDialog::getSaveFileName(tr("data.out"), QString::null,
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
      qmsg = "Reading Apollonius graph from file...";
      set_msg(qmsg);

      std::ifstream f(fileName);
      assert(f);

      CGAL::Timer timer;
      timer.start();
      f >> ag;
      timer.stop();

      int n_sites = static_cast<int>(ag.number_of_vertices()) +
	static_cast<int>(ag.number_of_hidden_sites());
      std::sprintf(msg,
			     "%d sites inserted. Insertion time: %f",
			     n_sites, timer.time());
      qmsg = qmsg + " done! " + msg;
      set_msg(qmsg);
      widget->redraw();
    }
  }

  void save_to_file()
  {
    set_msg("");
    QString fileName =
      QFileDialog::getSaveFileName(tr("data.out"), QString::null,
				   this, "Save as...");

    if ( !fileName.isNull() ) {
      qmsg = "Saving Apollonius graph to file...";
      set_msg(qmsg);
      std::ofstream f(fileName);
      assert(f);
      f << ag;
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
    ag.clear();
    widget->redraw();
  }


  void about()
  {
    QMessageBox::about( this, get_title(),
			"This is a demo for the 2D Apollonius graph\n\n"
			"Author: Menelaos Karavelas <mkaravel@tem.uoc.gr>\n\n"
			"Copyright(c) INRIA 2003,2004,2005");
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

#include "edit_vertex_layer.moc"
#include "qt_file_toolbar.moc"
#include "qt_layers_toolbar.moc"
#include "apollonius_graph_2.moc"

int
main(int argc, char* argv[])
{
  int size = 750;

  QApplication app( argc, argv );
  My_Window W(size,size);
  app.setMainWidget( &W );
#if !defined (__POWERPC__)
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  W.setIcon(cgal_icon);
#endif
  W.show();
  W.set_window(0,size,0,size);
  W.setCaption( W.get_title() );
  W.setMouseTracking(TRUE);
  return app.exec();
}


// moc_source_file: apollonius_graph_2.C
