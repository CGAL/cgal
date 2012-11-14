// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
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
// $URL$
// $Id$
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef MY_WINDOW_H
#define MY_WINDOW_H

#include <cassert>

#include <qlayout.h>
#include <qlabel.h>
#include <qradiobutton.h>

//#include <CGAL/IO/Color.h>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_get_circle.h>
#include <CGAL/IO/Qt_widget_get_point.h>
#include <CGAL/IO/pixmaps/demoicon.xpm>

#include "typedefs.h"
#include "qt_file_toolbar.h"
#include "qt_layers_toolbar.h"
#include "qt_layers.h"
#include "Virtual_Voronoi_diagram_2.h"

class Layout_widget
  : public QWidget
{
 public:
  Layout_widget(QWidget *parent, const char *name=0)
    : QWidget(parent, name)
  {
    QBoxLayout* topLayout = new QVBoxLayout(this, QBoxLayout::TopToBottom);

    // create/initialize the button group
    QString str = tr("Select Voronoi diagram type:");
    buttons = new QButtonGroup(3, Qt::Horizontal, str, this, "button group");
    buttons->setExclusive(true);
    QRadioButton* ag_button = new QRadioButton(buttons, "AG2");
    ag_button->setText("Apollonius");
    QRadioButton* dt_button = new QRadioButton(buttons, "DT2");
    dt_button->setText("Voronoi");
    QRadioButton* rt_button = new QRadioButton(buttons, "RT2");
    rt_button->setText("Power");
    buttons->insert(ag_button, 0);
    buttons->insert(dt_button, 1);
    buttons->insert(rt_button, 2);

    dt_button->setChecked(true);

    // create/initialize the label
    label = new QLabel(this, "label");
    label->setText("Voronoi diagram selected.");

    // create/initialize Qt_widget
    widget = new CGAL::Qt_widget(this);

    // add widgets to layout
    topLayout->addWidget(buttons, 0);
    topLayout->addWidget(widget, 1);
    topLayout->addWidget(label, 0);
  }

  ~Layout_widget(){}

  CGAL::Qt_widget* get_qt_widget() { return widget; }
  QButtonGroup*    get_button_group() { return buttons; }
  QButton*         get_button(int i) {
    switch (i) {
    case 1:
      return dt_button;
    case 2:
      return rt_button;
    default:
      return ag_button;
    }
  }
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
  QButtonGroup    *buttons;
  QButton         *ag_button;
  QButton         *dt_button;
  QButton         *rt_button;
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
class My_Window : public QMainWindow
{
  Q_OBJECT

  friend class Layers_toolbar;
 private:
  //  CGAL::Qt_widget *widget;
  Layout_widget  *widget;
  Layers_toolbar *layers_toolbar;
  File_toolbar   *file_toolbar;
  CGAL::Qt_widget_standard_toolbar *stoolbar;
  CGAL::Qt_widget_get_circle<Rep> get_circle;
  CGAL::Qt_widget_get_point<Rep> get_point;
  Input_mode input_mode;
  bool is_locate_mode;
  bool is_remove_mode;
  bool is_snap_mode;
  QString title_;
  char msg[400];
  CGAL::Concrete_Voronoi_diagram_2*     cvd;
  CGAL::Concrete_power_diagram_2*       cpd;
  CGAL::Concrete_Apollonius_diagram_2*  cad;
  CGAL::Virtual_Voronoi_diagram_2*      vvd;
  //  CGAL::Object locate_obj;
  //  ::Rep::Point_2 locate_q;

 public:
  My_Window(int x, int y)
  {
    cvd = new CGAL::Concrete_Voronoi_diagram_2();
    cpd = new CGAL::Concrete_power_diagram_2();
    cad = new CGAL::Concrete_Apollonius_diagram_2();
    vvd = cvd;

    //******************
    num_selected = 0;

    is_locate_mode = false;
    input_mode = VD_POINT;
    is_snap_mode = false;
    is_remove_mode = false;

    //    widget = new CGAL::Qt_widget(this);
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

    layers_toolbar = new Layers_toolbar(widget->get_qt_widget(), vvd,
					"Geometric Operations",
					this, this, FALSE,
					"Geometric Operations");

    this->addToolBar(stoolbar, Top, FALSE);
    this->addToolBar(file_toolbar, Top, FALSE);
    this->addToolBar(layers_toolbar, Top, FALSE);

    connect(widget->get_qt_widget(),
	    SIGNAL(new_cgal_object(CGAL::Object)), this,
	    SLOT(get_object(CGAL::Object)));

    connect(layers_toolbar, SIGNAL(inputModeChanged(Input_mode)), this,
    	    SLOT(get_input_mode(Input_mode)));

    connect(layers_toolbar, SIGNAL(insertModeChanged(bool)), this,
    	    SLOT(get_insert_mode(bool)));

    connect(layers_toolbar, SIGNAL(locateModeChanged(bool)), this,
	    SLOT(get_locate_mode(bool)));

    connect(file_toolbar, SIGNAL(fileToRead(const QString&)), this,
	    SLOT(read_from_file(const QString&)));

    connect(file_toolbar, SIGNAL(fileToWrite(const QString&)), this,
	    SLOT(write_to_file(const QString&)));

    connect(file_toolbar, SIGNAL(printScreen()), this,
	    SLOT(print_screen()));

    connect(file_toolbar, SIGNAL(clearAll()), this,
	    SLOT(remove_all()));

    connect(widget->get_button_group(),
	    SIGNAL(pressed(int)), this, SLOT(get_which_diagram(int)));

    setMouseTracking(true);
    widget->setMouseTracking(true);

    widget->attach(&get_point);
    widget->attach(&get_circle);

    get_point.activate();
    get_circle.deactivate();

    // file menu
    QPopupMenu* file = new QPopupMenu(this);
    menuBar()->insertItem("&File", file);
    file->insertItem("&Clear", this, SLOT(remove_all()), CTRL+Key_C);
    file->insertSeparator();
    file->insertItem("&Load Delaunay graph", this,
		     SLOT(open_from_file()), CTRL+Key_O);
    file->insertItem("&Save Delaunay graph", this,
		     SLOT(save_to_file()), CTRL+Key_S);
    file->insertSeparator();
    file->insertItem("&Read input data", this,
		     SLOT(read_input_from_file()), CTRL+Key_R);
    file->insertItem("&Save output data", this,
		     SLOT(write_output_to_file()), CTRL+Key_W);
    file->insertSeparator();
    file->insertItem("Print", this, SLOT(print_screen()), CTRL+Key_P);
    file->insertSeparator();
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

    title_ = tr("Voronoi diagram adaptor demo");
  }

 public:

  ~My_Window(){}

  const QString& title() const { return title_; }
  void set_title(const QString& title) { title_ = title; }

  void set_window(double xmin, double xmax,
		  double ymin, double ymax)
  {
    widget->set_window(xmin, xmax, ymin, ymax);
  }

 private:
  void get_object_locate_mode(CGAL::Object obj)
  {
    Rep::Point_2 q;

    if ( CGAL::assign(q, obj) ) {
      CGAL::Object o = vvd->locate(q);
      widget->redraw();
      *widget << CGAL::YELLOW;
      *widget << q;
      vvd->draw_feature(o, *widget->get_qt_widget());
    }
  }

private slots:
  void unfinished() {}

  void get_which_diagram(int id)
  {
    int selected_id = widget->get_button_group()->selectedId();
    if ( id == selected_id ) { return; }

    switch (id) {
    case 0:
      vvd = cad;
      widget->get_label()->setText("Apollonius diagram selected.");
      break;
    case 1:
      vvd = cvd;
      widget->get_label()->setText("Voronoi diagram selected.");
      break;
    default:
      vvd = cpd;
      widget->get_label()->setText("Power diagram selected.");
    }
    layers_toolbar->set_vvd(vvd);
    widget->redraw();
  }

  void get_object(CGAL::Object obj)
  {
    if ( is_locate_mode ) {
      get_object_locate_mode(obj);
      return;
    }

    if ( is_remove_mode ) {
#if 0
      // the true remove mode
      Rep::Point_2 q;

      if ( CGAL::assign(q, obj) ) {
	CGAL::Object o = vvd->locate(q);
	vvd->remove(o);
      }
      bool b = vvd->is_valid();
      if ( b ) {
	widget->get_label()->setText("Voronoi diagram is valid.");
      } else {
	widget->get_label()->setText("Voronoi diagram is NOT valid.");
      }
      widget->redraw();
      return;
#else
      // displaying the conflict region
      Rep::Point_2 q;
      Rep::Circle_2 c;
      CGAL::Object conflicts;

      widget->redraw();

      if ( CGAL::assign(q, obj) ) {
	conflicts = vvd->get_conflicts(q);
	*widget << CGAL::YELLOW;
	*widget << q;
	vvd->draw_conflicts(q, conflicts, *widget->get_qt_widget());
      } else if ( CGAL::assign(c, obj) ) {
	conflicts = vvd->get_conflicts(c);
	*widget << CGAL::YELLOW;
	*widget << c;
	vvd->draw_conflicts(c, conflicts, *widget->get_qt_widget());
      }
      return;
#endif
    }

    CGAL::Timer timer;

    if ( input_mode == VD_POINT ) {
      Rep::Point_2 p;
      if ( CGAL::assign(p, obj) ) {
	timer.start();
	vvd->insert(p);
	timer.stop();

	bool b = vvd->is_valid();
	QString msg_valid;
	if ( b ) {
	  msg_valid = "Voronoi diagram is valid.";
	} else {
	  msg_valid = "Voronoi diagram is NOT valid.";
	}

	std::sprintf(msg, "Insertion time: %f", timer.time());
	widget->get_label()->setText(msg_valid + " " + msg);
      }
    } else if ( input_mode == VD_CIRCLE ) {
      Circle_2 c;

      if( CGAL::assign(c, obj) ) {
	timer.start();
	vvd->insert(c);
	timer.stop();

	bool b = vvd->is_valid();
	QString msg_valid;
	if ( b ) {
	  msg_valid = "Voronoi diagram is valid.";
	} else {
	  msg_valid = "Voronoi diagram is NOT valid.";
	}

	std::sprintf(msg,	"Insertion time: %f", timer.time());
	widget->get_label()->setText(msg_valid + " " + msg);
      }
    }

    //    svd.is_valid(true,1);
    //    std::cout << std::endl;
    widget->redraw();
  }

  void get_locate_mode(bool b)
  {
    is_locate_mode = b;

    if ( is_locate_mode ) {
      get_point.activate();
      get_circle.deactivate();
    } else {
      if ( input_mode == VD_CIRCLE ) {
	get_point.deactivate();
	get_circle.activate();
      }
    }
  }

  void get_insert_mode(bool b)
  {
    is_remove_mode = b;
#if 0
    if ( is_remove_mode ) {
      get_point.activate();
      get_circle.deactivate();
    } else {
      if ( input_mode == VD_CIRCLE ) {
	get_point.deactivate();
	get_circle.activate();
      }
    }
#endif
  }

  void get_input_mode(Input_mode im) {
    input_mode = im;

    if ( input_mode == VD_POINT ) {
      get_point.activate();
      get_circle.deactivate();
    } else if ( input_mode == VD_CIRCLE ) {
      get_point.deactivate();
      get_circle.activate();
    }
  }

  void read_from_file(const QString& fileName)
  {
    //    typedef VVD2::Vertex_handle Vertex_handle;
    CGAL::Timer timer;
    vvd->clear();

    std::ifstream f(fileName.ascii());
    assert( f );

    int counter = 0;
    timer.start();

    Rep::Point_2 p;
    Rep::Circle_2 c;
    double r;
    if ( vvd == cad || vvd == cpd ) {
      while (f >> p >> r) {
	c = Rep::Circle_2(p,r);
	vvd->insert(c);
	counter++;

	if ( counter % 500 == 0 ) {
	  sprintf(msg, "%d sites have been inserted...", counter);
	  widget->get_label()->setText(msg);
	}
      } // endwhile
    } else {
      while (f >> p) {
	vvd->insert(p);
	counter++;

	if ( counter % 500 == 0 ) {
	  sprintf(msg, "%d sites have been inserted...", counter);
	  widget->get_label()->setText(msg);
	}
      } // endwhile
    }

    timer.stop();

    bool b = vvd->is_valid();
    QString msg_valid;
    if ( b ) {
      msg_valid = "Voronoi diagram is valid.";
    } else {
      msg_valid = "Voronoi diagram is NOT valid.";
    }

    std::sprintf(msg,
			   "%d sites inserted. Insertion time: %f",
			   counter, timer.time());
    widget->get_label()->setText(msg_valid + " " + msg);

    widget->redraw();
  }

  void write_to_file(const QString& /* fileName */)
  {
#if 0
    std::ofstream f(fileName);
    assert( f );
    f.precision(18);
    SVD_2::Input_sites_iterator sit;
    for (sit = svd.input_sites_begin();	sit != svd.input_sites_end(); ++sit) {
      f << (*sit) << std::endl;
    }
#endif
  }

  void read_input_from_file()
  {
#if 0
    QString fileName =
      QFileDialog::getOpenFileName(QString::null, QString::null,
				   this, "Open file...");

    if ( !fileName.isNull() ) {
      read_from_file(fileName);
    }
#endif
  }

  void write_output_to_file()
  {
#if 0
    QString fileName =
      QFileDialog::getSaveFileName(tr("data.out"), QString::null,
				   this, "Save as...");

    if ( !fileName.isNull() ) {
      write_to_file(fileName);
    }
#endif
  }

  void open_from_file()
  {
#if 0
    QString fileName =
      QFileDialog::getOpenFileName(QString::null, QString::null,
				   this, "Open file...");

    if ( !fileName.isNull() ) {
      std::ifstream f(fileName);
      assert(f);

      CGAL::Timer timer;
      timer.start();
      f >> svd;
      timer.stop();

      char msg[100];
      int n_sites = static_cast<int>(svd.number_of_input_sites());
      std::sprintf(msg,
			     "%d sites inserted. Insertion time: %f",
			     n_sites, timer.time());
      widget->get_label()->setText(msg);
      widget->redraw();
    }
#endif
  }

  void save_to_file()
  {
#if 0
    QString fileName =
      QFileDialog::getSaveFileName(tr("data.out"), QString::null,
				   this, "Save as...");

    if ( !fileName.isNull() ) {
      std::ofstream f(fileName);
      assert(f);
      f << svd;
    }
#endif
  }

  void print_screen()
  {
    widget->print_to_ps();
  }

  void remove_all()
    {
      //      sitelist.clear();
      num_selected = 0;
      vvd->clear();
      widget->redraw();
    }

  void about()
  {
    QMessageBox::about( this, title(),
			"This is a demo for the Voronoi diagram adaptor\n\n"
			"Author: Menelaos Karavelas "
			"<mkaravel@iacm.forth.gr>\n\n"
			"Copyright (c) CGAL 2006");
  }

  void aboutQt()
  {
    QMessageBox::aboutQt( this, title() );
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
