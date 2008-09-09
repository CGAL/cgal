// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>


#include <CGAL/basic.h>


#include <fstream>
#include <string>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/pixmaps/demoicon.xpm>
#include <CGAL/Random.h>

#include <qplatinumstyle.h>
#include <qapplication.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qmenudata.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qfiledialog.h>
#include <qtimer.h>
#include <qstatusbar.h>
#include <qstring.h>
#include <qiconset.h>
#include <qpixmap.h>
#include <qpainter.h>
#include <qpushbutton.h>

#include "typedefs.h"
#include "icons/vertices.xpm"
#include "icons/edges.xpm"
#include "icons/faces.xpm"

// The envelope diagram
Envelope_tri_diagram_2    tri_diag;
Envelope_sphere_diagram_2 sphere_diag;
Envelope_plane_diagram_2  plane_diag;
int current_state;
bool draw_v;
bool draw_e;
bool draw_f;

const QString my_title_string("Envelopes of 3D surfaces");

class Qt_layer_show_diag : public CGAL::Qt_widget_layer {
public:

    // default constructor
  Qt_layer_show_diag(){};

  // this method overrides the virtual method 'draw()' of Qt_widget_layer
  void draw()
  {
    widget->lock(); // widget have to be locked before drawing

    if(!tri_diag.is_empty())
      draw_arr(widget, tri_diag, draw_v, draw_e, draw_f);
    else if(!sphere_diag.is_empty())
      draw_arr(widget, sphere_diag, draw_v, draw_e, draw_f);
    else
      draw_arr(widget, plane_diag, draw_v, draw_e, draw_f);

    widget->unlock(); // widget have to be unlocked when finished drawing
  };

};//end class

class Qt_layer_show_statitics : public CGAL::Qt_widget_layer
{
public:
    // default constructor
  Qt_layer_show_statitics(){};

  void mousePressEvent(QMouseEvent* )
  {
    QString s("|V|=%1 |E|=%2 |F|=%3");
    unsigned int n_v = 0,
                 n_e = 0,
                 n_f = 0;
    if (!tri_diag.is_empty()) {
      n_v = tri_diag.number_of_vertices();
      n_e = tri_diag.number_of_edges();
      n_f = tri_diag.number_of_faces();
    }
    else if(!sphere_diag.is_empty()) {
      n_v = sphere_diag.number_of_vertices();
      n_e = sphere_diag.number_of_edges();
      n_f = sphere_diag.number_of_faces();
    }
    else if(!plane_diag.is_empty()) {
      n_v = plane_diag.number_of_vertices();
      n_e = plane_diag.number_of_edges();
      n_f = plane_diag.number_of_faces();
    }
    QMessageBox::information(widget, "Diagram size",
                             s.arg(n_v).arg(n_e).arg(n_f));
  }
};

/* The QMainWindow class provides a main application window,
 * with a menu bar, dock windows (e.g. for toolbars), and a status bar
 */
class MyWindow : public QMainWindow {
  Q_OBJECT

public:
  // constructor
  MyWindow(int w, int h)
  {
    // Constructs a widget which is a child of this window
    widget = new CGAL::Qt_widget(this);

    /* Sets the central widget for this main window to w.
     * The central widget is surrounded by the left, top, right and bottom
     * dock areas. The menu bar is above the top dock area
     */
    setCentralWidget(widget);

    curr_dir= QString::null;

    // create a timer for checking if somthing changed
    // constructs a timer whose parent is this window
    QTimer *timer = new QTimer(this);

    // connects the timer to the window
    connect(timer, SIGNAL(timeout()), this, SLOT(timer_done()));
    timer->start(200, FALSE);         // starts the timer with a timeout

    // file menu
    QPopupMenu * file = new QPopupMenu(this);
    menuBar()->insertItem("&File", file);
    file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
    file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
    file->insertSeparator();
    file->insertItem("&Open Triangles File", this, SLOT(open_triangles_file()),
                     CTRL+Key_O);
    file->insertSeparator();
    file->insertItem("&Open Spheres File", this, SLOT(open_spheres_file()),
                     CTRL+Key_S);
    file->insertSeparator();
    file->insertItem("&Open Planes File", this, SLOT(open_planes_file()),
                     CTRL+Key_H);
    file->insertSeparator();

    file->insertSeparator();
    file->insertItem("Print", widget, SLOT(print_to_ps()), CTRL+Key_P);
    file->insertSeparator();
    file->insertItem("&Close", this, SLOT(close()), CTRL+Key_X);
    file->insertItem("&Quit", qApp, SLOT(closeAllWindows()), CTRL+Key_Q);

    // help menu
    QPopupMenu * help = new QPopupMenu(this);
    menuBar()->insertItem("&Help", help);
    help->insertItem("How To", this, SLOT(howto()), Key_F1);
    help->insertSeparator();
    help->insertItem("&About", this, SLOT(about()), CTRL+Key_A);
    help->insertItem("About &Qt", this, SLOT(aboutQt()));

    //the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar(widget, this, "ST");

    show_toolbar = new QToolBar(this, "Show features");
    QIconSet set0(QPixmap((const char**)vertices_icon),
                  QPixmap((const char**)vertices_icon));

    QIconSet set1(QPixmap((const char**)edges), QPixmap((const char**)edges));
    QIconSet set2(QPixmap((const char**)faces), QPixmap((const char**)faces));

    show_v_button = new QToolButton(show_toolbar, "Show Vertices");
    show_v_button->setToggleButton(TRUE);
    show_v_button->setTextLabel("Show Vertices ");
    connect(show_v_button,SIGNAL(pressed()),
            this, SLOT(show_v_pressed()));
    show_v_button->setIconSet(set0);


    show_toolbar->addSeparator();
    show_e_button = new QToolButton(show_toolbar, "Show Edges");
    show_e_button->setToggleButton(TRUE);
    show_e_button->setTextLabel("Show Edges ");
    connect(show_e_button,SIGNAL(pressed()), this, SLOT(show_e_pressed()));
    show_e_button->setIconSet(set1);

    show_toolbar->addSeparator();
    show_f_button = new QToolButton(show_toolbar, "Show Faces");
    show_f_button->setToggleButton(TRUE);
    show_f_button->toggle();
    show_f_button->setTextLabel("Show Faces ");
    connect(show_f_button,SIGNAL(pressed()), this, SLOT(show_f_pressed()));
    show_f_button->setIconSet(set2);

    //layers
    widget->attach(&testlayer);
    widget->attach(&show_stat_layer);
    *widget <<CGAL::BackgroundColor (CGAL::BLACK);

    resize(w,h);
    widget->set_window(-1, 1, -1, 1);
    widget->setMouseTracking(TRUE);

    //application flag stuff
    old_state = 0;
  };

private:
  void something_changed(){current_state++;};

  template<class Arrangement>
  bool open_file(Arrangement &arr)
  {
    typedef typename Arrangement::Traits_3                      Traits_3;
    typedef typename Traits_3::Surface_3                        Surface_3;
    typedef typename Traits_3::Base_traits_3:: Surface_3        Base_surface_3;
    typedef typename Arrangement::Vertex_const_iterator
        Vertex_const_iterator;
    QString s = QFileDialog::getOpenFileName(curr_dir, QString::null, this,
                                             "open file dialog",
                                             "Choose a file");
    if (s == QString::null) return false;
    curr_dir = s;

    std::ifstream in_file(s);
    if (!in_file.is_open()) {
      QMessageBox::warning(widget,"Open","Can't open file");
      return false;
    }

    QCursor old = widget->cursor();
    widget->setCursor(Qt::WaitCursor);
    widget->lock();
    widget->clear_history();

    std::list<Surface_3> surfaces;
    int num_of_surfaces;
    in_file >> num_of_surfaces;
    CGAL::Random rand;
    for (int i=0 ; i<num_of_surfaces; i++) {
      int r = rand.get_int(128, 256);
      int g = rand.get_int(0, 256);
      int b = rand.get_int(0, 256);

      Base_surface_3 s;
      read_surface(in_file, s);
      surfaces.push_back(Surface_3(s, CGAL::Color(r, g, b)));
    }
    arr.clear();
    CGAL::lower_envelope_3(surfaces.begin(), surfaces.end(), arr);

    if (arr.number_of_vertices() != 0) {
      Vertex_const_iterator vit = arr.vertices_begin();
      double x_min = CGAL::to_double(vit->point().x());
      double x_max = CGAL::to_double(vit->point().x());
      double y_min = CGAL::to_double(vit->point().y());
      double y_max = CGAL::to_double(vit->point().y());

      for (++vit; vit != arr.vertices_end(); ++vit) {
        double curr_x = CGAL::to_double(vit->point().x());
        double curr_y = CGAL::to_double(vit->point().y());

        if (curr_x < x_min)
          x_min = curr_x;
        else if (curr_x > x_max)
          x_max = curr_x;

        if (curr_y < y_min)
          y_min = curr_y;
        else if (curr_y > y_max)
          y_max = curr_y;
      }
      double w = (x_max - x_min)/10;
      double h = (y_max - y_min)/10;

      // make sure the bbox is not degenerated
      if (w == 0.0)
        w+=10.0;
      if (h == 0.0)
        h+=10.0;
      widget->set_window(x_min - w, x_max + w, y_min - h, y_max + h);
    }
    widget->unlock();
    widget->setCursor(old);
    return true;
  }

  void read_surface(std::ifstream& is, Base_triangle_3& tri)
  {
    is >> tri;
  }

  void read_surface(std::ifstream& is, Base_sphere_3& s)
  {
    Rat_point_3 a;
    Rational sr;
    is >> a >> sr;
    s = Base_sphere_3(a, sr);
  }

  void read_surface(std::ifstream& is, Base_plane_3& p)
  {
    Coord_type a, b, c, d;
    is >> a >> b >> c >> d;
    p = Base_plane_3(Plane_3(a, b, c, d));
  }

public slots:

  void open_triangles_file()
  {
    if (open_file(tri_diag)) {
      sphere_diag.clear();
      plane_diag.clear();
      something_changed();
    }
  }

  void open_spheres_file()
  {
    if (open_file(sphere_diag)) {
      tri_diag.clear();
      plane_diag.clear();
      something_changed();
    }
  }

  void open_planes_file()
  {
    if (open_file(plane_diag)) {
      tri_diag.clear();
      sphere_diag.clear();
      something_changed();
    }
  }

  void new_instance()
  {
    widget->lock();
    clear_all_diags();
    widget->clear_history();
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
        // set the Visible Area to the Interval
    widget->unlock();

    something_changed();
  }

  void clear_all_diags()
  {
    tri_diag.clear();
    sphere_diag.clear();
    plane_diag.clear();
  }

  void show_v_pressed()
  {
    draw_v = !draw_v;
    something_changed();
  }

  void show_e_pressed()
  {
    draw_e = !draw_e;
    something_changed();
  }

  void show_f_pressed()
  {
    draw_f = !draw_f;
    something_changed();
  }

private slots:


  void about()
  {
    QMessageBox::about(this, my_title_string,
        "This is a demo for 3D envelopes of surfaces\n");

  }

  void aboutQt()
  {
    QMessageBox::aboutQt(this, my_title_string);
  }

  void howto()
  {
    QString home;
    home = "help/index.html";
    CGAL::Qt_help_window *help = new
      CGAL::Qt_help_window(home, ".", 0, "help viewer");
    help->resize(400, 400);
    help->setCaption("Demo HowTo");
    help->show();
  }

  void new_window(){
    MyWindow *ed = new MyWindow(500, 500);
    ed->setCaption("Layer");
    ed->widget->clear_history();
    ed->widget->set_window(widget->x_min(), widget->x_max(),
                           widget->y_min(), widget->y_max());
    ed->show();
    something_changed();
  }

  void timer_done()
  {
    if(old_state!=current_state)
    {
      widget->redraw();
      old_state = current_state;
    }
  }

private:
  CGAL::Qt_widget*                     widget;
  CGAL::Qt_widget_standard_toolbar*    stoolbar;
  QToolBar*                            show_toolbar;
  QToolButton*                         show_v_button;
  QToolButton*                         show_e_button;
  QToolButton*                         show_f_button;
  Qt_layer_show_diag                   testlayer;
  Qt_layer_show_statitics              show_stat_layer;
  int                                  old_state;
  QString                              curr_dir;
};

#include "envelope_3.moc"

int main(int argc, char **argv)
{
  QApplication app(argc, argv);
  MyWindow widget(600,400); // physical window size
  app.setMainWidget(&widget);
  widget.setCaption(my_title_string);
  widget.setMouseTracking(TRUE);
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  widget.setIcon(cgal_icon);
  widget.show();
  current_state = -1;
  draw_v = false;
  draw_e = false;
  draw_f = true;
  return app.exec();
}

