// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#include "envelope_3.h"
#include "icons/vertices.xpm"
#include "icons/edges.xpm"
#include "icons/faces.xpm"

// The envelope diagram
Envelope_tri_diagram_2    tri_diag;
Envelope_sphere_diagram_2 sphere_diag;
Envelope_plane_diagram_2  plane_diag;
bool draw_v;
bool draw_e;
bool draw_f;

const QString my_title_string("Envelopes of 3D surfaces");

// overrides the virtual method 'draw()' of Qt_widget_layer
void Qt_layer_show_diag::draw()
{
  widget->lock(); // widget have to be locked before drawing

  if (!tri_diag.is_empty())
    draw_arr(widget, tri_diag, draw_v, draw_e, draw_f);
  else if (!sphere_diag.is_empty())
    draw_arr(widget, sphere_diag, draw_v, draw_e, draw_f);
  else
    draw_arr(widget, plane_diag, draw_v, draw_e, draw_f);

  widget->unlock(); // widget have to be unlocked when finished drawing
}

void Qt_layer_show_statitics::mousePressEvent(QMouseEvent *)
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
  else if (!sphere_diag.is_empty()) {
    n_v = sphere_diag.number_of_vertices();
    n_e = sphere_diag.number_of_edges();
    n_f = sphere_diag.number_of_faces();
  }
  else if (!plane_diag.is_empty()) {
    n_v = plane_diag.number_of_vertices();
    n_e = plane_diag.number_of_edges();
    n_f = plane_diag.number_of_faces();
  }
  QMessageBox::information(widget, "Diagram size",
                           s.arg(n_v).arg(n_e).arg(n_f));
}

// constructor
MyWindow::MyWindow(int w, int h) :
  current_state(-1)
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
  QTimer * timer = new QTimer(this);

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
}

void MyWindow::read_surface(std::ifstream& is, Base_triangle_3& tri)
{
  is >> tri;
}

void MyWindow::read_surface(std::ifstream& is, Base_sphere_3& s)
{
  Rat_point_3 a;
  Rational sr;
  is >> a >> sr;
  s = Base_sphere_3(a, sr);
}

void MyWindow::read_surface(std::ifstream& is, Base_plane_3& p)
{
  Coord_type a, b, c, d;
  is >> a >> b >> c >> d;
  p = Base_plane_3(Plane_3(a, b, c, d));
}

void MyWindow::open_triangles_file()
{
  if (open_file(tri_diag)) {
    sphere_diag.clear();
    plane_diag.clear();
    something_changed();
  }
}

void MyWindow::open_spheres_file()
{
  if (open_file(sphere_diag)) {
    tri_diag.clear();
    plane_diag.clear();
    something_changed();
  }
}

void MyWindow::open_planes_file()
{
  if (open_file(plane_diag)) {
    tri_diag.clear();
    sphere_diag.clear();
    something_changed();
  }
}

void MyWindow::new_instance()
{
  widget->lock();
  clear_all_diags();
  widget->clear_history();
  widget->set_window(-1.1, 1.1, -1.1, 1.1);
  // set the Visible Area to the Interval
  widget->unlock();

  something_changed();
}

void MyWindow::clear_all_diags()
{
  tri_diag.clear();
  sphere_diag.clear();
  plane_diag.clear();
}

void MyWindow::show_v_pressed()
{
  draw_v = !draw_v;
  something_changed();
}

void MyWindow::show_e_pressed()
{
  draw_e = !draw_e;
  something_changed();
}

void MyWindow::show_f_pressed()
{
  draw_f = !draw_f;
  something_changed();
}

void MyWindow::about()
{
  QMessageBox::about(this, my_title_string,
                     "This is a demo for 3D envelopes of surfaces\n");
}

void MyWindow::aboutQt()
{
  QMessageBox::aboutQt(this, my_title_string);
}

void MyWindow::howto()
{
  QString home;
  home = "help/index.html";
  CGAL::Qt_help_window *help = new
    CGAL::Qt_help_window(home, ".", 0, "help viewer");
  help->resize(400, 400);
  help->setCaption("Demo HowTo");
  help->show();
}

void MyWindow::new_window()
{
  MyWindow * ed = new MyWindow(500, 500);
  ed->setCaption("Layer");
  ed->widget->clear_history();
  ed->widget->set_window(widget->x_min(), widget->x_max(),
                         widget->y_min(), widget->y_max());
  ed->show();
  something_changed();
}

void MyWindow::timer_done()
{
  if (old_state != current_state) {
    widget->redraw();
    old_state = current_state;
  }
}

#include "envelope_3.moc"

int main(int argc, char * argv[])
{
  QApplication app(argc, argv);
  MyWindow widget(600,400); // physical window size
  app.setMainWidget(&widget);
  widget.setCaption(my_title_string);
  widget.setMouseTracking(TRUE);
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  widget.setIcon(cgal_icon);
  widget.show();
  draw_v = false;
  draw_e = false;
  draw_f = true;
  return app.exec();
}

