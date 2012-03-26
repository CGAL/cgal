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
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>

#ifndef ENVELOPE_3_H
#define ENVELOPE_3_H

#include <fstream>
#include <string>

#include <CGAL/basic.h>
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

class Qt_layer_show_diag : public CGAL::Qt_widget_layer {
public:

  // default constructor
  Qt_layer_show_diag() {}

  // override the virtual method 'draw()' of Qt_widget_layer
  void draw();
};

class Qt_layer_show_statitics : public CGAL::Qt_widget_layer
{
public:
  // default constructor
  Qt_layer_show_statitics() {};

  void mousePressEvent(QMouseEvent *);
};

/* The QMainWindow class provides a main application window,
 * with a menu bar, dock windows (e.g. for toolbars), and a status bar
 */
class MyWindow : public QMainWindow {
  Q_OBJECT

public:
  // constructor
  MyWindow(int w, int h);

private:
  void something_changed() { current_state++; };

  template<class Arrangement>
  bool open_file(Arrangement & arr)
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

    std::ifstream in_file(s.ascii());
    if (!in_file.is_open()) {
      QMessageBox::warning(widget, "Open", "Can't open file");
      return false;
    }

    QCursor old = widget->cursor();
    widget->setCursor(Qt::WaitCursor);
    widget->lock();
    widget->clear_history();

    std::list<Surface_3> surfaces;
    unsigned int num_of_surfaces = 0;
    in_file >> num_of_surfaces;
    CGAL::Random rand;
    for (unsigned int i=0 ; i<num_of_surfaces; i++) {
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

  void read_surface(std::ifstream & is, Base_triangle_3 & tri);

  void read_surface(std::ifstream & is, Base_sphere_3 & s);

  void read_surface(std::ifstream & is, Base_plane_3 & p);

public slots:
  void open_triangles_file();

  void open_spheres_file();

  void open_planes_file();

  void new_instance();

  void clear_all_diags();

  void show_v_pressed();

  void show_e_pressed();

  void show_f_pressed();

private slots:
  void about();

  void aboutQt();

  void howto();

  void new_window();

  void timer_done();

private:
  CGAL::Qt_widget * widget;
  CGAL::Qt_widget_standard_toolbar * stoolbar;
  QToolBar * show_toolbar;
  QToolButton * show_v_button;
  QToolButton * show_e_button;
  QToolButton * show_f_button;
  Qt_layer_show_diag testlayer;
  Qt_layer_show_statitics show_stat_layer;
  int old_state;
  QString curr_dir;
  int current_state;
};

#endif
