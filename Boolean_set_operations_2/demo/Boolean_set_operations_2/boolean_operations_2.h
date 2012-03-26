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
// $URL: svn+ssh://guyzucke@scm.gforge.inria.fr/svn/cgal/trunk/Boolean_set_operations_2/demo/Boolean_set_operations_2/boolean_operations_2.cpp $
// $Id: boolean_operations_2.cpp 37292 2007-03-20 07:53:53Z afabri $
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>


// if QT is not installed, a message will be issued in runtime.
#ifndef BOOLEAN_SET_OPERATIONS_2_H
#define BOOLEAN_SET_OPERATIONS_2_H

#include <CGAL/basic.h>

#include <fstream>
#include <string>

#include "typedefs.h"

#include <CGAL/IO/Qt_widget.h>
#include "Qt_widget_circ_polygon.h"
#include <CGAL/Bbox_2.h>
#include <CGAL/iterator.h>


#include <CGAL/IO/Qt_widget_Polygon_2.h>

#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/pixmaps/demoicon.xpm>
#include <CGAL/IO/pixmaps/voronoi.xpm>
#include "boolean_operations_2_toolbar.h"
#include "Qt_widget_locate_layer.h"

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

#include <CGAL/IO/Dxf_bsop_reader.h>

class MyWindow;
class Tools_toolbar;
class Qt_layer_show_ch : public CGAL::Qt_widget_layer
{
 public:

  // default constructor
//make constructor with window parameter?  
  Qt_layer_show_ch(MyWindow* myWind);

  // this method overrides the virtual method 'draw()' of Qt_widget_layer
  void draw();
   private:
   MyWindow *myWin;
};//end class


/* The QMainWindow class provides a main application window,
 *  with a menu bar, dock windows (e.g. for toolbars), and a status bar
 */
class MyWindow : public QMainWindow
{
  Q_OBJECT
public:

  // constructor
  MyWindow(int w, int h);
  //initialize the window's testlayer. Must follow the MyWindow constructor
  void init_layer();
private: 
  void something_changed();

      	
  bool write_to_file(QString /* str */);
 
  template <class PgnItr>
    void draw_result(PgnItr begin, PgnItr end);
  
public slots:
  void open_dxf_file();
  void open_linear_polygon_file();
  void new_instance();
  void radio_selected(); //const
  void perform_intersection();
  void perform_union();
  void perform_diff();
  void perform_diff2();
  void perform_symm_diff();
  void perform_red_complement();
  void perform_blue_complement();
  void make_res_red();
  void make_res_blue(); 
  void perform_mink_sum();
  void mink_sum_warning();
  bool is_linear(const Polygon_2& pgn);
  bool is_disc(const Polygon_2& pgn);
  Coord_type get_radius(const Polygon_2& pgn);
  Linear_polygon_2 circ_2_linear(const Polygon_2& pgn);
  Polygon_with_holes linear_2_circ(const Linear_polygon_with_holes_2& pgn);
  Polygon_2 linear_2_circ(const Linear_polygon_2& pgn);
  void refresh();
  void delete_red_polygons();
  void delete_blue_polygons();
  
private slots:
   void get_new_object(CGAL::Object obj);
   void about();
   void aboutQt();
   void howto();
   void new_window();
   void timer_done();   
 
private:
  CGAL::Qt_widget*                     widget;
  CGAL::Qt_widget_standard_toolbar*    stoolbar;
	//refernce type instead of object  
  Qt_layer_show_ch*                     testlayer;
  QToolBar*                            radiotoolbar;
  QRadioButton*                        red_pgn;
  QRadioButton*                        blue_pgn;
  QVButtonGroup*                       radio_group;
  Tools_toolbar*                       newtoolbar;
  QToolBar*                            bops_toolbar;
  QToolButton *                         intersection_but;
  QToolButton *                         union_but;
  QToolButton *                         diff_but;
  QToolButton *                         diff_but2;
  QToolButton *                         symm_diff_but;
  QToolButton *                         mink_sum_but;
  QToolButton *                         red_complement_but;
  QToolButton *                         blue_complement_but;
  QToolButton *                         make_res_red_but;
  QToolButton *                         make_res_blue_but;
  QToolButton *                         refresh_but;
  QToolButton *                         delete_red_but;
  QToolButton *                         delete_blue_but;

  int                                  old_state;
  QString                              file_name;
 
public: 
 	int 											current_state;
  	bool 											red_active; //init to false 
  	Polygon_set&                   		red_set;
  	Polygon_set&                    	 	blue_set;
  	Polygon_set&                    	   res_set;
};

#endif

