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



#include <CGAL/basic.h>
 

#include <fstream>
#include <string>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/iterator.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/pixmaps/demoicon.xpm>
#include <CGAL/IO/pixmaps/voronoi.xpm>
 
#include "boolean_operations_2.h"
#include "typedefs.h"
#include "boolean_operations_2_toolbar.h"
#include "Qt_widget_circ_polygon.h"

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

 
#include "Qt_widget_X_monotone_circle_segment_2.h"
#include "icons/union.xpm"
#include "icons/intersection.xpm"
#include "icons/diff_PQ.xpm"
#include "icons/diff_QP.xpm"
#include "icons/symm_diff.xpm"
#include "icons/make_P.xpm"
#include "icons/make_Q.xpm"
#include "icons/comp_P.xpm"
#include "icons/comp_Q.xpm" 
#include "icons/del_P.xpm"
#include "icons/del_Q.xpm"
#include "icons/refresh.xpm"
#include "icons/mink_sum.xpm"

#include <CGAL/IO/Dxf_bsop_reader.h>
#include <qdeepcopy.h>
//#include<string>

const QString my_title_string("Boolean operations on polygons ");

//global variable to aid naming windows 
int winsOpened=2;

Qt_layer_show_ch::Qt_layer_show_ch(MyWindow* myWind) {
	myWin=myWind;
}
 
// this method overrides the virtual method 'draw()' of Qt_widget_layer
void Qt_layer_show_ch::draw()
{
  widget->lock(); // widget has to be locked before drawing
  RasterOp old_rasterop = widget->rasterOp();
  widget->get_painter().setRasterOp(XorROP);
  widget->setFilled (true);
  widget->setFillColor (CGAL::RED);  
  *widget <<  CGAL::BLACK; 
  std::list<Polygon_with_holes> red_pgns_list;
  myWin->red_set.polygons_with_holes(std::back_inserter(red_pgns_list));
  std::list<Polygon_with_holes>::iterator itpgn1 = red_pgns_list.begin();

  while (itpgn1 != red_pgns_list.end())
    {
      const Polygon_with_holes& pgn_with_hole = *itpgn1;
      const Polygon_2& outer_boundary = pgn_with_hole.outer_boundary();
      if (outer_boundary.is_empty()) 
	{
	  // no boundary -> unbounded polygon
	  Iso_rectangle rect(Point_2(widget->x_min(), widget->y_min()),
			     Point_2(widget->x_max(), widget->y_max()));
	  *widget << rect;
	}
      else
	*widget << outer_boundary;

      for(Hole_const_iterator hit = pgn_with_hole.holes_begin();
          hit != pgn_with_hole.holes_end();
          ++hit)
	{
	  *widget << *hit;
	}
      ++itpgn1;
    }


  widget->setFilled (true);
  widget->setFillColor (CGAL::BLUE);

  std::list<Polygon_with_holes> blue_pgns_list;
  myWin->blue_set.polygons_with_holes(std::back_inserter(blue_pgns_list));
  std::list<Polygon_with_holes>::iterator itpgn2 = blue_pgns_list.begin();

  while (itpgn2 != blue_pgns_list.end())
    {
      const Polygon_with_holes& pgn_with_hole = *itpgn2;
      const Polygon_2& outer_boundary = pgn_with_hole.outer_boundary();
      if (outer_boundary.is_empty())
	{
	  // no boundary -> unbounded polygon
	  Iso_rectangle rect(Point_2(widget->x_min(), widget->y_min()),
			     Point_2(widget->x_max(), widget->y_max()));
	  *widget << rect;
	}
      else
	{
	  *widget << outer_boundary;
	}

      for(Hole_const_iterator hit = pgn_with_hole.holes_begin();
          hit != pgn_with_hole.holes_end();
          ++hit)
	{
	  *widget << *hit;
	}
      ++itpgn2;
    }
  widget->get_painter().setRasterOp(old_rasterop);
  widget->setFilled (false);
  widget->unlock(); // widget have to be unlocked when finished drawing
}

/* The QMainWindow class provides a main application window,
 *  with a menu bar, dock windows (e.g. for toolbars), and a status bar
 */

// constructor
MyWindow::MyWindow(int w, int h) : 
  red_active(false),
  red_set(*(new Polygon_set)),
  blue_set(*(new Polygon_set)),
  res_set(*(new Polygon_set))										  
{
  widget = new CGAL::Qt_widget(this); //Constructs a widget which is a child of this window


  /* Sets the central widget for this main window to w.
   * The central widget is surrounded by the left, top, right and bottom dock areas.
   * The menu bar is above the top dock area
   */
  setCentralWidget(widget);

  file_name= QString::null;

  //create a timer for checking if somthing changed
  QTimer *timer = new QTimer( this ); // constructs a timer whose parent is this window

  connect( timer, SIGNAL(timeout()),
           this, SLOT(timer_done()) );  // connects the timer to the window
  timer->start( 200, FALSE ); // Starts the timer with a msec milliseconds timeout

  // file menu
  QPopupMenu * file = new QPopupMenu( this );
  menuBar()->insertItem( "&File", file );
  file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
  file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
  file->insertSeparator();
  file->insertItem("&Open Linear Polygon file", this, SLOT(open_linear_polygon_file()),CTRL+Key_O);
  file->insertItem("&Open DXF file", this, SLOT(open_dxf_file()),CTRL+Key_D);
  file->insertSeparator();
  //file->insertItem("&Save",this ,SLOT(save_file()),CTRL+Key_S);
  //file->insertItem("&Save as",this ,SLOT(save_file_as()));
  file->insertSeparator();
  file->insertItem("Print", widget, SLOT(print_to_ps()), CTRL+Key_P);
  file->insertSeparator();
  file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
  file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), CTRL+Key_Q );

  // help menu
  QPopupMenu * help = new QPopupMenu( this );
  menuBar()->insertItem( "&Help", help );
  help->insertItem("How To", this, SLOT(howto()), Key_F1);
  help->insertSeparator();
  help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
  help->insertItem("About &Qt", this, SLOT(aboutQt()) );

  //the standard toolbar
  stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this, "ST");

  radiotoolbar = new QToolBar(this, "polygon type");
  blue_pgn = new QRadioButton ("Blue", radiotoolbar);
  blue_pgn->toggle();
  red_pgn = new QRadioButton("Red", radiotoolbar);
  radio_group = new QVButtonGroup(this,"Radios");
  radio_group->insert(blue_pgn);
  radio_group->insert(red_pgn);
  radio_group->setRadioButtonExclusive(true);


  connect(blue_pgn, SIGNAL(toggled (bool)),
	  this, SLOT(radio_selected()));
  connect(red_pgn, SIGNAL(toggled (bool)),
	  this, SLOT(radio_selected()));


  //layers
  //widget->attach(&testlayer);

  //the new tools toolbar
  newtoolbar = new Tools_toolbar(widget, this);

  // voronoi toolbar
  bops_toolbar = new QToolBar(this, "Boolean operations");

  QIconSet set0(QPixmap( (const char**)intersection_xpm ),
		QPixmap( (const char**)intersection_xpm ));

  intersection_but = new QToolButton(bops_toolbar, "Boolean operations");
  intersection_but->setAutoRaise(TRUE);

  intersection_but->setIconSet(set0);
  intersection_but->setTextLabel("Intersection ");
  connect(intersection_but,SIGNAL(pressed()),
	  this, SLOT(perform_intersection()));

  QIconSet set1(QPixmap( (const char**)union_xpm ),
		QPixmap( (const char**)union_xpm ));

  bops_toolbar->addSeparator();
  union_but = new QToolButton(bops_toolbar, "Boolean operations");
  union_but->setAutoRaise(TRUE);

  union_but->setIconSet(set1);
  union_but->setTextLabel("Union ");
  connect(union_but,SIGNAL(pressed()),
	  this, SLOT(perform_union()));

  QIconSet set2(QPixmap( (const char**)diff_PQ_xpm ),
		QPixmap( (const char**)diff_PQ_xpm ));

  bops_toolbar->addSeparator();
  diff_but2 = new QToolButton(bops_toolbar, "Boolean operations");
  diff_but2->setAutoRaise(TRUE);

  diff_but2->setIconSet(set2);
  diff_but2->setTextLabel("Difference between Blue and Red");
  connect(diff_but2, SIGNAL(pressed()),
	  this, SLOT(perform_diff2()));

  QIconSet set3(QPixmap( (const char**)diff_QP_xpm ),
		QPixmap( (const char**)diff_QP_xpm ));

  bops_toolbar->addSeparator();
  diff_but = new QToolButton(bops_toolbar, "Boolean operations");
  diff_but->setAutoRaise(TRUE);

  diff_but->setIconSet(set3);
  diff_but->setTextLabel("Difference between Red and Blue");
  connect(diff_but, SIGNAL(pressed()),
	  this, SLOT(perform_diff()));

  QIconSet set4(QPixmap( (const char**)symm_diff_xpm ),
		QPixmap( (const char**)symm_diff_xpm ));
  bops_toolbar->addSeparator();

  symm_diff_but = new QToolButton(bops_toolbar, "Boolean operations");
  symm_diff_but->setAutoRaise(TRUE);

  symm_diff_but->setIconSet(set4);
  symm_diff_but->setTextLabel("Symmetric Difference ");
  connect(symm_diff_but, SIGNAL(pressed()),
	  this, SLOT(perform_symm_diff()));

  QIconSet set12(QPixmap( (const char**)mink_sum_xpm ),
		 QPixmap( (const char**)mink_sum_xpm ));
  bops_toolbar->addSeparator();
  mink_sum_but = new QToolButton(bops_toolbar, "Boolean operations");
  mink_sum_but->setAutoRaise(TRUE);
  mink_sum_but->setIconSet(set12);
  mink_sum_but->setTextLabel("Minkowski Sum ");
  connect(mink_sum_but, SIGNAL(pressed()),
	  this, SLOT(perform_mink_sum()));

  QIconSet set5(QPixmap( (const char**)comp_P_xpm ),
		QPixmap( (const char**)comp_P_xpm ));
  bops_toolbar->addSeparator();

  blue_complement_but = new QToolButton(bops_toolbar, "Boolean operations");
  blue_complement_but->setAutoRaise(TRUE);

  blue_complement_but->setIconSet(set5);
  blue_complement_but->setTextLabel("Blue Complement ");
  connect(blue_complement_but, SIGNAL(pressed()),
	  this, SLOT(perform_blue_complement()));

  QIconSet set6(QPixmap( (const char**)comp_Q_xpm ),
		QPixmap( (const char**)comp_Q_xpm ));
  bops_toolbar->addSeparator();

  red_complement_but = new QToolButton(bops_toolbar, "Boolean operations");
  red_complement_but->setAutoRaise(TRUE);

  red_complement_but->setIconSet(set6);
  red_complement_but->setTextLabel("Red Complement ");
  connect(red_complement_but, SIGNAL(pressed()),
	  this, SLOT(perform_red_complement()));


  QIconSet set7(QPixmap( (const char**)make_P_xpm ),
		QPixmap( (const char**)make_P_xpm ));
  bops_toolbar->addSeparator();
  make_res_blue_but = new QToolButton(bops_toolbar, "Boolean operations");
  make_res_blue_but->setAutoRaise(TRUE);


  make_res_blue_but->setIconSet(set7);
  make_res_blue_but->setTextLabel("Make Result Blue");
  connect(make_res_blue_but,SIGNAL(pressed()),
	  this, SLOT(make_res_blue()));

  QIconSet set8(QPixmap( (const char**)make_Q_xpm ),
		QPixmap( (const char**)make_Q_xpm ));
  bops_toolbar->addSeparator();
  make_res_red_but = new QToolButton(bops_toolbar, "Boolean operations");
  make_res_red_but->setAutoRaise(TRUE);


  make_res_red_but->setIconSet(set8);
  make_res_red_but->setTextLabel("Make Result Red");
  connect(make_res_red_but,SIGNAL(pressed()),
	  this, SLOT(make_res_red()));

  QIconSet set9(QPixmap( (const char**)refresh_xpm ),
		QPixmap( (const char**)refresh_xpm ));
  bops_toolbar->addSeparator();

  refresh_but = new QToolButton(bops_toolbar, "Boolean operations");
  refresh_but->setAutoRaise(TRUE);

  refresh_but->setIconSet(set9);
  refresh_but->setTextLabel("Refresh ");
  connect(refresh_but,SIGNAL(pressed()),
	  this, SLOT(refresh()));

  QIconSet set10(QPixmap( (const char**)del_P_xpm ),
		 QPixmap( (const char**)del_P_xpm ));
  bops_toolbar->addSeparator();

  delete_blue_but = new QToolButton(bops_toolbar, "Boolean operations");
  delete_blue_but->setAutoRaise(TRUE);

  delete_blue_but->setIconSet(set10);
  delete_blue_but->setTextLabel("Delete Blue Polygons");
  connect(delete_blue_but,SIGNAL(pressed()),
	  this, SLOT(delete_blue_polygons()));


  QIconSet set11(QPixmap( (const char**)del_Q_xpm ),
		 QPixmap( (const char**)del_Q_xpm ));
  bops_toolbar->addSeparator();

  delete_red_but = new QToolButton(bops_toolbar, "Boolean operations");
  delete_red_but->setAutoRaise(TRUE);

  delete_red_but->setIconSet(set11);
  delete_red_but->setTextLabel("Delete Red Polygons");
  connect(delete_red_but,SIGNAL(pressed()),
	  this, SLOT(delete_red_polygons()));




  *widget << CGAL::LineWidth(2) << CGAL::BackgroundColor (CGAL::BLACK);

  resize(w,h);
  widget->set_window(-1, 1, -1, 1);
  widget->setMouseTracking(TRUE);

  //connect the widget to the main function that receives the objects
  connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),
	  this, SLOT(get_new_object(CGAL::Object)));

  //application flag stuff
  old_state = 0;
  current_state = 1;
  red_active = false;
  red_set.clear();
  blue_set.clear();
  res_set.clear();
}


void MyWindow::something_changed(){current_state++;}
/*not necessary if rado toggle is not const
  void change_grp_color(int color) {
  if (color == 1)
  red_active=true;
  else
  red_active=false;   
  }*/ 


   
bool MyWindow::write_to_file(QString /* str */)
{
  /*std::ofstream out_file(str);
    if (!out_file.is_open())
    return false;

    std::list<Polygon_with_holes> red_pgns_list;
    red_set.polygons_with_holes(std::back_inserter(red_pgns_list));
    int num_of_red_pgns = red_pgns_list.size();
    out_file << num_of_red_pgns << std::endl;
    std::list<Polygon_with_holes>::iterator red_itr;
    for(red_itr = red_pgns_list.begin() ; red_itr != red_pgns_list.end() ; ++red_itr)
    {
    const Polygon_with_holes& pgn = *red_itr;
    out_file<< pgn << std::endl;
    }

    std::list<Polygon_with_holes> blue_pgns_list;
    int num_of_blue_pgns = blue_pgns_list.size();
    out_file << num_of_blue_pgns << std::endl;
    std::list<Polygon_with_holes>::iterator blue_itr;
    for(blue_itr = blue_pgns_list.begin() ; blue_itr != blue_pgns_list.end() ; ++blue_itr)
    {
    const Polygon_with_holes& pgn = *blue_itr;
    out_file<< pgn << std::endl;
    }
    out_file.close();
    return true;*/
  return true;
}

template <class PgnItr>
void MyWindow::draw_result(PgnItr begin, PgnItr end)
{
  widget->lock();
  widget->clear();
  widget->setFilled (true);
  widget->setColor(CGAL::GREEN);

  for(PgnItr itr = begin; itr != end; ++itr)
    {
      const Polygon_with_holes& pgn_with_hole = *itr;
      const Polygon_2& outer_boundary = pgn_with_hole.outer_boundary();
      widget->setFillColor (CGAL::ORANGE);
      if (outer_boundary.is_empty())
	{
	  // no boundary -> unbounded polygon
	  Iso_rectangle rect(Point_2(widget->x_min(), widget->y_min()),
			     Point_2(widget->x_max(), widget->y_max()));
	  *widget << rect;
	}
      else
	*widget << outer_boundary;

      widget->setFillColor (CGAL::BLACK) ;
      for(Hole_const_iterator hit = pgn_with_hole.holes_begin();
	  hit != pgn_with_hole.holes_end();
	  ++hit)
	{
	  *widget << *hit;
	}
    }

  widget->setFilled (false);
  widget->unlock();
}

void MyWindow::open_dxf_file()
{

  QString s = QFileDialog::getOpenFileName(file_name,
					   QString::null,
					   this,
					   "open file dialog",
					   "Choose a file" );
  if (s==QString::null)
    return;
  file_name=s;

  std::ifstream in_file(s.ascii());
  if (!in_file.is_open())
    {
      QMessageBox::warning( widget,"Open","Can't open file");
      return ;
    }
  bool are_simple;
  int answer = 0;
  answer =
    QMessageBox::question(this, QString("Open file"),
			  QString("Are all polygons simple and without holes?"),
			  QString("Yes"), QString("No"), QString::null, 0 , 0);
  if (answer == 0)
    are_simple = true;
  else
    are_simple = false;

  CGAL::Bbox_2 box = CGAL::Bbox_2 (widget->x_min(), widget->y_min(),
				   widget->x_max(), widget->y_max());
  QCursor old = widget->cursor();
  widget->setCursor(Qt::WaitCursor);
  widget->lock();
  widget->clear_history();

  CGAL::Dxf_bsop_reader<Kernel> reader;
  std::vector<Polygon_2>  circ_polygons;
  std::vector<Polygon_with_holes>  circ_polygons_with_holes;
  reader(in_file,
	 std::back_inserter(circ_polygons),
	 std::back_inserter(circ_polygons_with_holes),
	 !are_simple);

  if (red_active)
    {
      red_set.join(circ_polygons.begin(),
		   circ_polygons.end(),
		   circ_polygons_with_holes.begin(),
		   circ_polygons_with_holes.end());
    }
  else
    {
      blue_set.join(circ_polygons.begin(),
		    circ_polygons.end(),
		    circ_polygons_with_holes.begin(),
		    circ_polygons_with_holes.end());
    }

  if (!circ_polygons.empty())
    {
      box = circ_polygons.front().bbox();
    }
  else if (!circ_polygons_with_holes.empty())
    {
      box = circ_polygons_with_holes.front().outer_boundary().bbox();
    }

  std::vector<Polygon_2>::iterator itr1 = circ_polygons.begin();
  for(itr1 = circ_polygons.begin();
      itr1 != circ_polygons.end();
      ++itr1)
    {
      box = box + itr1->bbox();
    }
  std::vector<Polygon_with_holes>::iterator itr2;

  for (itr2 = circ_polygons_with_holes.begin();
       itr2 != circ_polygons_with_holes.end();
       ++itr2)
    {
      box = box + itr2->outer_boundary().bbox();
    }

  widget->set_window(box.xmin(),
		     box.xmax(),
		     box.ymin(),
		     box.ymax());
  widget->unlock();
  newtoolbar->reset();
  something_changed();
  widget->setCursor(old);
}

void MyWindow::open_linear_polygon_file()
{

  QString s = QFileDialog::getOpenFileName(file_name,
					   QString::null,
					   this,
					   "open file dialog",
					   "Choose a file" );
  if (s==QString::null)
    return;
  file_name=s;

  std::ifstream in_file(s.ascii());
  if (!in_file.is_open())
    {
      QMessageBox::warning( widget,"Open","Can't open file");
      return ;
    }

  CGAL::Bbox_2 box = CGAL::Bbox_2 (widget->x_min(), widget->y_min(),
				   widget->x_max(), widget->y_max());
  QCursor old = widget->cursor();
  widget->setCursor(Qt::WaitCursor);
  widget->lock();
  widget->clear_history();

  Linear_polygon_2 pgn;
  in_file >> pgn;
  if (pgn.is_empty())
    {
      widget->unlock();
      widget->setCursor(old);
      return;
    }
  if (pgn.orientation() != CGAL::COUNTERCLOCKWISE)
    pgn.reverse_orientation();

  const Polygon_2& circ_pgn = linear_2_circ(pgn);


  if (red_active)
    red_set.join(circ_pgn);
  else
    blue_set.join(circ_pgn);


  box = box + circ_pgn.bbox();
  widget->set_window(box.xmin(),
		     box.xmax(),
		     box.ymin(),
		     box.ymax());
  widget->unlock();
  newtoolbar->reset();
  something_changed();
  widget->setCursor(old);
}


/*void save_file()
  {
  if (file_name==QString::null)
  save_file_as();
  else
  {
  if (!write_to_file(file_name))
  QMessageBox::warning( widget,"Save","Can't write to file");
  }
  }*/


/*void save_file_as()
  {
  QString s = QFileDialog::getSaveFileName(
  "./",
  QString::null,
  this,
  "save file dialog",
  "Choose a filename to save under" );
  if (s==QString::null)
  return;
  if (!write_to_file(s))
  {
  QMessageBox::warning( widget,"Save","Can't write to file");
  return ;
  }
  file_name=s;
  }*/



void MyWindow::new_instance()
{
  newtoolbar->deactivate();
  widget->lock();
  file_name = QString::null;

  red_set.clear();
  blue_set.clear();
  widget->clear_history();
  widget->set_window(-1.1, 1.1, -1.1, 1.1);
  // set the Visible Area to the Interval
  widget->unlock();

  something_changed();
}

void MyWindow::radio_selected() //const
{
  if (red_pgn->isOn())
    red_active = true;
  else
    if (blue_pgn->isOn())
      red_active = false;
}



void MyWindow::perform_intersection()
{
  QCursor old = widget->cursor();
  widget->setCursor(Qt::WaitCursor);
  std::list<Polygon_with_holes> res_pgns;
  res_set = red_set;
  res_set.intersection(blue_set);
  res_set.polygons_with_holes(std::back_inserter(res_pgns));
  newtoolbar->reset();
  draw_result(res_pgns.begin(), res_pgns.end());
  widget->setCursor(old);
}

void MyWindow::perform_union()
{
  QCursor old = widget->cursor();
  widget->setCursor(Qt::WaitCursor);
  std::list<Polygon_with_holes> res_pgns;
  res_set = red_set;
  res_set.join(blue_set);
  res_set.polygons_with_holes(std::back_inserter(res_pgns));
  newtoolbar->reset();
  draw_result(res_pgns.begin(), res_pgns.end());
  widget->setCursor(old);
}

void MyWindow::perform_diff()
{
  QCursor old = widget->cursor();
  widget->setCursor(Qt::WaitCursor);
  std::list<Polygon_with_holes> res_pgns;
  res_set = red_set;
  res_set.difference(blue_set);
  res_set.polygons_with_holes(std::back_inserter(res_pgns));
  newtoolbar->reset();
  draw_result(res_pgns.begin(), res_pgns.end());
  widget->setCursor(old);
}

void MyWindow::perform_diff2()
{
  QCursor old = widget->cursor();
  widget->setCursor(Qt::WaitCursor);
  std::list<Polygon_with_holes> res_pgns;
  res_set = blue_set;
  res_set.difference(red_set);
  res_set.polygons_with_holes(std::back_inserter(res_pgns));
  newtoolbar->reset();

  draw_result(res_pgns.begin(), res_pgns.end());
  widget->setCursor(old);
}

void MyWindow::perform_symm_diff()
{
  QCursor old = widget->cursor();
  widget->setCursor(Qt::WaitCursor);
  std::list<Polygon_with_holes> res_pgns;
  res_set = red_set;
  res_set.symmetric_difference(blue_set);
  res_set.polygons_with_holes(std::back_inserter(res_pgns));
  newtoolbar->reset();

  draw_result(res_pgns.begin(), res_pgns.end());
  widget->setCursor(old);
}

void MyWindow::perform_red_complement()
{
  QCursor old = widget->cursor();
  widget->setCursor(Qt::WaitCursor);
  std::list<Polygon_with_holes> res_pgns;
  res_set = red_set;
  res_set.complement();
  res_set.polygons_with_holes(std::back_inserter(res_pgns));
  newtoolbar->reset();

  draw_result(res_pgns.begin(), res_pgns.end());
  widget->setCursor(old);
}

void MyWindow::perform_blue_complement()
{
  QCursor old = widget->cursor();
  widget->setCursor(Qt::WaitCursor);
  std::list<Polygon_with_holes> res_pgns;
  res_set = blue_set;
  res_set.complement();
  res_set.polygons_with_holes(std::back_inserter(res_pgns));
  newtoolbar->reset();

  draw_result(res_pgns.begin(), res_pgns.end());
  widget->setCursor(old);
}

void MyWindow::make_res_red()
{
  if (res_set.is_empty())
    {
      int answer = 0;
      answer = QMessageBox::warning(this, "Store result",
                                    QString( "Result is empty, all polygons will be deleted\n continue anyway?\n" ),
                                    "&Yes", "&No", QString::null, 1, 1 );
      if (answer == 1)
	{
	  // answer in 'no'
	  make_res_red_but->setDown(FALSE);
	  return;
	}
      make_res_red_but->setDown(FALSE);
    }

  red_set = res_set;
  res_set.clear();
  blue_set.clear();
  newtoolbar->reset();
  something_changed();
}

void MyWindow::make_res_blue()
{
  if (res_set.is_empty())
    {
      int answer = 0;
      answer = QMessageBox::warning(this, "Store result",
                                    QString( "Result is empty, all polygons will be deleted\n continue anyway?\n" ),
                                    "&Yes", "&No", QString::null, 1, 1 );
      if (answer == 1)
	{
	  // answer in 'no'
	  make_res_blue_but->setDown(FALSE);
	  return;
	}
      make_res_blue_but->setDown(FALSE);
    }
  blue_set = res_set;
  res_set.clear();
  red_set.clear();
  newtoolbar->reset();
  something_changed();
}

void MyWindow::perform_mink_sum()
{
  if (red_set.number_of_polygons_with_holes() > 1 ||
     blue_set.number_of_polygons_with_holes() > 1)
    {
      mink_sum_warning();
      return;
    }

  Polygon_with_holes red_p_wh;
  CGAL::Oneset_iterator<Polygon_with_holes> oi1(red_p_wh);
  red_set.polygons_with_holes(oi1);
  if (red_p_wh.has_holes() || red_p_wh.is_unbounded())
    {
      mink_sum_warning();
      return;
    }
  const Polygon_2& red_p = red_p_wh.outer_boundary();

  Polygon_with_holes blue_p_wh;
  CGAL::Oneset_iterator<Polygon_with_holes> oi2(blue_p_wh);
  blue_set.polygons_with_holes(oi2);
  if (blue_p_wh.has_holes() || blue_p_wh.is_unbounded())
    {
      mink_sum_warning();
      return;
    }
  QCursor old = widget->cursor();
  widget->setCursor(Qt::WaitCursor);
  const Polygon_2& blue_p = blue_p_wh.outer_boundary();
  if (is_linear(red_p) && is_linear(blue_p))
    {
      const Linear_polygon_2& linear_red_p  = circ_2_linear(red_p);
      const Linear_polygon_2& linear_blue_p = circ_2_linear(blue_p);

      const Linear_polygon_with_holes_2& res_p =
        CGAL::minkowski_sum_2(linear_red_p, linear_blue_p);

      Polygon_with_holes res_p_wh  = linear_2_circ(res_p);
      res_set.clear();
      res_set.insert(res_p_wh);
      newtoolbar->reset();
      std::list<Polygon_with_holes> res_pgns;
      res_pgns.push_back(res_p_wh);
      draw_result(res_pgns.begin(), res_pgns.end());
      widget->setCursor(old);
    }
  else if (is_disc(red_p) && is_linear(blue_p))
    {
      const Linear_polygon_2& linear_blue_p = circ_2_linear(blue_p);
      const Polygon_with_holes& res_p_wh =
        CGAL::approximated_offset_2 (linear_blue_p, get_radius(red_p), 0.00001);

      res_set.clear();
      res_set.insert(res_p_wh);
      newtoolbar->reset();
      std::list<Polygon_with_holes> res_pgns;
      res_pgns.push_back(res_p_wh);
      draw_result(res_pgns.begin(), res_pgns.end());
      widget->setCursor(old);
    }
  else if (is_disc(blue_p) && is_linear(red_p))
    {
      const Linear_polygon_2& linear_red_p = circ_2_linear(red_p);
      const Polygon_with_holes& res_p_wh =
        CGAL::approximated_offset_2 (linear_red_p, get_radius(blue_p), 0.00001);

      res_set.clear();
      res_set.insert(res_p_wh);
      newtoolbar->reset();
      std::list<Polygon_with_holes> res_pgns;
      res_pgns.push_back(res_p_wh);
      draw_result(res_pgns.begin(), res_pgns.end());
      widget->setCursor(old);
    } 
  else
    {
      mink_sum_warning();
      widget->setCursor(old);
      return;
    }
}

void MyWindow::mink_sum_warning()
{
  QMessageBox::warning(this,
		       "Minkowski Sum",
		       QString( "Minkowski sum can be performed on two linear polygons without holes\n\
 or on a linear polygon without holes and a disc\n" ),
		       "&Ok");
  mink_sum_but->setDown(FALSE);
}

bool MyWindow::is_linear(const Polygon_2& pgn)
{
  typedef Polygon_2::Curve_const_iterator    Curve_const_iterator;
  for(Curve_const_iterator i = pgn.curves_begin();
      i != pgn.curves_end();
      ++i)
    {
      if (i->is_circular())
        return false;
    }
  return true;
}

bool MyWindow::is_disc(const Polygon_2& pgn)
{
  if (pgn.size() != 2)
    return false;

  Polygon_2::Curve_const_iterator ci = pgn.curves_begin();

  if (!ci->is_circular())
    return false;

  const Circle& c1 = ci->supporting_circle();
  ++ci;
  if (!ci->is_circular())
    return false;

  const Circle& c2 = ci->supporting_circle();

  return ((c1.center() == c2.center()) &&
	  (c1.squared_radius() == c2.squared_radius()));
}

Coord_type MyWindow::get_radius(const Polygon_2& pgn)
{
  CGAL_assertion(is_disc(pgn));
  double r =
    CGAL::sqrt(CGAL::to_double(pgn.curves_begin()->supporting_circle().squared_radius()));
  return (Coord_type(r));
}

Linear_polygon_2 MyWindow::circ_2_linear(const Polygon_2& pgn)
{
  Linear_polygon_2 linear_pgn;
  typedef Polygon_2::Curve_const_iterator    Curve_const_iterator;
  for(Curve_const_iterator i = pgn.curves_begin();
      i != pgn.curves_end();
      ++i)
    {
      const Circular_point_2& circ_p = i->source();
      linear_pgn.push_back(Point_2(circ_p.x().alpha(),
                                   circ_p.y().alpha()));
    }

  return (linear_pgn);
}

Polygon_with_holes MyWindow::linear_2_circ(const Linear_polygon_with_holes_2& pgn)
{
  Polygon_with_holes p(linear_2_circ(pgn.outer_boundary()));

  typedef Linear_polygon_with_holes_2::Hole_const_iterator
    Hole_const_iterator;
  for(Hole_const_iterator hi = pgn.holes_begin();
      hi != pgn.holes_end();
      ++hi)
    {
      p.add_hole(linear_2_circ(*hi));
    }

  return (p);
}

Polygon_2 MyWindow::linear_2_circ(const Linear_polygon_2& pgn)
{
  Polygon_2 p;
  typedef Linear_polygon_2::Edge_const_iterator    Edge_const_iterator;
  for(Edge_const_iterator ei = pgn.edges_begin();
      ei != pgn.edges_end();
      ++ei)
    {
      XCurve cv(ei->source(), ei->target());
      p.push_back(cv);
    }

  return (p);
}

void MyWindow::refresh()
{
  newtoolbar->reset();
  something_changed();
}

void MyWindow::delete_red_polygons()
{
  red_set.clear();
  newtoolbar->reset();
  something_changed();
}

void MyWindow::delete_blue_polygons()
{
  blue_set.clear();
  newtoolbar->reset();
  something_changed();
}

void MyWindow::get_new_object(CGAL::Object obj)
{
  Polygon_2 pgn;
  if (CGAL::assign(pgn,obj))
    {
      if (pgn.orientation() == CGAL::CLOCKWISE)
        pgn.reverse_orientation();
      if (red_active)
        red_set.join(pgn);
      else
        blue_set.join(pgn);
      something_changed();
    }
  else
    {
      Circle circ;
      if (CGAL::assign(circ, obj))
	{
	  if (circ.is_degenerate()) // radius == 0
	    return;

	  if (circ.orientation() == CGAL::CLOCKWISE)
	    circ = circ.opposite();

	  std::vector<CGAL::Object> xcurves;
	  xcurves.reserve(2);
	  Traits tr;
	  Curve full_circ(circ);
	  tr.make_x_monotone_2_object()(full_circ, std::back_inserter(xcurves));

	  CGAL_assertion(xcurves.size() == 2);
	  XCurve half_circ1;
	  XCurve half_circ2;
	  CGAL::assign(half_circ1, xcurves[0]);
	  CGAL::assign(half_circ2, xcurves[1]);
	  pgn.push_back(half_circ1);
	  pgn.push_back(half_circ2);
	  if (red_active)
	    red_set.join(pgn);
	  else
	    blue_set.join(pgn);
	  something_changed();
	}
    }
}

void MyWindow::about()
{
  QMessageBox::about( this, my_title_string,
		      "This is a demo for boolean operations on polygons\n");

}

void MyWindow::aboutQt()
{
  QMessageBox::aboutQt( this, my_title_string );
}

void MyWindow::howto(){
  QString home;
  home = "help/index.html";
  CGAL::Qt_help_window *help = new
    CGAL::Qt_help_window(home, ".", 0, "help viewer");
  help->resize(400, 400);
  help->setCaption("Demo HowTo");
  help->show();
}
void MyWindow::init_layer() 
{ 
	testlayer = new Qt_layer_show_ch(this);	
	widget->attach(&(*testlayer));
}

void MyWindow::new_window(){
  MyWindow *ed = new MyWindow(500, 500);
  ed->init_layer();
  //give a number to new window
  //std::string str = "Window ";    
  QString new_title = "Window ";  
  //QDeepCopy<QString> new_title = my_title_string;  
  QString curnum;
  curnum.setNum(winsOpened,10);
  new_title.append(curnum);
  winsOpened++;
  ed->setCaption(new_title);
  ed->widget->clear_history();
  ed->widget->set_window(-1.1, 1.1, -1.1, 1.1);
  ed->show();
  something_changed();
}


void MyWindow::timer_done()
{
  if (old_state!=current_state) {
    widget->redraw();
    old_state = current_state;
  }
}


#include "boolean_operations_2.moc"

int main(int argc, char **argv)
{
  QApplication app( argc, argv );
  // physical window size  
  MyWindow widget(600,400);
  //initialize window's drawing layer  
  widget.init_layer(); 
  app.setMainWidget(&widget);
  widget.setCaption(my_title_string);
  widget.setMouseTracking(TRUE);
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  widget.setIcon(cgal_icon);
  widget.show();
  //current_state = -1;
  return app.exec();
}

