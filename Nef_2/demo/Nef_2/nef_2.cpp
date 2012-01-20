// Copyright (c) 2002  Max-Planck-Institute Saarbruecken (Germany)
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
// Author(s)     : Radu Ursu

#include <CGAL/basic.h>

#if !defined CGAL_USE_GMP
#include <iostream>
int main(int, char*)
{
  std::cout << "Sorry, this demo needs GMP...";
  std::cout << std::endl;
  return 0;
}
#else

#include <fstream>
#include <stack>
#include <set>
#include <list>

#include "cgal_types.h"
#include "nef_2.h"
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include "Qt_widget_toolbar.h"
#include "nef_2_layers.h"
#include <CGAL/IO/pixmaps/demoicon.xpm>

#include <qapplication.h>
#include <qfiledialog.h>
#include <qlayout.h>
#include <qlistbox.h>
#include <qmainwindow.h>
#include <qmenubar.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qsplitter.h>
#include <qsplitter.h>
#include <qstatusbar.h>
#include <qtimer.h>
#include <qtoolbar.h>
#include <qtoolbutton.h>
#include <qstring.h>

const QString my_title_string("Nef_2 Demo with"
			      " CGAL Qt_widget");

//global flags and variables
int current_state;

//This class contains a Nef_2 object and a text description
class Nef_description{
public:
  Nef_description(Nef_polyhedron n, QString name_p):
      N(n), name(name_p){}
  ~Nef_description(){}
  Nef_polyhedron N;
  QString name;
};

std::list<Nef_description>  nef_2_list;
bool is_the_first_widget = true; //true only before the creation of
                                 //the first widget

//This class is a generalisation of QListBox
class Nef_2_list_box : public QListBox{
  Q_OBJECT
public:
  Nef_2_list_box(QWidget *parent, const char *name):
    QListBox(parent, name)
  {
    setFocusPolicy(QWidget::StrongFocus);
  }
  ~Nef_2_list_box(){}

protected:
  void keyPressEvent(QKeyEvent* e){
    if ( e->key() == Qt::Key_Delete )
      if(isSelected(currentItem()))
        emit(delete_key(item(currentItem())));
  }
signals:
  void delete_key(QListBoxItem *);
};
int nef_index; //used to give names to polyhedrons

bool has_built_layout;


//This widget contains some widgets in a layout
class Layout_widget : public QWidget {
public:
  Layout_widget(QWidget *parent = 0, const char *name = 0)
    : QWidget(parent, name)
  {
    QBoxLayout *topLayout = new QVBoxLayout( this, 5 );
    QSplitter *top_splitter = new QSplitter(Qt::Vertical, this);
    topLayout->addWidget(top_splitter);
    QSplitter *splitter1 = new QSplitter(top_splitter);
    nef_list1 = new Nef_2_list_box(splitter1, "Nef_list_1");
    nef_list2 = new Nef_2_list_box(splitter1, "Nef_list_2");
    widget = new CGAL::Qt_widget(top_splitter);
  };
  CGAL::Qt_widget* get_qt_widget(){return widget;}
  Nef_2_list_box* get_nef_list1(){return nef_list1;}
  Nef_2_list_box* get_nef_list2(){return nef_list2;}
private:
  CGAL::Qt_widget *widget;
  Nef_2_list_box  *nef_list1, *nef_list2;
};

Layout_widget *cwidget;

class MyWindow : public QMainWindow
{
  Q_OBJECT
public:
  MyWindow(int w, int h){
    cwidget = new Layout_widget(this);
    widget = cwidget->get_qt_widget();
    list1 = cwidget->get_nef_list1();
    list2 = cwidget->get_nef_list2();
    setCentralWidget(cwidget);

    if(is_the_first_widget){
      Line l1(Point_2(0, 0), Point_2(0, 2));
      Nef_polyhedron N1(l1, Nef_polyhedron::INCLUDED);
      Nef_visible = N1;
      Line l2(Point_2(0, 0), Point_2(2, 0));
      Nef_polyhedron N2(l2, Nef_polyhedron::INCLUDED);
      Nef_visible2 = N2;
      insert_in_list(N1, QString("LeftHalf"));
      insert_in_list(N2, QString("TopHalf"));
    } else {
       Nef_visible = Nef_polyhedron(Nef_polyhedron::EMPTY);
       Nef_visible2 = Nef_polyhedron(Nef_polyhedron::EMPTY);
    }
    //create a timer for checking if somthing changed
    QTimer *timer = new QTimer( this );
    connect( timer, SIGNAL(timeout()),
           this, SLOT(sl_timer_done()) );
    timer->start( 200, FALSE );

    // file menu
    QPopupMenu * file = new QPopupMenu( this );
    menuBar()->insertItem( "&File", file );
    file->insertItem("&New", this, SLOT(sl_new_instance()), CTRL+Key_N);
    file->insertItem("New &Window", this, SLOT(sl_new_window()), CTRL+Key_W);
    file->insertSeparator();
    file->insertItem("Load Nef_2", this, SLOT(sl_load_nef()), CTRL+Key_L);
    file->insertItem("Save Nef_2", this, SLOT(sl_save_nef()), CTRL+Key_S);
    file->insertItem("Load Polygon_2", this,
                     SLOT(sl_load_polygon()), CTRL+Key_O);
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
    help->insertItem("&About", this, SLOT(sl_about()), CTRL+Key_A );
    help->insertItem("About &Qt", this, SLOT(sl_aboutQt()) );

    //the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this, "ST");
    //the new tools toolbar
    newtoolbar = new Tools_toolbar(widget, this);
    nef_layer1 = new Qt_layer_nef_blue<Nef_polyhedron>(Nef_visible);
    nef_layer2 = new Qt_layer_nef_gray<Nef_polyhedron>(Nef_visible2);
    widget->attach(nef_layer1);
    widget->attach(nef_layer2);
    //boolean operations toolbar
    QToolBar *optoolbar = new QToolBar("operations", this,
                                        QMainWindow::Right, true,"Operations");
    QToolButton *but_intersection = new QToolButton(
       QPixmap( (const char**)intersection_xpm ),
       "INTERSECTION", "BOOLEAN OPERATIONS", this,
       SLOT(sl_intersect()), optoolbar, "intersection");
    QToolButton *but_union = new QToolButton(
       QPixmap( (const char**)union_xpm ),
       "UNION", "BOOLEAN OPERATIONS", this,
       SLOT(sl_join()), optoolbar, "union");
    QToolButton *but_difference = new QToolButton(
       QPixmap( (const char**)difference_xpm ),
       "DIFFERENCE", "BOOLEAN OPERATIONS", this,
       SLOT(sl_difference()), optoolbar, "difference");
    QToolButton *but_symm_difference = new QToolButton(
       QPixmap( (const char**)symmetric_difference_xpm ),
       "SYMMETRIC DIFFERENCE", "BOOLEAN OPERATIONS", this,
       SLOT(sl_symm_difference()), optoolbar, "symm_difference");

    QToolButton *but_complement = new QToolButton(
       QPixmap( (const char**)complement_xpm ),
       "COMPLEMENT", "BOOLEAN OPERATIONS", this,
       SLOT(sl_complement()), optoolbar, "complement");
    QToolButton *but_interior = new QToolButton(
       QPixmap( (const char**)interior_xpm ),
       "INTERIOR", "BOOLEAN OPERATIONS", this,
       SLOT(sl_interior()), optoolbar, "interior");
    QToolButton *but_closure = new QToolButton(
       QPixmap( (const char**)closure_xpm ),
       "CLOSURE", "BOOLEAN OPERATIONS", this,
       SLOT(sl_closure()), optoolbar, "closure");
    QToolButton *but_boundary = new QToolButton(
       QPixmap( (const char**)boundary_xpm ),
       "BOUNDARY", "BOOLEAN OPERATIONS", this,
       SLOT(sl_boundary()), optoolbar, "boundary");
    QToolButton *but_regularization = new QToolButton(
       QPixmap( (const char**)regularization_xpm ),
       "REGULARIZATION", "BOOLEAN OPERATIONS", this,
       SLOT(sl_regularization()), optoolbar, "regularization");

    but_intersection->setOn(true);
    but_union->setOn(true);
    but_difference->setOn(true);
    but_symm_difference->setOn(true);
    but_complement->setOn(true);
    but_interior->setOn(true);
    but_closure->setOn(true);
    but_boundary->setOn(true);
    but_regularization->setOn(true);

    connect(list1, SIGNAL(delete_key(QListBoxItem*)),
            this, SLOT(sl_list1_delete_key_pressed(QListBoxItem*)));
    connect(list2, SIGNAL(delete_key(QListBoxItem*)),
            this, SLOT(sl_list2_delete_key_pressed(QListBoxItem*)));



    *widget << CGAL::BackgroundColor (CGAL::BLACK);

    widget->resize(w,h-100);
    list1->resize(w, 200);
    list2->resize(w, 200);

    resize(w, h + 100);
    widget->set_window(-1, 1, -1, 1);
    widget->setMouseTracking(TRUE);

    //connect the widget to the main function that receives the objects
    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),
    this, SLOT(sl_get_new_object(CGAL::Object)));

    //connect the signals doubleClicked from the lists to the specific slots
    connect(list1, SIGNAL(clicked(QListBoxItem*)),
      this, SLOT(sl_list1_selected_item(QListBoxItem*)));
    connect(list2, SIGNAL(clicked(QListBoxItem*)),
      this, SLOT(sl_list2_selected_item(QListBoxItem*)));

    if(is_the_first_widget){
      list1->setSelected(0, true);
      list2->setSelected(1, true);
    }
    //application flag stuff
    old_state = 0;
    is_the_first_widget =  false;
    has_no_selection_list1 = false;
    has_no_selection_list2 = false;
    widget->show();
  };


public slots:
  void sl_new_instance()
  {
    widget->lock();
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
			// set the Visible Area to the Interval
    Nef_polyhedron N_temp(Nef_polyhedron::EMPTY);
    Nef_visible = N_temp;
    Nef_visible2 = N_temp;
    list1->clear();
    list2->clear();
    nef_2_list.clear();
    has_no_selection_list1 = true;
    has_no_selection_list2 = true;
    nef_index = 0;
    widget->unlock();
    something_changed();
  }


  void  sl_load_nef()
  {
    QString s( QFileDialog::getOpenFileName( QString::null,
		    "CGAL files (*.cgal)", this ) );
    if ( s.isEmpty() )
      return;
    std::ifstream in(s.ascii());
    CGAL::set_ascii_mode(in);
    Nef_polyhedron N_temp(Nef_polyhedron::EMPTY);
    Nef_visible2 = N_temp;
    in >> Nef_visible;
    insert_in_list(Nef_visible, "loaded");
    list1->setSelected(list1->count()-1, true);
    list2->clearSelection();
    has_no_selection_list1 = false;
    has_no_selection_list2 = true;
    something_changed();
  }

  void  sl_save_nef()
  {
    QString fileName =
          QFileDialog::getSaveFileName( "nef_2.cgal",
                       "Cgal files (*.cgal)", this );
    if ( !fileName.isNull() ) {
      // got a file name
      std::ofstream out(fileName.ascii());
      CGAL::set_ascii_mode(out);
      out << Nef_visible << std::endl;
    }
  }//end save_nef()

  void  sl_load_polygon(){
    QString s( QFileDialog::getOpenFileName( QString::null,
		    "CGAL files (*.cgal)", this ) );
    if ( s.isEmpty() )
      return;
    std::ifstream in(s.ascii());
    CGAL::set_ascii_mode(in);

    Polygon_2_double poly;
    in >> poly;
      Vertex_iterator_double it = poly.vertices_begin();
      std::list<Point_2> l_of_p;
      while(it != poly.vertices_end()){
        CGAL::Gmpq p_q_x((*it).x());
        CGAL::Gmpq p_q_y((*it).y());
        RT wsx = p_q_x.numerator() * p_q_y.denominator();
        RT wsy = p_q_y.numerator() * p_q_x.denominator();
        RT wsh  = p_q_x.denominator() * p_q_y.denominator();
        Point_2 p1(wsx, wsy, wsh);
        l_of_p.push_back(p1);
        it++;
      }
      Nef_polyhedron Nt(l_of_p.begin(), l_of_p.end(),
                        Nef_polyhedron::INCLUDED);
      Nef_visible = Nt;
      QString tnr;
      tnr.setNum(poly.size());
      tnr.append("gon");
      insert_in_list(Nt, tnr);
      list1->setSelected(list1->count()-1, true);
      widget->set_window(poly.bbox().xmin(), poly.bbox().xmax(),
                         poly.bbox().ymin(), poly.bbox().ymax());
      Nef_polyhedron N_temp(Nef_polyhedron::EMPTY);
      Nef_visible2 = N_temp;
      list2->clearSelection();
      has_no_selection_list1 = false;
      has_no_selection_list2 = true;
      something_changed();
  }
private slots:
  void  sl_get_new_object(CGAL::Object obj)
  {
    Cartesian_point_2   p;
    Cartesian_polygon_2 poly;
    Cartesian_line_2    line;
    if(CGAL::assign(p, obj)) {
      /*
      CGAL::Quotient<RT> wsxq = double_to_quotient<RT>(p.x());
      CGAL::Quotient<RT> wsyq = double_to_quotient<RT>(p.y());
      RT wsx = wsxq.numerator() * wsyq.denominator();
      RT wsy = wsyq.numerator() * wsxq.denominator();
      RT wsh  = wsxq.denominator() * wsyq.denominator();
      */
      RT wsx = p.x().numerator() * p.y().denominator();
      RT wsy = p.y().numerator() * p.x().denominator();
      RT wsh  = p.x().denominator() * p.y().denominator();
      Point_2 p1(wsx, wsy, wsh);
      Point_2 pt[1] = {p1};
      Nef_polyhedron Nt(pt, pt+1);
      insert_in_list(Nt, "Point");
      Nef_visible = Nt;
      list1->setSelected(list1->count()-1, true);
      Nef_visible2 = Nef_polyhedron(Nef_polyhedron::EMPTY);
      list2->clearSelection();
      has_no_selection_list1 = false;
      has_no_selection_list2 = true;
      something_changed();
    } else if(CGAL::assign(poly, obj)){
      Vertex_iterator it = poly.vertices_begin();
      std::list<Point_2> l_of_p;
      while(it != poly.vertices_end()){
        //double xp = (*it).x();
        //double yp = (*it).y();
        //CGAL::Quotient<RT> wsxq = double_to_quotient<RT>(xp);
        //CGAL::Quotient<RT> wsyq = double_to_quotient<RT>(yp);
        RT wsx = (*it).x().numerator() * (*it).y().denominator();
        RT wsy = (*it).y().numerator() * (*it).x().denominator();
        RT wsh  = (*it).x().denominator() * (*it).y().denominator();
        Point_2 p1(wsx, wsy, wsh);
        l_of_p.push_back(p1);
        it++;
      }
      Nef_polyhedron Nt(l_of_p.begin(), l_of_p.end(),
                        Nef_polyhedron::INCLUDED);
      Nef_visible = Nt;
      QString tnr;
      tnr.setNum(poly.size());
      tnr.append("gon");
      insert_in_list(Nt, tnr);
      list1->setSelected(list1->count()-1, true);
      Nef_visible2 = Nef_polyhedron(Nef_polyhedron::EMPTY);
      list2->clearSelection();
      has_no_selection_list1 = false;
      has_no_selection_list2 = true;
      something_changed();
    } else if(CGAL::assign(line, obj)){
      /*
      CGAL::Quotient<RT> wsxq = double_to_quotient<RT>(line.point(0).x());
      CGAL::Quotient<RT> wsyq = double_to_quotient<RT>(line.point(0).y());
      RT wsx = wsxq.numerator() * wsyq.denominator();
      RT wsy = wsyq.numerator() * wsxq.denominator();
      RT wsh  = wsxq.denominator() * wsyq.denominator();
      */
      RT wsx = line.point(0).x().numerator()
               * line.point(0).y().denominator();
      RT wsy = line.point(0).y().numerator()
               * line.point(0).x().denominator();
      RT wsh  = line.point(0).x().denominator()
               * line.point(0).y().denominator();
      Point_2 p1(wsx, wsy, wsh);
      /*
      wsxq = double_to_quotient<RT>(line.point(1).x());
      wsyq = double_to_quotient<RT>(line.point(1).y());
      wsx = wsxq.numerator() * wsyq.denominator();
      wsy = wsyq.numerator() * wsxq.denominator();
      wsh  = wsxq.denominator() * wsyq.denominator();
      */
      wsx = line.point(1).x().numerator() * line.point(1).y().denominator();
      wsy = line.point(1).y().numerator() * line.point(1).x().denominator();
      wsh  = line.point(1).x().denominator()
             * line.point(1).y().denominator();
      Point_2 p2(wsx, wsy, wsh);

      Nef_polyhedron Nt(Line(p1, p2), Nef_polyhedron::INCLUDED);
      Nef_visible = Nt;
      insert_in_list(Nt, "HalfSpace");
      list1->setSelected(list1->count()-1, true);
      Nef_visible2 = Nef_polyhedron(Nef_polyhedron::EMPTY);
      list2->clearSelection();
      has_no_selection_list1 = false;
      has_no_selection_list2 = true;
      something_changed();
    }
  };

  void howto(){
    QString home;
    home = "help/index.html";
    CGAL::Qt_help_window *help = new
      CGAL::Qt_help_window(home, ".", 0, "help viewer");
    help->resize(400, 400);
    help->setCaption("Demo HowTo");
    help->show();
  }


  void  sl_about()
  {
    QMessageBox::about( this, my_title_string,
		"This is a demo for Nef_2,\n"
  		"Copyright CGAL @2002");
  };

  void  sl_aboutQt()
  {
    QMessageBox::aboutQt( this, my_title_string );
  }

  void  sl_new_window(){
    MyWindow *ed = new MyWindow(500, 500);
    ed->setCaption("Layer");
    ed->stoolbar->clear_history();
    ed->widget->set_window(-1.1, 1.1, -1.1, 1.1);
    ed->show();
    something_changed();
  }

  void  sl_timer_done()
  {
    if(old_state!=current_state){
      int temp_index_1 = list1->currentItem();
      int temp_index_2 = list2->currentItem();
      list1->clear();
      list2->clear();
      std::list<Nef_description>::iterator it = nef_2_list.begin();
      while(it!=nef_2_list.end()){
        list1->insertItem((*it).name);
        list2->insertItem((*it).name);
	it++;
      }
      if(!has_no_selection_list1)
        list1->setSelected(temp_index_1, true);
      if(!has_no_selection_list2)
        list2->setSelected(temp_index_2, true);
      widget->redraw();
      old_state = current_state;
    }
  }

  void  sl_list1_selected_item(QListBoxItem* item){
    if(item){
      if(Nef_visible != return_selected_nef(item->text())){
        Nef_visible = return_selected_nef(item->text());
	has_no_selection_list1 = false;
      } else {
        Nef_visible = Nef_polyhedron(Nef_polyhedron::EMPTY);
        list1->setSelected(item, false);
	has_no_selection_list1 = true;
      }
      widget->redraw();
    }
  }
  void  sl_list2_selected_item(QListBoxItem* item){
    if(item){
      if(Nef_visible2 != return_selected_nef(item->text())){
        Nef_visible2 = return_selected_nef(item->text());
	has_no_selection_list2 = false;
      } else {
        Nef_visible2 = Nef_polyhedron(Nef_polyhedron::EMPTY);
        list2->setSelected(item, false);
	has_no_selection_list2 = true;
      }
      widget->redraw();
    }
  }
  void  sl_list1_delete_key_pressed(QListBoxItem* item){
    if(item){
      int ci = list1->currentItem();
      std::list<Nef_description>::iterator it = nef_2_list.begin();
      for(int i=0; i<ci; i++)
	it++;
      nef_2_list.erase(it);
      Nef_visible = Nef_polyhedron(Nef_polyhedron::EMPTY);
      has_no_selection_list1 = true;
      if(list2->currentItem() == ci){
        Nef_visible2 = Nef_polyhedron(Nef_polyhedron::EMPTY);
	has_no_selection_list2 = true;
      }
      something_changed();
    }
  }
  void  sl_list2_delete_key_pressed(QListBoxItem* item){
    if(item){
      int ci = list2->currentItem();
      std::list<Nef_description>::iterator it = nef_2_list.begin();
      for(int i=0; i<ci; i++)
	it++;
      nef_2_list.erase(it);
      Nef_visible2 = Nef_polyhedron(Nef_polyhedron::EMPTY);
      has_no_selection_list2 = true;
      if(list1->currentItem() == ci){
        Nef_visible = Nef_polyhedron(Nef_polyhedron::EMPTY);
	has_no_selection_list1 = true;
      }
      something_changed();
    }
  }
  void  sl_intersect(){
    if( list1->isSelected(list1->currentItem())
        && list2->isSelected(list2->currentItem()) ){
      //if there is something selected
      Nef_polyhedron NT = Nef_visible.intersection(Nef_visible2);
      QString s = "AND(" + list1->currentText() + ", "
                  + list2->currentText() + ")";
      insert_in_list(NT, s);
      Nef_visible = NT;
      Nef_visible2 = Nef_polyhedron(Nef_polyhedron::EMPTY);
      list1->setSelected(list1->count()-1, true);
      list2->setSelected(list2->currentItem(), false);
      has_no_selection_list1 = false;
      has_no_selection_list2 = true;
      something_changed();
    }
  }
  void  sl_join(){
    if( list1->isSelected(list1->currentItem())
        && list2->isSelected(list2->currentItem()) ){
      //if there is something selected
      Nef_polyhedron NT = Nef_visible.join(Nef_visible2);
      QString s = "OR(" + list1->currentText() + ", "
                  + list2->currentText() + ")";
      insert_in_list(NT, s);
      Nef_visible = NT;
      Nef_visible2 = Nef_polyhedron(Nef_polyhedron::EMPTY);
      list1->setSelected(list1->count()-1, true);
      list2->setSelected(list2->currentItem(), false);
      has_no_selection_list1 = false;
      has_no_selection_list2 = true;
      something_changed();
    }
  }
  void  sl_difference(){
    if( list1->isSelected(list1->currentItem())
        && list2->isSelected(list2->currentItem()) ){
      //if there is something selected
      Nef_polyhedron NT = Nef_visible.difference(Nef_visible2);
      QString s = "Dif(" + list1->currentText() + ", "
                  + list2->currentText() + ")";
      insert_in_list(NT, s);
      Nef_visible = NT;
      Nef_visible2 = Nef_polyhedron(Nef_polyhedron::EMPTY);
      list1->setSelected(list1->count()-1, true);
      list2->setSelected(list2->currentItem(), false);
      has_no_selection_list1 = false;
      has_no_selection_list2 = true;
      something_changed();
    }
  }
  void  sl_symm_difference(){
    if( list1->isSelected(list1->currentItem())
        && list2->isSelected(list2->currentItem()) ){
      //if there is something selected
      Nef_polyhedron NT = Nef_visible.symmetric_difference(Nef_visible2);
      QString s = "SYM_DIF(" + list1->currentText() + ", "
                  + list2->currentText() + ")";
      insert_in_list(NT, s);
      Nef_visible = NT;
      Nef_visible2 = Nef_polyhedron(Nef_polyhedron::EMPTY);
      list1->setSelected(list1->count()-1, true);
      list2->setSelected(list2->currentItem(), false);
      has_no_selection_list1 = false;
      has_no_selection_list2 = true;
      something_changed();
    }
  }
  void  sl_complement(){
    if( list1->isSelected(list1->currentItem()) ){
      //if there is something from the 1st list selected
      Nef_polyhedron NT = Nef_visible.complement();
      QString s = "COMPLEMENT( " + list1->currentText() + " )";
      insert_in_list(NT, s);
      Nef_visible = NT;
      Nef_visible2 = Nef_polyhedron(Nef_polyhedron::EMPTY);
      list1->setSelected(list1->count()-1, true);
      has_no_selection_list1 = false;
      has_no_selection_list2 = true;
      something_changed();
    }
  }
  void  sl_interior(){
    if( list1->isSelected(list1->currentItem()) ){
      //if there is something from the 1st list selected
      Nef_polyhedron NT = Nef_visible.interior();
      QString s = "INTERIOR( " + list1->currentText() + " )";
      insert_in_list(NT, s);
      Nef_visible = NT;
      Nef_visible2 = Nef_polyhedron(Nef_polyhedron::EMPTY);
      list1->setSelected(list1->count()-1, true);
      has_no_selection_list1 = false;
      has_no_selection_list2 = true;
      something_changed();
    }
  }
  void  sl_closure(){
    if( list1->isSelected(list1->currentItem()) ){
      //if there is something from the 1st list selected
      Nef_polyhedron NT = Nef_visible.closure();
      QString s = "CLOSURE( " + list1->currentText() + " )";
      insert_in_list(NT, s);
      Nef_visible = NT;
      list1->setSelected(list1->count()-1, true);
      has_no_selection_list1 = false;
      has_no_selection_list2 = true;
      something_changed();
    }
  }
  void  sl_boundary(){
    if( list1->isSelected(list1->currentItem()) ){
      //if there is something from the 1st list selected
      Nef_polyhedron NT = Nef_visible.boundary();
      QString s = "BOUNDARY( " + list1->currentText() + " )";
      insert_in_list(NT, s);
      Nef_visible = NT;
      Nef_visible2 = Nef_polyhedron(Nef_polyhedron::EMPTY);
      list1->setSelected(list1->count()-1, true);
      has_no_selection_list1 = false;
      has_no_selection_list2 = true;
      something_changed();
    }
  }
  void  sl_regularization(){
    if( list1->isSelected(list1->currentItem()) ){
      //if there is something from the 1st list selected
      Nef_polyhedron NT = Nef_visible.regularization();
      QString s = "REGULARIZATION( " + list1->currentText() + " )";
      insert_in_list(NT, s);
      Nef_visible = NT;
      Nef_visible2 = Nef_polyhedron(Nef_polyhedron::EMPTY);
      list1->setSelected(list1->count()-1, true);
      has_no_selection_list1 = false;
      has_no_selection_list2 = true;
      something_changed();
    }
  }
private:
  void insert_in_list(Nef_polyhedron n, QString name)
  {
    QString tnr;
    tnr.setNum(nef_index++);
    tnr.append("gon");
    QString tname("N");
    tname.append(tnr);
    tname.append("= ");
    tname += name;
    Nef_description tempND(n, tname);
    nef_2_list.push_back(tempND);
    list1->insertItem(tname);
    list2->insertItem(tname);
  }

  Nef_polyhedron
  return_selected_nef(QString text){
    std::list<Nef_description>::const_iterator
      it = nef_2_list.begin();
    while(it != nef_2_list.end()) {
      if((*it).name == text)
        return (*it).N;
      ++it;
    }
    CGAL_error();
    return Nef_polyhedron(); // kill warning.
  }

  void  something_changed(){current_state+=2;};

  CGAL::Qt_widget             *widget;
  Tools_toolbar               *newtoolbar;
  CGAL::Qt_widget_standard_toolbar
                              *stoolbar;
  int                         old_state;
                              //used to refresh the current window
  Nef_2_list_box                     *list1, *list2;
  Qt_layer_nef_blue<Nef_polyhedron>  *nef_layer1;
  Qt_layer_nef_gray<Nef_polyhedron>  *nef_layer2;
  Nef_polyhedron Nef_visible;
  Nef_polyhedron Nef_visible2;

  bool has_no_selection_list2;
  bool has_no_selection_list1;

};

#include "nef_2.moc"


int
main(int argc, char **argv)
{
  QApplication app( argc, argv );
  has_built_layout = false;
  current_state = -1;
  nef_index = 0;
  MyWindow widget(500,500); // physical window size
  app.setMainWidget(&widget);
  widget.setCaption(my_title_string);
  widget.setMouseTracking(TRUE);
#if !defined (__POWERPC__)
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  widget.setIcon(cgal_icon);
#endif
  widget.show();
  return app.exec();
}

#endif // CGAL_USE_GMP
