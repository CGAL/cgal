// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// file          : regular_triangulation_2.C
// package       : Qt_widget
// author(s)     : Radu Ursu 
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================


// if QT is not installed, a message will be issued in runtime.
#ifndef CGAL_USE_QT
#include <iostream>
int main(int, char*)
{
  std::cout << "Sorry, this demo needs QT...";
  std::cout << std::endl;
  return 0;
}
#else

#include "regular_cgal_types.h"

//Qt_widget headers
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include "regular_triangulation_2_toolbar.h"
#include "regular_triangulation_2_toolbar_layers.h"
#include <CGAL/IO/pixmaps/demoicon.xpm>

//Qt headers
#include <qplatinumstyle.h>
#include <qapplication.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qfiledialog.h>
#include <qtimer.h>




const QString title_string("Regular_triangulation_2");

Regular_triangulation rt;
int		 current_state;
double xmin, ymin, xmax, ymax;

class Window : public QMainWindow
{
  Q_OBJECT
public:
  Window(int w, int h)
  { 
    widget = new CGAL::Qt_widget(this);
    setCentralWidget(widget);

    // file menu
    QPopupMenu * file = new QPopupMenu( this );
    menuBar()->insertItem( "&File", file );
    file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
    file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
    file->insertSeparator();
    file->insertItem("&Load Triangulation", this, 
		      SLOT(load_triangulation()), CTRL+Key_L);
    file->insertItem("&Save Triangulation", this,
		      SLOT(save_triangulation()), CTRL+Key_S);
    file->insertSeparator();
    file->insertItem("Print", widget, SLOT(print_to_ps()), CTRL+Key_P);
    file->insertSeparator();
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
    file->insertItem( "&Quit", qApp,
		      SLOT( closeAllWindows() ), CTRL+Key_Q );

    // help menu
    QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Help", help );
    help->insertItem("How To", this, SLOT(howto()), Key_F1);
    help->insertSeparator();
    help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );

    *widget << CGAL::BackgroundColor(CGAL::BLACK);
    *widget << CGAL::LineWidth(3);
    resize(w, h);
    widget->set_window(-1, 1, -1, 1);
    widget->setMouseTracking(TRUE);

    //the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this, "ST");
    //the input layers toolbar
    newtoolbar = new Tools_toolbar(widget, this, &rt);
    //the other layers toolbar
    ltoolbar = new Layers_toolbar(widget, this, &rt);


    //connect the widget to the main function that receives the objects
    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)), 
      this, SLOT(get_new_object(CGAL::Object)));
    
    //create a timer to check if something changed
    QTimer *timer = new QTimer( this );
    connect( timer, SIGNAL(timeout()),
    this, SLOT(timerDone()) );
    timer->start( 200, FALSE );
    old_state = 0; 

  }
  void 
  init_coordinates(){
    xmin = -1; xmax = 1;
    ymin = -1; ymax = 1;
  }
private slots:
  void 
  new_instance(){
    widget->lock();
    widget->clear();
    stoolbar->clear_history();
    rt.clear();
    // set the Visible Area to the Interval
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
    widget->unlock();
    something_changed();
  }
  void 
  new_window(){
    Window *ed = new Window(500, 500);
    ed->setCaption("Layer");    
    if(rt.number_of_vertices() > 1){
      Vertex_iterator it = rt.vertices_begin();
      xmin = xmax = (*it).point().x();
      ymin = ymax = (*it).point().y();
      while(it != rt.vertices_end()) {
        if(xmin > (*it).point().x())
          xmin = (*it).point().x();
        if(xmax < (*it).point().x())
          xmax = (*it).point().x();
        if(ymin > (*it).point().y())
          ymin = (*it).point().y();
        if(ymax < (*it).point().y())
          ymax = (*it).point().y();
        it++;
      }
    }
    ed->stoolbar->clear_history();
    ed->widget->set_window(xmin, xmax, ymin, ymax);
    ed->show();
    something_changed();
  }
  
  void 
  timerDone(){
    if(old_state!=current_state){
      widget->redraw();
      old_state = current_state;
    }
  }

  void 
  howto(){
    QString home;
    home = "help/rindex.html";
    CGAL::Qt_help_window *help = new 
      CGAL::Qt_help_window(home, ".", 0, "help viewer");
    help->resize(400, 400);
    help->setCaption("Demo HowTo");
    help->show();
  }

  void 
  about(){
    QMessageBox::about( this, title_string,
		"This is a demo for Regular_triangulation_2\n"
  		"Copyright CGAL @2003");
  }

  void 
  aboutQt(){
    QMessageBox::aboutQt( this, title_string );
  }

  void 
  get_new_object(CGAL::Object obj){
    Point p;
    Circle c;
    if(CGAL::assign(p,obj)) {
      rt.insert(p);
      something_changed();
    } else if (CGAL::assign(c, obj)){
      rt.insert(Gt::Weighted_point(Point(c.center()), c.squared_radius()));
      something_changed();
    }
  }


  void
  load_triangulation(){
    QString s( QFileDialog::getOpenFileName( QString::null,
			    "CGAL files (*.cgal)", this ) );
    if ( s.isEmpty() )
        return;
    rt.clear();

    //Gt::Weighted_point wp;
    
  }
  void
  save_triangulation(){
  }
private:
  inline  void something_changed(){current_state++;};
  CGAL::Qt_widget *widget;
  CGAL::Qt_widget_standard_toolbar  *stoolbar;
  Tools_toolbar                     *newtoolbar;
  Layers_toolbar                    *ltoolbar;
  int                               old_state;

};

#include "regular_triangulation_2.moc"

int
main(int argc, char **argv)
{
  QApplication app( argc, argv );
  Window W(500,500); // physical widget size
  app.setMainWidget(&W);
  W.setCaption(title_string);
  W.setMouseTracking(TRUE);
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  W.setIcon(cgal_icon);
  W.show();  
  W.init_coordinates();  
  current_state = -1;
  return app.exec();
}


#endif //CGAL_USE_QT
