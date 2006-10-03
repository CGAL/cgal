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


// if QT is not installed, a message will be issued in runtime.
#include <CGAL/basic.h>

#ifndef CGAL_USE_QT
#include <iostream>
int main(int, char*){
  std::cout << "Sorry, this demo needs QT...";
  std::cout << std::endl; return 0;
}
#else

#include <fstream>
#include <string>

#include "typedefs.h"

#include <CGAL/IO/Qt_widget.h>




#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/pixmaps/demoicon.xpm>

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


//The envelope diagram
Envelope_diagram_2 diag;
int current_state;

const QString my_title_string("envelopes of 3d triangles");


class Qt_layer_show_ch : public CGAL::Qt_widget_layer
{
public:
    
    // default constructor
  Qt_layer_show_ch(){};

  // this method overrides the virtual method 'draw()' of Qt_widget_layer
  void draw()
  {
    widget->lock(); // widget have to be locked before drawing 
    *widget <<  CGAL::BLUE; 
    for(Edge_const_iterator eit = diag.edges_begin();
        eit != diag.edges_end();
        ++eit)
    {
      *widget << eit->curve();
    }
   
    *widget <<  CGAL::RED; 
    for(Vertex_const_iterator vit = diag.vertices_begin();
        vit != diag.vertices_end();
        ++vit)
    {
      *widget << vit->point();
    }

    widget->unlock(); // widget have to be unlocked when finished drawing
  };    
  
};//end class 


/* The QMainWindow class provides a main application window, 
 *  with a menu bar, dock windows (e.g. for toolbars), and a status bar
 */
class MyWindow : public QMainWindow
{
  Q_OBJECT
public:

    // constructor
  MyWindow(int w, int h) 
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
    file->insertItem("&Open triangles file", this, SLOT(open_file()),CTRL+Key_O);
    file->insertSeparator();
    
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

   


    //layers
    widget->attach(&testlayer);  

   
  
    resize(w,h);
    widget->set_window(-1, 1, -1, 1);
    widget->setMouseTracking(FALSE);
    
    //application flag stuff
    old_state = 0;
  };

private:
  void something_changed(){current_state++;};

 

  

        
public slots:

    void open_file()
    {

      QString s = QFileDialog::getOpenFileName("./",
                                               QString::null,
                                               this,
                                               "open file dialog",
                                               "Choose a file" );
      if(s==QString::null)
        return;

      std::ifstream in_file(s);
      if(!in_file.is_open())
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

      std::list<Surface_3> triangles;
      int num_of_triangles;
      in_file >> num_of_triangles;
      for(int i=0 ; i<num_of_triangles; i++)
      {
        Triangle_3 tri;
        in_file >> tri;
        triangles.push_back(tri);
        box = box + bbox_2d(tri);
      }
      diag.clear();
      CGAL::lower_envelope_3(triangles.begin(), triangles.end(), diag);

      widget->set_window(box.xmin(),
                         box.xmax(),
                         box.ymin(),
                         box.ymax());
      widget->unlock();
      something_changed();
      widget->setCursor(old);
    }

 
  void new_instance()
  {
    widget->lock();
    file_name = QString::null;
    
    diag.clear();
    widget->clear_history();
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
        // set the Visible Area to the Interval
    widget->unlock();
    
    something_changed(); 
  }

private slots:
 

  void about()
  {
    QMessageBox::about( this, my_title_string,
        "This is a demo for 2d envelopes of triangles\n");
          
  }

  void aboutQt()
  {
    QMessageBox::aboutQt( this, my_title_string );
  }

  void howto(){
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
    ed->widget->set_window(-1.1, 1.1, -1.1, 1.1);
    ed->show();
    something_changed();
  }


  void timer_done()
  {
    if(old_state!=current_state){
      widget->redraw();
      old_state = current_state;
    }
  }    

  
private:
  CGAL::Qt_widget*                     widget;
  CGAL::Qt_widget_standard_toolbar*    stoolbar;
  Qt_layer_show_ch                     testlayer;
  

  int                                  old_state;
  QString                              file_name;
};

#include "envelope_3.moc"



int main(int argc, char **argv)
{
  QApplication app( argc, argv );
  MyWindow widget(600,400); // physical window size
  app.setMainWidget(&widget);
  widget.setCaption(my_title_string);
  widget.setMouseTracking(TRUE);
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  widget.setIcon(cgal_icon);
  widget.show();
  current_state = -1;
  return app.exec();   
}

#endif // CGAL_USE_QT
