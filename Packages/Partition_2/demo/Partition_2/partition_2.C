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

#include <qplatinumstyle.h>
#include <qapplication.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qfiledialog.h>
#include <qtimer.h>
#include <qthread.h>


#include "cgal_types.h"
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include "partition_2_toolbar.h"
#include "partition_2_toolbar_layers.h"
#include <CGAL/IO/pixmaps/demoicon.xpm> 

#include <fstream>
#include <iomanip>


const QString my_title_string("Polygon partition demo with"
			      " CGAL Qt_widget");


Cgal_Polygon polygon;

int current_state;

class MyWindow : public QMainWindow
{
  Q_OBJECT
public:
  MyWindow(int w, int h){
    widget = new CGAL::Qt_widget(this);
    setCentralWidget(widget);
    
    //create a timer for checking if somthing changed
    QTimer *timer = new QTimer( this );
    connect( timer, SIGNAL(timeout()),
           this, SLOT(timerDone()) );
    timer->start( 200, FALSE );

    // file menu
    QPopupMenu * file = new QPopupMenu( this );
    menuBar()->insertItem( "&File", file );
    file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
    file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
    file->insertSeparator();
    file->insertItem("&Load Polygon", this, SLOT(load_polygon()), CTRL+Key_L);
    file->insertItem("&Save Polygon", this, SLOT(save_polygon()), CTRL+Key_S);
    file->insertSeparator();
    file->insertItem("Print", widget, SLOT(print_to_ps()), CTRL+Key_P);
    file->insertSeparator();
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
    file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), CTRL+Key_Q );

    // drawing menu
    QPopupMenu * draw = new QPopupMenu( this );
    menuBar()->insertItem( "&Draw", draw );
    draw->insertItem("Generate Polygon", this, SLOT(gen_poly()), CTRL+Key_G);

    // help menu
    QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Help", help );
    help->insertItem("How To", this, SLOT(howto()), Key_F1);
    help->insertSeparator();
    help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );

    //the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this, "ST");
    //the new tools toolbar
    newtoolbar = new Tools_toolbar(widget, this);	
    //the new scenes toolbar
    vtoolbar = new Layers_toolbar(widget, this, &polygon);

    *widget << CGAL::BackgroundColor (CGAL::BLACK);
    resize(w,h);  
    widget->set_window(-1, 1, -1, 1);
    widget->setMouseTracking(true);
	
    //connect the widget to the main function that receives the objects
    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)), 
		    this, SLOT(get_new_object(CGAL::Object)));

    //application flag stuff
    old_state = 0;
  };

private:
  void something_changed(){current_state++;};
  
public slots:
  void new_instance()
  {
    widget->lock();
    widget->clear();
    polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
    // set the Visible Area to the Interval
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
    widget->unlock();
  }


private slots:
  void gen_poly(){
    stoolbar->clear_history();
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
    // set the Visible Area to the Interval
    polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
    CGAL::random_polygon_2(100,
			   std::back_inserter(polygon),
			   Point_generator(1));
    something_changed();
  }

  void get_new_object(CGAL::Object obj)
  {
    Cgal_Polygon poly;
    if (CGAL::assign(poly, obj))
    {
      polygon = poly;
      something_changed();
    }
    widget->redraw();
  };

	
  void about()
  {
    QMessageBox::about( this, my_title_string,
      "Polygon partition demo\n"
      "Copyright CGAL @2003");
  };

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
    ed->setCaption("View");
    ed->stoolbar->clear_history();
    ed->widget->set_window(-1.1, 1.1, -1.1, 1.1);
    ed->show();
    something_changed();
  }

  void timerDone()
  {
    if(old_state!=current_state){
      widget->redraw();
      old_state = current_state;
    }
  }	

	
  void save_polygon()
  {
    QString fileName = QFileDialog::getSaveFileName( 
		"polygon.cgal", "Cgal files (*.cgal)", this );
    if ( !fileName.isNull() ) {                 // got a file name
      std::ofstream out(fileName);
      //out << std::setprecision(15);
      CGAL::set_ascii_mode(out);
      out << polygon << std::endl;
    }
  }

  void load_polygon()
  {
    QString s( QFileDialog::getOpenFileName(
		QString::null, "CGAL files (*.cgal)", this ) );
    if ( s.isEmpty() )
        return;
    std::ifstream in(s);
    CGAL::set_ascii_mode(in);
    in >> polygon;
    something_changed();
  }

	

private:
  CGAL::Qt_widget        *widget;
  CGAL::Qt_widget_standard_toolbar
                         *stoolbar;
  Tools_toolbar          *newtoolbar;
  Layers_toolbar         *vtoolbar;
  int                    old_state;
};


#include "partition_2.moc"

int
main(int argc, char **argv)
{
  QApplication app( argc, argv );
  current_state = -1;
  
  MyWindow widget(500,500); // physical window size
  app.setMainWidget(&widget);
  widget.setCaption(my_title_string);
  widget.setMouseTracking(TRUE);
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  widget.setIcon(cgal_icon);
  widget.show();
  return app.exec();
  return 1;
}

#endif // CGAL_USE_QT
