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


#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/point_generators_2.h>

#include <fstream>
#include <iomanip>

typedef CGAL::MP_Float				    NT;
typedef CGAL::Cartesian<NT>                         K;
typedef CGAL::Partition_traits_2<K>                 Traits;
typedef Traits::Point_2                             Point_2;
typedef Traits::Polygon_2                           Polygon_2;
typedef CGAL::Random_points_in_square_2<Point_2>    Point_generator;


#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include "Qt_widget_toolbar.h"
#include "Qt_widget_toolbar_layers.h"

const QString my_title_string("Polygon partition demo with"
			      " CGAL Qt_widget");


Polygon_2 polygon;

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
    help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );

    //the new tools toolbar
    newtoolbar = new CGAL::Tools_toolbar(widget, this);	
    //the new scenes toolbar
    vtoolbar = new CGAL::Layers_toolbar(widget, this, &polygon);
    //the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this);
    this->addToolBar(stoolbar->toolbar(), Top, FALSE);
    this->addToolBar(newtoolbar->toolbar(), Top, FALSE);
    this->addToolBar(vtoolbar->toolbar(), Top, FALSE);

    *widget << CGAL::BackgroundColor (CGAL::BLACK);
    resize(w,h);
  
    widget->show();
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

  void set_window(double xmin, double xmax, double ymin, double ymax)
  {
    widget->set_window(xmin, xmax, ymin, ymax);
  }
  void new_instance()
  {
    widget->lock();
    widget->clear();
    widget->set_window(-1.1, 1.1, -1.1, 1.1); 
	// set the Visible Area to the Interval
    widget->unlock();
  }


private slots:
  void gen_poly(){
    widget->clear_history();
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
	// set the Visible Area to the Interval
    polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
    CGAL::random_polygon_2(100,
			   std::back_inserter(polygon),
			   Point_generator(1));
    widget->redraw();
  }

  void get_new_object(CGAL::Object obj)
  {
    Polygon_2 poly;
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
      "Polygon partition demo,\n"
      "Copyright CGAL @2001");
  };

  void aboutQt()
  {
    QMessageBox::aboutQt( this, my_title_string );
  }
  void new_window(){
    MyWindow *ed = new MyWindow(500, 500);
    ed->setCaption("View");
    ed->show();
    ed->set_window(-1.1, 1.1, -1.1, 1.1);
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
    widget->redraw();
  }

	

private:
  CGAL::Qt_widget	  *widget;	
  CGAL::Tools_toolbar	  *newtoolbar;
  CGAL::Layers_toolbar	  *vtoolbar;
  CGAL::Qt_widget_standard_toolbar  *stoolbar;
  int			old_state;	
};


#include "partition_2.moc"

int
main(int argc, char **argv)
{
  
  QApplication app( argc, argv );
  app.setStyle( new QPlatinumStyle );
  QPalette p( QColor( 250, 215, 100 ) );
  app.setPalette( p, TRUE );
  current_state = -1;
  
  MyWindow widget(500,500); // physical window size
  app.setMainWidget(&widget);
  widget.setCaption(my_title_string);
  widget.setMouseTracking(TRUE);
  widget.show();
  widget.set_window(-1, 1, -1, 1);
  // because Qt send resizeEvent only on show.
  
  return app.exec();
  return 1;
}

#endif // CGAL_USE_QT
