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

#include <fstream>
#include <stack>
#include <set>
#include <string>

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Constrained_triangulation_2.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_Constrained_triangulation_2.h>
#include <CGAL/IO/Qt_widget_get_point.h>
#include <qapplication.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qmenubar.h>

typedef double Coord_type;
typedef CGAL::Cartesian<Coord_type>  K;

typedef K::Point_2  Point;
typedef K::Segment_2  Segment;
typedef K::Triangle_2  Triangle;

typedef CGAL::Constrained_triangulation_2<K>  Constrained_triangulation;

typedef Constrained_triangulation::Constraint     Constraint;

typedef Constrained_triangulation::Face_handle    Face_handle;
typedef Constrained_triangulation::Vertex_handle  Vertex_handle;

const QString my_title_string("Contrained Triangulation Demo with"
			      " CGAL Qt_widget");


class MyWindow : public QMainWindow
{
  Q_OBJECT
public:
  MyWindow(int x, int y) {
    widget = new CGAL::Qt_widget(this);
    setCentralWidget(widget);
    widget->set_window(-1.1, 1.1, -1.1, 1.1, true);

    point_factory = new CGAL::Qt_widget_get_point<K>();
    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),
	    this, SLOT(new_point(CGAL::Object)) );
    widget->attach(point_factory);

    connect(widget, SIGNAL(custom_redraw()), this, SLOT(redrawWin()) );

    statusBar();
    
    // file menu
    QPopupMenu * file = new QPopupMenu( this );
    menuBar()->insertItem( "&File", file );
    file->insertItem("&Open", this, SLOT(open()), CTRL+Key_O );
    file->insertItem("&Quit", this, SLOT(close()), CTRL+Key_Q );

    // help menu
    QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Help", help );
    help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );

    widget->show();
    resize(x,y);
  };

  void draw_constraints()
  {
    *widget << CGAL::RED;
    widget->lock();
    std::list<Constraint>::iterator cit=lc.begin();
    for( ; cit != lc.end(); ++cit) {
      *widget << Segment((*cit).first,(*cit).second);
    }
    widget -> unlock();
    *widget << CGAL::BLUE;
  }

  void draw_connected_component(const Point&  p)
  {
    
    Face_handle fh = ct.locate(p);
    std::set<Face_handle> component; 
    std::list<Face_handle> st; 
    //std::list<Vertex_handle> stv; 
    // component includes the faces of the connected_component
    // stack includes the faces in component whose neighbors
    // have not yet been looked at
    
  
    st.push_back(fh);
    component.insert(fh);
    while (! st.empty()){
      fh = st.back();
      st.pop_back();
      for(int i = 0 ; i < 3 ; ++i){
	if ( (! fh->is_constrained(i)) && 
	     component.find(fh->neighbor(i)) == component.end() ) {
	  component.insert(fh->neighbor(i));
	  st.push_back(fh->neighbor(i));
	}
      }
    }

    // draw
    int width=widget->lineWidth();
    *widget << CGAL::FillColor(CGAL::GREEN) << CGAL::LineWidth(0);
    std::set<Face_handle>::iterator it;
    for ( it = component.begin(); it != component.end(); it++) {
      if (! ct.is_infinite( *it)) *widget << ct.triangle( *it);
      else *widget << ct.segment(*it, (*it)->index(ct.infinite_vertex()));
    }
    *widget << CGAL::LineWidth(width);

  };
  
  void init_paint()
  {
    widget->lock();
    load_file("data/fish");
    widget->unlock();
    statusBar()->message("Enter points with the left button");
  };

  void load_file(QString name)
  {
    std::ifstream is(name);
    lc.clear();

    int n;
    is >> n;
    qDebug("Reading %d constraints", n);
    Point p,q;
    for(; n > 0; n--) {
      is >> p >> q;
      lc.push_back(std::make_pair(p,q));
    }

    ct=lc;
    assert(ct.is_valid());
    redrawWin();
  };

  ~MyWindow()
  {
    delete(point_factory);
  };
public slots:

  void redrawWin()
  {
    widget->lock();
    widget->clear();
    *widget << CGAL::BLUE << ct;
    draw_constraints();
    widget->unlock();
  }

  void new_point(CGAL::Object obj)
  {
    Point p;
    if (CGAL::assign(p,obj))
      {
	widget->clear();
	widget->lock();
	*widget << CGAL::BLUE <<ct;
	draw_connected_component(p);
	draw_constraints();
	*widget << p ;
	widget->unlock();
      }
  };

  void about()
  {
    QMessageBox::about( this, my_title_string,
			"This is a demo from Mariette Yvinec courses,\n"
			"adapted to work with CGAL Qt_widget by\n"
			"Laurent Rineau ( rineau@clipper.ens.fr )");

  };

  void open()
  {
    QString fileName = QFileDialog::getOpenFileName( "data",
						     QString::null, this );
     if ( fileName.isEmpty() )
        return;
     load_file(fileName);
  };

private:
  CGAL::Qt_widget* widget;
  CGAL::Qt_widget_layer* point_factory;
  std::list<Constraint> lc;
  Constrained_triangulation ct;
};

// moc_source_file : constrained.C
#include "constrained.moc"

int
main(int argc, char **argv)
{
  QApplication app( argc, argv );
  MyWindow win(400,430); // physical window size
  app.setMainWidget(&win);
  win.setCaption(my_title_string);
  win.show();

  // initial paiting must be done after show()
  // because Qt send resizeEvent only on show.
  win.load_file("data/fish");
  win.statusBar()->message("Enter points with the left button");

  return app.exec();
}

#endif // CGAL_USE_QT
