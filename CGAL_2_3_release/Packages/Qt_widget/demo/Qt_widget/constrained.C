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

#include <CGAL/IO/Qt_Window.h>
#include <qapplication.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qmenubar.h>

typedef double Coord_type;
typedef CGAL::Cartesian<Coord_type>  Rep;

typedef CGAL::Point_2<Rep>  Point;
typedef CGAL::Segment_2<Rep>  Segment;
typedef CGAL::Triangle_2<Rep>  Triangle;

typedef CGAL::Triangulation_euclidean_traits_2<Rep> Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<Gt> Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Constrained_triangulation_2<Gt,Tds>  Constrained_triangulation;

typedef Constrained_triangulation::Constraint     Constraint;

typedef Constrained_triangulation::Face_handle    Face_handle;
typedef Constrained_triangulation::Vertex_handle  Vertex_handle;

typedef CGAL::Qt_widget Window_stream;

const QString my_title_string("Contrained Triangulation Demo with"
			      " CGAL Qt_widget");

void
draw_constraints(Window_stream &win, std::list<Constraint> & lc)
{
  win << CGAL::RED;
  win.lock();
  std::list<Constraint>::iterator cit=lc.begin();
  for( ; cit != lc.end(); ++cit) {
    win << Segment((*cit).first,(*cit).second);
  }
  win.unlock();
  win << CGAL::BLUE;
} 

void
input_constraints_from_file(std::list<Constraint> & list_contraintes,
			    std::ifstream& is)
{
  int n;
  is >> n;
  qDebug("Reading %d constraints", n);
  Point p,q;
  for(; n > 0; n--) {
    is >> p >> q;
    list_contraintes.push_back(std::make_pair(p,q));
  }
}

void
draw_connected_component(const Point&  p, 
			 const Constrained_triangulation& ct,
			 Window_stream& win)
{
  Face_handle fh = ct.locate(p);
  std::set<Face_handle> component; 
  std::stack<Face_handle> st; 
  // component includes the faces of the connected_component
  // stack includes the faces in component whose neighbors
  // have not yet been looked at

  st.push(fh);
  component.insert(fh);
  while (! st.empty()){
    fh = st.top();
    st.pop();
    for(int i = 0 ; i < 3 ; ++i){
      if ( (! fh->is_constrained(i)) && 
	   component.find(fh->neighbor(i)) == component.end() ) {
	component.insert(fh->neighbor(i));
	st.push(fh->neighbor(i));
      }
    }
  }

  // draw
  int width=win.lineWidth();
  win << CGAL::FillColor(CGAL::GREEN) << CGAL::LineWidth(0);
  std::set<Face_handle>::iterator it;
  for ( it = component.begin(); it != component.end(); it++) {
    if (! ct.is_infinite( *it)) win << ct.triangle( *it);
    else win << ct.segment(*it, (*it)->index(ct.infinite_vertex()));
  }
  win << CGAL::LineWidth(width);
  return;
}

class MyWindow : public QMainWindow
{
  Q_OBJECT
public:
  MyWindow(int x, int y): win(this) {
    setCentralWidget(&win);
    connect(&win, SIGNAL(mousePressed(QMouseEvent*)), this,
	    SLOT(mousePressedOnWin(QMouseEvent*)));
    connect(&win, SIGNAL(resized()), this, SLOT(redrawWin()));
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

    win.show();
    resize(x,y);
  };

  void init_paint()
  {
    win.lock();
    win.init(-1.1, 1.1, -1.1, 1.1); // we need to put it there because Qt
    // send resizeEvent only on show.
    load_file("data/fish");
    win.unlock();
    statusBar()->message("Enter points with the left button");
  };

  ~MyWindow() {};
public slots:

  void redrawWin()
  {
    win.lock();
    win.clear();
    win << CGAL::BLUE << ct;
    draw_constraints(win,lc);
    win.unlock();
  }

  void mousePressedOnWin(QMouseEvent* e)
  {
    statusBar()->message("Terminate with right button");
    if(e->button() == Qt::LeftButton)
      {
	Point p(win.x_real(e->x()),win.y_real(e->y()));
	win.clear();
	win.lock();
	win << CGAL::BLUE <<ct;
	draw_connected_component(p, ct, win);
	draw_constraints(win,lc);
	win << p ;
	win.unlock();
      }
    else if(e->button() == Qt::RightButton)
      {
	qApp->quit();
      }
  };
private slots:
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
  void load_file(QString name)
  {
    std::ifstream is(name);
    lc.clear();
    input_constraints_from_file(lc,is);
    ct=lc;
    assert(ct.is_valid());
    redrawWin();
  };

  CGAL::Qt_widget win;
  std::list<Constraint> lc;
  Constrained_triangulation ct;
};

#include "constrained.moc"

int
main(int argc, char **argv)
{
  QApplication app( argc, argv );
  MyWindow win(400,430); // physical window size
  app.setMainWidget(&win);
  win.setCaption(my_title_string);
  win.show();
  win.init_paint(); // initial paiting must be done after show()
  // because Qt send resizeEvent only on show.
  return app.exec();
}

#endif // CGAL_USE_QT
