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



#include "cgal_types.h"

//Qt_widget
#include "Qt_widget_toolbar_layers.h"
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_Alpha_shape_2.h>
#include "Qt_widget_toolbar.h"

//STL 
#include <fstream>
#include <stack>
#include <set>
#include <string>


//Qt
#include <qapplication.h>
#include <qfiledialog.h>
#include <qimage.h>
#include <qinputdialog.h>
#include <qlabel.h>
#include <qlayout.h>
#include <qmainwindow.h>
#include <qmessagebox.h>
#include <qmenubar.h>
#include <qpopupmenu.h>
#include <qplatinumstyle.h>
#include <qslider.h>
#include <qstatusbar.h>
#include <qtimer.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>

const QString my_title_string("Alpha_shapes_2 Demo with"
			      " CGAL Qt_widget");

//Global variables
Delaunay          tr1;
CGALPointlist     L;
std::list<Wpoint> LW;
Alpha_shape       A;
Alpha_shape_w     AW;
int               current_state;
double            alpha_index;
QImage            image;
Coord_type        xmin, ymin, xmax, ymax;


class Qt_layer_show_alpha_shape : public CGAL::Qt_widget_layer
{
public:
  Qt_layer_show_alpha_shape(){};
private:
  void draw(){
    A.set_alpha(alpha_index);
    A.set_mode(Alpha_shape::GENERAL);
    *widget << CGAL::LineWidth(2) << CGAL::GREEN;
    *widget << A;
  }
};

class Qt_layer_show_image : public CGAL::Qt_widget_layer{
public:
  Qt_layer_show_image(QImage* i) : image(i){}
  ~Qt_layer_show_image(){};
  void draw(){
    widget->get_painter().drawImage(0, 0, *image, 0, 0, widget->width(), widget->height(), 0);
  }
  QImage *image;
};
class Layout_widget : public QWidget{
public:
  Layout_widget(QWidget *parent, const char *name=0):
      QWidget(parent, name) {    
    QBoxLayout *topLayout = new QVBoxLayout( this, 5 );
    QBoxLayout *topLayout1 = new QHBoxLayout( this, 5 );
    slider = new QSlider(0, 10000, 1, 10, Qt::Vertical, this, "slider1");
    label = new QLabel(this, "label");
    label->setText("The current alpha value: 0.001");
    //QBoxLayout *bottomLayout = new QVBoxLayout( topLayout );
    widget = new CGAL::Qt_widget(this);
    topLayout1->addWidget(widget);
    topLayout1->addWidget(slider);
    topLayout->addWidget(label);
    topLayout->addLayout(topLayout1);
  }
  ~Layout_widget(){}
  CGAL::Qt_widget* get_qt_widget(){return widget;}
  QSlider*  get_slider(){return slider;}
  QLabel*   get_label(){return label;}
private:
	CGAL::Qt_widget *widget;
  QSlider         *slider;
  QLabel          *label;
};
class MyWindow : public QMainWindow
{
  Q_OBJECT
public:
  MyWindow(int w, int h) {
  Layout_widget *cwidget = new Layout_widget(this);
  widget = cwidget->get_qt_widget();
  slider = cwidget->get_slider();
  label  = cwidget->get_label();
  setCentralWidget(cwidget);

  //image = new QImage();
  image_layer = new Qt_layer_show_image(&image);
  widget->attach(image_layer);

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
  file->insertItem("&Load Triangulation", this, 
		   SLOT(load_triangulation()), CTRL+Key_L);
  file->insertItem("&Save Triangulation", 
		   this, SLOT(save_triangulation()), CTRL+Key_T);
  file->insertSeparator();
  file->insertItem("&Load Image", 
		   this, SLOT(load_image()), CTRL+Key_I);
  file->insertSeparator();
  file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
  file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), CTRL+Key_Q );

  // drawing menu
  QPopupMenu * edit = new QPopupMenu( this );
  menuBar()->insertItem( "&Edit", edit );
  edit->insertItem("&Generate Triangulation", this, 
		   SLOT(gen_tr()), CTRL+Key_G );
  edit->insertItem("&Change Alpha", this, 
		   SLOT(change_alpha()), CTRL+Key_C );
  edit->insertItem("C&rust", this, 
		   SLOT(show_crust()), CTRL+Key_R );
  edit->insertItem("&Weighted Alpha Shape", this, 
		   SLOT(show_weighted()), CTRL+Key_W );

  // help menu
  QPopupMenu * help = new QPopupMenu( this );
  menuBar()->insertItem( "&Help", help );
  help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
  help->insertItem("About &Qt", this, SLOT(aboutQt()) );

  //the new tools toolbar
  setUsesBigPixmaps(TRUE);
  newtoolbar = new CGAL::Tools_toolbar(widget, this, &tr1);	
  //the new scenes toolbar
  vtoolbar = new CGAL::Layers_toolbar(widget, this, &tr1);
  //the standard toolbar
  stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this);
  this->addToolBar(stoolbar->toolbar(), Top, FALSE);
  this->addToolBar(newtoolbar->toolbar(), Top, FALSE);
  this->addToolBar(vtoolbar->toolbar(), Top, FALSE);
  
  *widget << CGAL::LineWidth(2) << CGAL::BackgroundColor (CGAL::BLACK);

  resize(w,h);
  widget->show();
  widget->setMouseTracking(TRUE);  
	
  //connect the widget to the main function that receives the objects
  connect(widget, SIGNAL(new_cgal_object(CGAL::Object)), 
    this, SLOT(get_new_object(CGAL::Object)));
  connect(slider, SIGNAL(valueChanged(int)), this, SLOT(slider_changed(int)));

  //application flag stuff
  old_state = 0;
  //layers
  widget->attach(&alpha_shape_layer);
  alpha_index = 0.001;
  };

  ~MyWindow()
  {
  };

private:
  void something_changed(){current_state++;};
signals:
  void was_repainted();
  
public slots:
  void set_window(double xmin, double xmax, double ymin, double ymax)
  {
    widget->set_window(xmin, xmax, ymin, ymax);
  }
  void new_instance()
  {
    widget->lock();
    widget->clear();
    tr1.clear();
    A.clear();AW.clear();
    L.clear();LW.clear();
    widget->clear_history();
    widget->set_window(-1.1, 1.1, -1.1, 1.1); 
    // set the Visible Area to the Interval
    widget->unlock();
    something_changed();
  }

  void load_image(){
    QString s( QFileDialog::getOpenFileName( QString::null,
		    "Image files (*.bmp)", this ) );
    if ( s.isEmpty() )
        return;
    QImage img(s);
    image = img;
    something_changed();
  }
  void show_crust(){
    Delaunay T(tr1);
    Face_iterator f;
    typedef std::set<Point, Point_compare> point_set;
    point_set s;
    for (f=T.faces_begin(); f!=T.faces_end(); ++f)
      s.insert( T.dual(f) );
    widget->lock();
    point_set::iterator p;
    for (p=s.begin(); p!=s.end(); ++p)
      T.insert(*p);
    Edge_iterator e;
    for (e=T.edges_begin(); e!=T.edges_end(); ++e) {
      Face_handle f=(*e).first; int i=(*e).second;
      bool s1=s.find(f->vertex(f->ccw(i))->point())==s.end();
      bool s2=s.find(f->vertex(f->cw(i))->point())==s.end();
      if (s1&&s2)
        *widget << CGAL::YELLOW << T.segment(e);
    }
    widget->unlock();
  }
  void show_weighted(){
    AW.set_alpha(alpha_index);
    *widget << CGAL::LineWidth(2) << CGAL::GRAY;
    widget->lock();
    *widget << AW;
    widget->unlock();
  }
private slots:
  void get_new_object(CGAL::Object obj)
  {
    Point   p;
    Wpoint  wp;
    if(CGAL::assign(p,obj)) {
      wp = p;
      tr1.insert(p);
      LW.push_back(wp);
      L.push_back(p);
      A.make_alpha_shape(L.begin(), L.end());
      AW.make_alpha_shape(LW.begin(), LW.end());
      something_changed();
    } 
  };

  void about()
  {
    QMessageBox::about( this, my_title_string,
		"This is a demo for Triangulation,\n"
  		"Copyright CGAL @2001");
  };

  void aboutQt()
  {
    QMessageBox::aboutQt( this, my_title_string );
  }

  void new_window(){
    MyWindow *ed = new MyWindow(500, 500);
    ed->setCaption("Layer");
    ed->show();
    ed->set_window(xmin, xmax, ymin, ymax);
    something_changed();
  }

  void change_alpha() {
    bool ok = FALSE;
    double res = QInputDialog::getDouble( 
		 tr( "Please enter a decimal number" ),
		 "Between 0 and 1", alpha_index, 0, 1, 3, &ok, this );
    if ( ok ){
      alpha_index = res;
      slider->setValue(res*10000);
      widget->redraw();      
    }
  }

  void timerDone()
  {
    if(old_state!=current_state){
      widget->redraw();
      old_state = current_state;
    }
  }	

  void gen_tr()
  {
    tr1.clear();
    L.clear();
    LW.clear();
    A.clear();
    AW.clear();
    widget->clear_history();
    widget->lock();
    widget->set_window(-1.1, 1.1, -1.1, 1.1); 
    // set the Visible Area to the Interval
    widget->unlock();
    Wpoint pw;
    CGAL::Random_points_in_disc_2<Point> g(0.5);
    for(int count=0; count<200; count++) {
      tr1.insert(*g);
      pw = *g;
      L.push_back(*g++);
      LW.push_back(pw);
    }
    A.make_alpha_shape(L.begin(), L.end());
    AW.make_alpha_shape(LW.begin(), LW.end());
    something_changed();    
  }
  void save_triangulation()
  {
    QString fileName = 
      QFileDialog::getSaveFileName( "triangulation.cgal", 
				    "Cgal files (*.cgal)", this );
    if ( !fileName.isNull() ) {                 // got a file name
      std::ofstream out(fileName);
      CGAL::set_ascii_mode(out);
      out << tr1 << std::endl;    
    }
  }

  void load_triangulation()
  {
    QString s( QFileDialog::getOpenFileName( QString::null,
			    "CGAL files (*.cgal)", this ) );
    if ( s.isEmpty() )
        return;
    tr1.clear();
    A.clear();
    L.clear();
    std::ifstream in(s);
    CGAL::set_ascii_mode(in);
    in >> tr1;
    Vertex_iterator it = tr1.vertices_begin();
    while(it != tr1.vertices_end()) {
      L.push_back((*it).point());
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
    widget->set_window(xmin, xmax, ymin, ymax);
    A.make_alpha_shape(L.begin(), L.end());
    something_changed();
  }

  void slider_changed(int new_value)	{
    alpha_index = 0.0001 * new_value;
    label->setText(QString("The current alpha value: ") + QString::number(alpha_index));
    widget->redraw();   
  }

private:
  CGAL::Qt_widget         *widget;		
  CGAL::Tools_toolbar     *newtoolbar;
  CGAL::Layers_toolbar    *vtoolbar;
  CGAL::Qt_widget_standard_toolbar
                          *stoolbar;
  int                     old_state;
  Qt_layer_show_alpha_shape
                          alpha_shape_layer;
  QSlider                 *slider;
  QLabel                  *label;
  Qt_layer_show_image     *image_layer;
};

#include "alpha_shapes_2.moc"


int
main(int argc, char **argv)
{
  QApplication app( argc, argv );
    app.setStyle( new QPlatinumStyle );
    QPalette p( QColor( 250, 215, 100 ) );
    app.setPalette( p, TRUE );
  MyWindow win(600, 600); // physical window size
  app.setMainWidget(&win);
  win.setCaption(my_title_string);
  win.setMouseTracking(TRUE);
  win.show();
  // because Qt send resizeEvent only on show.
  win.set_window(-1, 1, -1, 1);
  current_state = -1;
  return app.exec();
}

#endif // CGAL_USE_QT
