
#include <CGAL/basic.h>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_handtool.h>
#include <CGAL/IO/pixmaps/zoom_out.xpm>
#include <CGAL/IO/pixmaps/zoom_in.xpm>
#include <CGAL/IO/pixmaps/hand.xpm>

#include <CGAL/IO/leda_window.h>
#include <CGAL/IO/Postscript_file_stream.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/IO/Pm_drawer.h>
#include <CGAL/IO/draw_pm.h>



#include <fstream>
#include <iostream>
#include <sstream>
#include <stack>
#include <set>
#include <string>
#include <list>
#include <vector>
#include <math.h>

#include <qaction.h>
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
#include <qtimer.h>
#include <qtabbar.h>
#include <qtabwidget.h>
#include <qstring.h>
#include <qdialog.h>
#include <qcombobox.h>
#include <qlabel.h> 
#include <qlayout.h> 
#include <qpushbutton.h> 
#include <qfile.h>
#include <qpainter.h>
#include <qprinter.h>
#include <qinputdialog.h> 
#include <qslider.h>
#include <qlcdnumber.h>
#include <qspinbox.h> 

#include "cgal_types1.h"

enum TraitsType { SEGMENT_TRAITS, POLYLINE_TRAITS , CONIC_TRAITS};
enum SnapMode   { NONE , GRID , POINT};
enum Mode       { INSERT , DELETE , POINT_LOCATION , DRAG };

////////////////////////////////////////////////////////////////////////

class OptionsForm : public QDialog
{
    Q_OBJECT
public:
    OptionsForm( QTabWidget * bar = 0 , QWidget* parent = 0 ,int number_of_tabs = 0 , const char* name = "options form",
		 bool modal = FALSE, WFlags f = 0  );
    ~OptionsForm() {}
  
    QLabel *textLabel1;
	QLabel *textLabel2;
    QComboBox *arrComboBox1;
	QComboBox *arrComboBox2;
    QPushButton *okPushButton;
    QPushButton *cancelPushButton;

protected:
    QVBoxLayout *optionsFormLayout;
    QHBoxLayout *arrLayout1;
	QHBoxLayout *arrLayout2;
    QHBoxLayout *buttonsLayout;

private:
	QTabWidget *myBar;
};
///////////////////////////////////////////////////////////////////////////////////
class PropertiesForm : public QDialog
{
    Q_OBJECT
public:
    PropertiesForm( QTabWidget * bar = 0 , QWidget* parent = 0 ,int number_of_tabs = 0 , const char* name = "options form",
		 bool modal = FALSE, WFlags f = 0  );
    ~PropertiesForm() {}
  
    QLabel *textLabel1;
	QLabel *textLabel2;
	QLabel *textLabel3;
	QSpinBox *box1;
	QSpinBox *box2;
	QSpinBox *box3;
    QPushButton *okPushButton;
    QPushButton *cancelPushButton;

protected:
    QVBoxLayout *optionsFormLayout;
    QHBoxLayout *arrLayout1;
	QHBoxLayout *arrLayout2;
	QHBoxLayout *arrLayout3;
    QHBoxLayout *buttonsLayout;

private:
	QTabWidget *myBar;
};
////////////////////////////////////////////////////////////////////////

class OverLay : public QDialog
{
    Q_OBJECT
public:
    OverLay( QTabWidget * bar = 0 , QWidget* parent = 0 ,int number_of_tabs = 0 , const char* name = "OverLay",
		 bool modal = FALSE, WFlags f = 0  );
    ~OverLay() {}
  
    QLabel *textLabel1;
	QLabel *textLabel2;
	QLabel *textLabel3;
    QComboBox *arrComboBox1;
	QComboBox *arrComboBox2;
	QComboBox *arrComboBox3;
    QPushButton *okPushButton;
    QPushButton *cancelPushButton;

protected:
    QVBoxLayout *optionsFormLayout;
    QHBoxLayout *arrLayout1;
	QHBoxLayout *arrLayout2;
	QHBoxLayout *arrLayout3;
    QHBoxLayout *buttonsLayout;

private:
	QTabWidget *myBar;
};
/////////////////////////////////////////////////////////////////////////////////////////////
class Qt_layer_show_ch : public CGAL::Qt_widget_layer
{
public:

    Qt_layer_show_ch( QTabWidget * );
	void draw();
 	
private:
	QTabWidget *myBar;
}; //end class 


////////////////////////////////////////////////////////////////////////////////////////////////
class Qt_widget_demo_tab;

class MyWindow : public QMainWindow
{
    Q_OBJECT
public:
    MyWindow(int w, int h);
	~MyWindow();

private:
    void something_changed();
	void skip_comments( std::ifstream& is, char* one_line );
	void ReadCurve(std::ifstream & is, Pm_base_conic_2 & cv);
	void init_widget(Qt_widget_demo_tab *widget);
	void load( const QString& filename );

private slots:
    void get_new_object(CGAL::Object obj);
    void about();
    void aboutQt();
    void howto();
	void add_segment_tab();
	void add_polyline_tab();
	void add_conic_tab();
	void remove_tab();
	void optionsSetOptions();
	void timer_done();
   	void updateTraitsType( QAction *action );
	void setTraits( TraitsType t );
	void updateSnapMode( QAction *action );
	void setSnapMode( SnapMode m );
	void updateMode( QAction *action );
	void setMode( Mode m );
    void new_instance();
	void update();
	void zoomin();
	void zoomout();
	void fileOpen();
	void fileOpenPolyline();
	void union_arr( int index1 , int index2 , TraitsType t );
	void make_overlay( int index1 , int index2 , int index3 , TraitsType t );
	void overlay();
    void fileSave();
	void fileSave_ps();
    void fileSaveAs();
	void print();
	void properties();

private:
	QTabWidget *myBar;
	int old_state;
	int tab_number;
	int number_of_tabs;
	bool overlay_flag;
	Qt_layer_show_ch *testlayer;
	QAction *setSegmentTraits;
	QAction *setPolylineTraits;
	QAction *setConicTraits;
	QAction *setNoneSnapMode;
	QAction *setGridSnapMode;
	QAction *setPointSnapMode;
	QAction *insertMode;
	QAction *deleteMode;
	QAction *pointLocationMode;
	QAction *dragMode;
	OptionsForm optionsForm;
	QString m_filename;
	int m_width, m_height; 
};
/////////////////////////////////////////////////////////////////////////////
class Qt_widget_demo_tab : public CGAL::Qt_widget
{
public:
	Qt_widget_demo_tab(QWidget *parent = 0, const char *name = 0, TraitsType t = SEGMENT_TRAITS):
	    CGAL::Qt_widget( parent , name ),
		snap_mode( NONE ),
		mode( INSERT ),
		traits_type(t)
	{}

	~Qt_widget_demo_tab() {}

	virtual void draw() = 0;
	virtual void mousePressEvent(QMouseEvent *e) = 0;
	virtual void mouseMoveEvent(QMouseEvent *e) = 0;
	virtual void leaveEvent(QEvent *e) = 0;
	virtual void remove_segment(QMouseEvent *e) = 0;
	void mousePressEvent_point_location(QMouseEvent *e);
	Coord_type getMid(Coord_type coord, int my_min, int my_max);
    bool is_pure(Qt::ButtonState s);
	Coord_segment convert(Pm_pol_point_2 & source , Pm_pol_point_2 & target);
	Coord_segment convert(const Pm_pol_point_2 & source , const Pm_pol_point_2 & target);
	Coord_point get_point(Coord_type x , Coord_type y);
	Coord_type dist(Coord_type x1, Coord_type y1, Coord_type x2, Coord_type y2);

    int         current_state;
	int         index;
	Coord_point pl_point;
	SnapMode    snap_mode;
	Mode        mode;
	TraitsType  traits_type;
	int         m_line_width;
	int m_xmin, m_xmax, m_ymin, m_ymax;

};
////////////////////////////////////////////////////////////////////////

class Qt_widget_segment_tab : public Qt_widget_demo_tab
{
public:

	Qt_widget_segment_tab(QWidget *parent = 0, const char *name = 0):
	    Qt_widget_demo_tab( parent , name , SEGMENT_TRAITS ),
		first_point_segment(false)
	{}

	~Qt_widget_segment_tab();

	virtual void draw();
	virtual void mousePressEvent(QMouseEvent *e);
	virtual void mouseMoveEvent(QMouseEvent *e);
	virtual void leaveEvent(QEvent *e);
	virtual void remove_segment(QMouseEvent *e);
    const Pm_base_seg_2* get_origin_curve(const Pm_xseg_2 & xseg );

	Pm_seg_list list_of_curves;
	Seg_arr segment_arr;
	
	bool first_point_segment; //true if the user left clicked once
    bool first_time_segment;  //true if the line is not drawn

	Coord_point m_p1,m_p2;
	 // from Qt_widget_get_segment:
	Coord_type m_x1; //the X of the first point
	Coord_type m_y1; //the Y of the first point
	Coord_type m_x2; //the old second point's X
    Coord_type m_y2; //the old second point's Y
};
////////////////////////////////////////////////////////////////////////

class Qt_widget_polyline_tab : public Qt_widget_demo_tab
{
public:

	Qt_widget_polyline_tab(QWidget *parent = 0, const char *name = 0) :
	    Qt_widget_demo_tab( parent , name , POLYLINE_TRAITS ),
		active(false),
		first_time_polyline(true)
	{}

	~Qt_widget_polyline_tab();
	
	virtual void draw();
	virtual void mousePressEvent(QMouseEvent *e);
	virtual void mouseMoveEvent(QMouseEvent *e);
	virtual void leaveEvent(QEvent *e);
	virtual void remove_segment(QMouseEvent *e);
	void get_polyline();
    const Pm_base_pol_2* get_origin_curve(const Pm_xpol_2 & xseg );
    bool compare(const Pm_base_pol_2 * org_seg , const Pm_pol_2 * closest);

 	Pm_pol_list list_of_curves;
	Pol_arr polyline_arr;

	bool active;              // true if the first point was inserted
    bool first_time_polyline; // true if it is the first time when draw the rubber band

	Coord_point rubber;				  // the new point of the rubber band
    Coord_point last_of_poly;		  // the last point of the polygon
    Coord_point rubber_old;			 // the old point of the rubber band
	std::vector<Pm_pol_point_2> points;

};

////////////////////////////////////////////////////////////////////////

class Qt_widget_conic_tab : public Qt_widget_demo_tab
{
public:
	Qt_widget_conic_tab(QWidget *parent = 0, const char *name = 0) :
	    Qt_widget_demo_tab( parent , name , CONIC_TRAITS )
	{}

	~Qt_widget_conic_tab();
	
	virtual void draw();
	virtual void mousePressEvent(QMouseEvent *e);
	virtual void mouseMoveEvent(QMouseEvent *e) {};
	virtual void leaveEvent(QEvent *e) {};
	virtual void remove_segment(QMouseEvent *e) {};
	Pm_base_conic_2* get_origin_curve( Pm_xconic_2 & xseg );
	void draw_curve(const Pm_xconic_2& c);

	Pm_xconic_list list_of_xcurves;
	Conic_arr conic_arr;
		
};


////////////////////////////////////////////////////////////////////////
OptionsForm::OptionsForm(  QTabWidget * bar, QWidget* parent ,int number_of_tabs , const char* name, bool modal, WFlags f  ): 
	QDialog( parent, name, modal, f ),
	myBar(bar)
{
    setCaption( "Union -- Options" );
    resize( 320, 290 );

    optionsFormLayout = new QVBoxLayout( this, 11, 6 );

    arrLayout1 = new QHBoxLayout( 0, 0, 6 );

    textLabel1 = new QLabel( "Arr 1", this );
    arrLayout1->addWidget( textLabel1 );

    arrComboBox1 = new QComboBox( FALSE, this );

	for (int i=0; i < number_of_tabs; i++)
		arrComboBox1->insertItem( myBar->label(i) );
		
	arrLayout1->addWidget( arrComboBox1 );
    optionsFormLayout->addLayout( arrLayout1 );

	arrLayout2 = new QHBoxLayout( 0, 0, 6 );

    textLabel2 = new QLabel( "Arr 2", this );
    arrLayout2->addWidget( textLabel2 );

    arrComboBox2 = new QComboBox( FALSE, this );
    
	for (int i=0; i < number_of_tabs; i++)
		arrComboBox2->insertItem( myBar->label(i) );
	
	arrLayout2->addWidget( arrComboBox2 );
    optionsFormLayout->addLayout( arrLayout2 );

	buttonsLayout = new QHBoxLayout( 0, 0, 6 );

    okPushButton = new QPushButton( "OK", this );
    okPushButton->setDefault( TRUE );
    buttonsLayout->addWidget( okPushButton );

    cancelPushButton = new QPushButton( "Cancel", this );
    buttonsLayout->addWidget( cancelPushButton );
    optionsFormLayout->addLayout( buttonsLayout );

    connect( okPushButton, SIGNAL( clicked() ), this, SLOT( accept() ) );
    connect( cancelPushButton, SIGNAL( clicked() ), this, SLOT( reject() ) );

    textLabel1->setBuddy( arrComboBox1 );
	textLabel2->setBuddy( arrComboBox2 );
   
}
////////////////////////////////////////////////////////////////////////
PropertiesForm::PropertiesForm(  QTabWidget * bar, QWidget* parent ,int number_of_tabs , const char* name, bool modal, WFlags f  ): 
	QDialog( parent, name, modal, f ),
	myBar(bar)
{
    setCaption( "Properties -- Options" );
    resize( 320, 290 );

    optionsFormLayout = new QVBoxLayout( this, 11, 6 );

    arrLayout1 = new QHBoxLayout( 0, 0, 6 );
    textLabel1 = new QLabel( "Width", this );
    arrLayout1->addWidget( textLabel1 );
	box1 = new QSpinBox( 300, 1000, 50, this, "box1" );
	box1->setValue(700);
	arrLayout1->addWidget( box1 );
    optionsFormLayout->addLayout( arrLayout1 );

	arrLayout2 = new QHBoxLayout( 0, 0, 6 );
    textLabel2 = new QLabel( "Hight", this );
    arrLayout2->addWidget( textLabel2 );
	box2 = new QSpinBox( 300, 1000, 50, this, "box2" );
	box2->setValue(700);
	arrLayout2->addWidget( box2 );
	optionsFormLayout->addLayout( arrLayout2 );

	arrLayout3 = new QHBoxLayout( 0, 0, 6 );
    textLabel3 = new QLabel( "Line Width", this );
    arrLayout3->addWidget( textLabel3 );
	box3 = new QSpinBox( 1, 5, 1, this, "box3" );
	box3->setValue(2);
	arrLayout3->addWidget( box3 );
	optionsFormLayout->addLayout( arrLayout3 );
    
	buttonsLayout = new QHBoxLayout( 0, 0, 6 );
    okPushButton = new QPushButton( "OK", this );
    okPushButton->setDefault( TRUE );
    buttonsLayout->addWidget( okPushButton );
    cancelPushButton = new QPushButton( "Cancel", this );
    buttonsLayout->addWidget( cancelPushButton );
    optionsFormLayout->addLayout( buttonsLayout );
    connect( okPushButton, SIGNAL( clicked() ), this, SLOT( accept() ) );
    connect( cancelPushButton, SIGNAL( clicked() ), this, SLOT( reject() ) );

    textLabel1->setBuddy( box1 );
	textLabel2->setBuddy( box2 );
	textLabel3->setBuddy( box3 );
   
}
////////////////////////////////////////////////////////////////////////
OverLay::OverLay(  QTabWidget * bar, QWidget* parent ,int number_of_tabs , const char* name, bool modal, WFlags f  ): 
	QDialog( parent, name, modal, f ),
	myBar(bar)
{
    setCaption( "Union -- Options" );
    resize( 320, 290 );

    optionsFormLayout = new QVBoxLayout( this, 11, 6 );

    arrLayout1 = new QHBoxLayout( 0, 0, 6 );

    textLabel1 = new QLabel( "Arr 1", this );
    arrLayout1->addWidget( textLabel1 );

    arrComboBox1 = new QComboBox( FALSE, this );

	for (int i=0; i < number_of_tabs; i++)
		arrComboBox1->insertItem( myBar->label(i) );
		
	arrLayout1->addWidget( arrComboBox1 );
    optionsFormLayout->addLayout( arrLayout1 );

	arrLayout2 = new QHBoxLayout( 0, 0, 6 );

    textLabel2 = new QLabel( "Arr 2", this );
    arrLayout2->addWidget( textLabel2 );

    arrComboBox2 = new QComboBox( FALSE, this );
    
	for (int i=0; i < number_of_tabs; i++)
		arrComboBox2->insertItem( myBar->label(i) );
	
	arrLayout2->addWidget( arrComboBox2 );
    optionsFormLayout->addLayout( arrLayout2 );

	arrLayout3 = new QHBoxLayout( 0, 0, 6 );

    textLabel3 = new QLabel( "Arr 3", this );
    arrLayout3->addWidget( textLabel3 );

    arrComboBox3 = new QComboBox( FALSE, this );

	for (int i=0; i < number_of_tabs; i++)
		arrComboBox3->insertItem( myBar->label(i) );
		
	arrLayout3->addWidget( arrComboBox3 );
    optionsFormLayout->addLayout( arrLayout3 );

	buttonsLayout = new QHBoxLayout( 0, 0, 6 );

    okPushButton = new QPushButton( "OK", this );
    okPushButton->setDefault( TRUE );
    buttonsLayout->addWidget( okPushButton );

    cancelPushButton = new QPushButton( "Cancel", this );
    buttonsLayout->addWidget( cancelPushButton );
    optionsFormLayout->addLayout( buttonsLayout );

    connect( okPushButton, SIGNAL( clicked() ), this, SLOT( accept() ) );
    connect( cancelPushButton, SIGNAL( clicked() ), this, SLOT( reject() ) );

    textLabel1->setBuddy( arrComboBox1 );
	textLabel2->setBuddy( arrComboBox2 );
	textLabel3->setBuddy( arrComboBox3 );
   
}
////////////////////////////////////////////////////////////////////////
#ifndef CGAL_USE_QT
#include <iostream>
int main(int, char*)
{

  std::cout << "Sorry, this demo needs QT...";
  std::cout << std::endl;

  return 0;
}

#else

//#include "demo1.h"


#include <CGAL/IO/pixmaps/movepoint.xpm>
#include <CGAL/IO/pixmaps/point.xpm>
#include <CGAL/IO/pixmaps/line.xpm>
#include <CGAL/IO/pixmaps/arrow.xpm>
#include <CGAL/IO/pixmaps/alpha_shape.xpm>
#include "polyline.xpm"
#include "insert.xpm"
#include "delete.xpm"
#include "grid.xpm"
#include "conic.xpm"
#include "none.xpm"
#include "po.xpm"

const QString my_title_string("Arrangement Demo with CGAL Qt_widget");

//////////////////////////////////////////////////////////////////////////////
	
Qt_layer_show_ch::Qt_layer_show_ch( QTabWidget * bar ) :
	myBar(bar)
{}

void Qt_layer_show_ch::draw()
{
	// We peform downcasting from QWigdet* to Qt_widget_demo_tab*, as we know that only
	// Qt_widget_demo_tab objects are stored in the tab pages.
	Qt_widget_demo_tab	*w_demo_p = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());

	w_demo_p->lock();

	if (w_demo_p->snap_mode == GRID)
	{
		//(*w_demo_p) << CGAL::GRAY;
		(*w_demo_p) << CGAL::DEEPBLUE;
		//(*w_demo_p) << CGAL::RED;
		//(*w_demo_p) << CGAL::ORANGE;
		//(*w_demo_p) << CGAL::PURPLE;

		(*w_demo_p) << CGAL::LineWidth(1);
		// get the edge coordinate
		int min_x = static_cast<int> (w_demo_p->x_min());
		int max_x = static_cast<int> (w_demo_p->x_max());
		int min_y = static_cast<int> (w_demo_p->y_min());
		int max_y = static_cast<int> (w_demo_p->y_max());

		// calculate cube size (minimum of 1)
		int cube_size_x = std::max(1, abs(max_x - min_x)/20);
		int cube_size_y = std::max(1, abs(max_y - min_y)/20);
		// draw the grid lines
		for (int i = min_x; i <= max_x; i += cube_size_x)
			(*w_demo_p) << Coord_segment(Coord_point( i , max_y + cube_size_y),Coord_point( i , min_y - cube_size_y));
		for (int i = min_y; i <= max_y; i += cube_size_y)
			(*w_demo_p) << Coord_segment(Coord_point( max_x + cube_size_x , i ),Coord_point( min_x - cube_size_x , i ));
	}
					
	w_demo_p->draw();

	w_demo_p->unlock();
}

  ////////////////////////////////////////////////////////////////////////////////////////////
Qt_widget_segment_tab::~Qt_widget_segment_tab()
{
	Pm_seg_iter itp;
	for (itp = list_of_curves.begin(); itp != list_of_curves.end(); ++itp)
		delete (*itp);
}

void Qt_widget_segment_tab::draw()
{
    (*this) << CGAL::GREEN;
    (*this) << CGAL::LineWidth(m_line_width);
	Pm_seg_const_iter itp;
	for (itp = list_of_curves.begin(); itp != list_of_curves.end(); ++itp)
	{
      (*this) << *(*itp);
    }

    if (mode == POINT_LOCATION && 
		! (segment_arr.halfedges_begin() == 
		   segment_arr.halfedges_end() ) ) 
    {
      (*this) << CGAL::LineWidth(3);
      (*this) << CGAL::YELLOW;

      Seg_locate_type lt;
      Pm_seg_point_2 temp_p (pl_point.x(), pl_point.y());
      Seg_halfedge_handle e = segment_arr.locate(temp_p, lt);
	
	  //std::cout << lt << std::endl;
      //color the face on the screen
      Seg_arr::Face_handle f = e->face();
	
      if (f->does_outer_ccb_exist()) // its an inside face
      {
		Seg_arr::Ccb_halfedge_circulator cc=f->outer_ccb();
		do {
			(*this) << cc->curve();
		} while (++cc != f->outer_ccb());
		
	  }

	  Seg_arr::Holes_iterator hit, eit = f->holes_end();
	  for (hit = f->holes_begin(); hit != eit; ++hit) 
	  {
		Seg_arr::Ccb_halfedge_circulator cc = *hit; 
		do 
		{
			(*this) << cc->curve();
			cc++;
		} while (cc != *hit);
	  }
	  (*this) << CGAL::LineWidth(m_line_width);
    }
  }

  
  void Qt_widget_segment_tab::mousePressEvent(QMouseEvent *e)
  {
	if (mode == POINT_LOCATION)
	{
		mousePressEvent_point_location( e );
		return;
	}
	if (mode == DELETE)
	{
		remove_segment( e );
		return;
	}
	if(e->button() == Qt::LeftButton 
       && !first_point_segment
       && is_pure(e->state()))
	{
        Coord_type x, y;
		x_real(e->x(), x);
		y_real(e->y(), y);
		m_p1 = get_point(x,y);
		m_p2 = m_p1;

		/*if (snap_mode == GRID)
		{
			m_x1 = getMid(x,static_cast<int> (x_min()),static_cast<int> (x_max()));
			m_y1 = getMid(y,static_cast<int> (y_min()),static_cast<int> (y_max()));
			m_x2 = getMid(x,static_cast<int> (x_min()),static_cast<int> (x_max()));
			m_y2 = getMid(y,static_cast<int> (y_min()),static_cast<int> (y_max()));

		}
		else
		{
			m_x1 = x;
			m_y1 = y;
			m_x2 = x;
			m_y2 = y;
		}*/
		first_point_segment = TRUE;
    } 
	else if((e->button() == Qt::LeftButton) && is_pure(e->state()))
	{
		Coord_type x, y;
		x_real(e->x(), x);
		y_real(e->y(), y);
		Coord_point p = get_point(x,y);
		
		/*if (snap_mode == GRID)
		{
			Coord_type x3, y3;
			x_real(e->x(), x3);
			y_real(e->y(), y3);
			x = getMid(x3, static_cast<int> (x_min()),static_cast<int> (x_max()));
			y = getMid(y3, static_cast<int> (y_min()),static_cast<int> (y_max()));
		}
		else
		{
			x_real(e->x(), x);
			y_real(e->y(), y);
		}*/
		if(m_p1.x() != p.x() || m_p1.y() != p.y()) 
		{
	        new_object( make_object( Coord_segment( m_p1 , p ) ) );
			//make_object(Coord_segment(Coord_point(m_x1, m_y1), Coord_point(x,y))));
			first_point_segment = FALSE;
		}    
    }
	return;
  }
 

void Qt_widget_segment_tab::remove_segment(QMouseEvent *e)
{
	if( list_of_curves.empty() )
		return;

	Coord_point p(x_real(e->x()) ,y_real(e->y()));

	bool is_first = true;
	Coord_type min_dist = 0;

	Pm_seg_iter itp;
	Pm_seg_iter it_seg = NULL;
	for (itp = list_of_curves.begin(); itp != list_of_curves.end(); ++itp)
	{
		const Pm_seg_2 & pm_seg = *(*itp);
		const Pm_seg_point_2 & source = pm_seg.source();
		const Pm_seg_point_2 & target = pm_seg.target();

		Coord_type x1 = CGAL::to_double(source.x());
		Coord_type y1 = CGAL::to_double(source.y());

		Coord_type x2 = CGAL::to_double(target.x());
		Coord_type y2 = CGAL::to_double(target.y());

		Coord_point coord_source(x1 , y1);
		Coord_point coord_target(x2 , y2);
		Coord_segment coord_seg(coord_source, coord_target);
		Coord_type dist = CGAL::squared_distance( p, coord_seg);

		if (is_first || dist < min_dist)
		{
			min_dist = dist;
			it_seg = itp;
			is_first = false;
		}
	}
	    
	// Remove curve from pmwx
	Seg_pm::Halfedge_iterator hei;

	typedef std::list<Seg_pm::Halfedge_iterator> Hafledge_list;
	
	Hafledge_list halfedge_list;
	Hafledge_list::iterator result;

	for (hei = segment_arr.halfedges_begin(); hei != segment_arr.halfedges_end(); ++hei) 
	{
        const Pm_xseg_2 & xseg = hei->curve();
		const Pm_base_seg_2 * org_seg = get_origin_curve(xseg);

		if ((*org_seg).source() == (**it_seg).source() && (*org_seg).target() == (**it_seg).target())
		{
		    result = std::find(halfedge_list.begin(), halfedge_list.end(), hei);
			if (result == halfedge_list.end())
			{
				halfedge_list.push_back(hei->twin());
				segment_arr.remove_edge(hei);
			}
		}				
	}

	list_of_curves.erase(it_seg);
	delete (*it_seg);

	redraw();
	
  } // remove_segment

const Pm_base_seg_2* Qt_widget_segment_tab::get_origin_curve(const Pm_xseg_2 & xseg )
{
	const Curve_data & curve_data = xseg.get_data();
	if (curve_data.m_type == Curve_data::LEAF) 
        return curve_data.m_ptr.m_curve;
    else // curve_data.m_type == Curve_data::INTERNAL
	{
		const Pm_xseg_2* xseg_p = static_cast<const Pm_xseg_2*>(curve_data.m_ptr.m_x_motonote_curve);
		return get_origin_curve( *xseg_p );
	}
}

  void Qt_widget_segment_tab::mouseMoveEvent(QMouseEvent *e)
  {
    if(first_point_segment)
    {
      Coord_type x, y;
      x_real(e->x(), x);
      y_real(e->y(), y);
      RasterOp old_raster = rasterOp();//save the initial raster mode
      QColor old_color = color();
      setRasterOp(XorROP);
      lock();
      *this << CGAL::GREEN;
      if(!first_time_segment)
     // *this << Coord_segment(Coord_point(m_x1, m_y1), Coord_point(m_x2, m_y2));
     // *this << Coord_segment(Coord_point(m_x1, m_y1), Coord_point(x, y));
      *this << Coord_segment( m_p1 , m_p2 );
      *this << Coord_segment( m_p1 , Coord_point(x, y));

	  unlock();
      setRasterOp(old_raster);
      setColor(old_color);

      //save the last coordinates to redraw the screen
	  m_p2 = Coord_point(x,y);
      //m_x2 = x;
      //m_y2 = y;
      first_time_segment = false;
    }
  }

    void Qt_widget_segment_tab::leaveEvent(QEvent *e)
  {
    if(first_point_segment)
    {
      RasterOp old_raster = rasterOp();//save the initial raster mode
      QColor old_color = color();
	  lock();
      setRasterOp(XorROP);
      *this << CGAL::GREEN;
      *this << Coord_segment(m_p1 , m_p2);
      setRasterOp(old_raster);
      setColor(old_color);
      unlock();
      first_time_segment = true;
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
Qt_widget_polyline_tab::~Qt_widget_polyline_tab()
{
	Pm_pol_iter itp;
	for (itp = list_of_curves.begin(); itp != list_of_curves.end(); ++itp)
		delete (*itp);
}

  void Qt_widget_polyline_tab::draw()
  {
	*this << CGAL::GREEN;
    *this << CGAL::LineWidth(m_line_width);
	Pm_pol_iter it = list_of_curves.begin();
    for (it  = list_of_curves.begin(); it != list_of_curves.end(); ++it)
	{
		Pm_pol_2 pol = **it;
        Pm_pol_2::const_iterator ps = pol.begin();
        Pm_pol_2::const_iterator pt = ps; pt++;

        while (pt != pol.end()) {
			const Pm_pol_point_2 & source = *ps;
			const Pm_pol_point_2 & target = *pt;
			Coord_segment coord_seg = convert(source , target);
			*this << coord_seg;
			ps++; pt++;
		}
	}

	// draw Point Location
    if(mode == POINT_LOCATION && !(polyline_arr.halfedges_begin() == polyline_arr.halfedges_end()) ) 
    {
	  *this << CGAL::LineWidth(3);
      *this << CGAL::YELLOW;

	  // PolyLine Arrangement
      Pol_locate_type ltp;
      Pm_pol_point_2 temp_pp(pl_point.x(), pl_point.y());
      Pol_halfedge_handle ep = polyline_arr.locate(temp_pp, ltp);
	
      //color the face on the screen
      Pol_face_handle fp = ep->face();
	
      if (fp->does_outer_ccb_exist()) 
      {
		Pol_ccb_halfedge_circulator cc = fp->outer_ccb();
		do {
			Pm_pol_2::const_iterator points_iter;
			for (points_iter = cc->curve().begin();
				points_iter != cc->curve().end(); ) 
				{
					Pm_pol_2::const_iterator source_point = points_iter;
					points_iter++;
    				Pm_pol_2::const_iterator target_point = points_iter;

					if (target_point == cc->curve().end())
						break;
    
					const Pm_pol_point_2 & source = *source_point;
					const Pm_pol_point_2 & target = *target_point;
					Coord_segment coord_seg = convert(source , target);
					*this << coord_seg;
				} // for
			} while (++cc != fp->outer_ccb());
	
	  }// if

      Pol_holes_iterator hitp = fp->holes_begin(), eitp = fp->holes_end();
      for (;hitp!=eitp; ++hitp) 
      {
		Pol_ccb_halfedge_circulator cc = *hitp; 
		do 
		{
			Pm_pol_2::const_iterator points_iter;
			for (points_iter = cc->curve().begin();
				   points_iter != cc->curve().end(); ) 
			{
				Pm_pol_2::const_iterator source_point = points_iter;
				points_iter++;
    			Pm_pol_2::const_iterator target_point = points_iter;
    
				if (target_point == cc->curve().end())
					break;
 
				const Pm_pol_point_2 & source = *source_point;
				const Pm_pol_point_2 & target = *target_point;
				Coord_segment coord_seg = convert(source , target);

				*this << coord_seg;
			}// for
		} while (++cc != *hitp);
      }		

	  *this << CGAL::LineWidth(m_line_width);
    }

  }

  void Qt_widget_polyline_tab::mousePressEvent(QMouseEvent *e)
  {
	if (mode == POINT_LOCATION)
	{
		mousePressEvent_point_location( e );
		return;
	}
	if (mode == DELETE)
	{
		remove_segment( e );
		return;
	}
	  if(e->button() == Qt::LeftButton && is_pure(e->state()))
      {
		Coord_type x, y;
		x_real(e->x(), x);
		y_real(e->y(), y);
		Coord_point p = get_point(x,y);

		//Pm_pol_point_2 p;
		//x_real(e->x(), x);
		//y_real(e->y(), y);
		//if (snap_mode == GRID)
		//	p = Pm_pol_point_2(getMid(x, static_cast<int> (x_min()), static_cast<int> (x_max())), 
		//					getMid(y, static_cast<int> (y_min()), static_cast<int> (y_max())));
		//else
		//	p = Pm_pol_point_2(x,y);

		if(!active)
		{
			active = true;
			last_of_poly = p;
			//poly.push_back(p);	
			points.push_back(Pm_pol_point_2(p.x(),p.y()));	
		} 
		else
		{
			if (last_of_poly == p) return;
			rubber_old = p;
		
			//poly.push_back(p);	
			points.push_back(Pm_pol_point_2(p.x(),p.y()));	
			//show the last rubber as edge of the polygon
			lock();
			RasterOp old_rasterop=rasterOp();
			get_painter().setRasterOp(XorROP);
			Coord_segment seg = Coord_segment(rubber, last_of_poly);
			 *this << CGAL::WHITE;
			 *this << seg;
			 *this << CGAL::GREEN;
			 *this << seg;
			setRasterOp(old_rasterop);
			unlock();
			last_of_poly = p;
		
		}
		return;
	  }
	  // finish polyline draw with right button click 
	  else if (active && e->button() == Qt::RightButton && is_pure(e->state())) {
        //new_object(make_object(poly));
		get_polyline();
        active = false;
        first_time_polyline = true;
        //poly.erase(poly.vertices_begin(), poly.vertices_end());
		points.clear();
        redraw();
      }
  }

  void Qt_widget_polyline_tab::get_polyline()
  {
      Pm_base_pol_2 * base_pol_p = new Pm_base_pol_2(points.begin(), points.end());
	  Curve_pol_data cd;
	  cd.m_type = Curve_pol_data::LEAF;
	  cd.m_index = index;
	  cd.m_ptr.m_curve = base_pol_p;
	  Pm_pol_2 * pol = new Pm_pol_2( *base_pol_p, cd );
	  //Pm_pol_2 *poly = new Pm_pol_2( Pm_base_pol_2( points.begin(), points.end()) , index );
	  //Pm_pol_2 *poly = new Pm_pol_2( points.begin(), points.end() );
	  list_of_curves.push_back(pol);
	  polyline_arr.insert( *pol );
	  redraw();
  }

void Qt_widget_polyline_tab::remove_segment(QMouseEvent *e)
{

      if( list_of_curves.empty() )
		return;

      Coord_type x=static_cast<Coord_type>(x_real(e->x()));
      Coord_type y=static_cast<Coord_type>(y_real(e->y()));

      Coord_point p(x,y);
      Coord_type min_dist=100000000;
     
	  Pm_pol_iter it_closest=NULL;
      Pm_pol_iter pit = list_of_curves.begin();

      Coord_segment closest_segment;

      while(pit!=list_of_curves.end())
      {
		Pm_pol_2 pol = **pit;
        Pm_pol_2::const_iterator ps = pol.begin();
        Pm_pol_2::const_iterator pt = ps; pt++;

        while (pt != pol.end()) {
			const Pm_pol_point_2 & source = *ps;
			const Pm_pol_point_2 & target = *pt;
			Coord_segment coord_seg = convert(source , target);
			Coord_type dist = CGAL::squared_distance( p, coord_seg);
			if( dist < min_dist)
			{
				min_dist = dist;
				it_closest = pit;
				closest_segment = coord_seg;
			}
			ps++; pt++;
		}
		pit++;
      }
    
	// Remove curve from pmwx
	Pol_pm::Halfedge_iterator hei;

	std::list<Pol_pm::Halfedge_iterator> halfedge_list;
	std::list<Pol_pm::Halfedge_iterator>::iterator result;

	const Pm_pol_2 * closest = *it_closest;

	for (hei = polyline_arr.halfedges_begin(); hei != polyline_arr.halfedges_end(); ++hei) 
	{
        const Pm_xpol_2 & xpol = hei->curve();
		const Pm_base_pol_2 * org_seg = get_origin_curve(xpol);
		if (compare( org_seg , closest ) )
		{
		    result = std::find(halfedge_list.begin(), halfedge_list.end(), hei);
			if (result == halfedge_list.end())
			{
				halfedge_list.push_back(hei->twin());
				polyline_arr.remove_edge(hei);
			}
		}
	}

    list_of_curves.erase( it_closest );
	delete (*it_closest);

    redraw();
}
  bool Qt_widget_polyline_tab::compare(const Pm_base_pol_2 * org_seg , const Pm_pol_2 * closest)
  {
	  Pm_base_pol_2::const_iterator ito = (*org_seg).begin();
	  Pm_pol_2::const_iterator itc = (*closest).begin();
	  while (ito != (*org_seg).end() && itc != (*closest).end())
	  {
		  if (*ito != *itc)
			  return false;
		  ito++;
		  itc++;
	  }
	  if (ito == (*org_seg).end() && itc == (*closest).end())
          return true;
	  else
		  return false;
  }
void Qt_widget_polyline_tab::mouseMoveEvent(QMouseEvent *e)
  {
    if (active)
    {
      Coord_type x, y;
      x_real(e->x(), x);
      y_real(e->y(), y);
	  rubber = get_point(x,y);

	  /*if (snap_mode == GRID)
          rubber = Pm_pol_point_2(getMid(x,static_cast<int> (x_min()),static_cast<int> (x_max())), 
							   getMid(y,static_cast<int> (y_min()),static_cast<int> (y_max())));
	  else
          rubber = Pm_pol_point_2(x, y);*/
      lock();
      RasterOp old_rasterop=rasterOp();
      get_painter().setRasterOp(XorROP);
      *this << CGAL::WHITE;    
      if(!first_time_polyline)
	  {
		  Coord_segment seg = Coord_segment(rubber_old, last_of_poly);
          *this << seg;
	  }
	  Coord_segment seg2 = Coord_segment(rubber, last_of_poly);
      *this << seg2;
      first_time_polyline = false;
      rubber_old = rubber;
      setRasterOp(old_rasterop);
      unlock();
    }
  }

  void Qt_widget_polyline_tab::leaveEvent(QEvent *)
  {
    if (active)
    {
      lock();
      RasterOp old_rasterop=rasterOp();
      get_painter().setRasterOp(XorROP);
      *this << CGAL::WHITE;
      Coord_segment seg = Coord_segment(rubber_old, last_of_poly);
      *this << seg;
      setRasterOp(old_rasterop);
      unlock();
      first_time_polyline = true;
    }
  }

const Pm_base_pol_2* Qt_widget_polyline_tab::get_origin_curve(const Pm_xpol_2 & xseg )
{
	const Curve_pol_data & curve_data = xseg.get_data();
	if (curve_data.m_type == Curve_pol_data::LEAF) 
		return curve_data.m_ptr.m_curve;
	else // curve_data.m_type == Curve_data::INTERNAL
	{
		const Pm_xpol_2* xseg_p = static_cast<const Pm_xpol_2*>(curve_data.m_ptr.m_x_motonote_curve);
		return get_origin_curve( *xseg_p );
	}
}
//////////////////////////////////////////////////////////////////////////////
Qt_widget_conic_tab::~Qt_widget_conic_tab()
{
	Pm_xconic_iter itp;
	for (itp = list_of_xcurves.begin(); itp != list_of_xcurves.end(); ++itp)
		delete (*itp);
}

void Qt_widget_conic_tab::draw_curve(const Pm_xconic_2& c)
{
	// Get the co-ordinates of the curve's source and target.
	double sx = CGAL::to_double(c.source().x()),
		sy = CGAL::to_double(c.source().y()),
		tx = CGAL::to_double(c.target().x()),
		ty = CGAL::to_double(c.target().y());

	if (c.is_segment())
	{
		Coord_point coord_source(sx , sy);
		Coord_point coord_target(tx , ty);
		Coord_segment coord_seg(coord_source, coord_target);
		
		*this << coord_seg;
	}
    else
	{
		// If the curve is monotone, than its source and its target has the
		// extreme x co-ordinates on this curve.
		if (c.is_x_monotone())
		{
			
			bool     is_source_left = (sx < tx);
			int      x_min = is_source_left ? (*this).x_pixel(sx) : 
											(*this).x_pixel(tx);
			int      x_max = is_source_left ? (*this).x_pixel(tx) :
												(*this).x_pixel(sx);
			double   prev_x = is_source_left ? sx : tx;
			double   prev_y = is_source_left ? sy : ty;
			double   end_x = is_source_left ? tx : sx;
			double   end_y = is_source_left ? ty : sy;
			double   curr_x, curr_y;
			int      x;

			Pm_conic_point_2 ps[2];
			int nps;

			for (x = x_min + 1; x < x_max; x++)
			{
				curr_x = (*this).x_real(x);
				nps = c.get_points_at_x(Pm_conic_point_2(curr_x, 0), ps);
				if (nps == 1)
				{
					curr_y = CGAL::to_double(ps[0].y());
					(*this) << Coord_segment( Coord_point(prev_x, prev_y) , Coord_point(curr_x, curr_y) );
					prev_x = curr_x;
					prev_y = curr_y;
						
				}
			}

			(*this) << Coord_segment( Coord_point(prev_x, prev_y) , Coord_point(end_x, end_y) );
		}
		else
		{
			// We should never reach here.
			CGAL_assertion(false);
		}
	}
}
		
void Qt_widget_conic_tab::draw()
{
    (*this) << CGAL::GREEN;
    (*this) << CGAL::LineWidth(m_line_width);
	Pm_xconic_const_iter itp;
	for (itp = list_of_xcurves.begin(); itp != list_of_xcurves.end(); ++itp)
	{
      draw_curve(**itp);
    }

    if (mode == POINT_LOCATION && 
		! (conic_arr.halfedges_begin() == 
		   conic_arr.halfedges_end() ) ) 
    {
      (*this) << CGAL::LineWidth(3);
      (*this) << CGAL::YELLOW;

      Conic_locate_type lt;
      Pm_conic_point_2 temp_p (pl_point.x(), pl_point.y());
      Conic_halfedge_handle e = conic_arr.locate(temp_p, lt);
	
	  //std::cout << lt << std::endl;
      //color the face on the screen
      Conic_arr::Face_handle f = e->face();
	
      if (f->does_outer_ccb_exist()) // its an inside face
      {
		Conic_arr::Ccb_halfedge_circulator cc=f->outer_ccb();
		do {
			draw_curve( cc->curve());
		} while (++cc != f->outer_ccb());
		
	  }

	  Conic_arr::Holes_iterator hit, eit = f->holes_end();
	  for (hit = f->holes_begin(); hit != eit; ++hit) 
	  {
		Conic_arr::Ccb_halfedge_circulator cc = *hit; 
		do 
		{
			draw_curve( cc->curve());
			cc++;
		} while (cc != *hit);
	  }
	  (*this) << CGAL::LineWidth(m_line_width);
    }
  }	

void Qt_widget_conic_tab::mousePressEvent(QMouseEvent *e)
{
	if (mode == POINT_LOCATION)
	{
		mousePressEvent_point_location( e );
		return;
	}
}

Pm_base_conic_2* Qt_widget_conic_tab::get_origin_curve(Pm_xconic_2 & xseg )
{
	const Curve_conic_data & curve_data = xseg.get_data();
	if (curve_data.m_type == Curve_conic_data::LEAF) 
        return curve_data.m_ptr.m_curve;
    else // curve_data.m_type == Curve_data::INTERNAL
	{
		Pm_xconic_2* xseg_p = static_cast<Pm_xconic_2*>(curve_data.m_ptr.m_x_motonote_curve);
		return get_origin_curve( *xseg_p );
	}
}
//////////////////////////////////////////////////////////////////////////////

  void Qt_widget_demo_tab::mousePressEvent_point_location(QMouseEvent *e)
  {
    if(e->button() == Qt::LeftButton
       && is_pure(e->state()))
    {
      Coord_type x, y;
      x_real(e->x(), x);
      y_real(e->y(), y);
      new_object(make_object(Coord_point(x, y)));
    }
  }

  bool Qt_widget_demo_tab::is_pure(Qt::ButtonState s){
    if((s & Qt::ControlButton) ||
       (s & Qt::ShiftButton) ||
       (s & Qt::AltButton))
      return 0;
    else
      return 1;
  } 

  Coord_point Qt_widget_demo_tab::get_point(Coord_type x, Coord_type y)
  {
	  int xmin = static_cast<int> (x_min());
	  int xmax = static_cast<int> (x_max());
	  int ymin = static_cast<int> (y_min());
	  int ymax = static_cast<int> (y_max());
	  Coord_type d = std::max(0.5 , (x_max() - x_min())/40);
	  switch ( snap_mode ) {
		case POINT:
			{
				bool first = true;
				Coord_type x1,y1,dt,min_dist = 0;
				Coord_point closest;

				if (traits_type == SEGMENT_TRAITS)
				{
					Qt_widget_segment_tab	*w_segment_p = static_cast<Qt_widget_segment_tab *> (this);
					// Go over all vertices and for each vertex check the indexes of the
					// half edges around, if you find 2 different indexes, print the vertic
					Seg_arr::Vertex_iterator   vit;
					for (vit = w_segment_p->segment_arr.vertices_begin(); vit != w_segment_p->segment_arr.vertices_end(); vit++)
					{
						const Pm_seg_point_2& p = (*vit).point();
						x1 = CGAL::to_double(p.x());
						y1 = CGAL::to_double(p.y());
						dt = dist(x1 , y1 , x , y);
						if ( dt < min_dist || first)
						{
							min_dist = dt;
							closest = Coord_point(x1 , y1);
							first = false;
						}
					}
				}
				else if (traits_type == POLYLINE_TRAITS)
				{
					Qt_widget_polyline_tab	*w_polyline_p = static_cast<Qt_widget_polyline_tab *> (this);
					// Go over all vertices and for each vertex check the indexes of the
					// half edges around, if you find 2 different indexes, print the vertic
					Pol_arr::Vertex_iterator   vit;
					for (vit = w_polyline_p->polyline_arr.vertices_begin(); vit != w_polyline_p->polyline_arr.vertices_end(); vit++)
					{
						const Pm_pol_point_2& p = (*vit).point();
						x1 = CGAL::to_double(p.x());
						y1 = CGAL::to_double(p.y());
						dt = dist(x1 , y1 , x , y);
						if ( dt < min_dist || first)
						{
							min_dist = dt;
							closest = Coord_point(x1 , y1);
							first = false;
						}
					}
				}
				if (min_dist <= d)
					return closest;
				else
					return Coord_point(x , y);

				break;
			}
		case GRID:
			return Coord_point(getMid(x, xmin, xmax), 
							   getMid(y, ymin, ymax) );
			break;
		default: // no snap
			return Coord_point(x,y);
		}
	  
  }
Coord_type Qt_widget_demo_tab::dist(Coord_type x1, Coord_type y1, Coord_type x2, Coord_type y2)
{
	return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}
   Coord_type Qt_widget_demo_tab::getMid(Coord_type coord, int my_min, int my_max)
  {
	  //int c = (my_max - my_min)/20;
	  int c = std::max(1, abs(my_max - my_min)/20);
	  Coord_type d = static_cast<Coord_type>(c)/2;
	  for (int i = my_min - c; i <= my_max; i += c)
	  {
		  Coord_type id = static_cast<Coord_type>(i); 
		  if (coord >= id - d && coord <= id + d)
		  {
			  Coord_type ans  = static_cast<Coord_type>(i);// + static_cast<Coord_type>(c);
			  return ans;
		  }
	  }
	  std::cout<< "Error: invalide frame size" << std::endl;
	  abort();
  }

Coord_segment Qt_widget_demo_tab::convert(Pm_pol_point_2 & source , Pm_pol_point_2 & target)
{
    Coord_type x1 = CGAL::to_double(source.x());
	Coord_type y1 = CGAL::to_double(source.y());
	Coord_point coord_source(x1, y1);

	Coord_type x2 = CGAL::to_double(target.x());
	Coord_type y2 = CGAL::to_double(target.y());
	Coord_point coord_target(x2, y2);

	return Coord_segment(coord_source, coord_target);
}
Coord_segment Qt_widget_demo_tab::convert(const Pm_pol_point_2 & source , const Pm_pol_point_2 & target)
{
    Coord_type x1 = CGAL::to_double(source.x());
	Coord_type y1 = CGAL::to_double(source.y());
	Coord_point coord_source(x1, y1);

	Coord_type x2 = CGAL::to_double(target.x());
	Coord_type y2 = CGAL::to_double(target.y());
	Coord_point coord_target(x2, y2);

	return Coord_segment(coord_source, coord_target);
}

//////////////////////////////////////////////////////////////////////////////

MyWindow::MyWindow(int w, int h){

	myBar = new QTabWidget(this);
    setCentralWidget(myBar);
	m_width = w;
	m_height = h;
	tab_number = 1;
	number_of_tabs = 0;
	overlay_flag = false; // flag for add_conic_tab
	testlayer = new Qt_layer_show_ch( myBar );

	// Traits Group
	QActionGroup *traitsGroup = new QActionGroup( this ); // Connected later
    traitsGroup->setExclusive( TRUE );

    setSegmentTraits = new QAction(
	    "Segment Traits", QPixmap( (const char**)line_xpm ),
	    "&Segment Traits", 0 ,traitsGroup, "Segment Traits" );
    setSegmentTraits->setToggleAction( TRUE );

	setPolylineTraits = new QAction(
	    "Polyline Traits", QPixmap( (const char**)polyline_xpm ),
	    "&Polyline Traits", 0 , traitsGroup, "Polyline Traits" );
    setPolylineTraits->setToggleAction( TRUE );

	setConicTraits = new QAction(
	    "Conic Traits", QPixmap( (const char**)conic_xpm ),
	    "&Conic Traits", 0 , traitsGroup, "Conic Traits" );
    setConicTraits->setToggleAction( TRUE );

	// Snap Mode Group
	QActionGroup *snapModeGroup = new QActionGroup( this ); // Connected later
    snapModeGroup->setExclusive( TRUE );

    setNoneSnapMode = new QAction(
	    "None Snap Mode", QPixmap( (const char**)none_xpm ),
	    "&None Snap Mode", 0 , snapModeGroup, "None Snap Mode" );
    setNoneSnapMode->setToggleAction( TRUE );

	setGridSnapMode = new QAction(
	    "Grid Snap Mode", QPixmap( (const char**)grid_xpm ),
	    "&Grid Snap Mode", 0 , snapModeGroup, "Grid Snap Mode" );
    setGridSnapMode->setToggleAction( TRUE );

	setPointSnapMode = new QAction(
	    "Point Snap Mode", QPixmap( (const char**)po_xpm ),
	    "&Point Snap Mode", 0 , snapModeGroup, "Point Snap Mode" );
    setPointSnapMode->setToggleAction( TRUE );

	// insert - delete - point_location Mode Group
	QActionGroup *modeGroup = new QActionGroup( this ); // Connected later
    modeGroup->setExclusive( TRUE );

    insertMode = new QAction(
	    "Insert", QPixmap( (const char**)insert_xpm ),
	    "&Insert", 0 , modeGroup, "Insert" );
    insertMode->setToggleAction( TRUE );

	deleteMode = new QAction(
	    "Delete", QPixmap( (const char**)delete_xpm ),
	    "&Delete", 0 , modeGroup, "Delete" );
    deleteMode->setToggleAction( TRUE );
   
	pointLocationMode = new QAction(
	    "PointLocation", QPixmap( (const char**)point_xpm ),
	    "&Point Location", 0 , modeGroup, "Point Location" );
    pointLocationMode->setToggleAction( TRUE );

	/*dragMode = new QAction(
    "Drag", QPixmap( (const char**)hand_xpm ),
    "&Drag", 0 , modeGroup, "Drag" );
    dragMode->setToggleAction( TRUE );*/

	// hand tool
    /*QToolButton* handtoolBt = new QToolButton(optionsTools, "handtool");
    handtoolBt->setPixmap(QPixmap(hand_xpm ));
    handtoolBt->setTextLabel("Pan tool");*/

    //connect(dragMode, SIGNAL(stateChanged(int)),
    //handtoollayer, SLOT(stateChanged(int)));


	//handtoollayer = new CGAL::Qt_widget_handtool(this, "handtoollayer");
 
	//create a timer for checking if somthing changed
    QTimer *timer = new QTimer( this );
    connect( timer, SIGNAL(timeout()),
           this, SLOT(timer_done()) );
    timer->start( 200, FALSE );

    // file menu
    QPopupMenu * file = new QPopupMenu( this );
    menuBar()->insertItem( "&File", file );
   	file->insertItem("&Open Conic File", this, SLOT(add_conic_tab()));
   	file->insertItem("&Open Polyline File", this, SLOT(fileOpenPolyline()));
   	file->insertItem("&Save", this, SLOT(fileSave()));
	file->insertItem("&Save As", this, SLOT(fileSaveAs()));
	file->insertItem("&Save to ps", this, SLOT(fileSave_ps()));
	file->insertSeparator();
	file->insertItem("&Print", this , SLOT(print()));
	file->insertSeparator();
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
    file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), CTRL+Key_Q );
	menuBar()->insertSeparator();

	// tab menu
	QPopupMenu * tab = new QPopupMenu( this );
    menuBar()->insertItem( "&Tab", tab );
    tab->insertItem("Add &Segment Tab", this, SLOT(add_segment_tab()));
    tab->insertItem("Add &Polyline Tab", this, SLOT(add_polyline_tab()));
	tab->insertItem("Add &Conic Tab", this, SLOT(add_conic_tab()));
	tab->insertSeparator();
	tab->insertItem("Remove &Tab", this, SLOT(remove_tab()));
	menuBar()->insertSeparator();
   
	// mode menu
	QPopupMenu * mode = new QPopupMenu( this );
    menuBar()->insertItem( "&Mode", mode );
	insertMode->addTo( mode );
    deleteMode->addTo( mode );
    pointLocationMode->addTo( mode );
	//dragMode->addTo( mode );
	menuBar()->insertSeparator();

	// snap mode menu
    QPopupMenu * snap_mode = new QPopupMenu( this );
    menuBar()->insertItem( "&Snap mode", snap_mode );
	setNoneSnapMode->addTo(snap_mode); 
	setGridSnapMode->addTo(snap_mode);
	setPointSnapMode->addTo(snap_mode);
	menuBar()->insertSeparator();

	 // options menu
    QPopupMenu * options = new QPopupMenu( this );
    menuBar()->insertItem( "&Options", options );
    setSegmentTraits->addTo(options); 
	setPolylineTraits->addTo(options); 
	setConicTraits->addTo(options); 
    options->insertSeparator();
	options->insertItem("Overlay 2 Arr", this, SLOT(optionsSetOptions()));
	menuBar()->insertSeparator();
	options->insertItem("Overlay 3 Arr", this, SLOT(overlay()));
	menuBar()->insertSeparator();
	options->insertItem("Properties", this, SLOT(properties()));

    // help menu
    QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Help", help );
    help->insertItem("How To", this, SLOT(howto()), Key_F1);
    help->insertSeparator();
    help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );

    //the new tools toolbar
    //newtoolbar = new Tools_toolbar(w_demo_p, this, &w_demo_p->list_of_curves);	
  
	// options toolbar

	QToolBar *optionsTools = new QToolBar( this, "options operations" );
    optionsTools->setLabel( "Options Operations" );
    insertMode->addTo( optionsTools );
    deleteMode->addTo( optionsTools );
	//dragMode->addTo( optionsTools );
    pointLocationMode->addTo( optionsTools );
    optionsTools->addSeparator();
    setNoneSnapMode->addTo( optionsTools );
    setGridSnapMode->addTo( optionsTools );
    setPointSnapMode->addTo( optionsTools );
    optionsTools->addSeparator();
	setSegmentTraits->addTo( optionsTools );
    setPolylineTraits->addTo( optionsTools );
    setConicTraits->addTo( optionsTools );
    optionsTools->addSeparator();
	optionsTools->addSeparator();
	
	// zoom button
	QToolButton* zoominBt =
      new QToolButton(QPixmap(zoomin_xpm ),
			       "Zoom in", 
			       0, 
			       this, 
			       SLOT(zoomin()), 
			       optionsTools, 
			       "Zoom in");
    zoominBt->setTextLabel("Scaling factor X2");

    QToolButton* zoomoutBt = 
      new QToolButton(QPixmap(zoomout_xpm ),
				"Zoom out", 
				0, 
				this, 
				SLOT(zoomout()), 
				optionsTools, 
				"Zoom out");
    zoomoutBt->setTextLabel("Scaling factor 1/2");


	//connect(dragMode, SIGNAL(stateChanged(int)),
    //handtoollayer, SLOT(stateChanged(int)));

    //connect(dragMode, SIGNAL(stateChanged(int) ), 
    //handtoollayer, SLOT(stateChanged(int)));

	// connect mode group

	connect( modeGroup, SIGNAL( selected(QAction*) ), 
		this, SLOT( updateMode(QAction*) ) );

	// connect Traits Group
	connect( traitsGroup, SIGNAL( selected(QAction*) ),
	     this, SLOT( updateTraitsType(QAction*) ) );

	// connect Snap Mode Group
	connect( snapModeGroup, SIGNAL( selected(QAction*) ),
	     this, SLOT( updateSnapMode(QAction*) ) );
	
	//application flag stuff
    old_state = 0;
   
	add_segment_tab();

	// connect the change of current tab
	connect( myBar, SIGNAL( currentChanged(QWidget * )  ),
	     this, SLOT( update() ) );
	
 }

MyWindow::~MyWindow()
{}

void MyWindow::something_changed()
{
	// We peform downcasting from QWigdet* to Qt_widget_demo_tab*, as we know that only
	// Qt_widget_demo_tab objects are stored in the tab pages.
	Qt_widget_demo_tab	*w_demo_p = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());

	w_demo_p->current_state++;
	//std::cout << myBar->currentPageIndex();
	//std::fflush( stdout );
}
   

void MyWindow::new_instance()
  {
	// We peform downcasting from QWigdet* to Qt_widget_demo_tab*, as we know that only
	// Qt_widget_demo_tab objects are stored in the tab pages.
	Qt_widget_demo_tab	*w_demo_p = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());

    w_demo_p->lock();
	if (w_demo_p->traits_type == SEGMENT_TRAITS)
	{
		Qt_widget_segment_tab	*w_segment_p = static_cast<Qt_widget_segment_tab *> (myBar->currentPage());
        w_segment_p->list_of_curves.clear();
		w_segment_p->segment_arr.clear();
	}
	else
	{
		Qt_widget_polyline_tab	*w_polyline_p = static_cast<Qt_widget_polyline_tab *> (myBar->currentPage());
        w_polyline_p->list_of_curves.clear();
		w_polyline_p->polyline_arr.clear();
	}
    w_demo_p->set_window(-100, 700, -100, 700);
			// set the Visible Area to the Interval
    w_demo_p->unlock();
    something_changed();
  }


void MyWindow::get_new_object(CGAL::Object obj)
  {
	// We peform downcasting from QWigdet* to Qt_widget_demo_tab*, as we know that only
	// Qt_widget_demo_tab objects are stored in the tab pages.
	Qt_widget_demo_tab	*w_demo_p = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());

	if (w_demo_p->traits_type == SEGMENT_TRAITS)
	{
		Qt_widget_segment_tab	*w_segment_p = static_cast<Qt_widget_segment_tab *> (myBar->currentPage());

		// Segment Object
		Coord_segment coord_seg;
		if(CGAL::assign(coord_seg, obj)) {
			// Convert inexact segment to exact segment
			const Coord_point & coord_source = coord_seg.source();
			const Coord_point & coord_target = coord_seg.target();
		    Pm_seg_point_2 source(coord_source.x(), coord_source.y());
			Pm_seg_point_2 target(coord_target.x(), coord_target.y());
			//Pm_seg_2 * seg = new Pm_seg_2( Pm_base_seg_2(source, target), w_segment_p->index);
			Pm_base_seg_2 * base_seg_p = new Pm_base_seg_2(source, target);
			Curve_data cd;
			cd.m_type = Curve_data::LEAF;
			cd.m_index = w_segment_p->index;
			cd.m_ptr.m_curve = base_seg_p;
			Pm_seg_2 * seg = new Pm_seg_2( *base_seg_p, cd );
			w_segment_p->list_of_curves.push_back(seg);
			w_segment_p->segment_arr.insert(*seg);
			something_changed();
		}
	}

	// point location
    Coord_point p;
    if(CGAL::assign(p,obj)) {

      w_demo_p->pl_point = p;
      //pl_valid = true;
      
      something_changed();
    }

  }

  void MyWindow::about()
  {
    QMessageBox::about( this, my_title_string,
		"This is a demo for the Arrangement package\n"
  		"Copyright CGAL @2003");
  }

  void MyWindow::aboutQt()
  {
    QMessageBox::aboutQt( this, my_title_string );
  }

  void MyWindow::howto(){
    QString home;
    home = "help/index.html";
    CGAL::Qt_help_window * help =
      new CGAL::Qt_help_window(home, ".", 0, "help viewer");
    help->resize(400, 400);
    help->setCaption("Demo HowTo");
    help->show();
  }

  void MyWindow::add_segment_tab()
  {
	Qt_widget_demo_tab *widget = new Qt_widget_segment_tab(this);
	// initialize the widget
	init_widget( widget );
	// add the new widget to myBar
	myBar->insertTab( widget, QString("Arr " + QString::number( widget->index ) ) , widget->index );
	myBar->setCurrentPage(myBar->indexOf(widget));
	something_changed();
  }

  void MyWindow::add_polyline_tab()
  {
	Qt_widget_demo_tab *widget = new Qt_widget_polyline_tab(this);
	// initialize the widget
	init_widget( widget );
	// add the new widget to myBar
	myBar->insertTab( widget, QString("Arr " + QString::number( widget->index ) ) , widget->index );
	myBar->setCurrentPage(myBar->indexOf(widget));
	something_changed();
  }
  
  void MyWindow::add_conic_tab()
  {
	Qt_widget_demo_tab *widget = new Qt_widget_conic_tab(this);
	// initialize the widget
	init_widget( widget );
	// add the new widget to myBar
	myBar->insertTab( widget, QString("Arr " + QString::number( widget->index ) ) , widget->index );
	myBar->setCurrentPage(myBar->indexOf(widget));
	if (!overlay_flag)
		fileOpen();
	something_changed();
  }

void MyWindow::init_widget(Qt_widget_demo_tab *widget)
{
// initialize the new tab widget
	*widget << CGAL::LineWidth(2) << CGAL::BackgroundColor (CGAL::WHITE);
	widget->set_window(-100, 700, -100, 700);
    widget->setMouseTracking(TRUE);
	connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),this, SLOT(get_new_object(CGAL::Object)));
    widget->attach(testlayer);
		
	widget->m_line_width = 2;
	widget->index = tab_number;
	tab_number++;
	number_of_tabs++;

	resize(m_width,m_height);
	
}

  void MyWindow::remove_tab()
  {
	  if (number_of_tabs > 1)
	  {
          myBar->removePage(myBar->currentPage());
		  number_of_tabs--;
	  }
	  else
	  {
		  QMessageBox::information( this, my_title_string,
			"Can not remove last tab");
	  }

  }


  void MyWindow::timer_done()
  {
	// We peform downcasting from QWigdet* to Qt_widget_demo_tab*, as we know that only
	// Qt_widget_demo_tab objects are stored in the tab pages.
	Qt_widget_demo_tab	*w_demo_p = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());

    if(old_state!=w_demo_p->current_state){
      w_demo_p->redraw();
      old_state = w_demo_p->current_state;
    }
  }	

	void MyWindow::updateTraitsType( QAction *action )
	{
		if (tab_number == 0) return;
		// We peform downcasting from QWigdet* to Qt_widget_demo_tab*, as we know that only
		// Qt_widget_demo_tab objects are stored in the tab pages.
		Qt_widget_demo_tab *old_widget = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());
		Qt_widget_demo_tab *widget;

		if (action == setSegmentTraits)
		{
			if (old_widget->traits_type == SEGMENT_TRAITS) return;
			widget = new Qt_widget_segment_tab(this);
		}
		else if (action == setPolylineTraits)
		{
			if (old_widget->traits_type == POLYLINE_TRAITS) return;
			widget = new Qt_widget_polyline_tab(this);
		}
		else if (action == setConicTraits)
		{
			if (old_widget->traits_type == CONIC_TRAITS) return;
			widget = new Qt_widget_conic_tab(this);
		}

		int old_index = old_widget->index;
		int index = myBar->currentPageIndex();
		QString label = myBar->label(index);
		myBar->removePage(myBar->currentPage());

			//initialize the new tab widget
		*widget << CGAL::LineWidth(2) << CGAL::BackgroundColor (CGAL::WHITE);
		widget->set_window(-100, 700, -100, 700);
		widget->setMouseTracking(TRUE);
		connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),this, SLOT(get_new_object(CGAL::Object)));
		widget->attach(testlayer);
		widget->m_line_width = 2;
		widget->index = old_index;

		// add the new widget to myBar
		myBar->insertTab( widget, label , index );

		myBar->setCurrentPage(index);
			
		resize(m_width,m_height);
		
		if (action == setConicTraits)
			fileOpen();

		something_changed();
	
	}

	void MyWindow::setTraits( TraitsType t )
	{
		//w_demo_p->traits_type = t;
		switch ( t ) {
		case SEGMENT_TRAITS:
			setSegmentTraits->setOn( TRUE );
			insertMode->setEnabled( TRUE );
			deleteMode->setEnabled( TRUE );
			break;
		case POLYLINE_TRAITS:
			setPolylineTraits->setOn( TRUE );
			insertMode->setEnabled( TRUE );
			deleteMode->setEnabled( TRUE );
			break;
		case CONIC_TRAITS:
			setConicTraits->setOn( TRUE );
			insertMode->setEnabled( FALSE );
			deleteMode->setEnabled( FALSE );
			break;
		}
	}

	void MyWindow::updateSnapMode( QAction *action )
	{
		// We peform downcasting from QWigdet* to Qt_widget_demo_tab*, as we know that only
		// Qt_widget_demo_tab objects are stored in the tab pages.
		Qt_widget_demo_tab	*w_demo_p = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());

		if ( action == setNoneSnapMode ) {
		w_demo_p->snap_mode = NONE;
		}
		else if ( action == setGridSnapMode ) {
		w_demo_p->snap_mode = GRID;
		}
		else if ( action == setPointSnapMode ) {
		w_demo_p->snap_mode = POINT;
		}

		something_changed();
	}

	void MyWindow::setSnapMode( SnapMode m )
	{
		// We peform downcasting from QWigdet* to Qt_widget_demo_tab*, as we know that only
		// Qt_widget_demo_tab objects are stored in the tab pages.
		Qt_widget_demo_tab	*w_demo_p = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());

		w_demo_p->snap_mode = m;
		switch ( m ) {
		case NONE:
			setNoneSnapMode->setOn( TRUE );
			something_changed();
			break;
		case GRID:
			setGridSnapMode->setOn( TRUE );
			something_changed();
			break;
		case POINT:
			setPointSnapMode->setOn( TRUE );
			something_changed();
			break;

		}
	}

	void MyWindow::updateMode( QAction *action )
	{
		// We peform downcasting from QWigdet* to Qt_widget_demo_tab*, as we know that only
		// Qt_widget_demo_tab objects are stored in the tab pages.
		Qt_widget_demo_tab	*w_demo_p = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());

		if ( action == insertMode ) 
		{
            w_demo_p->mode = INSERT;
			//handtoollayer->deactivate();
			something_changed();
		}
		else if ( action == deleteMode ) 
		{
			w_demo_p->mode = DELETE;
			//handtoollayer->deactivate();
			something_changed();
		}
		else if ( action == pointLocationMode ) 
		{
			w_demo_p->mode = POINT_LOCATION;
			//handtoollayer->deactivate();
			//something_changed();
		}
		/*else if ( action == dragMode ) 
		{
			w_demo_p->mode = DRAG;
			handtoollayer->activate();
			handtoollayer->stateChanged(w_demo_p->current_state);
			something_changed();
		}*/
	}

	void MyWindow::setMode( Mode m )
	{
		// We peform downcasting from QWigdet* to Qt_widget_demo_tab*, as we know that only
		// Qt_widget_demo_tab objects are stored in the tab pages.
		Qt_widget_demo_tab	*w_demo_p = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());

		w_demo_p->mode = m;
		switch ( m ) {
		case INSERT:
			insertMode->setOn( TRUE );
			break;
		case DELETE:
			deleteMode->setOn( TRUE );
			break;
		case POINT_LOCATION:
			pointLocationMode->setOn( TRUE );
			break;
		case DRAG:
			dragMode->setOn( TRUE );
			break;
		}
	}

	void MyWindow::update()
	{
		// We peform downcasting from QWigdet* to Qt_widget_demo_tab*, as we know that only
		// Qt_widget_demo_tab objects are stored in the tab pages.
		Qt_widget_demo_tab	*w_demo_p = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());
		setMode( w_demo_p->mode );
		setSnapMode( w_demo_p->snap_mode );
		setTraits( w_demo_p->traits_type );
				
	}

	void MyWindow::zoomin()
	{
        Qt_widget_demo_tab	*w_demo_p = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());
		w_demo_p->zoom(2);
	}

	void MyWindow::zoomout()
	{
        Qt_widget_demo_tab	*w_demo_p = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());
		w_demo_p->zoom(0.5);
	}

void MyWindow::optionsSetOptions()
{
    OptionsForm *optionsForm = new OptionsForm( myBar , this ,number_of_tabs );
    
    if ( optionsForm->exec() ) 
	{	
		if (optionsForm->arrComboBox1->currentItem() == optionsForm->arrComboBox2->currentItem())
			QMessageBox::information( this, my_title_string,"This is not a clever idea!!!");
			
		else
		{
			int index1 = optionsForm->arrComboBox1->currentItem();
			int index2 = optionsForm->arrComboBox2->currentItem();
            Qt_widget_demo_tab	*w_demo_p1 = static_cast<Qt_widget_demo_tab *> (myBar->page( index1 ));
			Qt_widget_demo_tab	*w_demo_p2 = static_cast<Qt_widget_demo_tab *> (myBar->page( index2 ));

			if (w_demo_p1->traits_type != w_demo_p2->traits_type)
				QMessageBox::information( this, my_title_string,"Can not union arr of different traits");
			else union_arr( index1 , index2 , w_demo_p2->traits_type);
		}
        //std::cout << optionsForm->arrComboBox1->currentItem() << std::endl;
    }
    delete optionsForm;
}

void MyWindow::properties()
{
    PropertiesForm *optionsForm = new PropertiesForm( myBar , this ,number_of_tabs );
    
    if ( optionsForm->exec() ) 
	{	
		Qt_widget_demo_tab	*w_demo_p1 = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());
		m_width = optionsForm->box1->value();
		m_height = optionsForm->box2->value();
		w_demo_p1->m_line_width = optionsForm->box3->value();
		resize(m_width,m_height);
		w_demo_p1->redraw();
		something_changed();
    }
    delete optionsForm;
}

void MyWindow::union_arr( int index1 , int index2 , TraitsType t)
{
	switch ( t ) {
		case SEGMENT_TRAITS:
		{
			add_segment_tab();
            Qt_widget_segment_tab *w_demo_p_new = static_cast<Qt_widget_segment_tab *> (myBar->currentPage());
            Qt_widget_segment_tab *w_demo_p1 = static_cast<Qt_widget_segment_tab *> (myBar->page( index1 ));
			Qt_widget_segment_tab *w_demo_p2 = static_cast<Qt_widget_segment_tab *> (myBar->page( index2 ));

			Pm_seg_const_iter itp;
			std::list<Pm_seg_2> seg_list;

			*w_demo_p_new << CGAL::RED;
			*w_demo_p_new << CGAL::LineWidth(3);
			for (itp = w_demo_p1->list_of_curves.begin(); itp != w_demo_p1->list_of_curves.end(); ++itp)
			{
				Pm_seg_2 * seg = new Pm_seg_2( **itp );
				w_demo_p_new->list_of_curves.push_back(seg);
				seg_list.push_back(*seg);
				//w_demo_p_new->segment_arr.insert(*seg);
				*w_demo_p_new << *(*itp);
			}
			
			*w_demo_p_new << CGAL::BLUE;
			for (itp = w_demo_p2->list_of_curves.begin(); itp != w_demo_p2->list_of_curves.end(); ++itp)
			{
				Pm_seg_2 * seg = new Pm_seg_2( **itp );
				w_demo_p_new->list_of_curves.push_back(seg);
				seg_list.push_back(*seg);
				//w_demo_p_new->segment_arr.insert(*seg);
				*w_demo_p_new << *(*itp);
			}

			w_demo_p_new->segment_arr.insert(seg_list.begin(),seg_list.end());

			*w_demo_p_new << CGAL::GREEN;
			*w_demo_p_new << CGAL::CROSS;
			// Go over all vertices and for each vertex check the indexes of the
			// half edges around, if you find 2 different indexes, print the vertic
			Seg_arr::Vertex_iterator   vit;
			for (vit = w_demo_p_new->segment_arr.vertices_begin(); vit != w_demo_p_new->segment_arr.vertices_end(); vit++)
			{
				Seg_arr::Halfedge_around_vertex_circulator 
				eit, first = (*vit).incident_halfedges();

				eit = first;

				int ind1;
				int ind2 = (*eit).curve().get_data().m_index;

				do 
				{
					ind1 = (*eit).curve().get_data().m_index;

					// Keep track of IDs we haven't seen before.
					if (ind1 != ind2)
					{
						const Pm_seg_point_2& p = (*vit).point();
						*w_demo_p_new << p;
						break;
					}

					eit++;

				} while (eit != first);
			}

			w_demo_p_new->current_state = old_state;

			//w_demo_p_new->m_xmin = std::min(w_demo_p1->m_xmin,w_demo_p2->m_xmin); 
			//w_demo_p_new->m_xmax = std::max(w_demo_p1->m_xmax,w_demo_p2->m_xmax);
			//w_demo_p_new->m_ymin = std::min(w_demo_p1->m_ymin,w_demo_p2->m_ymin);
			//w_demo_p_new->m_ymax = std::max(w_demo_p1->m_ymax,w_demo_p2->m_ymax);
			//w_demo_p_new->set_window(w_demo_p_new->m_xmin , w_demo_p_new->m_xmax , w_demo_p_new->m_ymin , w_demo_p_new->m_ymax);

			// update new planner map index
			Pm_seg_iter iter;
			for (iter = w_demo_p_new->list_of_curves.begin(); iter != w_demo_p_new->list_of_curves.end(); ++iter)
			{
				Curve_data cd = (**iter).get_data();
				cd.m_type = Curve_data::INTERNAL;
				cd.m_index = w_demo_p_new->index;
				(**iter).set_data( cd );
			}

			break;
		}
		case POLYLINE_TRAITS:
		{
			add_polyline_tab();
			Qt_widget_polyline_tab *w_demo_p_new = static_cast<Qt_widget_polyline_tab *> (myBar->currentPage());
            Qt_widget_polyline_tab *w_demo_p1 = static_cast<Qt_widget_polyline_tab *> (myBar->page( index1 ));
			Qt_widget_polyline_tab *w_demo_p2 = static_cast<Qt_widget_polyline_tab *> (myBar->page( index2 ));

			std::list<Pm_pol_2> pol_list;
			Pm_pol_const_iter itp;
			
			*w_demo_p_new << CGAL::RED;
			*w_demo_p_new << CGAL::LineWidth(3);
			for (itp = w_demo_p1->list_of_curves.begin(); itp != w_demo_p1->list_of_curves.end(); ++itp)
			{
				Pm_pol_2 * seg = new Pm_pol_2( **itp );
				w_demo_p_new->list_of_curves.push_back(seg);
				pol_list.push_back(*seg);
				//w_demo_p_new->polyline_arr.insert(*seg);
				*w_demo_p_new << *(*itp);
			}
			
			*w_demo_p_new << CGAL::BLUE;
			for (itp = w_demo_p2->list_of_curves.begin(); itp != w_demo_p2->list_of_curves.end(); ++itp)
			{
				Pm_pol_2 * seg = new Pm_pol_2( **itp );
				w_demo_p_new->list_of_curves.push_back(seg);
				pol_list.push_back(*seg);
				//w_demo_p_new->polyline_arr.insert(*seg);
				*w_demo_p_new << *(*itp);
			}

			w_demo_p_new->polyline_arr.insert(pol_list.begin(),pol_list.end());

			*w_demo_p_new << CGAL::GREEN;
			*w_demo_p_new << CGAL::CROSS;
			// Go over all vertices and for each vertex print the ID numbers of the
			// base curves that go through it.
			Pol_arr::Vertex_iterator   vit;
			for (vit = w_demo_p_new->polyline_arr.vertices_begin(); vit != w_demo_p_new->polyline_arr.vertices_end(); vit++)
			{
				Pol_arr::Halfedge_around_vertex_circulator 
				eit, first = (*vit).incident_halfedges();

				eit = first;

				int ind1;
				int ind2 = (*eit).curve().get_data().m_index;
				do 
				{
					ind1 = (*eit).curve().get_data().m_index;

					// Keep track of IDs we haven't seen before.
					if (ind1 != ind2)
					{
						const Pm_pol_point_2& p = (*vit).point();
						*w_demo_p_new << p;
						break;
					}

					eit++;

				} while (eit != first);
			}

			// allow the user to see the overlay
			w_demo_p_new->current_state = old_state;

		/*	w_demo_p_new->m_xmin = std::min(w_demo_p1->m_xmin,w_demo_p2->m_xmin); 
			w_demo_p_new->m_xmax = std::max(w_demo_p1->m_xmax,w_demo_p2->m_xmax);
			w_demo_p_new->m_ymin = std::min(w_demo_p1->m_ymin,w_demo_p2->m_ymin);
			w_demo_p_new->m_ymax = std::max(w_demo_p1->m_ymax,w_demo_p2->m_ymax);
			w_demo_p_new->set_window(w_demo_p_new->m_xmin , w_demo_p_new->m_xmax , w_demo_p_new->m_ymin , w_demo_p_new->m_ymax);
	*/
			// update new planner map index
			Pm_pol_iter iter;
			for (iter = w_demo_p_new->list_of_curves.begin(); iter != w_demo_p_new->list_of_curves.end(); ++iter)
			{
				Curve_pol_data cd = (**iter).get_data();
				cd.m_type = Curve_pol_data::INTERNAL;
				cd.m_index = w_demo_p_new->index;
				(**iter).set_data( cd );
			}

			break;
		}
		case CONIC_TRAITS:
			overlay_flag = true;
			add_conic_tab();
			overlay_flag = false;

            Qt_widget_conic_tab *w_demo_p_new = static_cast<Qt_widget_conic_tab *> (myBar->currentPage());
            Qt_widget_conic_tab *w_demo_p1 = static_cast<Qt_widget_conic_tab *> (myBar->page( index1 ));
			Qt_widget_conic_tab *w_demo_p2 = static_cast<Qt_widget_conic_tab *> (myBar->page( index2 ));

			Pm_xconic_const_iter itp;
			std::list<Pm_conic_2> curve_list;

			*w_demo_p_new << CGAL::RED;
			*w_demo_p_new << CGAL::LineWidth(3);
			for (itp = w_demo_p1->list_of_xcurves.begin(); itp != w_demo_p1->list_of_xcurves.end(); ++itp)
			{
				Pm_xconic_2 * seg = new Pm_xconic_2( **itp );
				w_demo_p_new->list_of_xcurves.push_back(seg);

				Pm_base_conic_2* org_conic = w_demo_p_new->get_origin_curve(*seg);
     			//curve_list.push_back(Pm_conic_2( *org_conic , seg->get_data() ));
	
				w_demo_p_new->conic_arr.insert(Pm_conic_2( *org_conic , seg->get_data() ));
				w_demo_p_new->draw_curve(**itp);
			}
			
			*w_demo_p_new << CGAL::BLUE;
			for (itp = w_demo_p2->list_of_xcurves.begin(); itp != w_demo_p2->list_of_xcurves.end(); ++itp)
			{
				Pm_xconic_2 * seg = new Pm_xconic_2( **itp );
				w_demo_p_new->list_of_xcurves.push_back(seg);

				Pm_base_conic_2* org_conic = w_demo_p_new->get_origin_curve(*seg);
				//curve_list.push_back(Pm_conic_2( *org_conic , seg->get_data() ));

				w_demo_p_new->conic_arr.insert(Pm_conic_2( *org_conic , seg->get_data() ));
				w_demo_p_new->draw_curve(**itp);
			}

			//w_demo_p_new->conic_arr.insert(curve_list.begin() , curve_list.end());

			*w_demo_p_new << CGAL::GREEN;
			*w_demo_p_new << CGAL::CROSS;
			// Go over all vertices and for each vertex print the ID numbers of the
			// base curves that go through it.
			Conic_arr::Vertex_iterator   vit;
			for (vit = w_demo_p_new->conic_arr.vertices_begin(); vit != w_demo_p_new->conic_arr.vertices_end(); vit++)
			{
				Conic_arr::Halfedge_around_vertex_circulator 
				eit, first = (*vit).incident_halfedges();

				eit = first;

				int ind1;
				int ind2 = (*eit).curve().get_data().m_index;

				do 
				{
					ind1 = (*eit).curve().get_data().m_index;

					// Keep track of IDs we haven't seen before.
					if (ind1 != ind2)
					{
						const Pm_conic_point_2& p = (*vit).point();
						*w_demo_p_new << p;
						break;
					}

					eit++;

				} while (eit != first);
			}

			w_demo_p_new->current_state = old_state;
	
	/*		w_demo_p_new->m_xmin = std::min(w_demo_p1->m_xmin,w_demo_p2->m_xmin); 
			w_demo_p_new->m_xmax = std::max(w_demo_p1->m_xmax,w_demo_p2->m_xmax);
			w_demo_p_new->m_ymin = std::min(w_demo_p1->m_ymin,w_demo_p2->m_ymin);
			w_demo_p_new->m_ymax = std::max(w_demo_p1->m_ymax,w_demo_p2->m_ymax);
			w_demo_p_new->set_window(w_demo_p_new->m_xmin , w_demo_p_new->m_xmax , w_demo_p_new->m_ymin , w_demo_p_new->m_ymax);*/

			// update new planner map index
			Pm_xconic_iter iter;
			for (iter = w_demo_p_new->list_of_xcurves.begin(); iter != w_demo_p_new->list_of_xcurves.end(); ++iter)
			{
				Curve_conic_data cd = (**iter).get_data();
				cd.m_type = Curve_conic_data::INTERNAL;
				cd.m_index = w_demo_p_new->index;
				(**iter).set_data( cd );
			}
			break;

	}
}

void MyWindow::fileOpenPolyline()
{
	add_polyline_tab();
	fileOpen();
}

void MyWindow::fileOpen()
{
    
    QString filename = QFileDialog::getOpenFileName(
			    QString::null, 0, this,
			    "file open", "Demo -- File Open" );
    if ( !filename.isEmpty() )
	load( filename );
    else
	statusBar()->message( "File Open abandoned", 2000 );
}

void MyWindow::load( const QString& filename )
{
		std::ifstream inputFile(filename);
        // Creates an ofstream object named inputFile
		if (! inputFile.is_open()) // Always test file open
		{
			std::cout << "Error opening input file" << std::endl;
			return;
		}

        Qt_widget_demo_tab	*w_demo = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());

		if (w_demo->traits_type == CONIC_TRAITS)
		{
			Qt_widget_conic_tab	*w_demo_p = static_cast<Qt_widget_conic_tab *> (myBar->currentPage());
			char dummy[256];
			Pm_base_conic_2* cv;
			int count;
			//std::list<Pm_conic_2> curve_list;
		
			inputFile >> count;
			inputFile.getline(dummy, sizeof(dummy));
			for (int i = 0; i < count; i++) 
			{
				cv = new Pm_base_conic_2();
				ReadCurve(inputFile, *cv);

				Curve_conic_data cd;
				cd.m_type = Curve_conic_data::LEAF;
				cd.m_index = w_demo_p->index;
				cd.m_ptr.m_curve = cv;

				//curve_list.push_back(Pm_conic_2( *cv , cd));
				w_demo_p->conic_arr.insert(Pm_conic_2( *cv , cd));
			}
			
			//w_demo_p->conic_arr.insert(curve_list.begin() , curve_list.end());
			// insert xcurve into xcurve list

			Conic_arr::Edge_iterator ei;

			for (ei = w_demo_p->conic_arr.edges_begin(); ei != w_demo_p->conic_arr.edges_end(); ++ei) 
			{
				Pm_xconic_2 *xseg = new Pm_xconic_2(ei->curve());
				w_demo_p->list_of_xcurves.push_back(xseg);
			}
		}

		else if (w_demo->traits_type == POLYLINE_TRAITS)
		{
			Qt_widget_polyline_tab	*w_demo_p = static_cast<Qt_widget_polyline_tab *> (myBar->currentPage());
			w_demo_p->list_of_curves.clear();
			w_demo_p->polyline_arr.clear();
			
			int num_polylines, num_segments;
			int ix, iy , xmin , xmax , ymin , ymax;
			std::vector<Pm_pol_point_2> points;
			int i, j;

			inputFile >> num_polylines;
			for (i = 0; i < num_polylines; i++) 
			{
				inputFile >> num_segments;
				points.clear();
				for (j = 0; j < num_segments; j++)
				{
					inputFile >> ix >> iy;
					
					points.push_back (Pm_pol_point_2(NT(ix),NT(iy)));
					if (j == 0 && i == 0) {
						xmin = xmax = ix;
						ymin = ymax = iy;
					} else {
						if (ix < xmin) xmin = ix;
						if (ix > xmax) xmax = ix;
						if (iy < ymin) ymin = iy;
						if (iy > ymax) ymax = iy;
					}
				}

				Pm_base_pol_2 *base_polyline = new Pm_base_pol_2(points.begin(), points.end());
				
				Curve_pol_data cd;
				cd.m_type = Curve_pol_data::LEAF;
				cd.m_index = w_demo_p->index;
				cd.m_ptr.m_curve = base_polyline;

				w_demo_p->polyline_arr.insert(Pm_pol_2( *base_polyline , cd));
				w_demo_p->list_of_curves.push_back(new Pm_pol_2( *base_polyline , cd));
		
			}
			//int px = (xmax - xmin)/10;
			//int py = (ymax - ymin)/10;
			//w_demo_p->m_xmin = xmin-px;
			//w_demo_p->m_xmax = xmax+px;
			//w_demo_p->m_ymin = ymin-py;
			//w_demo_p->m_ymax = ymax+py;

			//std::cout << xmin-px << " " << xmax+px << " " << ymin-py << " " << ymax+py;
			//std::fflush(stdout);

			w_demo_p->set_window(-100, 700, -100, 700);
			//w_demo_p->set_window(xmin-px , xmax+px , ymin-py , ymax+py);

		}

		inputFile.close();
		something_changed();
}


void MyWindow::ReadCurve(std::ifstream & is, Pm_base_conic_2 & cv)
{
      // Read a line from the input file.
      char one_line[128];
      
      skip_comments (is, one_line);
      std::string stringvalues(one_line);
      std::istringstream str_line (stringvalues, std::istringstream::in);
      
      // Get the arc type.
      char     type;
      bool     is_circle = false;              // Is this a circle.
      Pm_conic_circle_2 circle;
      CONIC_NT       r, s, t, u, v, w;               // The conic coefficients.
      
      str_line >> type;
      
      // An ellipse (full ellipse or a partial ellipse):
      if (type == 'f' || type == 'F' || type == 'e' || type == 'E')
      {  
          // Read the ellipse (using the format "a b x0 y0"):
          //
          //     x - x0   2      y - y0   2
          //  ( -------- )  + ( -------- )  = 1
          //       a               b
          //
          CONIC_NT     a, b, x0, y0;
          
          str_line >> a >> b >> x0 >> y0;
          
          CONIC_NT     a_sq = a*a;
          CONIC_NT     b_sq = b*b;
          
          if (a == b)
          {
              is_circle = true;
              circle = Pm_conic_circle_2 (Pm_conic_point_2 (x0, y0), a*b, CGAL::CLOCKWISE);
          }
          else
          {
              r = b_sq;
              s = a_sq;
              t = 0;
              u = -2*x0*b_sq;
              v = -2*y0*a_sq;
              w = x0*x0*b_sq + y0*y0*a_sq - a_sq*b_sq;
          }
          
          if (type == 'f' || type == 'F')
          {
              // Create a full ellipse (or circle).
              if (is_circle)
                  cv = Pm_base_conic_2 (circle);
              else
                  cv = Pm_base_conic_2 (r, s, t, u, v, w);
              
              return;
          }
      }
      else if (type == 'h' || type == 'H')
      {
          // Read the hyperbola (using the format "a b x0 y0"):
          //
          //     x - x0   2      y - y0   2
          //  ( -------- )  - ( -------- )  = 1
          //       a               b
          //
          CONIC_NT     a, b, x0, y0;
          
          str_line >> a >> b >> x0 >> y0;
          
          CONIC_NT     a_sq = a*a;
          CONIC_NT     b_sq = b*b;
          
          r = b_sq;
          s= -a_sq;
          t = 0;
          u = -2*x0*b_sq;
          v = 2*y0*a_sq;
          w = x0*x0*b_sq - y0*y0*a_sq - a_sq*b_sq;  
      }
      else if (type == 'p' || type == 'P')
      {
          // Read the parabola (using the format "c x0 y0"):
          //
          //                        2
          //  4c*(y - y0) = (x - x0)
          //
          CONIC_NT     c, x0, y0;
          
          str_line >> c >> x0 >> y0;
          
          r = 1;
          s = 0;
          t = 0;
          u = -2*x0;
          v = -4*c;
          w = x0*x0 + 4*c*y0;
      }
      else if (type == 'c' || type == 'C' || type == 'a' || type == 'A')
      {
          // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
          str_line >> r >> s >> t >> u >> v >> w;
          
          if (type == 'c' || type == 'C')
          {
              // Create a full conic (should work only for ellipses).
              cv = Pm_base_conic_2 (r, s, t, u, v, w);
              return;
          }
      }
      else if (type == 's' || type == 'S')
      {
          // Read a segment, given by its endpoints (x1,y1) and (x2,y2);
          CONIC_NT      x1, y1, x2, y2;
          
          str_line >> x1 >> y1 >> x2 >> y2;
          
          Pm_conic_point_2   source (x1, y1);
          Pm_conic_point_2   target (x2, y2);
          Pm_conic_segment_2 segment (source, target);
          
          // Create the segment.
          cv = Pm_base_conic_2(segment);
          return;
      }
      else if (type == 'i' || type == 'I')
      {
          // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
          str_line >> r >> s >> t >> u >> v >> w;
          
          // Read the approximated source, along with a general conic 
          // <r_1,s_1,t_1,u_1,v_1,w_1> whose intersection with <r,s,t,u,v,w>
          // defines the source.
          CONIC_NT     r1, s1, t1, u1, v1, w1;
          CONIC_NT     x1, y1;
          
          str_line >> x1 >> y1;
          str_line >> r1 >> s1 >> t1 >> u1 >> v1 >> w1;
          
          Pm_conic_point_2   app_source (x1, y1);
          
          // Read the approximated target, along with a general conic 
          // <r_2,s_2,t_2,u_2,v_2,w_2> whose intersection with <r,s,t,u,v,w>
          // defines the target.
          CONIC_NT     r2, s2, t2, u2, v2, w2;
          CONIC_NT     x2, y2;
          
          str_line >> x2 >> y2;
          str_line >> r2 >> s2 >> t2 >> u2 >> v2 >> w2;
          
          Pm_conic_point_2   app_target (x2, y2);
          
          // Create the conic arc.
          cv = Pm_base_conic_2 (r, s, t, u, v ,w,
                        app_source, r1, s1, t1, u1, v1, w1,
                        app_target, r2, s2, t2, u2, v2, w2);
          return;
      }
      else
      {
          std::cerr << "Illegal conic type specification: " << type << "."
                    << std::endl;
          return;
      }
      
      // Read the end points of the arc and create it.
      CONIC_NT    x1, y1, x2, y2;
      
      str_line >> x1 >> y1 >> x2 >> y2;
      
      Pm_conic_point_2 source (x1, y1);
      Pm_conic_point_2 target (x2, y2);
      
      // Create the conic (or circular) arc.
      if (is_circle)
      {
          cv = Pm_base_conic_2 (circle,
                        source, target);
      }
      else
      {
          cv = Pm_base_conic_2 (r, s, t, u, v, w,
                        source, target);
      }
      
      return;
}
    
void MyWindow::skip_comments( std::ifstream& is, char* one_line )
{
    while( !is.eof() )
	{
        is.getline( one_line, 128 );
        if( one_line[0] != '#' )
		{
            break;
        }
    }
}

void MyWindow::fileSaveAs()
{
    QString filename = QFileDialog::getSaveFileName(
			    QString::null, "Planar Map (*.pm)", this,
			    "file save as", "Planar Map -- File Save As" );
    if ( !filename.isEmpty() ) 
	{
		int answer = 0;
		if ( QFile::exists( filename ) )
			answer = QMessageBox::warning(
					this, "Overwrite File",
					QString( "Overwrite\n\'%1\'?" ).
					arg( filename ),
					"&Yes", "&No", QString::null, 1, 1 );
		if ( answer == 0 ) 
		{
			m_filename = filename;
			//updateRecentFiles( filename );
			fileSave();
			return;
		}
    }
    statusBar()->message( "Saving abandoned", 2000 );
}

void MyWindow::fileSave()
{
    if ( m_filename.isEmpty() ) {
	fileSaveAs();
	return;
    }

	std::ofstream outFile(m_filename);
    // Creates an ofstream object named outFile
	if (! outFile.is_open()) // Always test file open
	{
		std::cout << "Error opening input file" << std::endl;
		return;
	}

    Qt_widget_demo_tab	*w_demo_p1 = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());

	switch ( w_demo_p1->traits_type ) {
		case SEGMENT_TRAITS:
		{
			Qt_widget_segment_tab *w_demo_p = static_cast<Qt_widget_segment_tab *> (myBar->currentPage());
			outFile << w_demo_p->segment_arr;
			break;
		}
		case POLYLINE_TRAITS:
		{
			Qt_widget_polyline_tab *w_demo_p = static_cast<Qt_widget_polyline_tab *> (myBar->currentPage());
			outFile << w_demo_p->polyline_arr;
			break;
		}
		case CONIC_TRAITS:
		{
			Qt_widget_conic_tab *w_demo_p = static_cast<Qt_widget_conic_tab *> (myBar->currentPage());
			outFile << w_demo_p->conic_arr;
			break;
		}
	}  

	outFile.close();

    setCaption( QString( "Planar Map -- %1" ).arg( m_filename ) );
    statusBar()->message( QString( "Saved \'%1\'" ).arg( m_filename ), 2000 );
    //m_changed = FALSE;
}

void MyWindow::fileSave_ps()
{
#if 0
	CGAL::Postscript_file_stream ps_stream(m_width, m_height ,"pm.ps");

    Qt_widget_demo_tab	*w_demo_p1 = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());

	switch ( w_demo_p1->traits_type ) {
		case SEGMENT_TRAITS:
		{
			Qt_widget_segment_tab *w_demo_p = static_cast<Qt_widget_segment_tab *> (myBar->currentPage());
			
			//ps_stream.set_line_width(w_demo_p->m_line_width);
			CGAL::Pm_drawer<Seg_arr, CGAL::Postscript_file_stream> drawer(ps_stream);
			ps_stream << CGAL::BLUE;
			drawer.draw_halfedges(w_demo_p->segment_arr.halfedges_begin(), w_demo_p->segment_arr.halfedges_end());
			ps_stream << CGAL::RED;
			drawer.draw_vertices(w_demo_p->segment_arr.vertices_begin(), w_demo_p->segment_arr.vertices_end());
			break;
		}
		case POLYLINE_TRAITS:
		{
			//Qt_widget_polyline_tab *w_demo_p = static_cast<Qt_widget_polyline_tab *> (myBar->currentPage());
			//outFile << w_demo_p->polyline_arr;
			break;
		}
		case CONIC_TRAITS:
		{
			//Qt_widget_conic_tab *w_demo_p = static_cast<Qt_widget_conic_tab *> (myBar->currentPage());
			//outFile << w_demo_p->conic_arr;
			break;
		}
	}  

#endif 
}


void MyWindow::print()
{
    Qt_widget_demo_tab	*w_demo_p1 = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());
	w_demo_p1->print_to_ps();
}





void MyWindow::overlay()
{
    OverLay *optionsForm = new OverLay( myBar , this ,number_of_tabs );
    
    if ( optionsForm->exec() ) 
	{	
		if (optionsForm->arrComboBox1->currentItem() == optionsForm->arrComboBox2->currentItem()
			|| optionsForm->arrComboBox1->currentItem() == optionsForm->arrComboBox3->currentItem()
			|| optionsForm->arrComboBox2->currentItem() == optionsForm->arrComboBox3->currentItem())
			QMessageBox::information( this, my_title_string,"This is not a clever idea!!!");
			
		else
		{
			int index1 = optionsForm->arrComboBox1->currentItem();
			int index2 = optionsForm->arrComboBox2->currentItem();
			int index3 = optionsForm->arrComboBox3->currentItem();
            Qt_widget_demo_tab	*w_demo_p1 = static_cast<Qt_widget_demo_tab *> (myBar->page( index1 ));
			Qt_widget_demo_tab	*w_demo_p2 = static_cast<Qt_widget_demo_tab *> (myBar->page( index2 ));
			Qt_widget_demo_tab	*w_demo_p3 = static_cast<Qt_widget_demo_tab *> (myBar->page( index3 ));

			if (w_demo_p1->traits_type != w_demo_p2->traits_type || w_demo_p1->traits_type != w_demo_p3->traits_type)
				QMessageBox::information( this, my_title_string,"Can not union arr of different traits");
			else make_overlay( index1 , index2 , index3 , w_demo_p2->traits_type);
		}
        //std::cout << optionsForm->arrComboBox1->currentItem() << std::endl;
    }
    delete optionsForm;
}

void MyWindow::make_overlay( int index1 , int index2 , int index3 , TraitsType t)
{
	switch ( t ) 
	{
		case SEGMENT_TRAITS:
			return;
			break;
		case POLYLINE_TRAITS:
		{
			add_polyline_tab();
			Qt_widget_polyline_tab *w_demo_p_new = static_cast<Qt_widget_polyline_tab *> (myBar->currentPage());
            Qt_widget_polyline_tab *w_demo_p1 = static_cast<Qt_widget_polyline_tab *> (myBar->page( index1 ));
			Qt_widget_polyline_tab *w_demo_p2 = static_cast<Qt_widget_polyline_tab *> (myBar->page( index2 ));
			Qt_widget_polyline_tab *w_demo_p3 = static_cast<Qt_widget_polyline_tab *> (myBar->page( index3 ));

			std::list<Pm_pol_2> pol_list;
			Pm_pol_const_iter itp;
			
			*w_demo_p_new << CGAL::GREEN;
			*w_demo_p_new << CGAL::LineWidth(3);
			for (itp = w_demo_p1->list_of_curves.begin(); itp != w_demo_p1->list_of_curves.end(); ++itp)
			{
				Pm_pol_2 * seg = new Pm_pol_2( **itp );
				w_demo_p_new->list_of_curves.push_back(seg);
				pol_list.push_back(*seg);
				*w_demo_p_new << *(*itp);
			}
			
			*w_demo_p_new << CGAL::BLUE;
			for (itp = w_demo_p2->list_of_curves.begin(); itp != w_demo_p2->list_of_curves.end(); ++itp)
			{
				Pm_pol_2 * seg = new Pm_pol_2( **itp );
				w_demo_p_new->list_of_curves.push_back(seg);
				pol_list.push_back(*seg);
				*w_demo_p_new << *(*itp);
			}
			
			*w_demo_p_new << CGAL::BLACK;
			for (itp = w_demo_p3->list_of_curves.begin(); itp != w_demo_p3->list_of_curves.end(); ++itp)
			{
				Pm_pol_2 * seg = new Pm_pol_2( **itp );
				w_demo_p_new->list_of_curves.push_back(seg);
				pol_list.push_back(*seg);
				*w_demo_p_new << *(*itp);
			}

			w_demo_p_new->polyline_arr.insert(pol_list.begin(),pol_list.end());

			*w_demo_p_new << CGAL::RED;
			*w_demo_p_new << CGAL::DISC;
			// Go over all vertices and for each vertex print the ID numbers of the
			// base curves that go through it.
			Pol_arr::Vertex_iterator   vit;
			for (vit = w_demo_p_new->polyline_arr.vertices_begin(); vit != w_demo_p_new->polyline_arr.vertices_end(); vit++)
			{
				Pol_arr::Halfedge_around_vertex_circulator 
				eit, first = (*vit).incident_halfedges();

				eit = first;

				int ind1;
				int ind2 = (*eit).curve().get_data().m_index;
				do 
				{
					ind1 = (*eit).curve().get_data().m_index;

					// Keep track of IDs we haven't seen before.
					if (ind1 != ind2)
					{
						const Pm_pol_point_2& p = (*vit).point();
						*w_demo_p_new << p;
						break;
					}

					eit++;

				} while (eit != first);
			}

				// allow the user to see the overlay
			w_demo_p_new->current_state = old_state;

		/*	w_demo_p_new->m_xmin = std::min(w_demo_p1->m_xmin,std::min(w_demo_p2->m_xmin,w_demo_p3->m_xmin)); 
			w_demo_p_new->m_xmax = std::max(w_demo_p1->m_xmax,std::min(w_demo_p2->m_xmax,w_demo_p3->m_xmax));
			w_demo_p_new->m_ymin = std::min(w_demo_p1->m_ymin,std::min(w_demo_p2->m_ymin,w_demo_p3->m_ymin));
			w_demo_p_new->m_ymax = std::max(w_demo_p1->m_ymax,std::min(w_demo_p2->m_ymax,w_demo_p3->m_ymax));
			w_demo_p_new->set_window(w_demo_p_new->m_xmin , w_demo_p_new->m_xmax , w_demo_p_new->m_ymin , w_demo_p_new->m_ymax);*/
	
			// update new planner map index
			Pm_pol_iter iter;
			for (iter = w_demo_p_new->list_of_curves.begin(); iter != w_demo_p_new->list_of_curves.end(); ++iter)
			{
				Curve_pol_data cd = (**iter).get_data();
				cd.m_type = Curve_pol_data::INTERNAL;
				cd.m_index = w_demo_p_new->index;
				(**iter).set_data( cd );
			}

			break;
		}
		case CONIC_TRAITS:
			return;
			break;

		}


}







#include "demo1.moc"


int main(int argc, char **argv)
{
  QApplication app( argc, argv );
  MyWindow widget(700,700); // physical window size
  app.setMainWidget(&widget);
  widget.setCaption(my_title_string);
  widget.setMouseTracking(TRUE);
  widget.show();
  return app.exec();
}

#endif // CGAL_USE_QT

//std::cout << xseg << std::endl;
//std::fflush(stdout);



