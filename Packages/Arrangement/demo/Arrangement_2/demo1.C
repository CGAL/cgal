
#include <CGAL/basic.h>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>
//#include "Qt_widget_toolbar1.h"
//#include <CGAL/IO/Qt_widget_standard_toolbar.h>
//#include "Qt_widget_standard_toolbar.h"
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/pixmaps/zoom_out.xpm>
#include <CGAL/IO/pixmaps/zoom_in.xpm>


#include <fstream>
#include <iostream>
#include <stack>
#include <set>
#include <string>
#include <list>
#include <vector>

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
#include <qfiledialog.h>
#include <qtimer.h>
#include <qtabbar.h>
#include <qtabwidget.h>
#include <qstring.h>
#include <qdialog.h>
#include <qcombobox.h>
#include <qlabel.h> 
#include <qlayout.h> 
#include <qpushbutton.h> 

#include "cgal_types1.h"

enum TraitsType { SEGMENT_TRAITS, POLYLINE_TRAITS , CONIC_TRAITS};
enum SnapMode   { NONE , GRID };
enum Mode       { INSERT , DELETE , POINT_LOCATION };

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

class MyWindow : public QMainWindow
{
    Q_OBJECT
public:
    MyWindow(int w, int h);
	~MyWindow();

private:
    void something_changed();

private slots:
    void get_new_object(CGAL::Object obj);
    void about();
    void aboutQt();
    void howto();
	void add_segment_tab();
	void add_polyline_tab();
	void remove_tab();
	void optionsSetOptions();
	void timer_done();
    //void gen_segments();
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
	void load( const QString& filename );

private:
	QTabWidget *myBar;
	int old_state;
	int tab_number;
	int number_of_tabs;
	Qt_layer_show_ch *testlayer;
	QAction *setSegmentTraits;
	QAction *setPolylineTraits;
	QAction *setNoneSnapMode;
	QAction *setGridSnapMode;
	QAction *insertMode;
	QAction *deleteMode;
	QAction *pointLocationMode;
	OptionsForm optionsForm;
};

class Qt_widget_demo_tab : public CGAL::Qt_widget
{
public:
	Qt_widget_demo_tab(QWidget *parent = 0, const char *name = 0, TraitsType t = SEGMENT_TRAITS):
	    CGAL::Qt_widget( parent , name ),
		snap_mode( GRID ),
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

	
    int         current_state;
	Coord_point pl_point;
	SnapMode    snap_mode;
	Mode        mode;
	TraitsType  traits_type;

};
////////////////////////////////////////////////////////////////////////

class Qt_widget_segment_tab : public Qt_widget_demo_tab
{
public:
	typedef std::list<Pm_seg_2*>        Pm_seg_list;
	typedef Pm_seg_list::const_iterator Pm_seg_const_iter;
	typedef Pm_seg_list::iterator       Pm_seg_iter;

	Qt_widget_segment_tab(QWidget *parent = 0, const char *name = 0):
	    Qt_widget_demo_tab( parent , name , SEGMENT_TRAITS ),
		first_point_segment(false)
	{}

	~Qt_widget_segment_tab() {}

	void draw();
	void mousePressEvent(QMouseEvent *e);
	void mouseMoveEvent(QMouseEvent *e);
	void leaveEvent(QEvent *e);
	void remove_segment(QMouseEvent *e);

	Pm_seg_list list_of_segments;
	Seg_arr segment_arr;
	
	bool first_point_segment; //true if the user left clicked once
    bool first_time_segment;  //true if the line is not drawn

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
	typedef std::list<Pm_pol_2*>        Pm_pol_list;
	typedef Pm_pol_list::const_iterator Pm_pol_const_iter;
	typedef Pm_pol_list::iterator       Pm_pol_iter;

	Qt_widget_polyline_tab(QWidget *parent = 0, const char *name = 0) :
	    Qt_widget_demo_tab( parent , name , POLYLINE_TRAITS ),
		active(false),
		first_time_polyline(true)
	{}

	~Qt_widget_polyline_tab() {}
	
	void draw();
	void mousePressEvent(QMouseEvent *e);
	void mouseMoveEvent(QMouseEvent *e);
	void leaveEvent(QEvent *e);
	void remove_segment(QMouseEvent *e);
	void get_polyline();

 	Pm_pol_list list_of_polygons;
	Pol_arr polyline_arr;

	bool active;              // true if the first point was inserted
    bool first_time_polyline; // true if it is the first time when draw the rubber band

	Pm_pol_point_2 rubber;				  // the new point of the rubber band
    Pm_pol_point_2 last_of_poly;		  // the last point of the polygon
    Pm_pol_point_2 rubber_old;			 // the old point of the rubber band
	std::vector<Pm_pol_point_2> points;

};

////////////////////////////////////////////////////////////////////////

class Qt_widget_conic_tab : public Qt_widget_demo_tab
{
public:
	Qt_widget_conic_tab(QWidget *parent = 0, const char *name = 0) :
	    Qt_widget_demo_tab( parent , name , CONIC_TRAITS )
	{}

	~Qt_widget_conic_tab() {}
	
	void draw();
	void mousePressEvent(QMouseEvent *e) {}
	void mouseMoveEvent(QMouseEvent *e) {}
	void leaveEvent(QEvent *e) {}
	void remove_segment(QMouseEvent *e) {}
	
	Conic_arr conic_arr;

};

////////////////////////////////////////////////////////////////////////
OptionsForm::OptionsForm(  QTabWidget * bar, QWidget* parent ,int number_of_tabs , const char* name,
			  bool modal, WFlags f  ): 
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
		(*w_demo_p) << CGAL::BLUE;
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
void Qt_widget_segment_tab::draw()
{
    (*this) << CGAL::GREEN;
    (*this) << CGAL::LineWidth(2);
	Pm_seg_const_iter itp;
	for (itp = list_of_segments.begin(); itp != list_of_segments.end(); ++itp)
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
	  (*this) << CGAL::LineWidth(2);
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
		if (snap_mode == GRID)
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
		}
		first_point_segment = TRUE;
    } 
	else if((e->button() == Qt::LeftButton) && is_pure(e->state()))
	{
		Coord_type x, y;
		if (snap_mode == GRID)
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
		}
		if(m_x1 != x || m_y1 != y) 
		{
	        new_object(
			make_object(Coord_segment(Coord_point(m_x1, m_y1), Coord_point(x,y))));
	        first_point_segment = FALSE;
		}    
    }
	return;
  }
 

void Qt_widget_segment_tab::remove_segment(QMouseEvent *e)
{
	if( list_of_segments.empty() )
		return;

	Coord_point p(x_real(e->x()) ,y_real(e->y()));

	bool is_first = true;
	Coord_type min_dist = 0;

	Pm_seg_iter itp;
	Pm_seg_iter it_seg = NULL;
	for (itp = list_of_segments.begin(); itp != list_of_segments.end(); ++itp)
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

	std::list<Seg_pm::Halfedge_iterator> halfedge_list;
	std::list<Seg_pm::Halfedge_iterator>::iterator result;

	for (hei = segment_arr.halfedges_begin(); hei != segment_arr.halfedges_end(); ++hei) 
	{
        const Pm_xseg_2 & xseg = hei->curve();
		const Pm_seg_2 * org_seg = xseg.get_origin();
		if (org_seg == (*it_seg))
		{
		    result = std::find(halfedge_list.begin(), halfedge_list.end(), hei);
			if (result == halfedge_list.end())
			{
				halfedge_list.push_back(hei->twin());
				segment_arr.remove_edge(hei);
			}
		}

	}

	list_of_segments.erase(it_seg);
	delete (*it_seg);

	redraw();
	
  } // remove_segment

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
      *this << Coord_segment(Coord_point(m_x1, m_y1), Coord_point(m_x2, m_y2));
      *this << Coord_segment(Coord_point(m_x1, m_y1), Coord_point(x, y));
      unlock();
      setRasterOp(old_raster);
      setColor(old_color);

      //save the last coordinates to redraw the screen
      m_x2 = x;
      m_y2 = y;
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
      *this << Coord_segment(Coord_point(m_x1, m_y1), Coord_point(m_x2, m_y2));
      setRasterOp(old_raster);
      setColor(old_color);
      unlock();
      first_time_segment = true;
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////

  void Qt_widget_polyline_tab::draw()
  {
	*this << CGAL::GREEN;
    *this << CGAL::LineWidth(2);
	Pm_pol_iter it = list_of_polygons.begin();
    for (it  = list_of_polygons.begin(); it != list_of_polygons.end(); ++it)
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

	  *this << CGAL::LineWidth(2);
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
		Pm_pol_point_2 p;
		x_real(e->x(), x);
		y_real(e->y(), y);
		if (snap_mode == GRID)
			p = Pm_pol_point_2(getMid(x, static_cast<int> (x_min()), static_cast<int> (x_max())), 
							getMid(y, static_cast<int> (y_min()), static_cast<int> (y_max())));
		else
			p = Pm_pol_point_2(x,y);

		if(!active)
		{
			active = true;
			last_of_poly = p;
			//poly.push_back(p);	
			points.push_back(p);	
		} 
		else
		{
			if (last_of_poly == p) return;
			rubber_old = p;
		
			//poly.push_back(p);	
			points.push_back(p);
			//show the last rubber as edge of the polygon
			lock();
			RasterOp old_rasterop=rasterOp();
			get_painter().setRasterOp(XorROP);
			Coord_segment seg = convert(rubber, last_of_poly);
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
	  Pm_pol_2 *poly = new Pm_pol_2( points.begin(), points.end() );
	  list_of_polygons.push_back(poly);
	  polyline_arr.insert( *poly );
	  redraw();
  }

void Qt_widget_polyline_tab::remove_segment(QMouseEvent *e)
{

      if( list_of_polygons.empty() )
		return;

      Coord_type x=static_cast<Coord_type>(x_real(e->x()));
      Coord_type y=static_cast<Coord_type>(y_real(e->y()));

      Coord_point p(x,y);
      Coord_type min_dist=100000000;
     
	  Pm_pol_iter it_closest=NULL;
      Pm_pol_iter pit = list_of_polygons.begin();

      Coord_segment closest_segment;

      while(pit!=list_of_polygons.end())
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

	for (hei = polyline_arr.halfedges_begin(); hei != polyline_arr.halfedges_end(); ++hei) 
	{
        const Pm_xpol_2 & xpol = hei->curve();
		const Pm_pol_2 * org_pol = xpol.get_origin();
		if (org_pol == (*it_closest))
		{
		    result = std::find(halfedge_list.begin(), halfedge_list.end(), hei);
			if (result == halfedge_list.end())
			{
				halfedge_list.push_back(hei->twin());
				polyline_arr.remove_edge(hei);
			}
		}

	}

    list_of_polygons.erase( it_closest );
	delete (*it_closest);

    redraw();
}
  
  void Qt_widget_polyline_tab::mouseMoveEvent(QMouseEvent *e)
  {
    if (active)
    {
      Coord_type x, y;
      x_real(e->x(), x);
      y_real(e->y(), y);

	  if (snap_mode == GRID)
          rubber = Pm_pol_point_2(getMid(x,static_cast<int> (x_min()),static_cast<int> (x_max())), 
							   getMid(y,static_cast<int> (y_min()),static_cast<int> (y_max())));
	  else
          rubber = Pm_pol_point_2(x, y);
      lock();
      RasterOp old_rasterop=rasterOp();
      get_painter().setRasterOp(XorROP);
      *this << CGAL::WHITE;    
      if(!first_time_polyline)
	  {
		  Coord_segment seg = convert(rubber_old, last_of_poly);
          *this << seg;
	  }
	  Coord_segment seg2 = convert(rubber, last_of_poly);
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
      Coord_segment seg = convert(rubber_old, last_of_poly);
      *this << seg;
      setRasterOp(old_rasterop);
      unlock();
      first_time_polyline = true;
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

   Coord_type Qt_widget_demo_tab::getMid(Coord_type coord, int my_min, int my_max)
  {
	  //int c = (my_max - my_min)/20;
	  int c = std::max(1, abs(my_max - my_min)/20);
	  for (int i = my_min - c; i <= my_max; i += c)
	  {
		  if (coord >= i && coord <= i + c)
		  {
			  Coord_type ans  = static_cast<Coord_type>(i) + static_cast<Coord_type>(c)/2;
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

	tab_number = 1;
	number_of_tabs = 0;

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

	// Snap Mode Group
	QActionGroup *snapModeGroup = new QActionGroup( this ); // Connected later
    snapModeGroup->setExclusive( TRUE );

    setNoneSnapMode = new QAction(
	    "None Snap Mode", QPixmap( (const char**)alpha_shape_xpm ),
	    "&None Snap Mode", 0 , snapModeGroup, "None Snap Mode" );
    setNoneSnapMode->setToggleAction( TRUE );

	setGridSnapMode = new QAction(
	    "Grid Snap Mode", QPixmap( (const char**)grid_xpm ),
	    "&Grid Snap Mode", 0 , snapModeGroup, "Grid Snap Mode" );
    setGridSnapMode->setToggleAction( TRUE );

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

	//create a timer for checking if somthing changed
    QTimer *timer = new QTimer( this );
    connect( timer, SIGNAL(timeout()),
           this, SLOT(timer_done()) );
    timer->start( 200, FALSE );

    // file menu
    QPopupMenu * file = new QPopupMenu( this );
    menuBar()->insertItem( "&File", file );
    //file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
    file->insertItem("Add &Segment Tab", this, SLOT(add_segment_tab()));
    file->insertItem("Add &Polyline Tab", this, SLOT(add_polyline_tab()));
	file->insertItem("Remove &Tab", this, SLOT(remove_tab()));
    file->insertSeparator();
	//file->insertItem("Union 2 Arr", this, SLOT(optionsSetOptions()));
    //file->insertItem("Print", w_demo_p, SLOT(print_to_ps()), CTRL+Key_P);
    file->insertSeparator();
    //file->insertItem( "&Open", this, SLOT( fileOpen() ) );
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
    file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), CTRL+Key_Q );

    // drawing menu
    //QPopupMenu * draw = new QPopupMenu( this );
    //menuBar()->insertItem( "&Draw", draw );
    //draw->insertItem("&Generate segments", this,
	//			SLOT(gen_segments()), CTRL+Key_G );

	 // options menu
    QPopupMenu * options = new QPopupMenu( this );
    menuBar()->insertItem( "&Options", options );
	insertMode->addTo( options );
    deleteMode->addTo( options );
    pointLocationMode->addTo( options );
	options->insertSeparator();
	setNoneSnapMode->addTo(options); 
	setGridSnapMode->addTo(options); 
	options->insertSeparator();
    setSegmentTraits->addTo(options); 
	setPolylineTraits->addTo(options); 

    // help menu
    QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Help", help );
    help->insertItem("How To", this, SLOT(howto()), Key_F1);
    help->insertSeparator();
    help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );

    //the new tools toolbar
    //newtoolbar = new Tools_toolbar(w_demo_p, this, &w_demo_p->list_of_segments);	
  
	// options toolbar

	QToolBar *optionsTools = new QToolBar( this, "options operations" );
    optionsTools->setLabel( "Options Operations" );
    insertMode->addTo( optionsTools );
    deleteMode->addTo( optionsTools );
    pointLocationMode->addTo( optionsTools );
    optionsTools->addSeparator();
    setNoneSnapMode->addTo( optionsTools );
    setGridSnapMode->addTo( optionsTools );
    optionsTools->addSeparator();
	setSegmentTraits->addTo( optionsTools );
    setPolylineTraits->addTo( optionsTools );
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
        w_segment_p->list_of_segments.clear();
		w_segment_p->segment_arr.clear();
	}
	else
	{
		Qt_widget_polyline_tab	*w_polyline_p = static_cast<Qt_widget_polyline_tab *> (myBar->currentPage());
        w_polyline_p->list_of_polygons.clear();
		w_polyline_p->polyline_arr.clear();
	}
    w_demo_p->set_window(-10, 10, -10, 10); 
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
			Pm_seg_2 * seg = new Pm_seg_2(source, target);
			w_segment_p->list_of_segments.push_back(seg);
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
	
	// initialize the new tab widget
	*widget << CGAL::LineWidth(2) << CGAL::BackgroundColor (CGAL::WHITE);
	widget->set_window(-10, 10, -10, 10);
    widget->setMouseTracking(TRUE);
	connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),this, SLOT(get_new_object(CGAL::Object)));
    widget->attach(testlayer);
	
	// add the new widget to myBar
	//myBar->addTab( widget, QString("Arr " + QString::number( tab_number ) ) );
	myBar->insertTab( widget, QString("Arr " + QString::number( tab_number ) ) , tab_number );

	myBar->setCurrentPage(myBar->indexOf(widget));
	
	tab_number++;
	number_of_tabs++;

	resize(700,700);
	
	something_changed();
		
  }

  void MyWindow::add_polyline_tab()
  {
	Qt_widget_demo_tab *widget = new Qt_widget_polyline_tab(this);
	
	// initialize the new tab widget
	*widget << CGAL::LineWidth(2) << CGAL::BackgroundColor (CGAL::WHITE);
	widget->set_window(-10, 10, -10, 10);
    widget->setMouseTracking(TRUE);
	connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),this, SLOT(get_new_object(CGAL::Object)));
    widget->attach(testlayer);
	
	// add the new widget to myBar
	//myBar->addTab( widget, QString("Arr " + QString::number( tab_number ) ) );
	myBar->insertTab( widget, QString("Arr " + QString::number( tab_number ) ) , tab_number );

	myBar->setCurrentPage(myBar->indexOf(widget));
	
	tab_number++;
	number_of_tabs++;

	resize(700,700);
	
	something_changed();
			
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

 // void MyWindow::gen_segments()
 // {
	//// We peform downcasting from QWigdet* to Qt_widget_demo_tab*, as we know that only
	//// Qt_widget_demo_tab objects are stored in the tab pages.
	//Qt_widget_demo_tab	*w_demo_p = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());

 //   w_demo_p->set_window(-10, 10, -10, 10); 
	//	// set the Visible Area to the Interval

 //   // send resizeEvent only on show.
 //   CGAL::Random_points_in_square_2<Point> g(0.5);
 //   for(int count=0; count<25; count++) 
 //   {
 //     Point p1(*g++), p2(*g++);
 //     Coord_type scale(2);
 //     Segment s( Point(p1.x()*scale,p1.y()*scale)  ,
 //                Point(p2.x()*scale,p2.y()*scale) );
 //     w_demo_p->list_of_segments.push_back(s);
 //     w_demo_p->segment_arr.insert(s);
 //   }

 //   //pl_valid = false;

 //   something_changed();
 // }
	//
  
	void MyWindow::updateTraitsType( QAction *action )
	{
		// We peform downcasting from QWigdet* to Qt_widget_demo_tab*, as we know that only
		// Qt_widget_demo_tab objects are stored in the tab pages.
		Qt_widget_demo_tab	*w_demo_p = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());

		if (tab_number > 0 && 
			(action == setSegmentTraits && w_demo_p->traits_type != SEGMENT_TRAITS ||
			 action == setPolylineTraits && w_demo_p->traits_type != POLYLINE_TRAITS))
		{
			Qt_widget_demo_tab *widget;

			if ( action == setSegmentTraits && w_demo_p->traits_type != SEGMENT_TRAITS ) 
				widget = new Qt_widget_segment_tab(this);
			else if ( action == setPolylineTraits && w_demo_p->traits_type != POLYLINE_TRAITS) 
				widget = new Qt_widget_polyline_tab(this);

			int index = myBar->currentPageIndex();
			QString label = myBar->label(index);
			myBar->removePage(myBar->currentPage());

			 //initialize the new tab widget
			*widget << CGAL::LineWidth(2) << CGAL::BackgroundColor (CGAL::WHITE);
			widget->set_window(-10, 10, -10, 10);
			widget->setMouseTracking(TRUE);
			connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),this, SLOT(get_new_object(CGAL::Object)));
			widget->attach(testlayer);
			
			// add the new widget to myBar
			myBar->insertTab( widget, label , index );

			myBar->setCurrentPage(index);
				
			resize(700,700);
			
			something_changed();
		}
	}

	void MyWindow::setTraits( TraitsType t )
	{
		//w_demo_p->traits_type = t;
		switch ( t ) {
		case SEGMENT_TRAITS:
			setSegmentTraits->setOn( TRUE );
			break;
		case POLYLINE_TRAITS:
			setPolylineTraits->setOn( TRUE );
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
		}
	}

	void MyWindow::updateMode( QAction *action )
	{
		// We peform downcasting from QWigdet* to Qt_widget_demo_tab*, as we know that only
		// Qt_widget_demo_tab objects are stored in the tab pages.
		Qt_widget_demo_tab	*w_demo_p = static_cast<Qt_widget_demo_tab *> (myBar->currentPage());

		if ( action == insertMode ) {
		w_demo_p->mode = INSERT;
		something_changed();
		}
		else if ( action == deleteMode ) {
		w_demo_p->mode = DELETE;
		something_changed();
		}
		else if ( action == pointLocationMode ) {
		w_demo_p->mode = POINT_LOCATION;
		}
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
            Qt_widget_demo_tab	*w_demo_p1 = static_cast<Qt_widget_demo_tab *> (myBar->page( optionsForm->arrComboBox1->currentItem() ));
			Qt_widget_demo_tab	*w_demo_p2 = static_cast<Qt_widget_demo_tab *> (myBar->page( optionsForm->arrComboBox2->currentItem() ));

			if (w_demo_p1->traits_type != w_demo_p2->traits_type)
				QMessageBox::information( this, my_title_string,"Can not union arr of different traits");
			//else union_arr( w_demo_p1 , w_demo_p2);
		}
        //std::cout << optionsForm->arrComboBox1->currentItem() << std::endl;
    }
    delete optionsForm;
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
  //    	Qt_widget_demo_tab *widget = new Qt_widget_conic_tab(this);
	
		//// initialize the new tab widget
		//*widget << CGAL::LineWidth(2) << CGAL::BackgroundColor (CGAL::WHITE);
		//widget->set_window(-10, 10, -10, 10);
		//widget->setMouseTracking(TRUE);
		////connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),this, SLOT(get_new_object(CGAL::Object)));
		//widget->attach(testlayer);
		//// add the new widget to myBar
		////myBar->addTab( widget, QString("Arr " + QString::number( tab_number ) ) );
		//myBar->insertTab( widget, QString("Arr " + QString::number( tab_number ) ) , tab_number );
		//myBar->setCurrentPage(myBar->indexOf(widget));
		//tab_number++;
		//number_of_tabs++;
		//resize(700,700);
		//something_changed();

		//Qt_widget_conic_tab	*w_demo_p = static_cast<Qt_widget_conic_tab *> (myBar->currentPage());

		////init(); // Make sure we have colours
  //      //m_filename = filename;
  //      //QTextStream ts( &file );

		//std::ifstream inputFile(filename);
  //      // Creates an ofstream object named inputFile

		//if (! inputFile) // Always test file open
		//{
		//	std::cout << "Error opening input file" << std::endl;
		//	return;
		//}

		////Pm_writer writer(std::cout, pm);
		////Pm_scanner scanner(inputFile, w_demo_p->conic_arr);

		//inputFile >> w_demo_p->conic_arr;
		//w_demo_p->conic_arr.read( inputFile );
		//inputFile.close();

      
    }


#include "demo1.moc"


int
main(int argc, char **argv)
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



