/*! demo_tab.h contain the definetion and implementation of 
 *  the demo tab classes and the tabs traits classes         
 */
#include <qrect.h>
#include <qcursor.h>

#include "cgal_types1.h"
#include <CGAL/IO/pixmaps/hand.xpm>
#include <CGAL/IO/pixmaps/holddown.xpm>

/*! class Qt_widget_base_tab - inherits from CGAL::Qt_widget
 *  contain all data members that are not part of the traits 
 */
class Qt_widget_base_tab : public CGAL::Qt_widget
{
public:

	Qt_widget_base_tab(TraitsType  t , QWidget *parent = 0, int tab_number = 1):
		CGAL::Qt_widget( parent ),
		index(tab_number),
		snap_mode(NONE),
	    mode(INSERT),
		m_line_width(2),
		close_point(false),
		first_time(true),
		active(false),
		traits_type(t),
		bbox(CGAL::Bbox_2(-10, -10, 10, 10)),
		wasrepainted(true), 
		on_first(false)
	{
		*this << CGAL::LineWidth(2) << CGAL::BackgroundColor (CGAL::WHITE);
		set_window(-10, 10, -10, 10);
		setMouseTracking(TRUE);
	}

	/*! current_state - indecates when a tab state is changed */
	int current_state;
	/*! index - each tab has a uniqe index */
	int index;
	/*! pl_point - the point for point location */
	Coord_point pl_point;
	/*! snap_mode - current snap mode (none, grid or closest point)  */
	SnapMode snap_mode;
	/*! mode - current mode (insert,delete or point location) */
	Mode mode;
	/*! m_line_width - line width */
	int m_line_width;
	/*! close_point - boolean flag, true if found a close point in the
	 *                get_point function */
	bool close_point;
	/*! first_time - true when it is the first mouse click of the object */
	bool first_time;
	/*! active - true if the first point was inserted */
	bool active;
	/*! traits_type - the actual tab traits */
	TraitsType   traits_type;
	/*! bbox - bounding box */
    CGAL::Bbox_2 bbox;
	/*! oldcursor - when changing the cursor we want to save the old one */
	QCursor oldcursor;
	/*! for use of drag mode */
	int   first_x, first_y;
	int   x2, y2;
	bool	wasrepainted;
	bool	on_first;
};

/*! template class Qt_widget_demo_tab gets a Tab_traits class as 
 *  a template parameter. all the Tab_traits classes must support
 *  a set of demands.
 */
template <class Tab_traits>
class Qt_widget_demo_tab : public Qt_widget_base_tab
{
private:
	typedef typename Tab_traits::Curves_list Curves_list;
	typedef typename Tab_traits::Curves_arr Curves_arr;
	typedef typename Tab_traits::Curve Curve;
	typedef typename Tab_traits::Xcurve Xcurve;
	typedef typename Tab_traits::Base_curve Base_curve;
	typedef typename Tab_traits::Pm_curve_iter Pm_curve_iter;
	typedef typename Tab_traits::Pm_curve_const_iter Pm_curve_const_iter;
	typedef typename Tab_traits::Locate_type Locate_type;
	typedef typename Tab_traits::Pm_point_2 Pm_point_2;
	typedef typename Tab_traits::Halfedge_handle Halfedge_handle;
	typedef typename Tab_traits::Face_handle Face_handle;
	typedef typename Tab_traits::Ccb_halfedge_circulator Ccb_halfedge_circulator;
	typedef typename Tab_traits::Holes_iterator Holes_iterator;
	typedef typename Tab_traits::Halfedge_iterator Halfedge_iterator;
	typedef typename Tab_traits::Hafledge_list Hafledge_list;
	typedef typename Tab_traits::Hafledge_list_iterator Hafledge_list_iterator;
	typedef typename Tab_traits::Pm Pm;

public:
	/*! m_tab_traits - the traits object */
	Tab_traits m_tab_traits;
	/*! m_curves_list - list of inserted curves */
	Curves_list m_curves_list;
	/*! m_curves_arr - the tab planar map */
	Curves_arr m_curves_arr;

	/*! constractor */
	Qt_widget_demo_tab(TraitsType  t , QWidget *parent = 0 , int tab_number = 1):
	    Qt_widget_base_tab(t , parent, tab_number)
	{}
	
	/*! destructor - clean the Curves list */
	~Qt_widget_demo_tab() 
	{
		Pm_curve_iter itp;
		for (itp = m_curves_list.begin(); itp != m_curves_list.end(); ++itp)
			delete (*itp);
	}

	/*! draw - called everytime something changed, draw the PM and mark the 
	 *         point location if the mode is on. */
	void draw(const QColor c = Qt::green)
	{
		if (snap_mode == GRID)
			draw_grid();

		setColor(c);
		(*this) << CGAL::LineWidth(m_line_width);
		Pm_curve_const_iter itp;
		for (itp = m_curves_list.begin(); itp != m_curves_list.end(); ++itp)
		{
			m_tab_traits.draw_curve(this , **itp );
		}

		if (mode == POINT_LOCATION && 
			! (m_curves_arr.halfedges_begin() == m_curves_arr.halfedges_end() ) ) 
		{
			(*this) << CGAL::LineWidth(3);
			setColor(Qt::yellow);

			Locate_type lt;
			Pm_point_2 temp_p (pl_point.x(), pl_point.y());
			Halfedge_handle e = m_curves_arr.locate(temp_p, lt);
			
			//color the outer face 
			Face_handle f = e->face();
			if (f->does_outer_ccb_exist()) // its an inside face
			{
				Ccb_halfedge_circulator cc=f->outer_ccb();
				do {
					m_tab_traits.draw_xcurve(this , cc->curve() );
				} while (++cc != f->outer_ccb());
			}

			//color the holes
			Holes_iterator hit, eit = f->holes_end();
			for (hit = f->holes_begin(); hit != eit; ++hit) 
			{
				Ccb_halfedge_circulator cc = *hit; 
				do 
				{
					m_tab_traits.draw_xcurve(this , cc->curve() );
					cc++;
				} while (cc != *hit);
			}
			(*this) << CGAL::LineWidth(m_line_width);
		}

		if (mode == RAY_SHOOTING && 
			! (m_curves_arr.halfedges_begin() == m_curves_arr.halfedges_end() ) ) 
		{
			Locate_type lt;
			Pm_point_2 temp_p (pl_point.x(), pl_point.y());
			Halfedge_handle e = m_curves_arr.vertical_ray_shoot(temp_p, lt , true);
	
			setColor(Qt::blue);
			(*this) << CGAL::LineWidth(1);
			Coord_point up(pl_point.x() , y_max());
			(*this) << Coord_segment(pl_point , up );

			(*this) << CGAL::LineWidth(3);
			setColor(Qt::red);

			switch (lt) {
				case (Curves_arr::VERTEX):
					//std::cout << "VERTEX" << std::endl;
					*this << e->target()->point();
					break;
				case (Curves_arr::UNBOUNDED_VERTEX) :
					//std::cout << "UNBOUNDED_VERTEX" << std::endl;
					*this << e->target()->point();
					break;
				case (Curves_arr::EDGE):
					//std::cout << "EDGE" << std::endl;
					m_tab_traits.draw_xcurve(this , e->curve() );
					break;
	            case (Curves_arr::UNBOUNDED_EDGE) :
					//std::cout << "UNBOUNDED_EDGE" << std::endl;
					m_tab_traits.draw_xcurve(this , e->curve() );
					break;
				case (Curves_arr::FACE) :
					//std::cout << "FACE" << std::endl;
					break;
				case (Curves_arr::UNBOUNDED_FACE) :
					//std::cout << "UNBOUNDED_FACE" << std::endl;
					break;
	
			}				
						
			(*this) << CGAL::LineWidth(m_line_width);
		}

	}
	
	/*! draw_grid - draw the grid */
	void draw_grid()
	{
		//(*w_demo_p) << CGAL::GRAY;
		//setColor(Qt::deepblue);
		(*this) << CGAL::DEEPBLUE;
		//(*w_demo_p) << CGAL::RED;
		//(*w_demo_p) << CGAL::ORANGE;
		//(*w_demo_p) << CGAL::PURPLE;

		(*this) << CGAL::LineWidth(1);
		// get the edge coordinate
		int min_x = static_cast<int> (x_min());
		int max_x = static_cast<int> (x_max());
		int min_y = static_cast<int> (y_min());
		int max_y = static_cast<int> (y_max());

		// calculate cube size (minimum of 1)
		int cube_size_x = std::max(1, abs(max_x - min_x)/20);
		int cube_size_y = std::max(1, abs(max_y - min_y)/20);
		// draw the grid lines
		for (int i = min_x; i <= max_x; i += cube_size_x)
			(*this) << Coord_segment(Coord_point( i , max_y + cube_size_y),Coord_point( i , min_y - cube_size_y));
		for (int i = min_y; i <= max_y; i += cube_size_y)
			(*this) << Coord_segment(Coord_point( max_x + cube_size_x , i ),Coord_point( min_x - cube_size_x , i ));
	}

	/*! draw_curve - activate draw_curve of the traits object */
	void draw_curve( Curve c )
	{
		m_tab_traits.draw_curve(this , c );
	}

	/*! mousePressEvent - mouse click on the tab */
	void mousePressEvent(QMouseEvent *e)
	{
		if (mode == POINT_LOCATION || mode == RAY_SHOOTING)
		{
			mousePressEvent_point_location( e );
			return;
		}
		if (mode == DELETE)
		{
			remove_curve( e );
			return;
		}
		if (mode == INSERT)
		{
			Coord_type x, y;
			x_real(e->x(), x);
			y_real(e->y(), y);
			Coord_point p = get_point(x,y);

			lock();
			QColor old_color = color();
			RasterOp old_rasterop=rasterOp();
			get_painter().setRasterOp(XorROP);

			insert( e , p);
			
			setRasterOp(old_rasterop);
			setColor(old_color);
			unlock();

			return;
		}
		if (mode == DRAG)
		{
			mousePressEvent_drag(e);
			return;
		}
	}
	/*! insert - insert a curve to the planar map */
	void insert( QMouseEvent *e , Coord_point p)
	{
		if(e->button() == Qt::LeftButton && is_pure(e->state()))
		{
			if(!active)
			{
				active = true;
				m_tab_traits.first_point( p );
			} 
			else
			{
				//show the last rubber as edge of the polygon
				m_tab_traits.middle_point( p , this );
			}
		}
		// finish polyline draw with right button click 
		else if (active && e->button() == Qt::RightButton && is_pure(e->state())) {
			m_tab_traits.last_point( p , this );
		}
		
	}

	/*! mousePressEvent_point_location - creats the point location point */
	void mousePressEvent_point_location(QMouseEvent *e)
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
	
	/*! is_pure - insure no special button is pressed */
    bool is_pure(Qt::ButtonState s)
	{
        if((s & Qt::ControlButton) ||
           (s & Qt::ShiftButton) ||
           (s & Qt::AltButton))
		    return 0;
		else
			return 1;
    }

	/*! dist - the distance between 2 points */
	Coord_type dist(Coord_type x1, Coord_type y1, Coord_type x2, Coord_type y2)
	{
		return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
	}

	/*! getMid - return the closest grid point */
    Coord_type getMid(Coord_type coord, int my_min, int my_max)
    {
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

	/*! remove_curve - remove curve from the tab */
	void remove_curve(QMouseEvent *e)
	{
		if( m_curves_list.empty() )
			return;

		Coord_point p(x_real(e->x()) ,y_real(e->y()));

		bool is_first = true;
		Coord_type min_dist = 0;

		Pm_curve_iter itp;
		Pm_curve_iter closest_curve = NULL;
		for (itp = m_curves_list.begin(); itp != m_curves_list.end(); ++itp)
		{
			Coord_type dist = m_tab_traits.curve_point_distance( p, *itp);

			if (is_first || dist < min_dist)
			{
				min_dist = dist;
				closest_curve = itp;
				is_first = false;
			}
		}
		    
		// Remove curve from pmwx
		Halfedge_iterator hei;
		Hafledge_list halfedge_list;
		Hafledge_list_iterator result;

		for (hei = m_curves_arr.halfedges_begin(); hei != m_curves_arr.halfedges_end(); ++hei) 
		{
			const Xcurve & xcurve = hei->curve();
			const Base_curve * org_curve = m_tab_traits.get_origin_curve(xcurve);

			if (m_tab_traits.compare( org_curve , *closest_curve ) )
			{
				result = std::find(halfedge_list.begin(), halfedge_list.end(), hei);
				if (result == halfedge_list.end())
				{
					halfedge_list.push_back(hei->twin());
					m_curves_arr.remove_edge(hei);
				}
			}				
		}

		m_curves_list.erase(closest_curve);
		delete (*closest_curve);

		redraw();
		
	} // remove_curve

	/*! mouseMoveEvent - enable seeing the line to be drawen */
	void mouseMoveEvent(QMouseEvent *e)
	{
		if (mode == DRAG)
		{
			mouseMoveEvent_drag(e);
			return;
		}
		else if(active)
		{
			Coord_type x, y;
			x_real(e->x(), x);
			y_real(e->y(), y);
			Coord_point p = get_point(x,y);
			RasterOp old_raster = rasterOp();//save the initial raster mode
			setRasterOp(XorROP);
			lock();

			setColor(Qt::green);
			
			if(!first_time)
				m_tab_traits.draw_last_segment(this);

			m_tab_traits.draw_current_segment( p , this);

			unlock();
			setRasterOp(old_raster);
			first_time = false;
		}
	}

	/*! leaveEvent - hide the line if you leave the widget's area on the screen */
	void leaveEvent(QEvent *e)
	{
		if(active)
		{
			RasterOp old_raster = rasterOp();//save the initial raster mode
			QColor old_color = color();
			lock();
			setRasterOp(XorROP);
			setColor(Qt::green);
			m_tab_traits.draw_last_segment(this);
			setRasterOp(old_raster);
			setColor(old_color);
			unlock();
			first_time = true;
		}
	}

	/*! get_point - return a point according to the current snap mode and
	 *              recent points */
	Coord_point get_point(Coord_type x, Coord_type y)
    {
        int xmin = static_cast<int> (x_min());
		int xmax = static_cast<int> (x_max());
		int ymin = static_cast<int> (y_min());
		int ymax = static_cast<int> (y_max());
		Coord_type d = std::max(0.5 , (x_max() - x_min())/40);
		switch ( snap_mode ) {
			case POINT:
				{
					Coord_type min_dist = 0;
					Coord_point closest;

					if (m_curves_list.empty()) 
						return Coord_point(x , y);
						
					min_dist = m_tab_traits.closest_point(x,y,closest,this);				
					
					if (min_dist <= d)
					{
						close_point = true;
						return closest;
					}
					else
					{
						close_point = false;
						return Coord_point(x , y);
					}

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

	void mousePressEvent_drag(QMouseEvent *e)
	{
		if(e->button() == Qt::LeftButton 
		&& is_pure(e->state()))
		{
			setCursor(QCursor( QPixmap( (const char**)holddown_xpm)));
			if (!on_first){
				first_x = e->x();
				first_y = e->y();
				on_first = TRUE;
			}	
		}
	}

	void mouseReleaseEvent(QMouseEvent *e)
	{
		if(e->button() == Qt::LeftButton
			&& mode == DRAG
		&& is_pure(e->state()))
		{
		setCursor(QCursor( QPixmap( (const char**)hand_xpm)));
		double x, y, xfirst2, yfirst2;
		x_real(e->x(), x);
		y_real(e->y(), y);
		x_real(first_x, xfirst2);
		y_real(first_y, yfirst2);
				
		double	xmin, xmax, ymin, ymax, distx, disty;
		if(x < xfirst2) {xmin = x; xmax = xfirst2;}
		else {xmin = xfirst2; xmax = x;};
		if(y < yfirst2) {ymin = y; ymax = yfirst2;}
		else {ymin = yfirst2; ymax = y;};
		distx = xfirst2 - x;
		disty = yfirst2 - y;
		move_center(distx, disty);
		on_first = FALSE;
		}
	}

	void mouseMoveEvent_drag(QMouseEvent *e)
	{
		char tempc1[130], tempc2[40];
		double xcoord, ycoord;
		if(on_first)
		{
		int x = e->x();
		int y = e->y();
		//save the initial raster mode
		RasterOp old = rasterOp();	
		setRasterOp(XorROP);
		lock();
			setColor(Qt::gray);
		if(!wasrepainted) {
			x_real(x2 - first_x, xcoord);
			x_real(y2 - first_y, ycoord);
			CGAL_CLIB_STD::sprintf(tempc1, " dx=%20.6f", xcoord);
			CGAL_CLIB_STD::sprintf(tempc2, ", dy=%20.6f", ycoord);
			strcat(tempc1, tempc2);
			get_painter().drawLine(first_x, first_y, x2, y2);
			setColor(Qt::green);
			get_painter().drawText(x2, y2, tempc1, 49);
			setColor(Qt::gray);
		}
		x_real(x - first_x, xcoord);
		x_real(y - first_y, ycoord);
		CGAL_CLIB_STD::sprintf(tempc1, " dx=%20.6f", xcoord);
		CGAL_CLIB_STD::sprintf(tempc2, ", dy=%20.6f", ycoord);
		strcat(tempc1, tempc2);
		get_painter().drawLine(first_x, first_y, x, y);
		setColor(Qt::green);
		get_painter().drawText(x, y, tempc1, 49);
		unlock();
		setRasterOp(old);

      //save the last coordinates to redraw the screen
      x2 = x;
      y2 = y;
      wasrepainted = FALSE;
    }
  };


};
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

/*! class Segment_tab_traits defines the segment traits */
class Segment_tab_traits 
{
public:
	typedef Pm_seg_list Curves_list;
	typedef Seg_arr Curves_arr;
	typedef Pm_seg_iter Pm_curve_iter;
	typedef Pm_seg_const_iter Pm_curve_const_iter;
	typedef Seg_locate_type Locate_type;
	typedef Pm_seg_point_2 Pm_point_2;
	typedef Seg_halfedge_handle Halfedge_handle;
	typedef Pm_seg_2 Curve;
	typedef Pm_xseg_2 Xcurve;
    typedef Pm_base_seg_2 Base_curve;
	typedef Seg_face_handle Face_handle;
	typedef Seg_ccb_halfedge_circulator Ccb_halfedge_circulator;
	typedef Seg_holes_iterator Holes_iterator;
	typedef Seg_pm::Halfedge_iterator Halfedge_iterator;
	typedef std::list<Halfedge_iterator> Hafledge_list;
	typedef Hafledge_list::iterator Hafledge_list_iterator;
	typedef Curves_arr::Vertex_iterator Vertex_iterator;
	typedef Seg_pm Pm;

	/*! constructor */
	Segment_tab_traits()
	{}

	/*! distructor */
	~Segment_tab_traits()
	{}
    
	/*! draw_curve - use Qt_Widget operator to draw */
	void draw_curve(Qt_widget_demo_tab<Segment_tab_traits> * w , Curve c )
	{
		(*w) << c;
	}
	
	/*! draw_xcurve - use Qt_Widget operator to draw */
	void draw_xcurve(Qt_widget_demo_tab<Segment_tab_traits> * w , Xcurve c )
	{
		(*w) << c;
	}
	 
	/*! first_point - a first point of inserted sgment */
	void first_point( Coord_point p )
	{
		m_p1 = m_p2 = p;
	}

	/*! middle_point - the last point of a segment */
	void middle_point( Coord_point p , Qt_widget_demo_tab<Segment_tab_traits> * w)
	{
		if(m_p1.x() != p.x() || m_p1.y() != p.y()) 
		{
			get_segment( Coord_segment( m_p1 , p ) , w );
			w->active = false;
			//w->redraw();  // not working so I use new_object insted
			w->new_object(make_object(Coord_segment(m_p1 , p)));
		} 
	}

	/*! last_point - meaningless for segments */
	void last_point( Coord_point p , Qt_widget_demo_tab<Segment_tab_traits> * w )
	{
		return;
	}

	/*! get_segment - create a new segment, insert him into curves_list and planar map */
	void get_segment( Coord_segment coord_seg , Qt_widget_demo_tab<Segment_tab_traits> * w)
	{
		const Coord_point & coord_source = coord_seg.source();
		const Coord_point & coord_target = coord_seg.target();
		Pm_seg_point_2 source(coord_source.x(), coord_source.y());
		Pm_seg_point_2 target(coord_target.x(), coord_target.y());
		Pm_base_seg_2 * base_seg_p = new Pm_base_seg_2(source, target);
		Curve_data cd;
		cd.m_type = Curve_data::LEAF;
		cd.m_index = w->index;
		cd.m_ptr.m_curve = base_seg_p;
		Pm_seg_2 * seg = new Pm_seg_2( *base_seg_p, cd );
		w->m_curves_list.push_back(seg);
		w->m_curves_arr.insert(*seg);
	}

	/*! curve_point_distance - return the distance between a point and a segment */
	Coord_type curve_point_distance(Coord_point p, Curve * c)
	{
		const Pm_seg_point_2 & source = (*c).source();
		const Pm_seg_point_2 & target = (*c).target();

		Coord_type x1 = CGAL::to_double(source.x());
		Coord_type y1 = CGAL::to_double(source.y());

		Coord_type x2 = CGAL::to_double(target.x());
		Coord_type y2 = CGAL::to_double(target.y());

		Coord_point coord_source(x1 , y1);
		Coord_point coord_target(x2 , y2);
		Coord_segment coord_seg(coord_source, coord_target);
		return CGAL::squared_distance( p, coord_seg);
	}

	/*! compare - compare between base curve and curve by their points */
	bool compare(const Base_curve * org_curve , Curve * closest)
	{
		return ((*org_curve).source() == (*closest).source() && (*org_curve).target() == (*closest).target());
	}

	/*! get_origin_curve - return the origin base curve */
	const Pm_base_seg_2* get_origin_curve(const Pm_xseg_2 & xseg )
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

	/*! draw_last_segment - call from mouse move event */
	void draw_last_segment( Qt_widget_demo_tab<Segment_tab_traits> * w)
	{
		*w << Coord_segment( m_p1 , m_p2 );
	}
	
	/*! draw_current_segment - call from mouse move event */
    void draw_current_segment( Coord_point p , Qt_widget_demo_tab<Segment_tab_traits> * w)
	{
		*w << Coord_segment( m_p1 , p);
		m_p2 = p;
	}
	
	/*! closest_point - find the closest point in the planar map to a clicked point */
	Coord_type closest_point(Coord_type x, Coord_type y, Coord_point &closest, Qt_widget_demo_tab<Segment_tab_traits> * w)
	{
		bool first = true;
		Coord_type x1,y1,dt,min_dist = 0;
		Vertex_iterator vit;
		for (vit = w->m_curves_arr.vertices_begin(); vit != w->m_curves_arr.vertices_end(); vit++)
		{
			const Pm_point_2& p = (*vit).point();
			x1 = CGAL::to_double(p.x());
			y1 = CGAL::to_double(p.y());
			dt = w->dist(x1 , y1 , x , y);
			if ( dt < min_dist || first)
			{
				min_dist = dt;
				closest = Coord_point(x1 , y1);
				first = false;
			}
		}
		return min_dist;
	}

	/*! temporary points of the created segment */
	Coord_point m_p1,m_p2;
};
////////////////////////////////////////////////////////////////////////////////////////////

/*! class Polyline_tab_traits defines the polyline traits */
class Polyline_tab_traits  
{
public:
	typedef Pm_pol_list Curves_list;
	typedef Pol_arr Curves_arr;
	typedef Pm_pol_iter Pm_curve_iter;
	typedef Pm_pol_const_iter Pm_curve_const_iter;
	typedef Pol_locate_type Locate_type;
	typedef Pm_pol_point_2 Pm_point_2;
	typedef Pol_halfedge_handle Halfedge_handle;
	typedef Pm_pol_2 Curve;
	typedef Pm_xpol_2 Xcurve;
    typedef Pm_base_pol_2 Base_curve;
	typedef Pol_face_handle Face_handle;
	typedef Pol_ccb_halfedge_circulator Ccb_halfedge_circulator;
	typedef Pol_holes_iterator Holes_iterator;
	typedef Pol_pm::Halfedge_iterator Halfedge_iterator;
	typedef std::list<Halfedge_iterator> Hafledge_list;
	typedef Hafledge_list::iterator Hafledge_list_iterator;
	typedef std::vector<Pm_point_2>::iterator Point_vector_iterator;
	typedef Curve::const_iterator Curve_const_iterator;
	typedef Pol_pm Pm;

	/*! constructor */
	Polyline_tab_traits()
	{}

	/*! distructor */
	~Polyline_tab_traits()
	{}

	/*! draw_curve - go over the polyline parts and use Qt_Widget operator to draw */
	void draw_curve(Qt_widget_demo_tab<Polyline_tab_traits> * w , Curve pol )
	{
		Curve::const_iterator ps = pol.begin();
		Curve::const_iterator pt = ps; pt++;

		while (pt != pol.end()) {
			const Pm_point_2 & source = *ps;
			const Pm_point_2 & target = *pt;
			Coord_segment coord_seg = convert(source , target);
			*w << coord_seg;
			ps++; pt++;
		}
		
	}

	/*! draw_xcurve - go over the polyline parts and use Qt_Widget operator to draw */
	void draw_xcurve(Qt_widget_demo_tab<Polyline_tab_traits> * w , Xcurve pol )
	{
		Curve::const_iterator ps = pol.begin();
		Curve::const_iterator pt = ps; pt++;

		while (pt != pol.end()) {
			const Pm_point_2 & source = *ps;
			const Pm_point_2 & target = *pt;
			Coord_segment coord_seg = convert(source , target);
			*w << coord_seg;
			ps++; pt++;
		}
		
	}

	/*! first_point - a first point of inserted polyline */
	void first_point( Coord_point p )
	{
		last_of_poly = p;
		points.push_back(Pm_pol_point_2(p.x(),p.y()));	
	}

	/*! middle_point - a middle point of a polyline */
	void middle_point( Coord_point p , Qt_widget_demo_tab<Polyline_tab_traits> * w)
	{
		if (last_of_poly == p) return;
		rubber_old = p;

		points.push_back(Pm_pol_point_2(p.x(),p.y()));	
				
		*w << CGAL::WHITE;
		*w << Coord_segment(rubber, last_of_poly);
		*w << CGAL::GREEN;
		*w << Coord_segment(rubber, last_of_poly);

        last_of_poly = p;
	}

	/*! last_point - last point of the polyline, create new polyline and reset */
	void last_point( Coord_point p , Qt_widget_demo_tab<Polyline_tab_traits> * w )
	{
		get_polyline(w);
		points.clear();
		w->active = false;
		w->first_time = true;
		//w->redraw();  // not working so I use new_object insted
		w->new_object(make_object(Coord_segment(p , p)));
	}

	/*! curve_point_distance - return the distance between a point and a polyline */
	Coord_type curve_point_distance(Coord_point p, Curve * c)
	{
		Curve_const_iterator ps = c->begin();
		Curve_const_iterator pt = ps; pt++;
		bool first = true;
		Coord_type min_dist = 0;
		while (pt != c->end()) {
			const Pm_point_2 & source = *ps;
			const Pm_point_2 & target = *pt;
			Coord_segment coord_seg = convert(source , target);
			Coord_type dist = CGAL::squared_distance( p, coord_seg);
			if( dist < min_dist || first)
			{
				first = false;
				min_dist = dist;
			}
			ps++; pt++;
		}
		return min_dist;
	}

	/*! draw_last_segment - call from mouse move event */
	void draw_last_segment( Qt_widget_demo_tab<Polyline_tab_traits> * w)
	{
		*w << Coord_segment(rubber_old, last_of_poly);
	}

	/*! draw_current_segment - call from mouse move event */
    void draw_current_segment( Coord_point p , Qt_widget_demo_tab<Polyline_tab_traits> * w)
	{
		*w << Coord_segment(p, last_of_poly);
		rubber = p;
		rubber_old = p;
	}

	/*! closest_point - find the closest point in the planar map to a clicked point */
	Coord_type closest_point(Coord_type x, Coord_type y, Coord_point &closest, Qt_widget_demo_tab<Polyline_tab_traits> * w)
	{
		bool first = true;
        Coord_type x1,y1,dt,min_dist = 0;
		Halfedge_iterator heit;
		for (heit = w->m_curves_arr.halfedges_begin(); heit != w->m_curves_arr.halfedges_end(); heit++)
		{
			const Xcurve& curve = heit->curve();
			Curve_const_iterator cit;
			for (cit = curve.begin(); cit != curve.end(); cit++)
			{
				const Pm_point_2& p = *cit;
				x1 = CGAL::to_double(p.x());
				y1 = CGAL::to_double(p.y());
				dt = w->dist(x1 , y1 , x , y);
				if ( dt < min_dist || first)
				{
					min_dist = dt;
					closest = Coord_point(x1 , y1);
					first = false;
				}
			}
		}
		Point_vector_iterator it;
		for (it = points.begin(); it != points.end(); it++)
		{
			const Pm_pol_point_2& p = *it;
			x1 = CGAL::to_double(p.x());
			y1 = CGAL::to_double(p.y());
			dt = w->dist(x1 , y1 , x , y);
			if ( dt < min_dist || first)
			{
				min_dist = dt;
				closest = Coord_point(x1 , y1);
				first = false;
			}
		}
	return min_dist;
	}

	/*! compare - compare between base curve and curve by their points */
	bool compare(const Pm_base_pol_2 * org_seg , const Pm_pol_2 * closest)
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

	/*! get_origin_curve - return the origin base curve */
	const Pm_base_pol_2* get_origin_curve(const Pm_xpol_2 & xseg )
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

private:

	/*! get_polyline - create a new polyline */
	void get_polyline(Qt_widget_demo_tab<Polyline_tab_traits> * w)
	{
		Pm_base_pol_2 * base_pol_p = new Pm_base_pol_2(points.begin(), points.end());
		Curve_pol_data cd;
		cd.m_type = Curve_pol_data::LEAF;
		cd.m_index = w->index;
		cd.m_ptr.m_curve = base_pol_p;
		Pm_pol_2 * pol = new Pm_pol_2( *base_pol_p, cd );
		w->m_curves_list.push_back(pol);
		w->m_curves_arr.insert( *pol );
	}

	/*! convert - convert from Pm_pol_curve to Coord_segment */
	Coord_segment convert(Pm_pol_point_2 & source , Pm_pol_point_2 & target)
	{
		Coord_type x1 = CGAL::to_double(source.x());
		Coord_type y1 = CGAL::to_double(source.y());
		Coord_point coord_source(x1, y1);

		Coord_type x2 = CGAL::to_double(target.x());
		Coord_type y2 = CGAL::to_double(target.y());
		Coord_point coord_target(x2, y2);

		return Coord_segment(coord_source, coord_target);
	}

	/*! convert - convert from const Pm_pol_curve to Coord_segment */
	Coord_segment convert(const Pm_pol_point_2 & source , const Pm_pol_point_2 & target)
	{
		Coord_type x1 = CGAL::to_double(source.x());
		Coord_type y1 = CGAL::to_double(source.y());
		Coord_point coord_source(x1, y1);

		Coord_type x2 = CGAL::to_double(target.x());
		Coord_type y2 = CGAL::to_double(target.y());
		Coord_point coord_target(x2, y2);

		return Coord_segment(coord_source, coord_target);
	}

	/*! the new point of the rubber band */
	Coord_point rubber;
	/*! the last point of the polygon */
    Coord_point last_of_poly;
	/*! the old point of the rubber band */
    Coord_point rubber_old;
	/*! container to hold the point during polyline creation */
	std::vector<Pm_pol_point_2> points;

};

/////////////////////////////////////////////////////////////////////////////////////////
/*! */
class Conic_tab_traits
{
public:
	typedef Pm_xconic_list Curves_list;
	typedef Conic_arr Curves_arr;
	typedef Pm_xconic_iter Pm_curve_iter;
	typedef Pm_xconic_const_iter Pm_curve_const_iter;
	typedef Conic_locate_type Locate_type;
	typedef Pm_conic_point_2 Pm_point_2;
	typedef Conic_halfedge_handle Halfedge_handle;
	typedef Pm_xconic_2 Curve;
	typedef Pm_xconic_2 Xcurve;
	typedef Pm_base_conic_2 Base_curve;
	typedef Conic_face_handle Face_handle;
	typedef Conic_ccb_halfedge_circulator Ccb_halfedge_circulator;
	typedef Conic_holes_iterator Holes_iterator;
	typedef Conic_pm::Halfedge_iterator Halfedge_iterator;
	typedef std::list<Halfedge_iterator> Hafledge_list;
	typedef Hafledge_list::iterator Hafledge_list_iterator;
	typedef Conic_pm Pm;

	/*! constructor */
	Conic_tab_traits()
	{}

	/*! distructor */
	~Conic_tab_traits()
	{}

	/*! draw_curve - draw the conic curve point by point */
	void draw_curve(Qt_widget_demo_tab<Conic_tab_traits> * w , Curve c )
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
			
			*w << coord_seg;
		}
		else
		{
			// If the curve is monotone, than its source and its target has the
			// extreme x co-ordinates on this curve.
			if (c.is_x_monotone())
			{
				
				bool     is_source_left = (sx < tx);
				int      x_min = is_source_left ? (*w).x_pixel(sx) : 
												(*w).x_pixel(tx);
				int      x_max = is_source_left ? (*w).x_pixel(tx) :
													(*w).x_pixel(sx);
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
					curr_x = (*w).x_real(x);
					nps = c.get_points_at_x(Pm_conic_point_2(curr_x, 0), ps);
					if (nps == 1)
					{
						curr_y = CGAL::to_double(ps[0].y());
						(*w) << Coord_segment( Coord_point(prev_x, prev_y) , Coord_point(curr_x, curr_y) );
						prev_x = curr_x;
						prev_y = curr_y;
							
					}
				}

				(*w) << Coord_segment( Coord_point(prev_x, prev_y) , Coord_point(end_x, end_y) );
			}
			else
			{
				// We should never reach here.
				CGAL_assertion(false);
			}
		}
	}

	/*! draw_xcurve - same as draw_curve */
	void draw_xcurve(Qt_widget_demo_tab<Conic_tab_traits> * w , Xcurve c )
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
			
			*w << coord_seg;
		}
		else
		{
			// If the curve is monotone, than its source and its target has the
			// extreme x co-ordinates on this curve.
			if (c.is_x_monotone())
			{
				
				bool     is_source_left = (sx < tx);
				int      x_min = is_source_left ? (*w).x_pixel(sx) : 
												(*w).x_pixel(tx);
				int      x_max = is_source_left ? (*w).x_pixel(tx) :
													(*w).x_pixel(sx);
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
					curr_x = (*w).x_real(x);
					nps = c.get_points_at_x(Pm_conic_point_2(curr_x, 0), ps);
					if (nps == 1)
					{
						curr_y = CGAL::to_double(ps[0].y());
						(*w) << Coord_segment( Coord_point(prev_x, prev_y) , Coord_point(curr_x, curr_y) );
						prev_x = curr_x;
						prev_y = curr_y;
							
					}
				}

				(*w) << Coord_segment( Coord_point(prev_x, prev_y) , Coord_point(end_x, end_y) );
			}
			else
			{
				// We should never reach here.
				CGAL_assertion(false);
			}
		}
	}
		
	/*! get_origin_curve - return the origin base curve */
	Pm_base_conic_2* get_origin_curve( Pm_xconic_2 & xseg )
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
//////////////////////////////////////////////////////////////////////////////////////////////////////
// those functions will be implemented when support insertion and deletion of
// conic curves

	/*! first_point - a first point of inserted conic */
	void first_point( Coord_point p )
	{}

	/*! middle_point - a middle point of a conic */
	void middle_point( Coord_point p , Qt_widget_demo_tab<Conic_tab_traits> * w)
	{}

	/*! last_point - the last point of a conic */
	void last_point( Coord_point p , Qt_widget_demo_tab<Conic_tab_traits> * w )
	{}

	/*! curve_point_distance - return the distance between a point and a conic */
	Coord_type curve_point_distance(Coord_point p, Curve * c)
	{
		Coord_type a = CGAL::to_double(p.y());
		return a;
	}

    /*! compare - compare between base curve and curve by their points */
	bool compare(const Base_curve * org_curve , Curve * closest)
	{
		return false;
	}

	/*! get_origin_curve - return the origin base curve */
	const Base_curve* get_origin_curve(const Xcurve & xseg )
	{
		return 0;
	}

	/*! draw_last_segment - call from mouse move event */
	void draw_last_segment( Qt_widget_demo_tab<Conic_tab_traits> * w)
	{}

	/*! draw_current_segment - call from mouse move event */
    void draw_current_segment( Coord_point p , Qt_widget_demo_tab<Conic_tab_traits> * w)
	{}
	
	/*! closest_point - find the closest point in the planar map to a clicked point */
	Coord_type closest_point(Coord_type x, Coord_type y, Coord_point &closest, Qt_widget_demo_tab<Conic_tab_traits> * w)
	{
		return x;
	}



};

typedef Qt_widget_demo_tab<Segment_tab_traits> Qt_widget_segment_tab;
typedef Qt_widget_demo_tab<Polyline_tab_traits> Qt_widget_polyline_tab;
typedef Qt_widget_demo_tab<Conic_tab_traits> Qt_widget_conic_tab;
