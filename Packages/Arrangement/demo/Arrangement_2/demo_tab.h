/*! demo_tab.h contain the definetion and implementation of 
 *  the demo tab classes and the tabs traits classes.       
 *  all the possible shared code is in Qt_widget_demo_tab where
 *  the differences is in the traits classes.
 */
#include <math.h>

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
	m_vertex_width(3),
    close_point(false),
    first_time(true),
    active(false),
    traits_type(t),
    bbox(CGAL::Bbox_2(-10, -10, 10, 10)),
    wasrepainted(true), 
    on_first(false),
    snap(false),
    grid(false),
    conic_type(SEGMENT),
    cube_size(1),
	ray_shooting_direction(true),
	remove_org_curve(true),
	empty(true),
	first_time_merge(true)
  {
    *this << CGAL::LineWidth(2) << CGAL::BackgroundColor (CGAL::WHITE);
    set_window(-10, 10, -10, 10);
    setMouseTracking(TRUE);
    
    colors[1] = Qt::blue;
    colors[2] = Qt::gray;
    colors[3] = Qt::green;
    colors[4] = Qt::cyan;
    colors[5] = Qt::magenta;
    colors[6] = Qt::black;
    colors[7] = Qt::darkGreen;
    colors[8] = Qt::darkBlue;
    colors[9] = Qt::darkMagenta;
    colors[10] = Qt::darkCyan;
    colors[11] = Qt::yellow;
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
  /*! m_vertex_width - vertic width */
  int m_vertex_width;
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
  /*! for use of drag mode */
  int   first_x, first_y;
  int   x2, y2;
  bool    wasrepainted;
  bool    on_first;
  /*! array of colors */
  QColor colors[20];
  /*! snap flag */
  bool snap;
  /*! grid flag */
  bool grid;
  /*! conic insert type */
  ConicType conic_type;
  /*! grid cube size */
  int cube_size;
  /*! ray shooting direction */
  bool ray_shooting_direction; // true for up
  /*! remove all original curve or only a part */
  bool remove_org_curve;
  /*! true if pm is empty */
  bool empty;
  /*! true when it is the first time in merge mode */
  bool first_time_merge;
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
  typedef typename Tab_traits::Traits Traits;
  typedef typename Tab_traits::Curve Curve;
  typedef typename Tab_traits::Xcurve Xcurve;
  typedef typename Tab_traits::Base_curve Base_curve;
  typedef typename Tab_traits::Data Data;
  typedef typename Tab_traits::Pm_curve_iter Pm_curve_iter;
  typedef typename Tab_traits::Pm_curve_const_iter Pm_curve_const_iter;
  typedef typename Tab_traits::Locate_type Locate_type;
  typedef typename Tab_traits::Pm_point_2 Pm_point_2;
  typedef typename Tab_traits::Halfedge_handle Halfedge_handle;
  typedef typename Tab_traits::Face_handle Face_handle;
  typedef typename Tab_traits::Ccb_halfedge_circulator 
  Ccb_halfedge_circulator;
  typedef typename Tab_traits::Holes_iterator Holes_iterator;
  typedef typename Tab_traits::Halfedge_iterator Halfedge_iterator;
  typedef typename Tab_traits::Hafledge_list Hafledge_list;
  typedef typename Tab_traits::Hafledge_list_iterator 
  Hafledge_list_iterator;
  typedef typename Tab_traits::Pm Pm;
  typedef typename Tab_traits::Vertex_iterator Vertex_iterator;
  typedef typename Tab_traits::Halfedge_around_vertex_circulator
  Halfedge_around_vertex_circulator;
  typedef typename Tab_traits::Edge_iterator Edge_iterator;
  
public:
  /*! m_tab_traits - the traits object */
  Tab_traits m_tab_traits;
  /*! m_curves_arr - the tab planar map */
  Curves_arr m_curves_arr;
  /*! Original Traits */
  Traits m_traits;
  /*! the curve to be merged */
  Halfedge_iterator closest_curve;
  /*! the seconed curve to be merged */
  Halfedge_iterator seconed_curve;
  /*! the first point in the split curve */
  Pm_point_2 split_point;
  
  /*! constractor 
   *\ param t - widget traits type
   *\ param parent - widget parent window
   *\ param tab_number - widget program index
   */
  Qt_widget_demo_tab(TraitsType  t , QWidget *parent = 0 , int tab_number = 1):
    Qt_widget_base_tab(t , parent, tab_number)
  {}
  
  /*! destructor - clean the Curves list */
  ~Qt_widget_demo_tab() 
  {}
  
  /*! draw - called everytime something changed, draw the PM and mark the 
   *         point location if the mode is on. */
  void draw()
  {
    if (snap_mode == GRID || grid)
      draw_grid();
    
	(*this) << CGAL::LineWidth(m_line_width);
    Pm_curve_const_iter itp;
    
    for (Edge_iterator ei = m_curves_arr.edges_begin(); 
         ei != m_curves_arr.edges_end(); ++ei) 
    {
      int i = (ei->curve()).get_data().m_index;
      setColor(colors[i]);
      m_tab_traits.draw_xcurve(this , ei->curve() );
    }
    
    // Go over all vertices and for each vertex check the 
    // index numbers of the base curves that go through 
    // it and paint the point if they are different
    *this << CGAL::DISC;
    (*this) << CGAL::LineWidth(m_vertex_width);

    Vertex_iterator   vit;
    for (vit = m_curves_arr.vertices_begin(); 
         vit != m_curves_arr.vertices_end(); vit++)
    {
      Halfedge_around_vertex_circulator 
        eit, first = (*vit).incident_halfedges();
      
      eit = first;
	  setColor(Qt::red);
      int ind1;
      int ind2 = (*eit).curve().get_data().m_index;
      do 
      {
        ind1 = (*eit).curve().get_data().m_index;
        
        // Keep track of IDs we haven't seen before.
        if (ind1 != ind2)
        {
          const Pm_point_2& p = (*vit).point();
          *this << p;
          break;
        }
        
        eit++;
        
		if (eit == first)
		{
		  setColor(colors[ind1]);
		  const Pm_point_2& p = (*vit).point();
          *this << p;
		}
      } while (eit != first);
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
      Halfedge_handle e = 
		  m_curves_arr.vertical_ray_shoot(temp_p, lt ,ray_shooting_direction);
      
      setColor(Qt::black);
      (*this) << CGAL::LineWidth(1);

	  if (ray_shooting_direction)
	  {
	    if (lt == Curves_arr::UNBOUNDED_FACE)
		{
		  Coord_point up(pl_point.x() , y_max());
		  (*this) << Coord_segment(pl_point , up );
		}
		else // we shoot something
		{
		  Pm_point_2 p1c1(pl_point.x() , y_max());
		  Pm_point_2 p2c1(pl_point.x() , pl_point.y());
		  const Xcurve c1 = m_tab_traits.curve_make_x_monotone(p1c1 , p2c1);
		  const Xcurve c2 = e->curve();
		  Pm_point_2 p1;
		  Pm_point_2 p2;
      
		  m_traits.nearest_intersection_to_right(c1, c2, p2c1, p1, p2);
	      Coord_type x1 = CGAL::to_double(p1.x());
		  Coord_type y1 = CGAL::to_double(p1.y());
		  Coord_point up(x1,y1);
		  (*this) << Coord_segment(pl_point , up );
		}
	  }
	  else // down ray shooting
	  {
	    if (lt == Curves_arr::UNBOUNDED_FACE)
		{
		  Coord_point down(pl_point.x() , y_min());
		  (*this) << Coord_segment(pl_point , down );
		}
		else // we shoot something
		{
		  Pm_point_2 p1c1(pl_point.x() , y_min());
		  Pm_point_2 p2c1(pl_point.x() , pl_point.y());
		  const Xcurve c1 = m_tab_traits.curve_make_x_monotone(p1c1 , p2c1);
		  const Xcurve c2 = e->curve();
		  Pm_point_2 p1;
		  Pm_point_2 p2;
      
		  m_traits.nearest_intersection_to_left(c1, c2, p2c1, p1, p2);
	      Coord_type x1 = CGAL::to_double(p1.x());
		  Coord_type y1 = CGAL::to_double(p1.y());
		  Coord_point up(x1,y1);
		  (*this) << Coord_segment(pl_point , up );
		}
	  }
	    
      
      (*this) << CGAL::LineWidth(3);
      setColor(Qt::red);
      
      switch (lt) {
       case (Curves_arr::VERTEX):
        *this << e->target()->point();
        break;
       case (Curves_arr::UNBOUNDED_VERTEX) :
        *this << e->target()->point();
        break;
       case (Curves_arr::EDGE):
        m_tab_traits.draw_xcurve(this , e->curve() );
        break;
       case (Curves_arr::UNBOUNDED_EDGE) :
        m_tab_traits.draw_xcurve(this , e->curve() );
        break;
       case (Curves_arr::FACE) :
        break;
       case (Curves_arr::UNBOUNDED_FACE) :
        break;
        
      }                
      
      (*this) << CGAL::LineWidth(m_line_width);
    }
    
  }
  
  /*! draw_grid - draw the grid */
  void draw_grid()
  {
    (*this) << CGAL::DEEPBLUE;
    (*this) << CGAL::LineWidth(1);
    // get the edge coordinate
    int min_x = static_cast<int> (x_min());
    int max_x = static_cast<int> (x_max());
    int min_y = static_cast<int> (y_min());
    int max_y = static_cast<int> (y_max());
    
    // calculate cube size (minimum of 1)
    //int cube_size_x = std::max(1, abs(max_x - min_x)/20);
    //int cube_size_y = std::max(1, abs(max_y - min_y)/20);
    int cube_size_x = cube_size;
    int cube_size_y = cube_size;
    // draw the grid lines
    for (int i = min_x; i <= max_x; i += cube_size_x)
      (*this) << Coord_segment(Coord_point( i , max_y + cube_size_y),
                               Coord_point( i , min_y - cube_size_y));
    for (int i = min_y; i <= max_y; i += cube_size_y)
      (*this) << Coord_segment(Coord_point( max_x + cube_size_x , i ),
                               Coord_point( min_x - cube_size_x , i ));
  }
  
  /*! mousePressEvent - mouse click on the tab 
   *\ param e - mouse click event
   */
  void mousePressEvent(QMouseEvent *e)
  {
    QCursor old = cursor();
    setCursor(Qt::WaitCursor);
    if (mode == POINT_LOCATION || mode == RAY_SHOOTING)
    {
      mousePressEvent_point_location( e );
	  setCursor(old);
      return;
    }
    if (mode == DELETE)
    {
      remove_curve( e );
	  if( m_curves_arr.number_of_vertices() == 0 )
	    empty = true;
	  setCursor(old);
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
      if( m_curves_arr.number_of_vertices() > 0 )
	    empty = false;
	  setCursor(old);
      return;
    }
    if (mode == DRAG)
    {
      mousePressEvent_drag(e);
	  setCursor(old);
      return;
    }
	if (mode == MERGE)
    {
      mousePressEvent_merge(e);
	  setCursor(old);
      return;
    }
	if (mode == SPLIT)
    {
      Coord_type x, y;
      x_real(e->x(), x);
      y_real(e->y(), y);
      Coord_point p = get_point(x,y);
      
      lock();
      QColor old_color = color();
      RasterOp old_rasterop=rasterOp();
      get_painter().setRasterOp(XorROP);
      
      split( e , p);
      
      setRasterOp(old_rasterop);
      setColor(old_color);
      unlock();
	  first_time = true;
	  redraw();
	  setCursor(old);
      return;
    }
  }
  
  /*! insert - insert a curve to the planar map 
   *\ param e - mouse click event
   *\ param p - the pressed point
   */
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
    else if (active && e->button() == Qt::RightButton && is_pure(e->state())) 
	{
      m_tab_traits.last_point( p , this );
    }    
  }
  
  /*! split - split a xcurve in to 2 xcurves  
   *\ param e - mouse click event
   *\ param p - the pressed point
   */
  void split( QMouseEvent *e , Coord_point p)
  {
    if(e->button() == Qt::LeftButton && is_pure(e->state()))
    {
      if(!active)
      {
        active = true;
        m_tab_traits.first_point( p );
		split_point = Pm_point_2( p.x() , p.y() );
      } 
      else
      {	    
	    active = false;
		Pm_point_2 split_point2 = Pm_point_2( p.x() , p.y() );
        const Xcurve split_curve = 
			m_tab_traits.curve_make_x_monotone(split_point , split_point2);
	    Pm_point_2 p1;
	    Pm_point_2 p2;
		Pm_point_2 p_right;
		if (split_point.x() < split_point2.x())
		  p_right = split_point;
		else
		  p_right = split_point2;
		Halfedge_iterator hei;
		for (hei = m_curves_arr.halfedges_begin();
         hei != m_curves_arr.halfedges_end(); ++hei) 
        {
		  const Xcurve & xcurve = hei->curve();
		  m_tab_traits.draw_xcurve(this, xcurve);
		  if (m_traits.nearest_intersection_to_right(split_curve, 
			  xcurve, p_right, p1, p2))
			  break;

		}		
		m_tab_traits.draw_xcurve(this, hei->curve());
		if (hei != m_curves_arr.halfedges_end())
		  m_tab_traits.split_edge(hei , p1 , this);
		
	  }
    }    
  }

  /*! mousePressEvent_point_location - creats the point location point 
   *\ param e - mouse click event
   */
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
  
  /*! is_pure - insure no special button is pressed 
   *\ param s - keyboard modifier flags that existed 
   *  immediately before the event occurred.
   *\ return true if one of them existed, false otherway.
   */
  bool is_pure(Qt::ButtonState s)
  {
    if((s & Qt::ControlButton) ||
       (s & Qt::ShiftButton) ||
       (s & Qt::AltButton))
      return 0;
    else
      return 1;
  }
  
  /*! dist  
   *\ param x1,y1,x2,y2 - points coordinates
   *\ return the distance between 2 points 
   */
  Coord_type dist(Coord_type x1, Coord_type y1, Coord_type x2, Coord_type y2)
  {
    return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
  }
  
  /*! getMid 
   *\ param coord - coord vulue apon the grid
   *\ return the closest grid point 
   */
  Coord_type getMid(Coord_type coord, int my_min, int my_max)
  {
    //int cube_size = std::max(1, abs(my_max - my_min)/20);
    Coord_type d = static_cast<Coord_type>(cube_size)/2;
    for (int i = my_min - cube_size; i <= my_max; i += cube_size)
    {
      Coord_type id = static_cast<Coord_type>(i); 
      if (coord >= id - d && coord <= id + d)
      {
        Coord_type ans  = static_cast<Coord_type>(i);
        return ans;
      }
    }
    return 0;
  }
  
  /*! remove_curve - remove curve from the tab 
   *\ param e - mouse click event
   */
  void remove_curve(QMouseEvent *e)
  {
    if( m_curves_arr.number_of_vertices() == 0 )
      return;
    
    Coord_point p(x_real(e->x()) ,y_real(e->y()));
    
    bool is_first = true;
    Coord_type min_dist = 0;
    Halfedge_iterator hei;
    Halfedge_iterator closest_curve;
    
    for (hei = m_curves_arr.halfedges_begin();
         hei != m_curves_arr.halfedges_end(); ++hei) 
    {
      Xcurve & xcurve = hei->curve();
      Coord_type dist = m_tab_traits.xcurve_point_distance( p, xcurve , this);
      if (is_first || dist < min_dist)
      {
        min_dist = dist;
        closest_curve = hei;
        is_first = false;
      }    
    }

	if (remove_org_curve)
	{
      const Base_curve * org_curve = m_tab_traits.get_origin_curve(closest_curve->curve());
	  Hafledge_list li;
	  Hafledge_list_iterator result;
      for (hei = m_curves_arr.halfedges_begin();
           hei != m_curves_arr.halfedges_end(); ++hei) 
	  {
        const Base_curve * curve = m_tab_traits.get_origin_curve(hei->curve());
		result = std::find(li.begin(), li.end(), hei);
	    if (curve == org_curve && result == li.end())
		{
		  li.push_back(hei->twin());
		  m_curves_arr.remove_edge(hei);
		}
	  }
	}

	else
	  m_curves_arr.remove_edge(closest_curve);
    
    redraw();
    
  } // remove_curve
  
    /*! mouseMoveEvent - enable seeing the line to be drawen 
     *\ param e - mouse click event
     */
  void mouseMoveEvent(QMouseEvent *e)
  {
    (*this) << CGAL::LineWidth(m_line_width);
    if (mode == DRAG)
    {
      mouseMoveEvent_drag(e);
      return;
    }
	if (mode == MERGE && !first_time_merge)
    {
	  if (seconed_curve != m_curves_arr.halfedges_end())
	  {
        int i = (seconed_curve->curve()).get_data().m_index;
        setColor(colors[i]);
	    m_tab_traits.draw_xcurve(this,seconed_curve->curve());
	  }
	  Coord_type x, y;
      x_real(e->x(), x);
      y_real(e->y(), y);
      Coord_point p(x,y);
      seconed_curve = m_curves_arr.halfedges_end();
	  m_tab_traits.find_close_curve(closest_curve ,seconed_curve ,p ,this ,true);
	  setColor(Qt::red);
	  m_tab_traits.draw_xcurve(this,closest_curve->curve());
	  if (seconed_curve != m_curves_arr.halfedges_end())
	  {
	    setColor(Qt::green);
	    m_tab_traits.draw_xcurve(this,seconed_curve->curve());
	  }
	  else
	  {
	    first_time_merge = true;
	   	redraw();
	  }
      return;
    }
	
    if(active)
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
  
  /*! leaveEvent - hide the line if you leave the widget's area on the screen 
   *\ param e - mouse click event
   */
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
  
  /*! get_point 
   *\ params x,y - the mouse clicked point coordinates
   *\    return a point according to the current snap mode and
   *  recent points. 
   */
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
       
       if( m_curves_arr.number_of_vertices() == 0 ) 
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
  
  /*! mousePressEvent_drag - change the Cursor on the drag mode
   *  mouse pressed event
   *\ param e - mouse click event
   */
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
  
  /*! mouseReleaseEvent - change the Cursor on the drag mode
   *  mouse pressed event and move the widget center according
   *  to the drag distance.
   *\ param e - mouse release event
   */
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
      
      double    xmin, xmax, ymin, ymax, distx, disty;
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
  
  /*! mouseMoveEvent_drag - calculate new widget position
   *\ param e - mouse release event
   */
  void mouseMoveEvent_drag(QMouseEvent *e)
  {
    if(on_first)
    {        
      int x = e->x();
      int y = e->y();
      //save the last coordinates to redraw the screen
      x2 = x;
      y2 = y;
      wasrepainted = FALSE;
    }
  };
  
  /*! mousePressEvent_merge -  merge mode
   *  mouse pressed event
   *\ param e - mouse click event
   */
  void mousePressEvent_merge(QMouseEvent *e)
  {
    if(e->button() == Qt::LeftButton 
       && is_pure(e->state()))
    {
	  if( m_curves_arr.number_of_vertices() == 0 )
      return;

      setColor(Qt::red);
      Coord_point p(x_real(e->x()) ,y_real(e->y()));
	  bool is_first = true;
      Coord_type min_dist = 0;
	  
	  if (first_time_merge)
	  {
		first_time_merge = false;
        Halfedge_iterator hei;
		//closest_curve = m_curves_arr.halfedges_end();
    
        for (hei = m_curves_arr.halfedges_begin();
             hei != m_curves_arr.halfedges_end(); ++hei) 
        {
		  Vertex_iterator   vis = hei->source();
		  Vertex_iterator   vit = hei->target();
          if (vis->degree() != 2 && vit->degree() != 2)
		    continue;
          Xcurve & xcurve = hei->curve();
		  Coord_type dist = m_tab_traits.xcurve_point_distance( p, xcurve , this);
          if (is_first || dist < min_dist)
          {		   
            min_dist = dist;
            closest_curve = hei;
            is_first = false;
          }    
        }
		if (is_first) // we didn't find any "good" curve
		{
		  first_time_merge = true;
		  return;
		}
	    m_tab_traits.draw_xcurve(this , closest_curve->curve() );
		seconed_curve = m_curves_arr.halfedges_end();
	  }
	  else
	  {
		first_time_merge = true;
	    m_tab_traits.find_close_curve(closest_curve ,seconed_curve ,p ,this ,false);
		redraw();	
	  }	  
	}
  }


};
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/*! class Segment_tab_traits defines the segment traits */
class Segment_tab_traits 
{
public:
  typedef Pm_seg_list Curves_list;
  typedef Seg_arr Curves_arr;
  typedef Seg_traits Traits;
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
  typedef Curves_arr::Vertex_iterator Vertex_iterator;
  typedef Curves_arr::Halfedge_around_vertex_circulator
  Halfedge_around_vertex_circulator;
  typedef Curves_arr::Edge_iterator Edge_iterator;
  typedef Curve_data Data;
  
  /*! constructor */
  Segment_tab_traits()
  {}
  
  /*! distructor */
  ~Segment_tab_traits()
  {}
  
  /*! draw_xcurve - use Qt_Widget operator to draw 
   *\ param w - the demo widget
   *\ c - curve to be drawen
   */
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
  void middle_point( Coord_point p , 
	  Qt_widget_demo_tab<Segment_tab_traits> * w)
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
  void last_point( Coord_point p , 
	  Qt_widget_demo_tab<Segment_tab_traits> * w )
  {
    return;
  }
  
  /*! get_segment - create a new segment, insert him into curves_list
      and planar map
   */
  void get_segment( Coord_segment coord_seg ,
                    Qt_widget_demo_tab<Segment_tab_traits> * w)
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
    w->m_curves_arr.insert(*seg);
  }
  
  /*! curve_point_distance - return the distance between a point 
      and a segment
	  */
  Coord_type curve_point_distance(Coord_point p, Curve * c ,
                                  Qt_widget_demo_tab<Segment_tab_traits> * w)
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
  
  /*! xcurve_point_distance - return the distance between a point
      and a xsegment */
  Coord_type xcurve_point_distance(Coord_point p, Xcurve & c ,
                                   Qt_widget_demo_tab<Segment_tab_traits> * w)
  {
    const Pm_seg_point_2 & source = c.source();
    const Pm_seg_point_2 & target = c.target();
    
    Coord_type x1 = CGAL::to_double(source.x());
    Coord_type y1 = CGAL::to_double(source.y());
    
    Coord_type x2 = CGAL::to_double(target.x());
    Coord_type y2 = CGAL::to_double(target.y());
    
    Coord_point coord_source(x1 , y1);
    Coord_point coord_target(x2 , y2);
    Coord_segment coord_seg(coord_source, coord_target);
    return CGAL::squared_distance( p, coord_seg);
  }
  
  /*! get_origin_curve - return the origin base curve */
  const Pm_base_seg_2* get_origin_curve(const Pm_xseg_2 & xseg )
  {
    const Curve_data & curve_data = xseg.get_data();
    if (curve_data.m_type == Curve_data::LEAF) 
      return curve_data.m_ptr.m_curve;
    else // curve_data.m_type == Curve_data::INTERNAL
    {
      const Pm_xseg_2* xseg_p =
        static_cast<const Pm_xseg_2*>(curve_data.m_ptr.m_x_motonote_curve);
      return get_origin_curve( *xseg_p );
    }
  }
  
  /*! draw_last_segment - call from mouse move event */
  void draw_last_segment( Qt_widget_demo_tab<Segment_tab_traits> * w)
  {
    *w << Coord_segment( m_p1 , m_p2 );
  }
  
  /*! draw_current_segment - call from mouse move event */
  void draw_current_segment( Coord_point p ,
                             Qt_widget_demo_tab<Segment_tab_traits> * w)
  {
    *w << Coord_segment( m_p1 , p);
    m_p2 = p;
  }
  
  /*! closest_point - find the closest point in the planar map
      to a clicked point
   */
  Coord_type closest_point(Coord_type x, Coord_type y, Coord_point &closest,
                           Qt_widget_demo_tab<Segment_tab_traits> * w)
  {
    bool first = true;
    Coord_type x1,y1,dt,min_dist = 0;
    Vertex_iterator vit;
    for (vit = w->m_curves_arr.vertices_begin();
         vit != w->m_curves_arr.vertices_end(); vit++)
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

  /*! curve_make_x_monotone */
  const Xcurve curve_make_x_monotone(Pm_point_2 p1 , Pm_point_2 p2)
  {
    Data d;
    const Curve cv(Base_curve(p1 , p2) , d);
	std::list<Xcurve> xcurve_list;
	m_traits.curve_make_x_monotone(cv, std::back_inserter(xcurve_list));
    const Xcurve c1 = xcurve_list.front();
	return c1;
  }

  void find_close_curve(Halfedge_iterator &closest_curve ,Halfedge_iterator 
	&seconed_curve ,Coord_point &p ,Qt_widget_demo_tab<Segment_tab_traits> *w,
	bool move_event)
  {
    bool is_first = true;
    Coord_type min_dist = 0;
	  
	Halfedge_iterator hei;
	Pm_point_2 s = closest_curve->source()->point();
	Pm_point_2 t = closest_curve->target()->point();
	Vertex_iterator   vis = closest_curve->source();
	Vertex_iterator   vit = closest_curve->target();
    for (hei = w->m_curves_arr.halfedges_begin();
          hei != w->m_curves_arr.halfedges_end(); ++hei) 
    {	
	  Pm_point_2 s1 = hei->source()->point();
	  Pm_point_2 t1 = hei->target()->point();
	  if ((s1 == s && t1 == t) || (s1 == t && t1 == s))
		continue;      // same curve as closest_curve
	  if (((s == s1 || s == t1) && vis->degree() != 2) ||
		  ((t == s1 || t == t1) && vit->degree() != 2))
		continue;      // vertex degree > 2
	  if (((s1 == s || t1 == t) &&   // source or target points are connected
	      ((s.y()-t.y())/(s.x()-t.x()) == 
		  (t1.y()-s1.y())/(t1.x()-s1.x()))) || // same "shipua"
		  ((s1 == t || t1 == s) &&   // points are connected source to target
	      ((s.y()-t.y())/(s.x()-t.x()) == 
		  (s1.y()-t1.y())/(s1.x()-t1.x())))) // same "shipua"
	  { 		  
        Xcurve & xcurve = hei->curve();
        Coord_type dist = xcurve_point_distance( p, xcurve , w);
        if (is_first || dist < min_dist)
        {
	      min_dist = dist;
          seconed_curve = hei;
          is_first = false;
        }    
      }
	}
	if (is_first) // didn't find any "good" curve
	{
	  return;
	}
	else if (!move_event)	  
	{	  
	  Pm_point_2 s1 = seconed_curve->source()->point();
	  Pm_point_2 t1 = seconed_curve->target()->point();
	  Base_curve *base;
	  
	  if ( t == s1 )
	    base = new Base_curve(s, t1);
	  else if ( t == t1 )
	    base = new Base_curve(s, s1);
	  else if ( s == s1 )
	    base = new Base_curve(t, t1);
	  else if ( s == t1 )
	  	base = new Base_curve(t, s1);

      Curve_data cd;
      cd.m_type = Curve_data::LEAF;
      cd.m_index = w->index;
      cd.m_ptr.m_curve = base;
      Curve *seg = new Curve( *base, cd );
      std::list<Xcurve> xcurve_list;
	  m_traits.curve_make_x_monotone(*seg, std::back_inserter(xcurve_list));
      const Xcurve c = xcurve_list.front();
	  w->m_curves_arr.merge_edge( closest_curve , seconed_curve , c); 
	}
  }

  void split_edge(Halfedge_iterator &hei ,Pm_point_2 &p ,
	  Qt_widget_demo_tab<Segment_tab_traits> *w)
  {
    std::list<Xcurve> xcurve_list;
    Pm_point_2 s = hei->source()->point();
	Pm_point_2 t = hei->target()->point();

	Base_curve *base1 = new Base_curve(s, p);
	Curve_data cd1;
    cd1.m_type = Curve_data::LEAF;
    cd1.m_index = hei->curve().get_data().m_index;
    cd1.m_ptr.m_curve = base1;
    Curve *seg1 = new Curve( *base1, cd1 );
    m_traits.curve_make_x_monotone(*seg1, std::back_inserter(xcurve_list));
    const Xcurve c1 = xcurve_list.front();
	
	xcurve_list.clear();

	Base_curve *base2 = new Base_curve(p, t);
	Curve_data cd2;
    cd2.m_type = Curve_data::LEAF;
    cd2.m_index = hei->curve().get_data().m_index;
    cd2.m_ptr.m_curve = base2;
    Curve *seg2 = new Curve( *base2, cd2 );
    m_traits.curve_make_x_monotone(*seg2, std::back_inserter(xcurve_list));
    const Xcurve c2 = xcurve_list.front();

	w->m_curves_arr.split_edge(hei , c1 , c2);
  }

  /*! temporary points of the created segment */
  Traits m_traits;
  Coord_point m_p1,m_p2;
};

//////////////////////////////////////////////////////////////////////////////

/*! class Polyline_tab_traits defines the polyline traits */
class Polyline_tab_traits  
{
public:
  typedef Pm_pol_list Curves_list;
  typedef Pol_arr Curves_arr;
  typedef Pol_traits Traits;
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
  typedef Curves_arr::Vertex_iterator Vertex_iterator;
  typedef Curves_arr::Halfedge_around_vertex_circulator
  Halfedge_around_vertex_circulator;
  typedef Curves_arr::Edge_iterator Edge_iterator;
  typedef Curve_pol_data Data;
  
  /*! constructor */
  Polyline_tab_traits()
  {}
  
  /*! distructor */
  ~Polyline_tab_traits()
  {}
  
  /*! draw_xcurve - go over the polyline parts and use Qt_Widget operator
      to draw
   */
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
  void middle_point( Coord_point p ,
	  Qt_widget_demo_tab<Polyline_tab_traits> * w)
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
  
  /*! last_point - last point of the polyline, create new 
      polyline and reset 
   */
  void last_point( Coord_point p ,Qt_widget_demo_tab<Polyline_tab_traits> * w)
  {
    get_polyline(w);
    points.clear();
    w->active = false;
    w->first_time = true;
    //w->redraw();  // not working so I use new_object insted
    w->new_object(make_object(Coord_segment(p , p)));
  }
  
  /*! curve_point_distance - return the distance between a point
      and a polyline
   */
  Coord_type curve_point_distance(Coord_point p, Curve * c ,
                                  Qt_widget_demo_tab<Polyline_tab_traits> * w)
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
  
  /*! xcurve_point_distance - return the distance between a point
      and a polyline
   */
  Coord_type xcurve_point_distance(Coord_point p, Xcurve & c ,
                                  Qt_widget_demo_tab<Polyline_tab_traits> * w)
  {
    Curve_const_iterator ps = c.begin();
    Curve_const_iterator pt = ps; pt++;
    bool first = true;
    Coord_type min_dist = 0;
    while (pt != c.end()) {
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
  void draw_current_segment( Coord_point p ,
                             Qt_widget_demo_tab<Polyline_tab_traits> * w)
  {
    *w << Coord_segment(p, last_of_poly);
    rubber = p;
    rubber_old = p;
  }
  
  /*! closest_point - find the closest point in the planar map
      to a clicked point
  */
  Coord_type closest_point(Coord_type x, Coord_type y, Coord_point &closest,
                           Qt_widget_demo_tab<Polyline_tab_traits> * w)
  {
    bool first = true;
    Coord_type x1,y1,dt,min_dist = 0;
    Halfedge_iterator heit;
    for (heit = w->m_curves_arr.halfedges_begin();
         heit != w->m_curves_arr.halfedges_end(); heit++)
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
  
  /*! get_origin_curve - return the origin base curve */
  const Pm_base_pol_2* get_origin_curve(const Pm_xpol_2 & xseg )
  {
    const Curve_pol_data & curve_data = xseg.get_data();
    if (curve_data.m_type == Curve_pol_data::LEAF) 
      return curve_data.m_ptr.m_curve;
    else // curve_data.m_type == Curve_data::INTERNAL
    {
      const Pm_xpol_2* xseg_p =
        static_cast<const Pm_xpol_2*>(curve_data.m_ptr.m_x_motonote_curve);
      return get_origin_curve( *xseg_p );
    }
  }

   /*! curve_make_x_monotone */
  const Xcurve curve_make_x_monotone(Pm_point_2 p1 , Pm_point_2 p2)
  {
	std::vector<Pm_pol_point_2> temp_points;
	temp_points.push_back(p1);
	temp_points.push_back(p2);
    Data d;
    Curve cv(Base_curve(temp_points.begin(), temp_points.end() ) , d);
	std::list<Xcurve> xcurve_list;
	m_traits.curve_make_x_monotone(cv, std::back_inserter(xcurve_list));
    const Xcurve c1 = xcurve_list.front();
	return c1;
  }
 
  void find_close_curve(Halfedge_iterator &closest_curve ,Halfedge_iterator 
	&seconed_curve ,Coord_point &p ,Qt_widget_demo_tab<Polyline_tab_traits> *w
	,bool move_event)
  {
    bool is_first = true;
    Coord_type min_dist = 0;
	  
	Halfedge_iterator hei;
	Pm_point_2 s;
	Pm_point_2 t;
    Vertex_iterator   vis;
    Vertex_iterator   vit;

	if (closest_curve->source()->point() == *(closest_curve->curve().begin()))
	{
	  vis = closest_curve->source();
      vit = closest_curve->target();
	}
	else
	{
	  vit = closest_curve->source();
      vis = closest_curve->target();
	}
	s = *(closest_curve->curve().begin());
	t = *(closest_curve->curve().rbegin());

    for (hei = w->m_curves_arr.halfedges_begin();
          hei != w->m_curves_arr.halfedges_end(); ++hei) 
    {
	  Pm_point_2 s1 = *(hei->curve().begin());
	  Pm_point_2 t1 = *(hei->curve().rbegin());
	  
	  if (m_traits.curve_equal(closest_curve->curve(), hei->curve()))
	    continue;  // same curve as closest_curve
	  if (!(s1 == s || t1 == t || s1 == t || t1 == s))
	    continue; // no shared points
	  if (((s == s1 || s == t1) && vis->degree() != 2) ||
		  ((t == s1 || t == t1) && vit->degree() != 2))
	  {
		continue;      // vertex degree > 2
	  }
      if ((s == s1 && ((t.x() < s.x() && t1.x() > s.x()) ||
		               (t.x() > s.x() && t1.x() < s.x()))) ||
	      (s == t1 && ((t.x() < s.x() && s1.x() > s.x()) ||
		               (t.x() > s.x() && s1.x() < s.x()))) ||
	      (t == s1 && ((t.x() < s.x() && t1.x() < t.x()) ||
		               (t.x() > s.x() && t1.x() > t.x()))) ||
	      (t == t1 && ((t.x() < s.x() && s1.x() < t.x()) ||
		               (t.x() > s.x() && s1.x() > t.x()))))
	  {
	    Xcurve & xcurve = hei->curve();
        Coord_type dist = xcurve_point_distance( p, xcurve , w);
        if (is_first || dist < min_dist)
        {
	      min_dist = dist;
          seconed_curve = hei;
          is_first = false;
        }    
      }
	}
	if (is_first) // didn't find any "good" curve
	{
	std::cout << "can't find a curve to merge" << std::endl;
      std::fflush(stdout);
	  return;
	}
	else if (!move_event)
	{
	  Xcurve & c = closest_curve->curve();
	  Xcurve & c1 = seconed_curve->curve();

	  Pm_point_2 s1 = *(seconed_curve->curve().begin());
	  Pm_point_2 t1 = *(seconed_curve->curve().rbegin());
	  std::vector<Pm_pol_point_2> temp_points;
	  Curve_const_iterator cit;
	 
	  if (t == s1 || t == t1)
	  {
        for (cit = c.begin(); cit != c.end(); cit++)
	      temp_points.push_back(*cit);
	  }
	  else
	  {
		for (cit = c.rbegin(); cit != c.rend(); cit++)
		  temp_points.push_back(*cit);
		  
	  }

	  if (s1 == s || s1 == t)
	  {
		cit = c1.begin(), cit++;
		for (; cit != c1.end(); cit++)
		  temp_points.push_back(*cit);		  
	  }
	  else
	  {
		cit = c1.rbegin(), cit++;
		for (; cit != c1.rend(); cit++)
		  temp_points.push_back(*cit);		
	  }

	  Base_curve *base = new Base_curve(temp_points.begin(), temp_points.end());
      Curve_pol_data cd;
      cd.m_type = Curve_pol_data::LEAF;
      cd.m_index = w->index;
      cd.m_ptr.m_curve = base;
      Curve *curve = new Curve( *base, cd );
	  std::list<Xcurve> xcurve_list;
	  m_traits.curve_make_x_monotone(*curve, std::back_inserter(xcurve_list));
      const Xcurve cc = xcurve_list.front();
	  w->m_curves_arr.merge_edge( closest_curve , seconed_curve , cc); 
	}
  }

  void split_edge(Halfedge_iterator &hei ,Pm_point_2 &p ,
	  Qt_widget_demo_tab<Polyline_tab_traits> *w)
  {
    Xcurve & c = hei->curve();
    Pm_point_2 s = *(c.begin());
	Pm_point_2 t = *(c.rbegin());
	std::list<Xcurve> xcurve_list;
	std::vector<Pm_pol_point_2> temp_points1;
	std::vector<Pm_pol_point_2> temp_points2;
	Curve_const_iterator cit;
	
	temp_points2.push_back(p);
	if (s.x() < p.x())
	{
	  for (cit = c.begin(); (*cit).x() < p.x(); cit++)
        temp_points1.push_back(*cit);
      for (; cit != c.end(); cit++)
        temp_points2.push_back(*cit);
	}
	else
	{
	  for (cit = c.begin(); (*cit).x() > p.x(); cit++)
        temp_points1.push_back(*cit);
	   for (; cit != c.end(); cit++)
          temp_points2.push_back(*cit);
	}
   
	temp_points1.push_back(p);

	Base_curve *base1 = new Base_curve(temp_points1.begin(), temp_points1.end());
	Curve_pol_data cd1;
    cd1.m_type = Curve_pol_data::LEAF;
    cd1.m_index = hei->curve().get_data().m_index;
    cd1.m_ptr.m_curve = base1;
    Curve *seg1 = new Curve( *base1, cd1 );
    m_traits.curve_make_x_monotone(*seg1, std::back_inserter(xcurve_list));
    const Xcurve c1 = xcurve_list.front();
	
	xcurve_list.clear();

	Base_curve *base2 = new Base_curve(temp_points2.begin(), temp_points2.end());
	Curve_pol_data cd2;
    cd2.m_type = Curve_pol_data::LEAF;
    cd2.m_index = hei->curve().get_data().m_index;
    cd2.m_ptr.m_curve = base2;
    Curve *seg2 = new Curve( *base2, cd2 );
    m_traits.curve_make_x_monotone(*seg2, std::back_inserter(xcurve_list));
    const Xcurve c2 = xcurve_list.front();
	w->m_curves_arr.split_edge(hei , c1 , c2);

	points.clear();
  }

private:
  
  /*! get_polyline - create a new polyline */
  void get_polyline(Qt_widget_demo_tab<Polyline_tab_traits> * w)
  {
    Pm_base_pol_2 * base_pol_p =
      new Pm_base_pol_2(points.begin(), points.end());
    Curve_pol_data cd;
    cd.m_type = Curve_pol_data::LEAF;
    cd.m_index = w->index;
    cd.m_ptr.m_curve = base_pol_p;
    Pm_pol_2 * pol = new Pm_pol_2( *base_pol_p, cd );
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
  Coord_segment convert(const Pm_pol_point_2 & source ,
                        const Pm_pol_point_2 & target)
  {
    Coord_type x1 = CGAL::to_double(source.x());
    Coord_type y1 = CGAL::to_double(source.y());
    Coord_point coord_source(x1, y1);
    
    Coord_type x2 = CGAL::to_double(target.x());
    Coord_type y2 = CGAL::to_double(target.y());
    Coord_point coord_target(x2, y2);
    
    return Coord_segment(coord_source, coord_target);
  }
  
  Traits m_traits;
  /*! the new point of the rubber band */
  Coord_point rubber;
  /*! the last point of the polygon */
  Coord_point last_of_poly;
  /*! the old point of the rubber band */
  Coord_point rubber_old;
  /*! container to hold the point during polyline creation */
  std::vector<Pm_pol_point_2> points;
  
};

//////////////////////////////////////////////////////////////////////////////
/*! */
class Conic_tab_traits
{
public:
  typedef Pm_xconic_list Curves_list;
  typedef Conic_arr Curves_arr;
  typedef Conic_traits Traits;
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
  typedef Curves_arr::Vertex_iterator Vertex_iterator;
  typedef Curves_arr::Halfedge_around_vertex_circulator
  Halfedge_around_vertex_circulator;
  typedef Curves_arr::Edge_iterator Edge_iterator;
  typedef Curve_conic_data Data;
  
  /*! constructor */
  Conic_tab_traits()
  {}
  
  /*! distructor */
  ~Conic_tab_traits()
  {}
  
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
            (*w) << Coord_segment( Coord_point(prev_x, prev_y) ,
                                   Coord_point(curr_x, curr_y) );
            prev_x = curr_x;
            prev_y = curr_y;
            
          }
        }
        
        (*w) << Coord_segment( Coord_point(prev_x, prev_y) ,
                               Coord_point(end_x, end_y) );
      }
      else
      {
        // We should never reach here.
        CGAL_assertion(false);
      }
    }
  }
  
  /*! get_origin_curve - return the origin base curve */
  const Pm_base_conic_2* get_origin_curve(const Pm_xconic_2 & xseg )
  {
    const Curve_conic_data & curve_data = xseg.get_data();
    if (curve_data.m_type == Curve_conic_data::LEAF) 
      return curve_data.m_ptr.m_curve;
    else // curve_data.m_type == Curve_data::INTERNAL
    {
      Pm_xconic_2* xseg_p =
        static_cast<Pm_xconic_2*>(curve_data.m_ptr.m_x_motonote_curve);
      return get_origin_curve( *xseg_p );
    }
  }
  
  ////////////////////////////////////////////////////////////////////////////
  
  /*! first_point - a first point of inserted sgment */
  void first_point( Coord_point p )
  {
    m_p1 = m_p2 = p;
	first_time_hyperbula = true;
  }
  
  /*! middle_point - the last point of a segment */
  void middle_point( Coord_point p , Qt_widget_demo_tab<Conic_tab_traits> * w)
  {
    if(m_p1.x() != p.x() || m_p1.y() != p.y()) 
    {
      CONIC_NT r, s, t, u, v, ww;  // The conic coefficients.
      CONIC_NT a, b, c, a_sq, b_sq, x, y, x1, y1, x0, y0, temp;
      x = m_p1.x();
      y = m_p1.y();
      x1 = p.x();
      y1 = p.y();
      Pm_base_conic_2 *cv;
	  Coord_type dist;
	  switch (w->conic_type)
	  {
	  case CIRCLE:
	    dist = pow(m_p1.x() - p.x(), 2) + pow(m_p1.y() - p.y(), 2);
	    cv = new Pm_base_conic_2(Pm_conic_circle_2 
                             (Pm_conic_point_2(x,y ), dist, CGAL::CLOCKWISE));
		break;
	  case SEGMENT:
	    cv = new Pm_base_conic_2(Pm_conic_segment_2(Pm_conic_point_2(x,y),
                                                    Pm_conic_point_2(x1,y1)));
		break;
	  case ELLIPSE:
		a = abs(x1 - x)/2;
		b = abs(y1 - y)/2;
        a_sq = a*a;
		b_sq = b*b;
		x0 = (x + x1)/2;
		y0 = (y + y1)/2;

		r = b_sq;
        s = a_sq;
        t = 0;
        u = -2*x0*b_sq;
        v = -2*y0*a_sq;
        ww = x0*x0*b_sq + y0*y0*a_sq - a_sq*b_sq;

		cv = new Pm_base_conic_2(r, s, t, u, v, ww);
	    break;
	  case PARABOLA:
	    if (x > x1)
		{
		  temp = x1;
		  x1 = x;
		  x = temp;
		}		  
	    x0 = (x + x1)/2;
		y0 = y;
		a = (y1 - y0)/((x1 - x0)*(x1 - x0));
		b = -2*x0*a;
		c = y0 + x0*x0*a;
        
		r = a;
        s = 0;
        t = 0;
        u = b;
        v = -1;
        ww = c;

		cv = new Pm_base_conic_2(r, s, t, u, v, ww, Pm_conic_point_2(x1,y1), 
			Pm_conic_point_2(x,y1));
	    break;
	  case HYPERBOLA:
	    if (first_time_hyperbula)
		{
		  *w << CGAL::RED;
		  *w << Coord_segment( Coord_point(m_p1.x(),(m_p1.y() + p.y())/2) ,
			                   Coord_point(p.x()   ,(m_p1.y() + p.y())/2) );   
		  first_time_hyperbula = false;
		  m_p3 = m_p2;

		  *w << CGAL::DISC;
		  *w << m_p1;
		  *w << m_p2;
		  return;
		}
		else
		{                        /*  x,y             */
		  x1 = m_p3.x();         /*     \    /       */
		  y1 = m_p3.y();         /*      \  /        */
		                         /*  x0,y0\/         */
		  x0 = (x + x1)/2;       /*       /\         */
		  y0 = (y + y1)/2;       /*      /  \        */
                                 /*     /    \       */
		  a = abs(x0 - p.x());   /*           x1,y1  */
		  b = ((y - y1)/(x - x1))*a;
		 
		  if (p.x() == x0 || (p.x() < x && x < x1 ) || (p.x() > x && x > x1 )
			           || (p.x() < x1 && x1 < x ) || (p.x() > x1 && x1 > x ))
		    return; // p is out of rectangle bounds
		  r = b*b;
		  s = -1*a*a;
		  t = 0;
		  u = -2*b*b*x0;
		  v = 2*a*a*y0;
		  ww = b*b*x0*x0 - a*a*y0*y0 - a*a*b*b;

		  CONIC_NT y2, y3 ,root;

		  root = CGAL::sqrt(v*v - 4*s*(ww + u*x1 + r*x1*x1));
		  y2 = (-1*v + root)/(2*s);
		  y3 = (-1*v - root)/(2*s);

		  if (x1 > x && p.x() > x0)
		    cv = new Pm_base_conic_2(r, s, t, u, v, ww, 
			         Pm_conic_point_2(x1,y3), Pm_conic_point_2(x1,y2));	
		  else if (x1 < x && p.x() < x0)
		    cv = new Pm_base_conic_2(r, s, t, u, v, ww, 
			         Pm_conic_point_2(x1,y2), Pm_conic_point_2(x1,y3));	
		  else if (x1 < x && p.x() > x0)
		    cv = new Pm_base_conic_2(r, s, t, u, v, ww, 
			         Pm_conic_point_2(x,y3), Pm_conic_point_2(x,y2));	
		  else if (x1 > x && p.x() < x0)
		    cv = new Pm_base_conic_2(r, s, t, u, v, ww, 
			         Pm_conic_point_2(x,y2), Pm_conic_point_2(x,y3));	
		
		}		
		break;	
	  }

      Curve_conic_data cd;
      cd.m_type = Curve_conic_data::LEAF;
      cd.m_index = w->index;
      cd.m_ptr.m_curve = cv;
      w->m_curves_arr.insert(Pm_conic_2( *cv , cd));
        
      w->active = false;
      //w->redraw();  // not working so I use new_object insted
      w->new_object(make_object(Coord_segment(m_p1 , p)));
    } 
  }
  
  /*! last_point - meaningless for conics - at least for now */
  void last_point( Coord_point p , Qt_widget_demo_tab<Conic_tab_traits> * w )
  {   
      return;
  }
  /*! draw_last_segment - call from mouse move event */
  void draw_last_segment( Qt_widget_demo_tab<Conic_tab_traits> * w)
  {
    if (w->mode == SPLIT)
	  *w << Coord_segment( m_p1 , m_p2 );
	else
	{
      switch (w->conic_type)
	  {
	    case CIRCLE:
	    *w << Coord_circle(m_p1,
		    pow(m_p1.x() - m_p2.x(), 2) + pow(m_p1.y() - m_p2.y(),2));
	    break;
	    case SEGMENT:
    	  *w << Coord_segment( m_p1 , m_p2 );
	    break;
        case ELLIPSE:
	    {
 		  *w << Coord_segment( Coord_point(m_p1.x(),m_p2.y()) , m_p2 );
		  *w << Coord_segment( Coord_point(m_p1.x(),m_p2.y()) , m_p1 );
		  *w << Coord_segment( Coord_point(m_p2.x(),m_p1.y()) , m_p2 );
		  *w << Coord_segment( Coord_point(m_p2.x(),m_p1.y()) , m_p1 );
		  break;
	    }		
	    case PARABOLA:
	    {
 	  	  *w << Coord_segment( Coord_point(m_p1.x(),m_p2.y()) , m_p2 );
		  *w << Coord_segment( Coord_point(m_p1.x(),m_p2.y()) , m_p1 );
		  *w << Coord_segment( Coord_point(m_p2.x(),m_p1.y()) , m_p2 );
		  *w << Coord_segment( Coord_point(m_p2.x(),m_p1.y()) , m_p1 );
		  break;
	    }		
	    case HYPERBOLA:
	    {
	      if (first_time_hyperbula)
	 	  {
 		    *w << Coord_segment( Coord_point(m_p1.x(),m_p2.y()) , m_p2 );
		    *w << Coord_segment( Coord_point(m_p1.x(),m_p2.y()) , m_p1 );
		    *w << Coord_segment( Coord_point(m_p2.x(),m_p1.y()) , m_p2 );
		    *w << Coord_segment( Coord_point(m_p2.x(),m_p1.y()) , m_p1 );
		    *w << Coord_segment( m_p1 , m_p2 );
		    *w << Coord_segment( Coord_point(m_p2.x(),m_p1.y()) , Coord_point(m_p1.x(),m_p2.y()) );
		  }
		  break;
	    }		
	  }
	}
  }
  
  /*! draw_current_segment - call from mouse move event */
  void draw_current_segment( Coord_point p ,
                             Qt_widget_demo_tab<Conic_tab_traits> * w)
  {
    m_p2 = p;
	draw_last_segment(w);
  }
	
  /*! closest_point - find the closest point in the planar map
      to a clicked point
  */
  Coord_type closest_point(Coord_type x, Coord_type y, Coord_point &closest,
                           Qt_widget_demo_tab<Conic_tab_traits> * w)
  {
    bool first = true;
    Coord_type x1,y1,dt,min_dist = 0;
    Vertex_iterator vit;
    for (vit = w->m_curves_arr.vertices_begin();
         vit != w->m_curves_arr.vertices_end(); vit++)
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
  
  /*! curve_point_distance - return the distance between 
      a point and a conic */
  Coord_type curve_point_distance(Coord_point p, Curve * c ,
                                  Qt_widget_demo_tab<Conic_tab_traits> * w)
  {
    // Get the co-ordinates of the curve's source and target.
    double sx = CGAL::to_double(c->source().x()),
      sy = CGAL::to_double(c->source().y()),
      tx = CGAL::to_double(c->target().x()),
      ty = CGAL::to_double(c->target().y());
    
    if (c->is_segment())
    {
      Coord_point coord_source(sx , sy);
      Coord_point coord_target(tx , ty);
      Coord_segment coord_seg(coord_source, coord_target);
      
      return CGAL::squared_distance( p, coord_seg);
    }
    else
    {
      // If the curve is monotone, than its source and its target has the
      // extreme x co-ordinates on this curve.
      if (c->is_x_monotone())
      {
        
        bool     is_source_left = (sx < tx);
        int      x_min = is_source_left ? (*w).x_pixel(sx) : 
          (*w).x_pixel(tx);
        int      x_max = is_source_left ? (*w).x_pixel(tx) :
          (*w).x_pixel(sx);
        double   prev_x = is_source_left ? sx : tx;
        double   prev_y = is_source_left ? sy : ty;
        //double   end_x = is_source_left ? tx : sx;
        //double   end_y = is_source_left ? ty : sy;
        double   curr_x, curr_y;
        int      x;
        
        Pm_conic_point_2 ps[2];
        int nps;
        
        bool first = true;
        Coord_type min_dist = 0;
        
        
        for (x = x_min + 1; x < x_max; x++)
        {
          curr_x = (*w).x_real(x);
          nps = c->get_points_at_x(Pm_conic_point_2(curr_x, 0), ps);
          if (nps == 1)
          {
            curr_y = CGAL::to_double(ps[0].y());
            
            Coord_segment coord_seg( Coord_point(prev_x, prev_y) ,
                                     Coord_point(curr_x, curr_y) );
            Coord_type dist = CGAL::squared_distance( p, coord_seg);
            if( dist < min_dist || first)
            {
              first = false;
              min_dist = dist;
            }
            
            
            prev_x = curr_x;
            prev_y = curr_y;
            
          }
        }
        return min_dist;
      }// if is_x_monotone
      else
      {
        // We should never reach here.
        CGAL_assertion(false);
        return 0;
      }
    } // else 
  } // 
  
    /*! xcurve_point_distance - return the distance between
        a point and a conic
     */
  Coord_type xcurve_point_distance(Coord_point p, Xcurve & c ,
                                   Qt_widget_demo_tab<Conic_tab_traits> * w)
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
      
      return CGAL::squared_distance( p, coord_seg);
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
        double   curr_x, curr_y;
        int      x;
        
        Pm_conic_point_2 ps[2];
        int nps;
        
        bool first = true;
        Coord_type min_dist = 100000000;
        
        
        for (x = x_min + 1; x < x_max; x++)
        {
          curr_x = (*w).x_real(x);
          nps = c.get_points_at_x(Pm_conic_point_2(curr_x, 0), ps);
          if (nps == 1)
          {
            curr_y = CGAL::to_double(ps[0].y());
            
            Coord_segment coord_seg( Coord_point(prev_x, prev_y) ,
                                     Coord_point(curr_x, curr_y) );
            Coord_type dist = CGAL::squared_distance( p, coord_seg);
            if( dist < min_dist || first)
            {
              first = false;
              min_dist = dist;
            }
            
            
            prev_x = curr_x;
            prev_y = curr_y;
            
          }
        }
        return min_dist;
      }// if is_x_monotone
      else
      {
        // We should never reach here.
        CGAL_assertion(false);
        return 0;
      }
    } // else 
  } // 
  
  /*! curve_make_x_monotone */
  const Xcurve curve_make_x_monotone(Pm_point_2 p1 , Pm_point_2 p2)
  {
    Data d;
	Pm_base_conic_2 base(Pm_conic_segment_2(p1 , p2));
    const Pm_conic_2 cv( base , d);
	std::list<Xcurve> xcurve_list;
	m_traits.curve_make_x_monotone(cv, std::back_inserter(xcurve_list));
    const Xcurve c1 = xcurve_list.front();
	return c1;
  }
  
   void find_close_curve(Halfedge_iterator &closest_curve ,Halfedge_iterator 
	&seconed_curve ,Coord_point &p ,Qt_widget_demo_tab<Conic_tab_traits> *w,
	bool move_event)
  {
    bool is_first = true;
    Coord_type min_dist = 0;
	const Xcurve first = closest_curve->curve();
	Halfedge_iterator hei;
	Pm_point_2 s = closest_curve->curve().source();
	Pm_point_2 t = closest_curve->curve().target();
	Vertex_iterator   vis;
	Vertex_iterator   vit;
	if ( closest_curve->curve().source() == closest_curve->source()->point())
	{
	  vis = closest_curve->source();
	  vit = closest_curve->target();
	}
	else
	{
	  vit = closest_curve->source();
	  vis = closest_curve->target();
	}
    for (hei = w->m_curves_arr.halfedges_begin();
          hei != w->m_curves_arr.halfedges_end(); ++hei) 
    {
	  if (first.has_same_base_conic(hei->curve()))
	  { 		  
	    Pm_point_2 s1 = hei->curve().source();
	    Pm_point_2 t1 = hei->curve().target();
	    if ((s1 == s && t1 == t) || (s1 == t && t1 == s))
		  continue;      // same curve as closest_curve
	    if (((s == s1 || s == t1) && vis->degree() != 2) ||
		    ((t == s1 || t == t1) && vit->degree() != 2))
		  continue;      // vertex degree > 2
	 	if ((s == s1 && ((t.x() < s.x() && t1.x() > s.x()) ||
		               (t.x() > s.x() && t1.x() < s.x()))) ||
	      (s == t1 && ((t.x() < s.x() && s1.x() > s.x()) ||
		               (t.x() > s.x() && s1.x() < s.x()))) ||
	      (t == s1 && ((t.x() < s.x() && t1.x() < t.x()) ||
		               (t.x() > s.x() && t1.x() > t.x()))) ||
	      (t == t1 && ((t.x() < s.x() && s1.x() < t.x()) ||
		               (t.x() > s.x() && s1.x() > t.x()))))
	    { // the connected curve will be xmonnotone	  
          Xcurve & xcurve = hei->curve();
          Coord_type dist = xcurve_point_distance( p, xcurve , w);
          if (is_first || dist < min_dist)
          {
	        min_dist = dist;
            seconed_curve = hei;
            is_first = false;
           }    
		}
      }
	}
	if (is_first) // didn't find any "good" curve
	{
	  return;
	}
	else if (!move_event)	  
	{	  
	  Xcurve curve;
	  const Xcurve seconed = seconed_curve->curve();

	  m_traits.curve_merge(curve , first , seconed );
	
	  Curve_conic_data cd1;
      cd1.m_type = Curve_conic_data::LEAF;
      cd1.m_index = first.get_data().m_index;
      cd1.m_ptr.m_curve = first.get_data().m_ptr.m_curve;
	
	  w->m_curves_arr.remove_edge(closest_curve);
	  w->m_curves_arr.remove_edge(seconed_curve);
	  w->m_curves_arr.insert(Pm_conic_2( curve , cd1));	
	}
  }

  /*!
   */
  void split_edge(Halfedge_iterator &hei ,Pm_point_2 &p ,
	  Qt_widget_demo_tab<Conic_tab_traits> *w)
  {    
    
	const Xcurve split_curve = hei->curve();
	Xcurve sbc1;
	Xcurve sbc2;

	m_traits.curve_split(split_curve , sbc1 , sbc2 , p);
	
	Curve_conic_data cd1;
    cd1.m_type = Curve_conic_data::LEAF;
    cd1.m_index = hei->curve().get_data().m_index;
    cd1.m_ptr.m_curve = hei->curve().get_data().m_ptr.m_curve;

	Pm_conic_2   sub_curve1 (sbc1, cd1);

	Curve_conic_data cd2;
    cd2.m_type = Curve_conic_data::LEAF;
    cd2.m_index = hei->curve().get_data().m_index;
    cd2.m_ptr.m_curve =hei->curve().get_data().m_ptr.m_curve;

	Pm_conic_2   sub_curve2 (sbc2, cd2);
 
	w->m_curves_arr.remove_edge(hei);
	w->m_curves_arr.insert(sub_curve1);
	w->m_curves_arr.insert(sub_curve2);
  }

  Traits m_traits;
  /*! temporary points of the created conic */
  Coord_point m_p1,m_p2,m_p3;
  /*! bool flag for hyperbula insertion */
  bool first_time_hyperbula;
};

typedef Qt_widget_demo_tab<Segment_tab_traits> Qt_widget_segment_tab;
typedef Qt_widget_demo_tab<Polyline_tab_traits> Qt_widget_polyline_tab;
typedef Qt_widget_demo_tab<Conic_tab_traits> Qt_widget_conic_tab;

