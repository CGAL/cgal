#ifndef DEMO_TAB_H
#define DEMO_TAB_H

/*! demo_tab.h contain the definetion and implementation of 
 *  the demo tab classes and the tabs traits classes.       
 *  all the possible shared code is in Qt_widget_demo_tab where
 *  the differences is in the traits classes.
 */
#include <math.h>

#include <qrect.h>
#include <qcursor.h>
#include <qmessagebox.h> 
#include <qcolor.h>
#include <qpainter.h> 
#include <qpen.h>

#include "cgal_types.h"
#include "seg_notif.h"
#include "pol_notif.h"
#include "conic_notif.h"
#include <CGAL/IO/pixmaps/hand.xpm>
#include <CGAL/IO/pixmaps/holddown.xpm>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>

#include<vector>

/*! class Qt_widget_base_tab - inherits from CGAL::Qt_widget
 *  contain all data members that are not part of the traits 
 */
class Qt_widget_base_tab : public CGAL::Qt_widget
{
public:

  /*!
   */
  Qt_widget_base_tab(TraitsType  t , Strategy _straregy, QWidget *parent = 0,
                     int tab_number = 1 ) :
    CGAL::Qt_widget( parent ),
    current_state(0),
    
    index(tab_number),
    snap_mode(NONE),
    mode(INSERT),
    m_line_width(2),
    m_vertex_width(3),
    first_time(true),
    active(false),
    traits_type(t),
    bbox(CGAL::Bbox_2(-10, -10, 10, 10)),
    wasrepainted(true), 
    on_first(false),
    change_pm_color(false),
    snap(false),
    grid(false),
    conic_type(SEGMENT),
    cube_size(1),
    ray_shooting_direction(true),
    remove_org_curve(true),
    read_from_file(false),
    empty(true),
    first_time_merge(true),
    draw_vertex(true),
    fill_face_color(def_bg_color),
    strategy(_straregy)
  {
    static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(2) <<
      CGAL::BackgroundColor (CGAL::BLACK);
    set_window(0, 700, 0, 700);
   
    setMouseTracking(TRUE);
    
    colors[1] = Qt::blue;
    colors[2] = Qt::gray;
    colors[3] = Qt::green;
    colors[4] = Qt::cyan;
    colors[5] = Qt::magenta;
    colors[6] = Qt::darkRed;
    colors[7] = Qt::darkGreen;
    colors[8] = Qt::darkBlue;
    colors[9] = Qt::darkMagenta;
    colors[10] = Qt::darkCyan;
    colors[11] = Qt::yellow;

    pm_color = colors[index];
  }

  /*! Destructor */
  virtual ~Qt_widget_base_tab(){}
  
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

  /*! m_vertex_width - vertex width */
  int m_vertex_width;

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

  /*! planar map color */
  QColor pm_color;

  /*! flag to know that pm color has changed */
  bool change_pm_color;

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

  /*! pm read from file - need a spcial treatment */
  bool read_from_file;

  /*! true if pm is empty */
  bool empty;

  /*! true when it is the first time in merge mode */
  bool first_time_merge;

  /*! true when you want to draw all vertex, false if
   * only the intersecting vertex
   */
  bool draw_vertex;

  /*! the color for filling faces ( obtained by fill button) */
  QColor fill_face_color;

  /*! point location strategy */
  Strategy strategy;

  /*! get the color of the unbounded face (its the same as the background
   * color of the tab)
   */
  QColor unbounded_face_color() { return this->backgroundColor(); } 

  /*! set the colo of the unbounded face (its the same as the background
   * color of the tab)
   */
  void set_unbounded_face_color(QColor c) { this->setBackgroundColor(c); } 

  /*! increment current_state to inidicate that something has changed
   */
  void something_changed(){ current_state++ ; }
};


/*! check if  two segemnts are mergeable
 */
template<class Pm_point_2>
bool merge_segments_precondition(Pm_point_2 s, Pm_point_2 t, Pm_point_2 s1,
                                 Pm_point_2 t1 ,int s_degree , int t_degree)
{
  return
    //  cannot be the same curve as closest_curve
    (!((s1 == s && t1 == t) || (s1 == t && t1 == s)) && 
     // must be one shared point
     (s1 == s || t1 == t || s1 == t || t1 == s) && 
     (!(((s == s1 || s == t1) && s_degree != 2) ||
        // vertex degree must be exactly 2
        ((t == s1 || t == t1) && t_degree != 2))) && 
     ((s == s1 && ((t.x() < s.x() && t1.x() > s.x()) ||
                   (t.x() > s.x() && t1.x() < s.x()))) ||
      (s == t1 && ((t.x() < s.x() && s1.x() > s.x()) ||
                   (t.x() > s.x() && s1.x() < s.x()))) ||
      (t == s1 && ((t.x() < s.x() && t1.x() < t.x()) ||
                   (t.x() > s.x() && t1.x() > t.x()))) ||
      (t == t1 && ((t.x() < s.x() && s1.x() < t.x()) ||
                   // checking x-monotone
                   (t.x() > s.x() && s1.x() > t.x()))) || 
      // dealing with vertical curves
      ((s.x() == s1.x())&&( t.x() == t1.x()) && (s.x()==t.x()))));  
}

/*!
 */
template<class ColorsIterator>
QColor blend_colors(ColorsIterator begin_colors , ColorsIterator  end_colors)
{
  ColorsIterator clr_itr;
  int red = 0,
    green = 0,
    blue = 0,
    counter = 0; 

  QColor avg_color = def_bg_color; // the average color

  for(clr_itr = begin_colors ; clr_itr!= end_colors ; ++clr_itr)
  {
    QColor curr_color = *clr_itr;
    if( curr_color.isValid() && curr_color != def_bg_color)
    {
      red += curr_color.red();
      green += curr_color.green();
      blue += curr_color.blue();
      counter++;
    }
  }

  if( counter)  // make sure counter is not equal zero
    avg_color = QColor(red / counter , green / counter , blue / counter);

  return avg_color ; 
}


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
  typedef typename Tab_traits::Halfedge      Halfedge;
  typedef typename Tab_traits::Face_iterator  Face_iterator;
  typedef typename Tab_traits::Trap_point_location Trap_point_location;
  typedef typename Tab_traits::Naive_point_location Naive_point_location;
  typedef typename Tab_traits::Simple_point_location Simple_point_location;
  typedef typename Tab_traits::Walk_point_location Walk_point_location;

private:
  // function object - FillFace
  class FillFace
  {
    /*! */
    Qt_widget_demo_tab<Tab_traits> *ptr;

  public:
    /*! constructor */
    FillFace(Qt_widget_demo_tab<Tab_traits>* tab) : ptr(tab){}

    /*!
     */
    void
    operator()(typename Qt_widget_demo_tab<Tab_traits>::Face_handle& face) 
    {
      ptr->m_tab_traits.fill_face(ptr,face);
    }
  };
  
  /*!
   */
  class OverLay
  {
    Qt_widget_demo_tab<Tab_traits> *ptr;
    std::vector<QColor> overlay_colors ;

  public:
    /*! constructor */
    OverLay(Qt_widget_demo_tab<Tab_traits> * tab,
            std::vector<QColor> colors) : ptr(tab),overlay_colors(colors){}

    /*!
     */
    void
    operator()(typename Qt_widget_demo_tab<Tab_traits>::Face_handle& face) 
    {
      ptr->overlay(face,overlay_colors);
    }
  };

public:
  /*! m_tab_traits - the traits object */
  Tab_traits m_tab_traits;

  /*! m_curves_arr - pointer for the tab planar map */
  Curves_arr *m_curves_arr;

  /*! prev_curves_arr - for undo operation */
  Curves_arr m_prev_curves_arr;

  /*! Original Traits */
  Traits m_traits;

  /*! the curve to be merged */
  Halfedge_iterator closest_curve;

  /*! the second curve to be merged */
  Halfedge_iterator second_curve;

  /*! the first point in the split curve */
  Pm_point_2 split_point;

  /*! list of removable halfedge iterators (created by move event when
   * DELETE option is on
   */
  Hafledge_list removable_halfedges;

  /*! pointer to the point location strategy */
  void * m_strategy;

  
  /*! constructor 
   *\ param t - widget traits type
   *\ param parent - widget parent window
   *\ param tab_number - widget program index
   */
  Qt_widget_demo_tab(TraitsType t, Strategy _strategy, QWidget * parent = 0,
                     int tab_number = 1):
    Qt_widget_base_tab(t ,_strategy, parent, tab_number)
  {
    switch(strategy)
    {
     case NAIVE:
      {
       m_strategy = new Naive_point_location();
       m_curves_arr = new Curves_arr((Naive_point_location*)m_strategy);
       break;
      }

     case TRAP:
      {
       m_strategy = new Trap_point_location();
       m_curves_arr = new Curves_arr((Trap_point_location*)m_strategy);
       break;
      }

     case SIMPLE:
      {
       m_strategy = new Simple_point_location();
       m_curves_arr = new Curves_arr((Simple_point_location*)m_strategy);
       break;
      }

     case WALK:
      {
       m_strategy = new Walk_point_location();
       m_curves_arr = new Curves_arr((Walk_point_location*)m_strategy);
       break;
      }
    }

    // set the unbounded face initial color
    m_curves_arr->unbounded_face()->set_info(def_bg_color);   
  }
  
  /*! destructor - delete the planar map and the point location pointer
   */
  virtual ~Qt_widget_demo_tab() 
  {
    delete m_curves_arr;
    // delete the point location strategy object
     switch(strategy)
    {
     case NAIVE:
      delete (Naive_point_location*)m_strategy;
      break;

     case TRAP:
      delete (Trap_point_location*)m_strategy;
      break;

     case SIMPLE:
      delete (Simple_point_location*)m_strategy;
      break;

     case WALK:
      delete (Walk_point_location*)m_strategy;
      break;
    }
  }
  
  ///*! save previous pm */
  //void save_prev_pm()
  //{
  //  m_prev_curves_arr = m_curves_arr;
  //}


  /*! draw - called everytime something changed, draw the PM and mark the 
   *         point location if the mode is on.
   */
  void draw()
  {
    QCursor old = cursor();
    setCursor(Qt::WaitCursor);

    if ( mode == FILLFACE ) 
    {
      Locate_type lt;
      Face_handle f;
      Pm_point_2 temp_p (pl_point.x(), pl_point.y());
      Halfedge_handle e = m_curves_arr->locate(temp_p, lt);
      if(lt == Curves_arr::UNBOUNDED_FACE )
        f = m_curves_arr->unbounded_face();
      else
        f = e->face();
      
      set_face_color(f,fill_face_color);
    }

    // draw all faces (fill them with their color)
    visit_faces(FillFace(this));  

    if (snap_mode == GRID || grid)
      draw_grid();

    for (Edge_iterator ei = m_curves_arr->edges_begin(); 
         ei != m_curves_arr->edges_end(); ++ei) 
    {
      if (change_pm_color)
        setColor(pm_color);
      else
      {
        int i = (ei->curve()).get_data().m_index;
        setColor(colors[i]);
      }
      m_tab_traits.draw_xcurve(this , ei->curve() );
    }
    
    // Go over all vertices and for each vertex check the 
    // index numbers of the base curves that go through 
    // it and paint the point if they are different (beacuse ew want to
    // color red the intersection opints between two different planar maps
    // which are overlayed
    *this << CGAL::DISC;
    static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(m_vertex_width);

    Vertex_iterator   vit;
    for (vit = m_curves_arr->vertices_begin(); 
         vit != m_curves_arr->vertices_end(); vit++)
    {
      Halfedge_around_vertex_circulator 
        eit,eit2, first = (*vit).incident_halfedges();
      
      eit = first;
      setColor(Qt::red);
      int ind1, ind2;
      
      eit2 = first;
      do
      { // find (if exist) a different index than the tab index
        ind2 = (*eit2).curve().get_data().m_index;
        if (ind2 != index)
          break;
        eit2++;
      } while (eit2 != first);

      do 
      {
        ind1 = (*eit).curve().get_data().m_index;
        
        // Keep track of IDs we haven't seen before.
        if (ind1 != ind2 && ind1 != index && ind2 != index)
        {
          //const Pm_point_2& p = (*vit).point();
          Coord_point p(CGAL::to_double((*vit).point().x()) /
                        m_tab_traits.COORD_SCALE, 
                        CGAL::to_double((*vit).point().y()) /
                        m_tab_traits.COORD_SCALE); 
          static_cast<CGAL::Qt_widget&>(*this) << p;
          break;
        }
        
        eit++;
        
        // draw all vertexes of the planar map is 'draw_vertex' is true
        // draw_vertex is a flag that indicates if we draw the vertexes
        if (eit == first && draw_vertex)
        {
          if (change_pm_color)
            setColor(pm_color);
          else
            setColor(colors[ind2]);
          //const Pm_point_2& p = (*vit).point();
          Coord_point p(CGAL::to_double((*vit).point().x()) /
                        m_tab_traits.COORD_SCALE, 
                        CGAL::to_double((*vit).point().y()) /
                        m_tab_traits.COORD_SCALE); 
          static_cast<CGAL::Qt_widget&>(*this) << p;
        }
      } while (eit != first);
    }
    
    if (mode == POINT_LOCATION && 
        ! (m_curves_arr->halfedges_begin() == m_curves_arr->halfedges_end()))
    {
      static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(3);
      setColor(Qt::yellow);
      
      Locate_type lt;
      Pm_point_2 temp_p (pl_point.x(), pl_point.y());
      Halfedge_handle e = m_curves_arr->locate(temp_p, lt);
      
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
      static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(m_line_width);
    }
    
    if (mode == RAY_SHOOTING && 
        !(m_curves_arr->halfedges_begin() == m_curves_arr->halfedges_end()))
    {
      Coord_point up;
      Locate_type lt;
      Pm_point_2 temp_p (pl_point.x(), pl_point.y());
      Halfedge_handle e = 
        m_curves_arr->vertical_ray_shoot(temp_p, lt ,ray_shooting_direction);
      
      setColor(Qt::black);
      static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(1);

      Coord_point pl_draw(pl_point.x() / m_tab_traits.COORD_SCALE , 
                          pl_point.y() / m_tab_traits.COORD_SCALE);

      if (ray_shooting_direction)
      {
        if (lt == Curves_arr::UNBOUNDED_FACE)
        {
          up = Coord_point(pl_draw.x() , y_max());
          static_cast<CGAL::Qt_widget&>(*this) << Coord_segment(pl_draw, up);
        }
        else // we shoot something
        {
          Pm_point_2 p1c1(pl_point.x() , y_max() * m_tab_traits.COORD_SCALE);
          Pm_point_2 p2c1(pl_point.x() , pl_point.y());
          const Xcurve c1 = m_tab_traits.curve_make_x_monotone(p1c1 , p2c1);
          const Xcurve c2 = e->curve();
          Pm_point_2 most_left(x_min() * m_tab_traits.COORD_SCALE,
                               pl_point.y());
          CGAL::Object res =
            m_traits.nearest_intersection_to_right(c1, c2, most_left);
          Pm_point_2 p1;
          if (CGAL::assign(p1, res)) {
            Coord_type y1 =
              CGAL::to_double(p1.y()) / m_tab_traits.COORD_SCALE;
            up = Coord_point(pl_draw.x(), y1);
            static_cast<CGAL::Qt_widget&>(*this) << Coord_segment(pl_draw,
                                                                  up);
          }
          //! \todo what if the intersection is empty, or a subcurve?
        }
  
        // draw an arrow that points to 'up' point
        int x = this->x_pixel(CGAL::to_double(up.x()));
        int y = this->y_pixel(CGAL::to_double(up.y()));

        this->get_painter().drawLine(x-7 , y+7 , x , y);
        this->get_painter().drawLine(x+7 , y+7 , x , y);
      }
      else // down ray shooting
      {
        Coord_point down;
        if (lt == Curves_arr::UNBOUNDED_FACE)
        {
          down = Coord_point(pl_draw.x() , y_min());
          static_cast<CGAL::Qt_widget&>(*this) << Coord_segment(pl_draw ,
                                                                down);
        }
        else // we shoot something
        {
          Pm_point_2 p1c1(pl_point.x() , y_min() * m_tab_traits.COORD_SCALE);
          Pm_point_2 p2c1(pl_point.x() , pl_point.y());
          const Xcurve c1 = m_tab_traits.curve_make_x_monotone(p1c1 , p2c1);
          const Xcurve c2 = e->curve();
          Pm_point_2 most_right(x_max() * m_tab_traits.COORD_SCALE,
                                pl_point.y());

          CGAL::Object res =
            m_traits.nearest_intersection_to_left(c1, c2, most_right);
          Pm_point_2 p1;
          if (CGAL::assign(p1, res)) {
            Coord_type y1 =
              CGAL::to_double(p1.y()) / m_tab_traits.COORD_SCALE;
            down = Coord_point(pl_draw.x(),y1);
            static_cast<CGAL::Qt_widget&>(*this) << Coord_segment(pl_draw,
                                                                  down);
          }
          //! \todo what if the intersection is empty, or a subcurve?
        }
      
        // draw an arrow that points to 'down' point
        int x = this->x_pixel(CGAL::to_double(down.x()));
        int y = this->y_pixel(CGAL::to_double(down.y()));

        this->get_painter().drawLine(x-7 , y-7 , x , y);
        this->get_painter().drawLine(x+7 , y-7 , x , y);
      }
      
      static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(3);
      setColor(Qt::red);
      
      switch (lt) {
       case (Curves_arr::VERTEX):
        {
         Coord_point p(CGAL::to_double(e->target()->point().x()) /
                       m_tab_traits.COORD_SCALE, 
                       CGAL::to_double(e->target()->point().y()) /
                       m_tab_traits.COORD_SCALE); 
         static_cast<CGAL::Qt_widget&>(*this) << p;
         break;
        }
       case (Curves_arr::UNBOUNDED_VERTEX) :
        {
          Coord_point p(CGAL::to_double(e->target()->point().x()) /
                        m_tab_traits.COORD_SCALE, 
                        CGAL::to_double(e->target()->point().y()) /
                        m_tab_traits.COORD_SCALE); 
          static_cast<CGAL::Qt_widget&>(*this) << p;
        break;
        }
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
      static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(m_line_width);
    }
    setCursor(old);
  }

  /*!
   */
  void update_colors(std::vector<QColor> colors)
  {
    OverLay overlay_ob(this, colors);
    visit_faces(overlay_ob);
  }
  
  /*!
   */
  void set_face_color(Face_handle &f ,QColor& c)
  {  
    f->set_info(c);
    if( ! f->does_outer_ccb_exist() || empty)
      set_unbounded_face_color(c); 
  }

  /*!
   */
  template <class Function>
  void visit_faces(Function func)
  {
    Face_iterator  fi=m_curves_arr->faces_begin();
    for( ; fi != m_curves_arr->faces_end() ; ++fi )
      fi->set_visited(false);
    Face_handle ub = m_curves_arr->unbounded_face();
    visit_face_rec (ub,func) ;
}


  /*! antenna - return true if the halfedge and its
   *  twin point to the same face.
   */
  bool antenna(Halfedge & h)
  {
    Halfedge twin = *(h.opposite());
    return (twin.face() == h.face());
  }

  /*! draw a face and all its holes recursively
   */
  template<class Function>
  void visit_face_rec( Face_handle &f, Function func )
  {
    if (! f->visited())
    {
      Holes_iterator hit; // holes iterator
      func(f);
      f->set_visited(true);
      for(hit= f->holes_begin() ; hit!=f->holes_end() ; ++hit)
      {
        Ccb_halfedge_circulator cc = *hit;
        do {
          Halfedge he = *cc; 
          Halfedge he2 = *(he.opposite());
          Face_handle inner_face = he2.face(); // not sure about that
          if (antenna(he))
            continue;
          
          // move on to next hole
          visit_ccb_faces(inner_face , func); 
        }while (++cc != *hit);
      }// for
    }
  }// visit_face_rec

  template <class Function>
  void visit_ccb_faces(Face_handle & fh, Function func)
  {
    visit_face_rec(fh,func);
    Ccb_halfedge_circulator cc=fh->outer_ccb();
    do {
      Halfedge he = *cc;
      if(! he.opposite()->face()->visited())
      {
        Face_handle nei = (Face_handle) he.opposite()->face();
        visit_ccb_faces( nei ,func );
      }
      //created from the outer boundary of the face
    } while (++cc != fh->outer_ccb());
  }

  /*!
   */
  void overlay (Face_handle f , std::vector<QColor> & colors)
  {
    f->assign_overlay_info(colors); 
    if (f->does_outer_ccb_exist())  // f is not the unbounded face
    {
      /* running  around the outer of the face  */
      Ccb_halfedge_circulator cc = f->outer_ccb();
      do {
        Halfedge he = *cc;  // get the halfedge
        Data d = he.curve().get_data();  // get the data of the halfedge
        // get the original halfedge (from the original PM)
        Halfedge_handle original_he = d.halfedge_handle, temp_he;

        // getting the index of the original PM of the curve 
        int index = d.m_index;   
        if(m_tab_traits.curve_has_same_direction(cc))
          temp_he = original_he;
        else
          temp_he=(Halfedge_handle)original_he->opposite();

        // the original color of the face        
        QColor original_face_color = temp_he->face()->info(); 
        // update the colors vector in the face
        f->set_overlay_info( index , original_face_color);  
        // we want to update the colors vector to pass
        colors[index] = original_face_color; 
        // the information down the recurtion process.
        
      } while (++cc != f->outer_ccb());
      // blend the colors
      QColor c = blend_colors(f->OverlayInfoBegin() , f-> OverlayInfoEnd());
      f->set_info(c);
    }  
    else // f is the unbounded face
    {  
      // blend the colors
      QColor c = blend_colors(colors.begin() , colors.end() );
      f->set_info(c);
      set_unbounded_face_color(c);
    }
  }

  /*! draw_grid - draw the grid
   */
  void draw_grid()
  {
    setColor(Qt::white);
    static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(1);
    // get the edge coordinate
    int min_x = static_cast<int> (x_min());
    int max_x = static_cast<int> (x_max());
    int min_y = static_cast<int> (y_min());
    int max_y = static_cast<int> (y_max());
    
    // calculate cube size (minimum of 1)
    //int cube_size_x = std::max(1, abs(max_x - min_x)/20);
    //int cube_size_y = std::max(1, abs(max_y - min_y)/20);
    if (cube_size < std::abs(max_x - min_x)/40 ||
        cube_size < std::abs(max_y - min_y)/40)
      cube_size = std::max(std::max(1, std::abs(max_x - min_x)/20),
                           std::max(1, std::abs(max_y - min_y)/20));
    
    int cube_size_x = cube_size;
    int cube_size_y = cube_size;
    // draw the grid lines
    for (int i = min_x; i <= max_x; i += cube_size_x)
      static_cast<CGAL::Qt_widget&>(*this) <<
        Coord_segment(Coord_point(i, max_y + cube_size_y),
                      Coord_point( i , min_y - cube_size_y));
    for (int i = min_y; i <= max_y; i += cube_size_y)
      static_cast<CGAL::Qt_widget&>(*this) <<
        Coord_segment(Coord_point( max_x + cube_size_x , i ),
                      Coord_point( min_x - cube_size_x , i ));
  }
  
  /*! mousePressEvent - mouse click on the tab 
   *\ param e - mouse click event
   */
  void mousePressEvent(QMouseEvent *e)
  {
    QCursor old = cursor();
    setCursor(Qt::WaitCursor);

    if (mode == POINT_LOCATION || mode == RAY_SHOOTING || mode == FILLFACE)
    {
      mousePressEvent_point_location( e );
      setCursor(old);
      return;
    }

    if (mode == DELETE)
    {
      Hafledge_list_iterator itr;
      for (itr = removable_halfedges.begin();
           itr != removable_halfedges.end() ; ++itr)
        m_curves_arr->remove_edge(*itr);

      removable_halfedges.clear(); // clear the list which is now irrelevant
      redraw();

      if( m_curves_arr->number_of_vertices() == 0 )
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
      if( m_curves_arr->number_of_vertices() > 0 )
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
        m_tab_traits.first_point( p , mode );
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
        m_tab_traits.first_point( p , mode);
        split_point = Pm_point_2( p.x() * m_tab_traits.COORD_SCALE ,
                                  p.y() * m_tab_traits.COORD_SCALE);
      } 
      else
      {
        active = false;
        Pm_point_2 split_point2 =
          Pm_point_2(p.x() * m_tab_traits.COORD_SCALE,
                     p.y() * m_tab_traits.COORD_SCALE);
        const Xcurve split_curve = 
          m_tab_traits.curve_make_x_monotone(split_point , split_point2);
        Pm_point_2 p1;
        Pm_point_2 p_right;
        if (split_point.x() < split_point2.x())
          p_right = split_point;
        else
          p_right = split_point2;
        Halfedge_iterator hei;
        for (hei = m_curves_arr->halfedges_begin();
             hei != m_curves_arr->halfedges_end(); ++hei) 
        {
          const Xcurve & xcurve = hei->curve();
          m_tab_traits.draw_xcurve(this, xcurve);
          CGAL::Object res =
            m_traits.nearest_intersection_to_right(split_curve,
                                                   xcurve,p_right);
          if (CGAL::assign(p1, res))
            break;
        }

        Vertex_iterator   vit;
        bool is_already_vertex = false;
        // we dont want to split an already existed vertex...
        for (vit = m_curves_arr->vertices_begin(); 
             vit != m_curves_arr->vertices_end(); vit++)
        {
          Pm_point_2 temp = (*vit).point();
          if (p1 == temp)
            is_already_vertex = true;
        }
        m_tab_traits.draw_xcurve(this, hei->curve());
        if (hei != m_curves_arr->halfedges_end() && !is_already_vertex)
          m_tab_traits.split_edge(hei , p1 , this);
        
      }// else
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

      new_object(make_object(Coord_point(x * m_tab_traits.COORD_SCALE, 
                                         y * m_tab_traits.COORD_SCALE)));
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
  
  /*! find_removeable_halfedges - find removable curve in the tab 
   *\ param e - mouse click event
   */
  void find_removable_halfedges(QMouseEvent *e)
  {
    //  if the PM is empty do nothing
    if( m_curves_arr->number_of_vertices() == 0)
      return;
   
    if(!removable_halfedges.empty())
    {
      int i = (removable_halfedges.front()->curve()).get_data().m_index;
      setColor(colors[i]);
      Hafledge_list_iterator itr;
      for (itr = removable_halfedges.begin();
          itr != removable_halfedges.end() ; ++itr)
        m_tab_traits.draw_xcurve(this,(*itr)->curve());

      // clear the previous list which is now irrelevant
      removable_halfedges.clear(); 
    }

    // get the point of the mouse
    Coord_point p(x_real(e->x()) * m_tab_traits.COORD_SCALE ,
                  y_real(e->y()) * m_tab_traits.COORD_SCALE);
    
    bool is_first = true;       
    Coord_type min_dist = 0;
    Halfedge_iterator hei;
    Halfedge_iterator closest_hei;
    
    for (hei = m_curves_arr->halfedges_begin();
         hei != m_curves_arr->halfedges_end(); ++hei) 
    {
      Xcurve & xcurve = hei->curve();
      Coord_type dist = m_tab_traits.xcurve_point_distance(p, xcurve , this);
      if (is_first || dist < min_dist)
      {
        min_dist = dist;
        closest_hei = hei;
        is_first = false;
      }    
    }
    // now 'closest_hei' holds the cloeset halfedge to the point of the mouse

    if (remove_org_curve && !read_from_file)
    {
      const Base_curve * org_curve =
        m_tab_traits.get_origin_curve(closest_hei->curve());
      Hafledge_list li;
      Hafledge_list_iterator result;
      for (hei = m_curves_arr->halfedges_begin();
           hei != m_curves_arr->halfedges_end(); ++hei) 
      {
        const Base_curve * curve =
          m_tab_traits.get_origin_curve(hei->curve());
        result = std::find(li.begin(), li.end(), hei);
        if (curve == org_curve && result == li.end())
        {
          li.push_back(hei->twin());
          //m_curves_arr->remove_edge(hei);
          setColor(Qt::red);  // highlight the removable edge with red color
          m_tab_traits.draw_xcurve(this,hei->curve());
          removable_halfedges.push_back(hei);
        }
      }
    }
    else
    {
      // m_curves_arr->remove_edge(closest_hei);
      setColor(Qt::red);  // highlight the removable edge with red color
      m_tab_traits.draw_xcurve(this,closest_hei->curve());
      removable_halfedges.push_back(closest_hei);
    }
  }

  /*! mouseMoveEvent - enable seeing the line to be drawen 
   *\ param e - mouse click event
   */
  void mouseMoveEvent(QMouseEvent *e)
  {
    static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(m_line_width);
    if (mode == DELETE)
      // find removable edges , store them in the list 
      find_removable_halfedges(e);
    //'removable_halfedges' and highlight them

    if (mode == DRAG)
    {
      mouseMoveEvent_drag(e);
      return;
    }
    if (mode == MERGE && !first_time_merge)
    {
      if (second_curve != m_curves_arr->halfedges_end())
      {
        int i = (second_curve->curve()).get_data().m_index;
        setColor(colors[i]);
        m_tab_traits.draw_xcurve(this,second_curve->curve());
      }
      Coord_type x, y;
      x_real(e->x(), x);
      y_real(e->y(), y);
      Coord_point p(x * m_tab_traits.COORD_SCALE,
                    y * m_tab_traits.COORD_SCALE);
      second_curve = m_curves_arr->halfedges_end();
      m_tab_traits.find_close_curve(closest_curve, second_curve, p, this,
                                    true);
      setColor(Qt::red);
      m_tab_traits.draw_xcurve(this,closest_curve->curve());
      if (second_curve != m_curves_arr->halfedges_end())
      {
        setColor(Qt::green);
        m_tab_traits.draw_xcurve(this,second_curve->curve());
      }
      else
      {
        first_time_merge = true;
        redraw();
      }
      return;
    }// merge

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
       
       if( m_curves_arr->number_of_vertices() == 0 ) 
         return Coord_point(x , y);
       
       min_dist = m_tab_traits.closest_point(x,y,closest,this);                
       
       if (min_dist <= d)
         return closest;
       else
         return Coord_point(x , y);
       
       break;
      }
     case GRID:
      return Coord_point(getMid(x, xmin, xmax), 
                         getMid(y, ymin, ymax) );

     case NONE: break;
    }
    return Coord_point(x,y);
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
    if(e->button() == Qt::LeftButton && is_pure(e->state()))
    {
      if( m_curves_arr->number_of_vertices() == 0 )
      return;

      setColor(Qt::red);
      Coord_point p(x_real(e->x()) * m_tab_traits.COORD_SCALE ,
                    y_real(e->y()) * m_tab_traits.COORD_SCALE);
      bool is_first = true;
      Coord_type min_dist = 0;

      if (first_time_merge)
      {
        first_time_merge = false;
        Halfedge_iterator hei;
        //closest_curve = m_curves_arr->halfedges_end();
    
        for (hei = m_curves_arr->halfedges_begin();
             hei != m_curves_arr->halfedges_end(); ++hei) 
        {
          Vertex_iterator   vis = hei->source();
          Vertex_iterator   vit = hei->target();
          if (vis->degree() != 2 && vit->degree() != 2)
            continue;
          Xcurve & xcurve = hei->curve();
          Coord_type dist = m_tab_traits.xcurve_point_distance(p, xcurve,
                                                               this);
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
        second_curve = m_curves_arr->halfedges_end();
      }
      else
      {
        first_time_merge = true;
        m_tab_traits.find_close_curve(closest_curve, second_curve, p, this,
                                      false);
        redraw();
      }
    }
  }

  /*!
   */
  void update_curves_date_after_split(Halfedge_handle & e1, Xcurve & c1,
                                      Xcurve & c2)
  {
    Xcurve my_c1 , my_c2;  //temp variable

    // check that e1 is indeed the halfede of c1
    if(m_tab_traits.is_curve_and_halfedge_same_direction(e1,c1)) 
    {
      my_c1 = c1;
      my_c2 = c2;
    }
    else
    {
      my_c1 = c2;
      my_c2 = c1;
    }

    // update my_c1
    Data d1 = my_c1.get_data();

    if(m_tab_traits.is_curve_and_halfedge_same_direction(e1 , my_c1))
      d1.halfedge_handle = e1;
    else
      d1.halfedge_handle = (Halfedge_handle)e1->opposite();

    e1->curve().set_data( d1 );
    e1->opposite()->curve().set_data( d1 );

    // update my_c2
    Data d2 = my_c2.get_data();
    if (m_tab_traits.is_curve_and_halfedge_same_direction
        ((Halfedge_handle)e1->next() , my_c2))
      d2.halfedge_handle = (Halfedge_handle)e1->next();
    else
      d2.halfedge_handle = (Halfedge_handle) e1->next()->opposite();

    e1->next()->curve().set_data( d2 );
    e1->next()->opposite()->curve().set_data( d2 );
  }

  /*!
   */
  void update_curve_data_after_merge(Halfedge_handle &h ,Xcurve &c)
  { 
    Data d = c.get_data();
    if(m_tab_traits.is_curve_and_halfedge_same_direction(h , c))
      d.halfedge_handle = h;
    else
      d.halfedge_handle = (Halfedge_handle)h->opposite();

    h->curve().set_data( d );
    h->opposite()->curve().set_data( d );
  }
};

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/*! class Segment_tab_traits defines the segment traits
 */
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
  typedef Seg_halfedge           Halfedge;
  typedef  Seg_face_iterator    Face_iterator;
  
  //point location
  typedef Seg_trap_point_location      Trap_point_location;
  typedef Seg_naive_point_location     Naive_point_location;
  typedef Seg_simple_point_location    Simple_point_location;
  typedef Seg_walk_point_location      Walk_point_location;
  //typedef Seg_lenmarks_point_location  Lenmarks_point_location;
 
 
public:
 
  /*! coordinate scale - used in conics*/
  int COORD_SCALE;

  /*! constructor */
  Segment_tab_traits():
  COORD_SCALE(1)
  {}
  
  /*! distructor */
  ~Segment_tab_traits()
  {}

  /*! curve_has_same_direction - return true if the halfegde and
   *  its curve has the same direction
   */
  bool curve_has_same_direction( Ccb_halfedge_circulator &cc)
  {
    return (cc->curve().source() == cc->source()->point());
  }

  /*! check if curve and its halfedge are at the same direction */
  bool is_curve_and_halfedge_same_direction (const Halfedge_handle & he,
                                             const Xcurve & cv)
  {
    return (he->source()->point() == cv.source());
  }

  /*! fill_face - fill a face with its color (which is stored at the face)
   * it creates a polyong from the outer boundary of the face and 
   * uses CGAL opertaor << of polygons 
   */
  void fill_face(Qt_widget_demo_tab<Segment_tab_traits> * w , Face_handle f)
  {
    if (f->does_outer_ccb_exist())  // f is not the unbounded face
    {
      std::list< Coord_point > pts; // holds the points of the polygon         

      /* running with around the outer of the face and generate from it
       * polygon
       */
      Ccb_halfedge_circulator cc=f->outer_ccb();
      do {
        Coord_type x = CGAL::to_double(cc->source()->point().x());
        Coord_type y = CGAL::to_double(cc->source()->point().y());
        Coord_point coord_source(x , y);
        pts.push_back(coord_source );
        //created from the outer boundary of the face
      } while (++cc != f->outer_ccb());

      // make polygon from the outer ccb of the face 'f'
      Polygon pgn (pts.begin() , pts.end()); 

      w->setFilled(true);

      // fill the face according to its color (stored at any of her
      // incidents curves)
      if (! f->info().isValid())
        w->setFillColor(def_bg_color);
      else
        w->setFillColor(f->info());

      QPen old_penstyle = w->get_painter().pen();
      w->get_painter().setPen(Qt::NoPen);
      (*w) << pgn ;  // draw the polyong 
      w->setFilled(false);       
      w->get_painter().setPen(old_penstyle);
    }
    else
    {
      Coord_point points[4];
      points[0] = (Coord_point(w->x_min(),w->y_min()));
      points[1] = (Coord_point(w->x_min(),w->y_max()));
      points[2] = (Coord_point(w->x_max(),w->y_max()));
      points[3] = (Coord_point(w->x_max(),w->y_min()));
     
      w->setFilled(true);
      w->setFillColor(w->unbounded_face_color());
      
      QPen old_penstyle = w->get_painter().pen();
      w->get_painter().setPen(Qt::NoPen);
      (*w)<<Polygon(points , points +4 );
      w->setFilled(false);
       w->get_painter().setPen(old_penstyle);
    }
  }


  /*! draw_xcurve - use Qt_Widget operator to draw 
   *\ param w - the demo widget
   *\ c - curve to be drawen
   */
  void draw_xcurve(Qt_widget_demo_tab<Segment_tab_traits> * w , Xcurve c )
  {
    (*w) << c;
  }
  
  /*! first_point - a first point of inserted sgment
   */
  void first_point( Coord_point p , Mode m )
  {
    m_p1 = m_p2 = p;
  }
  
  /*! middle_point - the last point of a segment
   */
  void middle_point(Coord_point p , 
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
  
  /*! last_point - meaningless for segments
   */
  void last_point( Coord_point p , 
                   Qt_widget_demo_tab<Segment_tab_traits> * w )
  {
    return;
  }
  
  /*! get_segment - create a new segment, insert him into curves_list
   * and planar map
   */
  void get_segment( Coord_segment coord_seg ,
                    Qt_widget_demo_tab<Segment_tab_traits> * w)
  {
    Seg_notification  seg_notif;
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
    Halfedge_handle hh = w->m_curves_arr->insert(*seg, &seg_notif);
    CGAL::Bbox_2 curve_bbox = seg->bbox();
    w->bbox = w->bbox + curve_bbox;

    //Curve_data d = hh->curve().get_data();  // get the data of the halfedge
    //   Halfedge_handle original_he = d.halfedge_handle;
  }
  
  /*! curve_point_distance - return the distance between a point 
   * and a segment
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
   * and a xsegment
   */
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
  
  /*! get_origin_curve - return the origin base curve
   */
  const Pm_base_seg_2* get_origin_curve(const Pm_xseg_2 & xseg )
  {
    const Curve_data & curve_data = xseg.get_data();
    if (curve_data.m_type == Curve_data::LEAF) 
      return curve_data.m_ptr.m_curve;
    else // curve_data.m_type == Curve_data::INTERNAL
    {
      const Pm_xseg_2* xseg_p = curve_data.m_ptr.m_x_motonote_curve;      
      return get_origin_curve( *xseg_p );
    }
  }
  
  /*! draw_last_segment - call from mouse move event
   */
  void draw_last_segment( Qt_widget_demo_tab<Segment_tab_traits> * w)
  {
    *w << Coord_segment( m_p1 , m_p2 );
  }
  
  /*! draw_current_segment - call from mouse move event
   */
  void draw_current_segment( Coord_point p ,
                             Qt_widget_demo_tab<Segment_tab_traits> * w)
  {
    *w << Coord_segment( m_p1 , p);
    m_p2 = p;
  }
  
  /*! closest_point - find the closest point in the planar map
   * to a clicked point
   */
  Coord_type closest_point(Coord_type x, Coord_type y, Coord_point &closest,
                           Qt_widget_demo_tab<Segment_tab_traits> * w)
  {
    bool first = true;
    Coord_type x1,y1,dt,min_dist = 0;
    Vertex_iterator vit;
    for (vit = w->m_curves_arr->vertices_begin();
         vit != w->m_curves_arr->vertices_end(); vit++)
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

  /*! curve_make_x_monotone
   */
  const Xcurve curve_make_x_monotone(Pm_point_2 p1 , Pm_point_2 p2)
  {
    Data d;
    const Curve cv(Base_curve(p1 , p2) , d);
    std::list<Xcurve> xcurve_list;
    m_traits.curve_make_x_monotone(cv, std::back_inserter(xcurve_list));
    const Xcurve c1 = xcurve_list.front();
    return c1;
  }

  /*!
   */
  void find_close_curve(Halfedge_iterator &closest_curve ,Halfedge_iterator 
                        &second_curve ,Coord_point &p ,
                        Qt_widget_demo_tab<Segment_tab_traits> *w,
                        bool move_event)
  {
    bool is_first = true;
    Coord_type min_dist = 0;
    Halfedge_iterator hei;
    Pm_point_2 s = closest_curve->source()->point();
    Pm_point_2 t = closest_curve->target()->point();
    Vertex_iterator   vis = closest_curve->source();
    Vertex_iterator   vit = closest_curve->target();
    Kernel            ker;
    for (hei = w->m_curves_arr->halfedges_begin();
         hei != w->m_curves_arr->halfedges_end(); ++hei) 
    {
      Pm_point_2 s1 = hei->source()->point();
      Pm_point_2 t1 = hei->target()->point();

      if (merge_segments_precondition<Pm_point_2>(s, t, s1, t1,
                                                  vis->degree(),
                                                  vit->degree()) && 
          ker.compare_slope_2_object()(Kernel::Segment_2(s, t),
                                       Kernel::Segment_2(s1, t1)) ==
          CGAL::EQUAL)
      {
        Xcurve & xcurve = hei->curve();
        Coord_type dist = xcurve_point_distance( p, xcurve , w);
        if (is_first || dist < min_dist)
        {
          min_dist = dist;
          second_curve = hei;
          is_first = false;
        }    
      }
    }
    if (is_first) // didn't find any "good" curve
      return;

    else if (!move_event)
    {
      Pm_point_2 s1 = second_curve->source()->point();
      Pm_point_2 t1 = second_curve->target()->point();
      Base_curve * base = 0;

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
      Xcurve c = xcurve_list.front();
      Halfedge_handle h = w->m_curves_arr->merge_edge( closest_curve ,
                                                       second_curve , c); 
      
      w->update_curve_data_after_merge(h , c); 
      //// update c
      //
      // Curve_data d = c.get_data();
      //if(is_curve_and_halfedge_same_direction(h , c))
      //  d.halfedge_handle = h;
      //else
      //  d.halfedge_handle = (Halfedge_handle)h->opposite();
      
      //h->curve().set_data( d );
      //h->opposite()->curve().set_data( d );
      
    }
  }

  /*!
   */
  void split_edge(Halfedge_iterator &hei ,Pm_point_2 &p ,
                  Qt_widget_demo_tab<Segment_tab_traits> *w)
  {
    std::list<Xcurve> xcurve_list;

    Base_curve *base1 = new Base_curve;
    Base_curve *base2 = new Base_curve;

    Base_seg_traits  base_seg_traits;

    base_seg_traits.curve_split (hei->curve(),
                                 *base1, *base2,
                                 p);

    Curve_data cd1;
    cd1.m_type = Curve_data::LEAF;
    cd1.m_index = hei->curve().get_data().m_index;
    cd1.m_ptr.m_curve = base1;
    Curve *seg1 = new Curve( *base1, cd1 );

    m_traits.curve_make_x_monotone(*seg1, std::back_inserter(xcurve_list));
    Xcurve c1 = xcurve_list.front();

    Curve_data cd2;
    cd2.m_type = Curve_data::LEAF;
    cd2.m_index = hei->curve().get_data().m_index;
    cd2.m_ptr.m_curve = base2;
    Curve *seg2 = new Curve( *base2, cd2 );
    xcurve_list.clear();
    m_traits.curve_make_x_monotone(*seg2, std::back_inserter(xcurve_list));
    Xcurve c2 = xcurve_list.front();

    // split hei into curves : c1 ,c2 and store halfedge e1 (halfedge of c1)
    Halfedge_handle e1 = w->m_curves_arr->split_edge(hei , c1 , c2);

    // update the halfedge_handle field in the Curve_Data of c1 ,c2
    // Curve_Data

    w->update_curves_date_after_split(e1 , c1, c2);
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
  typedef Pol_halfedge           Halfedge;
  typedef  Pol_face_iterator    Face_iterator;

  //point location
  typedef Pol_trap_point_location      Trap_point_location;
  typedef Pol_naive_point_location     Naive_point_location;
  typedef Pol_simple_point_location    Simple_point_location;
  typedef Pol_walk_point_location      Walk_point_location;
  //typedef Pol_lenmarks_point_location  Lenmarks_point_location;
 

  /*! coordinate scale - used in conics*/
  int COORD_SCALE;
  
  /*! constructor */
  Polyline_tab_traits():
  COORD_SCALE(1)
  {}
  
  /*! distructor */
  ~Polyline_tab_traits()
  {}

 
  /*! curve_has_same_direction - return true if the curve and
   *  the halfedge has the same direction
   */
  bool curve_has_same_direction( Ccb_halfedge_circulator &cc)
  {
    return ( *(cc->curve().begin()) == cc->source()->point());
  }

  /*!
   */
  bool is_curve_and_halfedge_same_direction(const Halfedge_handle & he,
                                            const Xcurve & cv)
  {
    return (he->source()->point() == *(cv.begin()));
  }

  /*! fill_face - fill a face with its color (which is stored at the curves)
   * it creates a polyong from the outer boundary of the face and 
   * uses CGAL opertaor << of polygons 
   */
  void fill_face(Qt_widget_demo_tab<Polyline_tab_traits> * w , Face_handle f)
  {
    if (f->does_outer_ccb_exist())  // f is not the unbounded face
    {
      std::list< Coord_point > pts; // holds the points of the polygon
      Xcurve::const_iterator   pt_itr;
      Xcurve cv;

      /* running with around the outer of the face and generate from it
       * polygon
       */
      Ccb_halfedge_circulator cc=f->outer_ccb();
      do {
        cv = cc->curve();
        if( curve_has_same_direction (cc) )
        {
          for( pt_itr = cv.begin() , ++pt_itr ; pt_itr != cv.end(); ++pt_itr)
          {
            Coord_type x = CGAL::to_double((*pt_itr).x());
            Coord_type y = CGAL::to_double((*pt_itr).y());
            Coord_point coord_source(x , y);
            pts.push_back(coord_source );
          }
        }
        else
        {
          for( pt_itr = cv.rbegin() , ++pt_itr ; pt_itr != cv.rend();
               ++pt_itr)
          {
            Coord_type x = CGAL::to_double((*pt_itr).x());
            Coord_type y = CGAL::to_double((*pt_itr).y());
            Coord_point coord_source(x , y);
            pts.push_back(coord_source );
          }
        }
        //created from the outer boundary of the face
      } while (++cc != f->outer_ccb());

      // make polygon from the outer ccb of the face 'f'
      Polygon pgn (pts.begin() , pts.end());

      w->setFilled(true);

      // fill the face according to its color (stored at any of her
      // incidents curves)
      if (! f->info().isValid())
        w->setFillColor(def_bg_color);
      else
        w->setFillColor(f->info());
      QPen old_penstyle = w->get_painter().pen();
      w->get_painter().setPen(Qt::NoPen);

      (*w) << pgn ;  // draw the polyong 
      w->setFilled(false);  
      w->get_painter().setPen(old_penstyle);
    }
    else
    {
      Coord_point points[4];
      points[0] = (Coord_point(w->x_min(),w->y_min()));
      points[1] = (Coord_point(w->x_min(),w->y_max()));
      points[2] = (Coord_point(w->x_max(),w->y_max()));
      points[3] = (Coord_point(w->x_max(),w->y_min()));
     
      w->setFilled(true);
      w->setFillColor(w->unbounded_face_color());
      
      QPen old_penstyle = w->get_painter().pen();
      w->get_painter().setPen(Qt::NoPen);
      (*w)<<Polygon(points , points +4 );
      w->setFilled(false);
       w->get_painter().setPen(old_penstyle);
    }
  }
  
  /*! draw_xcurve - go over the polyline parts and use Qt_Widget operator
   * to draw
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
  
  /*! first_point - a first point of inserted polyline or a splitter
   */
  void first_point( Coord_point p , Mode m)
  {
    last_of_poly = p;
    if (m == INSERT)
      points.push_back(Pm_pol_point_2(p.x(),p.y()));    
  }
  
  /*! middle_point - a middle point of a polyline
   */
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
   * polyline and reset 
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
   * and a polyline
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
   * and a polyline
   */
  Coord_type xcurve_point_distance(Coord_point p, Xcurve & c,
                                   Qt_widget_demo_tab<Polyline_tab_traits> *
                                   w)
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
  
  /*! draw_last_segment - call from mouse move event
   */
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
   * to a clicked point
   */
  Coord_type closest_point(Coord_type x, Coord_type y, Coord_point &closest,
                           Qt_widget_demo_tab<Polyline_tab_traits> * w)
  {
    bool first = true;
    Coord_type x1,y1,dt,min_dist = 0;
    Halfedge_iterator heit;
    for (heit = w->m_curves_arr->halfedges_begin();
         heit != w->m_curves_arr->halfedges_end(); heit++)
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
  
  /*! get_origin_curve - return the origin base curve
   */
  const Pm_base_pol_2* get_origin_curve(const Pm_xpol_2 & xseg )
  {
    const Curve_pol_data & curve_data = xseg.get_data();
    if (curve_data.m_type == Curve_pol_data::LEAF) 
      return curve_data.m_ptr.m_curve;
    else // curve_data.m_type == Curve_data::INTERNAL
    {
      const Pm_xpol_2* xseg_p = curve_data.m_ptr.m_x_motonote_curve;        
      return get_origin_curve( *xseg_p );
    }
  }

  /*! curve_make_x_monotone
   */
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
 
  /*!
   */
  void find_close_curve(Halfedge_iterator & closest_curve,
                        Halfedge_iterator & second_curve, Coord_point & p,
                        Qt_widget_demo_tab<Polyline_tab_traits> * w,
                        bool move_event)
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

    for (hei = w->m_curves_arr->halfedges_begin();
         hei != w->m_curves_arr->halfedges_end(); ++hei) 
    {
      Pm_point_2 s1 = *(hei->curve().begin());
      Pm_point_2 t1 = *(hei->curve().rbegin());

      if (merge_segments_precondition<Pm_point_2>(s, t, s1, t1,
                                                  vis->degree(),
                                                  vit->degree()))
      {
        Xcurve & xcurve = hei->curve();
        Coord_type dist = xcurve_point_distance( p, xcurve , w);
        if (is_first || dist < min_dist)
        {
          min_dist = dist;
          second_curve = hei;
          is_first = false;
        }    
      }
    }
    if (is_first) // didn't find any "good" curve
      return;
    else if (!move_event)
    {
      Xcurve & c = closest_curve->curve();
      Xcurve & c1 = second_curve->curve();

      Pm_point_2 s1 = *(second_curve->curve().begin());
      Pm_point_2 t1 = *(second_curve->curve().rbegin());
      std::vector<Pm_pol_point_2> temp_points;
      Curve_const_iterator cit;
 
      if (t == s1 || t == t1)
      {
        for (cit = c.begin(); cit != c.end(); cit++)
        {
          temp_points.push_back(*cit);
        }
      }
      else
      {
        for (cit = c.rbegin(); cit != c.rend(); cit++)
        {
          temp_points.push_back(*cit);
        }
      }
      Kernel ker;
      if (ker.compare_slope_2_object()(Kernel::Segment_2 (s, t),
                                       Kernel::Segment_2 (s1, t1)) ==
          CGAL::EQUAL)
      {
        temp_points.pop_back();
      }

      if (s1 == s || s1 == t)
      {
        cit = c1.begin(), cit++;
        for ( ; cit != c1.end(); cit++)
        {
          temp_points.push_back(*cit);
        }
      }
      else
      {
        cit = c1.rbegin(), cit++;
        for (; cit != c1.rend(); cit++)
        {
          temp_points.push_back(*cit);
        }
      }

      Base_curve *base = new Base_curve(temp_points.begin(),
                                        temp_points.end());
      Curve_pol_data cd;
      cd.m_type = Curve_pol_data::LEAF;
      cd.m_index = w->index;
      cd.m_ptr.m_curve = base;
      Curve *curve = new Curve( *base, cd );
      std::list<Xcurve> xcurve_list;
      m_traits.curve_make_x_monotone(*curve, std::back_inserter(xcurve_list));
      Xcurve cc = xcurve_list.front();
      Halfedge_handle h = w->m_curves_arr->merge_edge(closest_curve ,
                                                      second_curve, cc); 
      w->update_curve_data_after_merge(h , cc);
    }
  }

  /*!
   */
  void split_edge(Halfedge_iterator &hei ,Pm_point_2 &p ,
                  Qt_widget_demo_tab<Polyline_tab_traits> *w)
  {
    Base_curve *base1 = new Base_curve;
    Base_curve *base2 = new Base_curve;

    Base_pol_traits  base_pol_traits;

    base_pol_traits.curve_split (hei->curve(), *base1, *base2, p);

    std::list<Xcurve> xcurve_list;

    //Base_curve *base1 = new Base_curve(temp_points1.begin(),
    // temp_points1.end());
    Curve_pol_data cd1;
    cd1.m_type = Curve_pol_data::LEAF;
    cd1.m_index = hei->curve().get_data().m_index;
    cd1.m_ptr.m_curve = base1;
    Curve *seg1 = new Curve( *base1, cd1 );
    m_traits.curve_make_x_monotone(*seg1, std::back_inserter(xcurve_list));
    Xcurve c1 = xcurve_list.front();

    xcurve_list.clear();

    // Base_curve *base2 = new Base_curve(temp_points2.begin(),
    // temp_points2.end());
    Curve_pol_data cd2;
    cd2.m_type = Curve_pol_data::LEAF;
    cd2.m_index = hei->curve().get_data().m_index;
    cd2.m_ptr.m_curve = base2;
    Curve *seg2 = new Curve( *base2, cd2 );
    m_traits.curve_make_x_monotone(*seg2, std::back_inserter(xcurve_list));
    Xcurve c2 = xcurve_list.front();
    Halfedge_handle e1 = w->m_curves_arr->split_edge(hei , c1 , c2);

    w->update_curves_date_after_split(e1 , c1, c2);
  }

private:
  
  /*! get_polyline - create a new polyline
   */
  void get_polyline(Qt_widget_demo_tab<Polyline_tab_traits> * w)
  {
    Pol_notification pol_notif;
    Pm_base_pol_2 * base_pol_p =
      new Pm_base_pol_2(points.begin(), points.end());
    Curve_pol_data cd;
    cd.m_type = Curve_pol_data::LEAF;
    cd.m_index = w->index;
    cd.m_ptr.m_curve = base_pol_p;
    Curve * pol = new Curve( *base_pol_p, cd );
    w->m_curves_arr->insert( *pol , &pol_notif);
    CGAL::Bbox_2 curve_bbox = pol->bbox();
    w->bbox = w->bbox + curve_bbox;
  }
  
  /*! convert - convert from Pm_pol_curve to Coord_segment
   */
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
  
  /*! convert - convert from const Pm_pol_curve to Coord_segment
   */
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
/*!
 */
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
  typedef Conic_halfedge           Halfedge;
  typedef  Conic_face_iterator    Face_iterator;
  
  //point location
  typedef Conic_trap_point_location      Trap_point_location;
  typedef Conic_naive_point_location     Naive_point_location;
  typedef Conic_simple_point_location    Simple_point_location;
  typedef Conic_walk_point_location      Walk_point_location;
  //typedef Conic_lenmarks_point_location  Lenmarks_point_location;
 
  /*! coordinate scale - used in conics*/
  int COORD_SCALE;
  int DRAW_FACTOR;

  /*! constructor */
  Conic_tab_traits():
  COORD_SCALE(1),
  DRAW_FACTOR(5)
  {}
  
  /*! distructor */
  ~Conic_tab_traits()
  {}

  /*! curve_has_same_direction - return true if the curve and
   *  the halfedge has the same direction
   */
  bool curve_has_same_direction( Ccb_halfedge_circulator &cc)
  {
    return (cc->curve().source() == cc->source()->point());
  }

  /*! check if curve and its halfedge are at the same direction
   */
  bool is_curve_and_halfedge_same_direction(const Halfedge_handle & he,
                                            const Xcurve & cv)
  {
    return (he->source()->point() == cv.source());
  }

  /*!
   */
  void fill_face(Qt_widget_demo_tab<Conic_tab_traits> * w , Face_handle f)
  {
    if (f->does_outer_ccb_exist())  // f is not the unbounded face
    {
      std::list< Coord_point > pts; // holds the points of the polygon         
      /* running with around the outer of the face and generate from it
       * polygon
       */
      Ccb_halfedge_circulator cc=f->outer_ccb();
      do {
        if(w->antenna(*cc))
          continue;
        Xcurve c = cc->curve(); 
        // Get the co-ordinates of the curve's source and target.
        double sx = CGAL::to_double(cc->source()->point().x()),
               sy = CGAL::to_double(cc->source()->point().y()),
               tx = CGAL::to_double(cc->target()->point().x()),
               ty = CGAL::to_double(cc->target()->point().y());
        
        Coord_point coord_source(sx / COORD_SCALE, sy / COORD_SCALE);
        Coord_point coord_target(tx / COORD_SCALE, ty / COORD_SCALE);

        if (c.is_segment())
            pts.push_back(coord_source ); 
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
            double   curr_x, curr_y;
            int      x;
  
            Pm_conic_point_2 ps[2];
            int nps;

            pts.push_back(coord_source ); 

            if (is_source_left)
            {
              for (x = x_min + DRAW_FACTOR; x < x_max; x+=DRAW_FACTOR)
                //= COORD_SCALE)
              {
                curr_x = (*w).x_real(x);
                nps = c.get_points_at_x(Pm_conic_point_2(curr_x, 0), ps);
                if (nps == 1)
                {
                  curr_y = CGAL::to_double(ps[0].y());
                  pts.push_back(Coord_point(curr_x / COORD_SCALE,
                                            curr_y / COORD_SCALE));
                }// if
              }// for            
            }
            else
            {              
              for (x = x_max; x > x_min; x-=DRAW_FACTOR)
              {
                curr_x = (*w).x_real(x);
                nps = c.get_points_at_x(Pm_conic_point_2(curr_x, 0), ps);
                if (nps == 1)
                {
                  curr_y = CGAL::to_double(ps[0].y());
                  pts.push_back(Coord_point(curr_x / COORD_SCALE,
                                            curr_y / COORD_SCALE));
                }// if
              }// for       
            }// else
            pts.push_back(coord_target ); 
          }// if
          else
          {
            // We should never reach here.
            CGAL_assertion(false);
          }// else
        }// else
        //created from the outer boundary of the face
      } while (++cc != f->outer_ccb());

       // make polygon from the outer ccb of the face 'f'
      Polygon pgn (pts.begin() , pts.end());
      QPen old_penstyle = w->get_painter().pen();
      w->get_painter().setPen(Qt::NoPen);
      w->setFilled(true);

      // fill the face according to its color (stored at any of her incidents
      // curves)
      if (! f->info().isValid())
        w->setFillColor(def_bg_color);
      else
        w->setFillColor(f->info());

      (*w) << pgn ;  // draw the polyong 
      w->get_painter().setPen(old_penstyle);
      w->setFilled(false);
    }
    else
    {
      Coord_point points[4];
      points[0] = (Coord_point(w->x_min(),w->y_min()));
      points[1] = (Coord_point(w->x_min(),w->y_max()));
      points[2] = (Coord_point(w->x_max(),w->y_max()));
      points[3] = (Coord_point(w->x_max(),w->y_min()));
     
      w->setFilled(true);
      w->setFillColor(w->unbounded_face_color());
      
      QPen old_penstyle = w->get_painter().pen();
      w->get_painter().setPen(Qt::NoPen);
      (*w)<<Polygon(points , points +4 );
      w->setFilled(false);
       w->get_painter().setPen(old_penstyle);
    }
  }
   
  /*! draw_xcurve - same as draw_curve
   */
  void draw_xcurve(Qt_widget_demo_tab<Conic_tab_traits> * w , Xcurve c )
  {
    // Get the co-ordinates of the curve's source and target.
    double sx = CGAL::to_double(c.source().x()),
      sy = CGAL::to_double(c.source().y()),
      tx = CGAL::to_double(c.target().x()),
      ty = CGAL::to_double(c.target().y());
    
    if (c.is_segment())
    {
      Coord_point coord_source(sx / COORD_SCALE, sy / COORD_SCALE);
      Coord_point coord_target(tx / COORD_SCALE, ty / COORD_SCALE);
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
        
        for (x = x_min + DRAW_FACTOR; x < x_max; x+=DRAW_FACTOR)
        {
          curr_x = (*w).x_real(x);
          nps = c.get_points_at_x(Pm_conic_point_2(curr_x, 0), ps);
          if (nps == 1)
          {
            curr_y = CGAL::to_double(ps[0].y());
            (*w) << Coord_segment (Coord_point(prev_x / COORD_SCALE,
                                               prev_y / COORD_SCALE),
                                   Coord_point(curr_x / COORD_SCALE,
                                               curr_y / COORD_SCALE));
            prev_x = curr_x;
            prev_y = curr_y;
            
          }
        }
        
        (*w) << Coord_segment (Coord_point(prev_x / COORD_SCALE,
                                           prev_y / COORD_SCALE),
                               Coord_point(end_x / COORD_SCALE,
                                           end_y / COORD_SCALE));
      }
      else
      {
        // We should never reach here.
        CGAL_assertion(false);
      }
    }
  }
  
  /*! get_origin_curve - return the origin base curve
   */
  const Pm_base_conic_2* get_origin_curve(const Pm_xconic_2 & xseg )
  {
    const Curve_conic_data & curve_data = xseg.get_data();
    if (curve_data.m_type == Curve_conic_data::LEAF) 
      return curve_data.m_ptr.m_curve;
    else // curve_data.m_type == Curve_data::INTERNAL
    {
      Pm_xconic_2* xseg_p = curve_data.m_ptr.m_x_motonote_curve;        
      return get_origin_curve( *xseg_p );
    }
  }
  
  ////////////////////////////////////////////////////////////////////////////
  
  /*! first_point - a first point of inserted sgment
   */
  void first_point( Coord_point p , Mode m)
  {
    m_p_old = m_p1 = m_p2 = m_p3 = m_p4 = p;
    num_points = 1;
    first_time = true;
  }
  
  /*! middle_point - the last point of a segment
   */
  void middle_point( Coord_point p , Qt_widget_demo_tab<Conic_tab_traits> * w)
  {
    CfNT r, s, t, u, v, ww;  // The conic coefficients.
    CfNT a, b, c, a_sq, b_sq;
    CfNT x, y, x1, y1, x0, y0, temp;
    CfNT sq_rad;

    x1 = CfNT(static_cast<int>(COORD_SCALE * m_p1.x()));
    y1 = CfNT(static_cast<int>(COORD_SCALE * m_p1.y()));
    x = CfNT(static_cast<int>(COORD_SCALE * p.x()));
    y = CfNT(static_cast<int>(COORD_SCALE * p.y()));

    if(x != x1 || y != y1) 
    {      
      Pm_base_conic_2 *cv = NULL;

      switch (w->conic_type)
      {
       case CIRCLE:
        sq_rad = CGAL::square(x - x1) + CGAL::square(y - y1);
        cv = new Pm_base_conic_2(Int_circle_2 (Int_point_2(x1, y1), sq_rad)); 
        break;

       case SEGMENT:
        cv = new Pm_base_conic_2(Int_segment_2 (Int_point_2(x,y),
                                                Int_point_2(x1,y1)));
        break;

       case ELLIPSE:
        if (y == y1 || x == x1)
        {
          QMessageBox::information( w, "Insert Ellipse", "Invalid Ellipse");
          w->active = false;
          return;
        }

        a = CORE::abs(x1 - x)/2;
        b = CORE::abs(y1 - y)/2;
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

        // RWRW: Do nothing ...
       case PARABOLA:
        *w << CGAL::LineWidth(3);
        if (num_points == 1)
        {
          m_p2 = p;
          num_points++;
          *w << m_p2;
          return;
        }
        if (num_points == 2)
        {
          CfNT x2 = CfNT(static_cast<int>(COORD_SCALE * m_p2.x()));
          CfNT y2 = CfNT(static_cast<int>(COORD_SCALE * m_p2.y()));
          Int_kernel ker;
          // the three points of the parabola cannot be collinear (
          if (ker.collinear_2_object()(Int_point_2(x1,y1),
                                       Int_point_2(x,y),Int_point_2(x2,y2)))
          {
            QMessageBox::information( w, "Insert Conic", "Invalid Conic");
            w->active = false;
            *w << m_p1 << m_p2 ;
            return; 
          }
          cv = new Pm_base_conic_2 (Int_point_2(x1,y1),Int_point_2(x2,y2),
                                    Int_point_2(x,y));
        }
        break;

       case HYPERBOLA:
        *w << CGAL::LineWidth(3);
        if (num_points == 1)
        {
          m_p2 = p;
          num_points++;
          *w << m_p2;
          return;
        }
        if (num_points == 2)
        {
          m_p3 = p;
          num_points++;
          *w << m_p3;
          return;
        }
        if (num_points == 3)
        {
          m_p4 = p;
          num_points++;
          *w << m_p4;
          return;
        }
        if (num_points == 4)
        { 
          *w << p;
          CfNT x2 = CfNT(static_cast<int>(COORD_SCALE * m_p2.x()));
          CfNT y2 = CfNT(static_cast<int>(COORD_SCALE * m_p2.y()));
          CfNT x3 = CfNT(static_cast<int>(COORD_SCALE * m_p3.x()));
          CfNT y3 = CfNT(static_cast<int>(COORD_SCALE * m_p3.y()));
          CfNT x4 = CfNT(static_cast<int>(COORD_SCALE * m_p4.x()));
          CfNT y4 = CfNT(static_cast<int>(COORD_SCALE * m_p4.y()));
          cv = new Pm_base_conic_2 (Int_point_2(x1,y1),Int_point_2(x2,y2),
                                    Int_point_2(x3,y3),Int_point_2(x4,y4),
                                    Int_point_2(x,y));
          if (! cv->is_valid())
          {
            QMessageBox::information( w, "Insert Conic", "Invalid Conic");
            w->active = false;
            *w << m_p1 << m_p2 << m_p3 << m_p4 << p;
            return;
          }
        }
        break;
      }
      Curve_conic_data cd;
      cd.m_type = Curve_conic_data::LEAF;
      cd.m_index = w->index;
      cd.m_ptr.m_curve = cv;
      Conic_notification conic_notif;
      w->m_curves_arr->insert(Pm_conic_2( *cv , cd), & conic_notif);
      CGAL::Bbox_2 curve_bbox = cv->bbox();
      w->bbox = w->bbox + curve_bbox; 
      w->active = false;
      //w->redraw();  // not working so I use new_object insted
      w->new_object(make_object(Coord_segment(m_p1 , p)));
    } 
  }
  
  /*! last_point - meaningless for conics - at least for now
   */
  void last_point( Coord_point p , Qt_widget_demo_tab<Conic_tab_traits> * w )
  {   
    return;
  }

  /*! draw_last_segment - call from mouse move event
   */
  void draw_last_segment( Qt_widget_demo_tab<Conic_tab_traits> * w)
  {
    if (w->mode == SPLIT)
      *w << Coord_segment( m_p1 , m_p_old );
    else
    {
      switch (w->conic_type)
      {
       case CIRCLE:
        *w << Coord_circle(m_p1,
                           pow(m_p1.x() - m_p_old.x(), 2) +
                           pow(m_p1.y() - m_p_old.y(),2));
        break;
       case SEGMENT:
        *w << Coord_segment( m_p1 , m_p_old );
        break;
       case ELLIPSE:
        {
         *w << Coord_segment( Coord_point(m_p1.x(),m_p_old.y()) , m_p_old );
         *w << Coord_segment( Coord_point(m_p1.x(),m_p_old.y()) , m_p1 );
         *w << Coord_segment( Coord_point(m_p_old.x(),m_p1.y()) , m_p_old );
         *w << Coord_segment( Coord_point(m_p_old.x(),m_p1.y()) , m_p1 );
         break;
        }
       case PARABOLA:
        if (first_time)
        {
          *w << CGAL::LineWidth(3);
          *w << m_p1;
          first_time = false;
        }
        break;
       case HYPERBOLA:
        if (first_time)
        {
          *w << CGAL::LineWidth(3);
          *w << m_p1;
          first_time = false;
        }
        break;
      }
    }
  }
  
  /*! draw_current_segment - call from mouse move event
   */
  void draw_current_segment( Coord_point p ,
                             Qt_widget_demo_tab<Conic_tab_traits> * w)
  {
    m_p_old = p;
    draw_last_segment(w);
  }

  /*! closest_point - find the closest point in the planar map
   * to a clicked point
   */
  Coord_type closest_point(Coord_type x, Coord_type y, Coord_point &closest,
                           Qt_widget_demo_tab<Conic_tab_traits> * w)
  {
    bool first = true;
    Coord_type x1,y1,dt,min_dist = 0;
    //Coord_type x_scale = COORD_SCALE * x;
    //Coord_type y_scale = COORD_SCALE * y;
    Vertex_iterator vit;
    for (vit = w->m_curves_arr->vertices_begin();
         vit != w->m_curves_arr->vertices_end(); vit++)
    {
      const Pm_point_2& p = (*vit).point();
      x1 = CGAL::to_double(p.x()) / COORD_SCALE;
      y1 = CGAL::to_double(p.y()) / COORD_SCALE;
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
   * a point and a conic
   */
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
   * a point and a conic
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
  
  /*! curve_make_x_monotone
   */
  const Xcurve curve_make_x_monotone(Pm_point_2 p1 , Pm_point_2 p2)
  {
    /*! */
    Data d;

    // RWRW:
    Int_point_2          my_p1 (static_cast<int> (CGAL::to_double(p1.x())),
                                static_cast<int> (CGAL::to_double(p1.y())));
    Int_point_2          my_p2 (static_cast<int> (CGAL::to_double(p2.x())),
                                static_cast<int> (CGAL::to_double(p2.y())));
    Pm_base_conic_2 base (Int_segment_2 (my_p1, my_p2));
    // Pm_base_conic_2 base(Pm_conic_segment_2(p1 , p2));

    const Pm_conic_2 cv( base , d);
    std::list<Xcurve> xcurve_list;
    m_traits.curve_make_x_monotone(cv, std::back_inserter(xcurve_list));
    const Xcurve c1 = xcurve_list.front();
    return c1;
  }

  /*!
   */
  void find_close_curve(Halfedge_iterator & closest_curve,
                        Halfedge_iterator & second_curve, Coord_point & p,
                        Qt_widget_demo_tab<Conic_tab_traits> * w,
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
    for (hei = w->m_curves_arr->halfedges_begin();
         hei != w->m_curves_arr->halfedges_end(); ++hei) 
    {
      if (first.has_same_base_conic(hei->curve()))
      {
        Pm_point_2 s1 = hei->curve().source();
        Pm_point_2 t1 = hei->curve().target();

        if(  merge_segments_precondition<Pm_point_2>(s,t, s1,t1,
                                                     vis->degree(),
                                                     vit->degree()))
        { 
          Xcurve & xcurve = hei->curve();
          Coord_type dist = xcurve_point_distance( p, xcurve , w);
          if (is_first || dist < min_dist)
          {
            min_dist = dist;
            second_curve = hei;
            is_first = false;
          }    
        }
      }
    }
    if (is_first) // didn't find any "good" curve
      return;
    else if (!move_event)
    {
      Xcurve curve;
      const Xcurve second = second_curve->curve();

      m_traits.curve_merge (first, second, curve);

      Curve_conic_data cd1;
      cd1.m_type = Curve_conic_data::LEAF;
      cd1.m_index = first.get_data().m_index;
      cd1.m_ptr.m_curve = first.get_data().m_ptr.m_curve;

      Conic_notification conic_notif;

      // save the colors of the incident faces (to restore them later)
      QColor color1 = closest_curve->face()->info();
      QColor color2 = closest_curve->opposite()->face()->info();

      Pm_point_2 closest_curve_source = closest_curve->source()->point();
      Pm_point_2 closest_curve_target = closest_curve->target()->point();

      w->m_curves_arr->remove_edge(closest_curve);
      w->m_curves_arr->remove_edge(second_curve);

      Pm_conic_2 c( curve , cd1);
      Halfedge_handle h = w->m_curves_arr->insert(c );

      //restore the color of the incident faces of 'h'
      //check that 'closest_curev' and 'h' are at the dame direction
      if(closest_curve_source == h->source()->point() || 
         closest_curve_target == h->target()->point()    )
      {
        h->face()->set_info(color1);
        h->opposite()->face()->set_info(color2);
      }
      else  
      {
        h->face()->set_info(color2);
        h->opposite()->face()->set_info(color1);
      }

      Data d = c.get_data();
      if(is_curve_and_halfedge_same_direction(h , curve))
        d.halfedge_handle = h;
      else
        d.halfedge_handle = (Halfedge_handle)h->opposite();

      h->curve().set_data( d );
      h->opposite()->curve().set_data( d );
      
      //w->update_curve_data_after_merge(h , c);
      //// update c
      //
      // Curve_conic_data d = c.get_data();
      //d.halfedge_handle = h;
      
      //h->curve().set_data( d );
      //h->opposite()->curve().set_data( d );
    }
  }

  /*!
   */
  void split_edge(Halfedge_iterator &hei ,Pm_point_2 & p,
                  Qt_widget_demo_tab<Conic_tab_traits> * w)
  {    
    const Xcurve split_curve = hei->curve();
    Xcurve sbc1;
    Xcurve sbc2;
    
    m_traits.curve_split(split_curve , sbc1 , sbc2 , p);

    Curve_conic_data cd1;
    cd1.m_type = Curve_conic_data::LEAF;
    cd1.m_index = hei->curve().get_data().m_index;
    cd1.m_ptr.m_curve = hei->curve().get_data().m_ptr.m_curve;

    Pm_conic_2    sub_curve1 (sbc1, cd1);

    Curve_conic_data cd2;
    cd2.m_type = Curve_conic_data::LEAF;
    cd2.m_index = hei->curve().get_data().m_index;
    cd2.m_ptr.m_curve =hei->curve().get_data().m_ptr.m_curve;

    Pm_conic_2    sub_curve2 (sbc2, cd2);
    Conic_notification conic_notif;
    
    //  w->m_curves_arr->remove_edge(hei);
    Halfedge_handle h1 = w->m_curves_arr->insert(sub_curve1 , & conic_notif);
    Halfedge_handle h2 = w->m_curves_arr->insert(sub_curve2 , & conic_notif);

    // 
    // w->update_curves_date_after_split(hei , sub_curve1, sub_curve2);
    
    //// update data
    //
    // Curve_conic_data d1 = sub_curve1.get_data();
    //d1.halfedge_handle = h1;

    //h1->curve().set_data( d1 );
    //h1->opposite()->curve().set_data( d1 );

    // Curve_conic_data d2 = sub_curve2.get_data();
    //d2.halfedge_handle = h2;

    //h2->curve().set_data( d2 );
    //h2->opposite()->curve().set_data( d2 );
  }

  Traits m_traits;
  /*! temporary points of the created conic */
  Coord_point m_p_old,m_p1,m_p2,m_p3,m_p4;
  /*! bool flag for hyperbola insertion */
  bool first_time;
  /*! counter for the number of points */
  int num_points;
};

typedef Qt_widget_demo_tab<Segment_tab_traits> Qt_widget_segment_tab;
typedef Qt_widget_demo_tab<Polyline_tab_traits> Qt_widget_polyline_tab;
typedef Qt_widget_demo_tab<Conic_tab_traits> Qt_widget_conic_tab;

#endif //DEMO_TAB_H
