// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/features/gsoc2012-Arrangement_on_surface_2-demo-atsui/Arrangement_on_surface_2/demo/Arrangement_on_surface_2/demo_tab.h $
// $Id: demo_tab.h 67117 2012-01-13 18:14:48Z lrineau $
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Efi Fogel       <efif@post.tau.ac.il>

#ifndef DEMO_TAB_H
#define DEMO_TAB_H

/*! demo_tab.h contain the definetion and implementation of
 *  the demo tab classes and the tabs traits classes.
 *  all the possible shared code is in Qt_widget_demo_tab where
 *  the differences is in the traits classes.
 */

#include <stddef.h>
#include <math.h>
#include <limits>

#include <qrect.h>
#include <qcursor.h>
#include <qmessagebox.h>
#include <qcolor.h>
#include <qpainter.h>
#include <qpen.h>

#include "cgal_types.h"

#include <CGAL/IO/pixmaps/hand.xpm>
#include <CGAL/IO/pixmaps/holddown.xpm>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>
#include <CGAL/iterator.h> //for CGAL::Oneset_iterator<T>
#include <CGAL/Object.h>
#include <CGAL/envelope_2.h>
#include <CGAL/Envelope_diagram_1.h>

#include <vector>

//global variables defined at MyWindow_operations.cpp
extern bool lower_env;
extern bool upper_env;
/*! class Qt_widget_base_tab - inherits from CGAL::Qt_widget
 *  contain all data members that are not part of the traits
 */
class Qt_widget_base_tab : public CGAL::Qt_widget
{
private:
  //! Default color scheme ordering initialized in demo_tab.cpp 
  static QColor s_color_order[];
	
public:

  /*!
   */
  Qt_widget_base_tab(TraitsType  t ,  QWidget *parent ,
                     int tab_number , QColor color ) :
    CGAL::Qt_widget( parent ),
    change_background_flag(FALSE),    
    current_state(0),
    index(tab_number),
    snap_mode(SNAP_NONE),
    mode(MODE_INSERT),
    m_line_width(2),
    m_vertex_width(3),
    first_time(true),
    active(false),
    traits_type(t),
    bbox(CGAL::Bbox_2(-10, -10, 10, 10)),
    wasrepainted(true),
    on_first(false),
    edge_color(color),
    change_edge_color(FALSE),
    vertex_color(color),
    change_vertex_color(FALSE),
    snap(false),
    grid(false),
    conic_type(SEGMENT),
    cube_size(1),
    ray_shooting_direction(true),
    remove_org_curve(true),
    read_from_file(false),
    first_time_merge(true),
    draw_vertex(true),
    fill_face_color(def_bg_color)
  {
    static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(2) <<
      CGAL::BackgroundColor (CGAL::BLACK);
    set_window(-10, 10, -10, 10);

    setMouseTracking(TRUE);
  }

  /*! Destructor */
  virtual ~Qt_widget_base_tab(){}

  /*! current_state - indecates when a tab state is changed */
  bool change_background_flag;

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

  /*! edge color */
  QColor edge_color;

  /*! flag to know that edge color has changed */
  bool change_edge_color;

  /*! vertices color */
  QColor vertex_color;
  
  /*! flag to know that vertices color has changed */  
  bool change_vertex_color;

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

  /*! true when it is the first time in merge mode */
  bool first_time_merge;

  /*! true when you want to draw all vertex, false if
   * only the intersecting vertex
   */
  bool draw_vertex;

  /*! the color for filling faces ( obtained by fill button) */
  QColor fill_face_color;

  /*! get the color of the unbounded face (its the same as the background
   * color of the tab)
   */
  QColor unbounded_face_color() { return this->backgroundColor(); }

  /*! increment current_state to inidicate that something has changed
   */
  void something_changed(){ current_state++ ; }

  virtual void change_strategy(Strategy /* s */){}

  virtual bool is_empty(){return true;}
  
  /*function that aids selecting a color that will be visible for drawing.
   *par_color_order is an array of (at least 4) different colors sorted
   *according to preferance. Default argument defined in demo_tab.cpp.
   *bad_color is an optional color that we wish to avoid (usually the color of the
   *face the point is in.
   *Since 3 comparisons are done for each color order_color must have more than 3
   *colors to guarrantee success*/
  void setCorrectColor(QColor bad_color=Qt::black, 
                       QColor *par_color_order = s_color_order,
                       unsigned int order_size = 4)
  {
		unsigned int i=0;
		while (i<order_size) {				
		if ((unbounded_face_color() != par_color_order[i]) && 
			(edge_color != par_color_order[i]) &&
			(bad_color != par_color_order[i])) 
			{
				setColor(par_color_order[i]);
				return;
			}
			++i;
		}
  		return;
  }  
};


/*! template class Qt_widget_demo_tab gets a Tab_traits class as
 *  a template parameter. all the Tab_traits classes must support
 *  a set of demands.
 */
template <class Tab_traits>
class Qt_widget_demo_tab : public Qt_widget_base_tab
{
private:
  typedef typename Tab_traits::Curves_list             Curves_list;
  typedef typename Tab_traits::Arrangement_2           Arrangement_2;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Halfedge_around_vertex_const_circulator
    Halfedge_around_vertex_const_circulator;
  typedef typename Tab_traits::Traits                  Traits;
  typedef typename Tab_traits::Curve_2                 Curve_2;
  typedef typename Tab_traits::X_monotone_curve_2      X_monotone_curve_2;
  typedef typename Tab_traits::Arr_curve_iter          Arr_curve_iter;
  typedef typename Tab_traits::Arr_curve_const_iter    Arr_curve_const_iter;
  typedef typename Tab_traits::Point_2                 Point_2;
  typedef typename Tab_traits::Halfedge_handle         Halfedge_handle;
  typedef typename Tab_traits::Face_handle             Face_handle;

  typedef typename Tab_traits::Face_const_handle       Face_const_handle;
  typedef typename Tab_traits::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Tab_traits::Vertex_const_handle     Vertex_const_handle;

  typedef typename Tab_traits::Ccb_halfedge_circulator
    Ccb_halfedge_circulator;
  typedef typename Tab_traits::Holes_iterator          Holes_iterator;
  typedef typename Arrangement_2::Hole_const_iterator  Holes_const_iterator;
  typedef typename Arrangement_2::Isolated_vertex_const_iterator
                                               Isolated_vertex_const_iterator;
  typedef typename Tab_traits::Halfedge_iterator       Halfedge_iterator;
  typedef typename Tab_traits::Hafledge_list           Hafledge_list;
  typedef typename Tab_traits::Hafledge_list_iterator
    Hafledge_list_iterator;
  typedef typename Tab_traits::Vertex_iterator Vertex_iterator;
  typedef typename Tab_traits::Halfedge_around_vertex_circulator
    Halfedge_around_vertex_circulator;
  typedef typename Tab_traits::Edge_iterator           Edge_iterator;
  typedef typename Tab_traits::Halfedge                Halfedge;
  typedef typename Tab_traits::Face_iterator           Face_iterator;
  typedef typename Tab_traits::Trap_point_location     Trap_point_location;
  typedef typename Tab_traits::Simple_point_location   Simple_point_location;
  typedef typename Tab_traits::Walk_point_location     Walk_point_location;
  typedef typename Tab_traits::Lanmarks_point_location Lanmarks_point_location;

  typedef My_observer<Arrangement_2>                    Observer;
  typedef typename Arrangement_2::Originating_curve_iterator
                                                    Originating_curve_iterator;
  typedef typename Arrangement_2::Induced_edge_iterator Induced_edge_iterator;
  typedef typename Arrangement_2::Curve_handle          Curve_handle;
  typedef CGAL::Envelope_diagram_1<Traits>              Diagram_1;
  
  typedef typename Traits::Multiplicity                 Multiplicity; 

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



public:
  /*! m_tab_traits - the traits object */
  Tab_traits       m_tab_traits;

  /*! m_curves_arr - pointer for the tab planar map */
  Arrangement_2*   m_curves_arr;

  Observer         m_observer;

  CGAL::Object     m_point_location;

  /*! Original Traits */
  Traits           m_traits;

  /*! the curve to be merged */
  Halfedge_iterator closest_curve;

  /*! the second curve to be merged */
  Halfedge_iterator second_curve;

  /*! the first point in the split curve */
  Point_2 split_point;

  /*! a removable halfedge iterators (created by move event when
   * DELETE option is on
   */
  Halfedge_iterator removable_halfedge;



  /*! constructor
   *\ param t - widget traits type
   *\ param parent - widget parent window
   *\ param tab_number - widget program index
   */
  Qt_widget_demo_tab(TraitsType t, QWidget * parent , int tab_number, QColor c):
    Qt_widget_base_tab(t , parent, tab_number, c),
    m_curves_arr (new Arrangement_2()),
    m_observer(*m_curves_arr),
    removable_halfedge()
  {
    // set the unbounded face initial color
    m_curves_arr->unbounded_face()->set_color(def_bg_color);
    m_point_location =
      CGAL::make_object(new Walk_point_location(*m_curves_arr));
  }

  /*! destructor - delete the planar map and the point location pointer
   */
  virtual ~Qt_widget_demo_tab()
  {
    m_observer.detach ();
    delete m_curves_arr;
  }


  /*! draw - called everytime something changed, draw the PM and mark the
   *         point location if the mode is on.
   */
  void draw()
  {
    QCursor old = cursor();
    setCursor(Qt::WaitCursor);

    if ( (mode == MODE_FILLFACE)  && (!change_background_flag) 
           && (!change_edge_color) && (!change_vertex_color) )
    {
      Point_2 temp_p (pl_point.x(), pl_point.y());
      CGAL::Object obj = locate(temp_p);
      
      Face_const_handle f;
      if (CGAL::assign (f, obj))
      {
        Face_handle ncf = m_curves_arr->non_const_handle(f);
        set_face_color(ncf, fill_face_color);
      }
    }
    
    // draw all faces (fill them with their color)
    visit_faces(FillFace(this));
    
    //reset change_background,edge, and vertex flags to FALSE
    //in order to allow future fill operations. 
	 change_background_flag=FALSE;
	 change_edge_color=FALSE;
	 change_vertex_color=FALSE;	
    if (snap_mode == SNAP_GRID || grid)
      draw_grid();

    for (Edge_iterator ei = m_curves_arr->edges_begin();
         ei != m_curves_arr->edges_end(); ++ei)
    {
      setColor(edge_color);
      m_tab_traits.draw_xcurve(this , ei->curve() );
    }
    // Go over all vertices and for each vertex check the
    // index numbers of the base curves that go through
    // it and paint the point if they are different (beacuse we want to
    // color red the intersection opints between two different planar maps
    // which are overlayed
    *this << CGAL::DISC;
    static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(m_vertex_width);

    Vertex_iterator   vit;
    for (vit = m_curves_arr->vertices_begin();
         vit != m_curves_arr->vertices_end(); vit++)
    {
        // draw all vertexes of the planar map if 'draw_vertex' is true
        // draw_vertex is a flag that indicates if we draw the vertexes

          setColor(vertex_color);
          Coord_point p(CGAL::to_double((*vit).point().x()) /
                        m_tab_traits.COORD_SCALE,
                        CGAL::to_double((*vit).point().y()) /
                        m_tab_traits.COORD_SCALE);
          static_cast<CGAL::Qt_widget&>(*this) << p;

    }
    if (mode == MODE_POINT_LOCATION)
    {
      static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(3);


      Point_2 temp_p (pl_point.x(), pl_point.y());
      CGAL::Object obj = locate(temp_p);

      Face_const_handle f = get_face(obj);
		
		/* more prudent color selection that selects the drawing color
		according to my_prefrance. replaced setColor(Qt::yellow)*/
		QColor my_preferance[4]= {Qt::yellow,Qt::green,Qt::red,Qt::blue};		
		setCorrectColor(f->color(),my_preferance, 4);		
		
      if (!f->is_unbounded()) // its an inside face
      {
        Ccb_halfedge_const_circulator cc = f->outer_ccb();
        do
        {
          m_tab_traits.draw_xcurve(this , cc->curve() );
        }
        while (++cc != f->outer_ccb());
      }


      //color the holes of the located face
      Holes_const_iterator hit, eit = f->holes_end();
      for (hit = f->holes_begin(); hit != eit; ++hit)
      {
        Ccb_halfedge_const_circulator cc = *hit;
        do
        {
          m_tab_traits.draw_xcurve(this , cc->curve() );
          cc++;
        }
        while (cc != *hit);
      }

      //color isolated vertices
      Isolated_vertex_const_iterator ivit = f->isolated_vertices_begin();
      for (; ivit != f->isolated_vertices_end(); ++ivit)
      {
        static_cast<CGAL::Qt_widget&>(*this) << ivit->point();
      }

      static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(m_line_width);
    }

    if (mode == MODE_RAY_SHOOTING_UP)
    {
		//relevant_face_color will keep the fill color of the face we pass through		
		QColor relevant_face_color=unbounded_face_color();      
      Coord_point up;
      Point_2 temp_p (pl_point.x(), pl_point.y());
      Coord_point pl_draw(pl_point.x() / m_tab_traits.COORD_SCALE ,
                          pl_point.y() / m_tab_traits.COORD_SCALE);
      CGAL::Object    obj = ray_shoot_up (temp_p);
      if (!obj.is_empty())
      {
        Face_const_handle ubf;
        if (CGAL::assign(ubf, obj))
        {
          CGAL_assertion(ubf->is_unbounded());
			 //relevant_face_color = unbounded_face_color() as initialized          
          up = Coord_point(pl_draw.x() , y_max());
          static_cast<CGAL::Qt_widget&>(*this) << Coord_segment(pl_draw, up);
        }
        // we shoot something
        else
        {
          Halfedge_const_handle he;
          if (CGAL::assign(he, obj))
          {
            Point_2 p1c1(pl_point.x() , y_max() * m_tab_traits.COORD_SCALE);
            Point_2 p2c1(pl_point.x() , pl_point.y());
            const X_monotone_curve_2 c1 =
              m_tab_traits.curve_make_x_monotone(p1c1 , p2c1);
            const X_monotone_curve_2 c2 = he->curve();

            CGAL::Object             res;
            CGAL::Oneset_iterator<CGAL::Object> oi(res);

            m_traits.intersect_2_object()(c1, c2, oi);
            std::pair<Point_2,Multiplicity> p1;
            if (CGAL::assign(p1, res))
            {
              Coord_type y1 =
                CGAL::to_double(p1.first.y())/ m_tab_traits.COORD_SCALE;
              up = Coord_point(pl_draw.x(), y1);
            }
            else
            {
              up = pl_draw;
            }
				relevant_face_color = he->face()->color();				
				/*choose color to mark the edge that differs from the current 
				edge_color, the background, and the relevant face color*/				
				setCorrectColor(relevant_face_color);				 
            m_tab_traits.draw_xcurve(this , he->curve() );
          }
          else
          {
            Vertex_const_handle v;
            CGAL_assertion(CGAL::assign(v, obj));
            CGAL::assign(v, obj);
            up = Coord_point(CGAL::to_double(v->point().x()) /
                        m_tab_traits.COORD_SCALE,
                            CGAL::to_double(v->point().y()) /
                        m_tab_traits.COORD_SCALE);
            
				//locate face that arrow will be drawn in, and retrieve its color 
      		CGAL::Object obj1 = locate(temp_p);
      		Face_const_handle f1 = get_face(obj1);
      	   relevant_face_color=f1->color();
            
            /*choose color to mark the vertice so that it differs from the 
            edge_color, the background, and the relevant_face_color*/				
				setCorrectColor(relevant_face_color);				     
            static_cast<CGAL::Qt_widget&>(*this) << up;
          }
        }

        //select arrow color that differs from the color of the face it is in
		  setCorrectColor(relevant_face_color);        
        
        static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(2);
        static_cast<CGAL::Qt_widget&>(*this) << Coord_segment(pl_draw,up);

        // draw an arrow that points to 'up' point
        int x = this->x_pixel(CGAL::to_double(up.x()));
        int y = this->y_pixel(CGAL::to_double(up.y()));

        this->get_painter().drawLine(x-7 , y+7 , x , y);
        this->get_painter().drawLine(x+7 , y+7 , x , y);
        static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(m_line_width);
      }
    }
    if (mode == MODE_RAY_SHOOTING_DOWN)
    {
		//relevant_face_color will keep the fill color of the face we pass through		
		QColor relevant_face_color=unbounded_face_color();      
      Coord_point up;
      Point_2 temp_p (pl_point.x(), pl_point.y());
      Coord_point pl_draw(pl_point.x() / m_tab_traits.COORD_SCALE ,
                          pl_point.y() / m_tab_traits.COORD_SCALE);
      CGAL::Object    obj = ray_shoot_down (temp_p);
      
      
      if (!obj.is_empty())
      {
        Coord_point down;
        Face_const_handle ubf;
        if (CGAL::assign(ubf, obj))
        {
			 //relevant_face_color = unbounded_face_color() as initialized           
          down = Coord_point(pl_draw.x() , y_min());
          static_cast<CGAL::Qt_widget&>(*this) << Coord_segment(pl_draw, down);
        }
        // we shoot something
        else
        {
          Halfedge_const_handle he;
          if (CGAL::assign(he, obj))
          {
            Point_2 p1c1(pl_point.x() , y_min() * m_tab_traits.COORD_SCALE);
            Point_2 p2c1(pl_point.x() , pl_point.y());
            const X_monotone_curve_2 c1 =
              m_tab_traits.curve_make_x_monotone(p1c1 , p2c1);
            const X_monotone_curve_2 c2 = he->curve();

            CGAL::Object             res;
            CGAL::Oneset_iterator<CGAL::Object> oi(res);

            m_traits.intersect_2_object()(c1, c2, oi);
            std::pair<Point_2,Multiplicity> p1;
            if (CGAL::assign(p1, res))
            {
              Coord_type y1 =
                CGAL::to_double(p1.first.y()) / m_tab_traits.COORD_SCALE;
              down = Coord_point(pl_draw.x(),y1);
            }
            else
            {
              down = pl_draw;
            }
				relevant_face_color = he->face()->color();				
				/*choose color to mark the edge that differs from the edge_color
				the background, and the relevant face color*/				
				setCorrectColor(relevant_face_color);
            m_tab_traits.draw_xcurve(this , he->curve() );
          }
          else
          {
            Vertex_const_handle v;
            CGAL_assertion(CGAL::assign(v, obj));
            CGAL::assign(v, obj);
            down = Coord_point(CGAL::to_double(v->point().x()) /
                          m_tab_traits.COORD_SCALE,
                          CGAL::to_double(v->point().y()) /
                          m_tab_traits.COORD_SCALE);
                          
				//locate face that arrow will be drawn in, and retrieve its color 
      		CGAL::Object obj1 = locate(temp_p);
      		Face_const_handle f1 = get_face(obj1);
      	   relevant_face_color=f1->color();
            
            /*choose color to mark the vertice so that it differs from the 
            edge_color, the background, and the relevant face color*/				
				setCorrectColor(relevant_face_color);				       
            static_cast<CGAL::Qt_widget&>(*this) << down;
          }
        }

        //select arrow color that differs from the color of the face it is in
		  setCorrectColor(relevant_face_color);   
        
        static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(2);
        static_cast<CGAL::Qt_widget&>(*this) << Coord_segment(pl_draw,down);
        // draw an arrow that points to 'down' point
        int x = this->x_pixel(CGAL::to_double(down.x()));
        int y = this->y_pixel(CGAL::to_double(down.y()));

        this->get_painter().drawLine(x-7 , y-7 , x , y);
        this->get_painter().drawLine(x+7 , y-7 , x , y);
        static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(m_line_width);
      }
    }

    if (lower_env || upper_env)
    {
      std::list<X_monotone_curve_2> xcurves;
      for (Edge_iterator ei = m_curves_arr->edges_begin();
           ei != m_curves_arr->edges_end();
           ++ei)
      {
        xcurves.push_back(ei->curve());
      }

      Diagram_1 diag;
      if (lower_env)
      {
        CGAL::lower_envelope_x_monotone_2(xcurves.begin(), xcurves.end(), diag);
        print_diag_1(diag);
      }
      diag.clear();
      if (upper_env)
      {
        CGAL::upper_envelope_x_monotone_2(xcurves.begin(), xcurves.end(), diag);
        print_diag_1(diag);
      }
    }
    setCursor(old);
  }

  void print_diag_1(const Diagram_1 & diag)
  {
    // Print the minimization diagram.
    typename Diagram_1::Edge_const_handle     e = diag.leftmost();
    typename Diagram_1::Vertex_const_handle   v;
    typename Diagram_1::Curve_const_iterator  cit;

    setCorrectColor();
    while (e != diag.rightmost())
    {
      if (! e->is_empty())
      {
        // The edge is not empty: draw a representative curve.
        // Note that the we only draw the portion of the curve
        // that overlaps the x-range defined by the two vertices
        // that are incident to this edge.
        m_tab_traits.draw_xcurve_segment(this , e->curve(),
                                         e->left()->point(),
                                         e->right()->point());
      }

      v = e->right();

      // Draw the point associated with the current vertex.
      Coord_point p(CGAL::to_double(v->point().x()) /
                    m_tab_traits.COORD_SCALE,
                    CGAL::to_double(v->point().y()) /
                    m_tab_traits.COORD_SCALE);
      static_cast<CGAL::Qt_widget&>(*this) << p;
      
      e = v->right();
    }
  }



  /*!
   */
  void set_face_color(Face_handle f ,QColor& c)
  {
    f->set_color(c);
    if ( f->is_unbounded())
      this->setBackgroundColor(c);;
  }

  /*!
   */
  template <class Function>
  void visit_faces(Function func)
  {
    Face_iterator  fi = m_curves_arr->faces_begin();
    for( ; fi != m_curves_arr->faces_end() ; ++fi )
      fi->set_visited(false);
    Face_handle ub = m_curves_arr->unbounded_face();    
    visit_face_rec (ub,func) ;
}


  /*! antenna - return true if the halfedge and its
   *  twin point to the same face.
   */
  bool antenna(Halfedge_handle h)
  {
    Halfedge_handle twin = h->twin();
    return (twin->face() == h->face());
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
          Halfedge_handle he = cc;
          Halfedge_handle he2 = he->twin();
          Face_handle inner_face = he2->face();
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
      if (! he.twin()->face()->visited())
      {
        Face_handle nei = (Face_handle) he.twin()->face();
        visit_ccb_faces( nei ,func );
      }
      //created from the outer boundary of the face
    } while (++cc != fh->outer_ccb());
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
    //int cube_size_x = (CGAL::max)(1, abs(max_x - min_x)/20);
    //int cube_size_y = (CGAL::max)(1, abs(max_y - min_y)/20);
    if (cube_size < std::abs(max_x - min_x)/40 ||
        cube_size < std::abs(max_y - min_y)/40)
      cube_size = (CGAL::max)((CGAL::max)(1, std::abs(max_x - min_x)/20),
                              (CGAL::max)(1, std::abs(max_y - min_y)/20));

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

    if (mode == MODE_POINT_LOCATION || mode == MODE_RAY_SHOOTING_UP ||
        mode == MODE_RAY_SHOOTING_DOWN || mode == MODE_FILLFACE)
    {
      mousePressEvent_point_location( e );
      setCursor(old);
      return;
    }

    if (mode == MODE_DELETE)
    {
      if (removable_halfedge == Halfedge_handle())
      {
        setCursor(old);
        return;
      }
      if (remove_org_curve)
      {
        Originating_curve_iterator  ocit, temp,
          ocit_end = m_curves_arr->originating_curves_end (removable_halfedge);
        Curve_handle                ch;
        ocit  = m_curves_arr->originating_curves_begin (removable_halfedge);
        while (ocit != ocit_end)
        {
          temp = ocit;
          ++temp;
          ch = ocit;
          CGAL::remove_curve(*m_curves_arr, ocit);
          ocit = temp;
        }
      }
      else
        m_curves_arr->remove_edge(removable_halfedge);

      removable_halfedge = Halfedge_handle();
      redraw();

      setCursor(old);
      return;
    }

    if (mode == MODE_INSERT)
    {
      Coord_type x, y;
      x_real(e->x(), x);
      y_real(e->y(), y);
      Coord_point p = point(x,y);

      lock();
      QColor old_color = color();
      RasterOp old_rasterop=rasterOp();
      get_painter().setRasterOp(XorROP);

      insert( e , p);

      setRasterOp(old_rasterop);
      setColor(old_color);
      unlock();

      setCursor(old);
      return;
    }
    if (mode == MODE_DRAG)
    {
      mousePressEvent_drag(e);
      setCursor(old);
      return;
    }
    if (mode == MODE_MERGE)
    {
      mousePressEvent_merge(e);
      setCursor(old);
      removable_halfedge = Halfedge_handle();
      return;
    }
    if (mode == MODE_SPLIT)
    {
      Coord_type x, y;
      x_real(e->x(), x);
      y_real(e->y(), y);
      Coord_point p = point(x,y);

      lock();
		//retain old color and rasterOp t preserve way things
		// are written to the paint device      
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
    if (e->button() == Qt::LeftButton && is_pure(e->state()))
    {
      if (!active)
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

  /*! split - split a xcurve into 2 xcurves. If several xcurves intersect
  	* with the inserted curve they are all split at the intersection point.
   *\ param e - mouse click event
   *\ param p - the pressed point
   */
  void split( QMouseEvent *e , Coord_point p)
  {
    if (e->button() == Qt::LeftButton && is_pure(e->state()))
    {
      if (!active)
      {
        active = true;
        m_tab_traits.first_point( p , mode);
        split_point = Point_2( p.x() * m_tab_traits.COORD_SCALE ,
                               p.y() * m_tab_traits.COORD_SCALE);
      }
      else
      {
        active = false;
        Point_2 split_point2 =
          Point_2(p.x() * m_tab_traits.COORD_SCALE,
                  p.y() * m_tab_traits.COORD_SCALE);
        const X_monotone_curve_2 split_curve =
          m_tab_traits.curve_make_x_monotone(split_point , split_point2);
        std::pair<Point_2,Multiplicity> p1;
        Point_2 p_right;
        if (split_point.x() < split_point2.x())
          p_right = split_point;
        else
          p_right = split_point2;
        Halfedge_iterator hei;
        for (hei = m_curves_arr->halfedges_begin();
             hei != m_curves_arr->halfedges_end(); ++hei)
        {
          const X_monotone_curve_2 & xcurve = hei->curve();
          //m_tab_traits.draw_xcurve(this, xcurve); removed, not necessary 
          CGAL::Object             res;
          CGAL::Oneset_iterator<CGAL::Object> oi(res);

          m_traits.intersect_2_object()(split_curve, xcurve, oi);

          if (CGAL::assign(p1, res)) {
          	if (hei == m_curves_arr->halfedges_end())
      			return; 
          	// we dont want to split an already existed vertex...
       		if (m_traits.equal_2_object()(hei->source()->point(), p1.first) ||
          		m_traits.equal_2_object()(hei->target()->point(), p1.first))
           		continue;
          	//m_tab_traits.draw_xcurve(this, hei->curve());
				//split the desired half edge at the intersection stored by p1        		
        		m_curves_arr->split_edge(hei , p1.first);
          }
        } //for loop
      }// else
    }
  }

  /* mousePressEvent_point_location - creats the point location point
   param e - mouse click event*/
  void mousePressEvent_point_location(QMouseEvent *e)
  {
    if (e->button() == Qt::LeftButton && is_pure(e->state()))
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
    if ((s & Qt::ControlButton) || (s & Qt::ShiftButton) || (s & Qt::AltButton))
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
    //int cube_size = (CGAL::max)(1, abs(my_max - my_min)/20);
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
    //  if the arrangement is empty do nothing
    if ( m_curves_arr->number_of_edges() == 0)
      return;

    // get the point of the mouse
    if (removable_halfedge != Halfedge_handle())
    {
      setColor(edge_color);
      if (remove_org_curve)
      {

        Originating_curve_iterator  ocit, temp;
        ocit  = m_curves_arr->originating_curves_begin (removable_halfedge);
        while (ocit !=
               m_curves_arr->originating_curves_end (removable_halfedge))
        {
          temp = ocit;
          ++temp;
          Curve_handle          ch = ocit;
          Induced_edge_iterator itr;
          for(itr = m_curves_arr->induced_edges_begin(ch);
              itr != m_curves_arr->induced_edges_end(ch);
              ++itr)
          {
            m_tab_traits.draw_xcurve(this,(*itr)->curve());
          }
          ocit = temp;
        }
      }
      else
        m_tab_traits.draw_xcurve(this, removable_halfedge->curve());
    }


    Coord_point p(x_real(e->x()) * m_tab_traits.COORD_SCALE ,
                  y_real(e->y()) * m_tab_traits.COORD_SCALE);

    bool is_first = true;
    Coord_type min_dist = 0;
    Halfedge_iterator hei;
    Halfedge_iterator closest_hei;

    for (hei = m_curves_arr->halfedges_begin();
         hei != m_curves_arr->halfedges_end();
         ++hei)
    {
      X_monotone_curve_2 & xcurve = hei->curve();
      Coord_type dist = m_tab_traits.xcurve_point_distance(p, xcurve , this);
      if (is_first || dist < min_dist)
      {
        min_dist = dist;
        closest_hei = hei;
        is_first = false;
      }
    }
    // now 'closest_hei' holds the cloeset halfedge to the point of the mouse

    removable_halfedge = closest_hei;
    if (remove_org_curve)
    {
      setColor(Qt::red);  // highlight the removable edge with red color

       Originating_curve_iterator  ocit, temp;
       ocit  = m_curves_arr->originating_curves_begin (removable_halfedge);
       while(ocit != m_curves_arr->originating_curves_end (removable_halfedge))
       {
         temp = ocit;
         ++temp;

         Curve_handle          ch = ocit;
         Induced_edge_iterator itr;
         for(itr = m_curves_arr->induced_edges_begin(ch);
             itr != m_curves_arr->induced_edges_end(ch);
             ++itr)
         {
           m_tab_traits.draw_xcurve(this,(*itr)->curve());
         }
        ocit = temp;
      }
    }
    else
    {
      setColor(Qt::red);  // highlight the removable edge with red color
      m_tab_traits.draw_xcurve(this,closest_hei->curve());
    }
  }

  /*! mouseMoveEvent - enable seeing the line to be drawn
   *\ param e - mouse click event
   */
  void mouseMoveEvent(QMouseEvent *e)
  {
    static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(m_line_width);
    if (mode == MODE_DELETE)
      // find removable edges , store them in the list
      find_removable_halfedges(e);
    //'removable_halfedges' and highlight them

    if (mode == MODE_DRAG)
    {
      mouseMoveEvent_drag(e);
      return;
    }
    if (mode == MODE_MERGE && !first_time_merge)
    {//after closest_edge was selected, highlight second curve according to 
     //the mouse movement
      if (second_curve != m_curves_arr->halfedges_end())
      {//case a second curve exists recolor it to edge_color before searching
      //for a new second_curve using the new mouse position  
        setColor(edge_color);
        m_tab_traits.draw_xcurve(this,second_curve->curve());
      }
      Coord_type x, y;
      x_real(e->x(), x);
      y_real(e->y(), y);
      Coord_point p(x * m_tab_traits.COORD_SCALE,
                    y * m_tab_traits.COORD_SCALE);
      second_curve = m_curves_arr->halfedges_end();
      //search for a new second_curve 
      find_close_curve(closest_curve, second_curve, p, true);
		//color the halfedges that are about to be merged      
      setColor(Qt::red);
      m_tab_traits.draw_xcurve(this,closest_curve->curve());
      if (second_curve != m_curves_arr->halfedges_end())
      {
        setColor(Qt::green);
        m_tab_traits.draw_xcurve(this,second_curve->curve());
      }
      else
      { //did not find mergable half edges 
        first_time_merge = true;
        redraw();
      }
      return;
    }// merge

    if (active) //case for split action 
    {
      Coord_type x, y;
      x_real(e->x(), x);
      y_real(e->y(), y);
      Coord_point p = point(x,y);
      RasterOp old_raster = rasterOp();//save the initial raster mode
      setRasterOp(XorROP);
      lock();

      setColor(Qt::green);

      if (!first_time)
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
  void leaveEvent(QEvent * /* e */)
  {
    if (active)
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

  /*! point
   *\ params x,y - the mouse clicked point coordinates
   *\    return a point according to the current snap mode and
   *  recent points.
   */
  Coord_point point(Coord_type x, Coord_type y)
  {
    int xmin = static_cast<int> (x_min());
    int xmax = static_cast<int> (x_max());
    int ymin = static_cast<int> (y_min());
    int ymax = static_cast<int> (y_max());
    Coord_type d = (CGAL::max)(0.5 , (x_max() - x_min())/40);
    switch ( snap_mode ) {
     case SNAP_POINT:
      {
       Coord_type min_dist = 0;
       Coord_point closest;

       if ( m_curves_arr->number_of_vertices() == 0 )
         return Coord_point(x , y);

       min_dist = m_tab_traits.closest_point(x,y,closest,this);

       if (min_dist <= d)
         return closest;
       else
         return Coord_point(x , y);

       break;
      }
     case SNAP_GRID:
      return Coord_point(getMid(x, xmin, xmax),
                         getMid(y, ymin, ymax) );

     case SNAP_NONE: break;
    }
    return Coord_point(x,y);
  }

  /*! mousePressEvent_drag - change the Cursor on the drag mode
   *  mouse pressed event
   *\ param e - mouse click event
   */
  void mousePressEvent_drag(QMouseEvent *e)
  {
    if (e->button() == Qt::LeftButton && is_pure(e->state()))
    {
      setCursor(QCursor( QPixmap( (const char**)holddown_xpm)));
      if (!on_first) {
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
    if (e->button() == Qt::LeftButton
       && mode == MODE_DRAG
       && is_pure(e->state()))
    {
      setCursor(QCursor( QPixmap( (const char**)hand_xpm)));
      double x, y, xfirst2, yfirst2;
      x_real(e->x(), x);
      y_real(e->y(), y);
      x_real(first_x, xfirst2);
      y_real(first_y, yfirst2);

      // double    xmin, xmax, ymin, ymax;
      // if (x < xfirst2) { xmin = x; xmax = xfirst2; }
      // else { xmin = xfirst2; xmax = x; }
      // if (y < yfirst2) { ymin = y; ymax = yfirst2; }
      // else { ymin = yfirst2; ymax = y; }
      double distx = xfirst2 - x;
      double disty = yfirst2 - y;
      move_center(distx, disty);
      on_first = FALSE;
    }
  }

  /*! mouseMoveEvent_drag - calculate new widget position
   *\ param e - mouse release event
   */
  void mouseMoveEvent_drag(QMouseEvent *e)
  {
    if (on_first)
    {
      int x = e->x();
      int y = e->y();
      //save the last coordinates to redraw the screen
      x2 = x;
      y2 = y;
      wasrepainted = FALSE;
    }
  }

  /*! mousePressEvent_merge -  merge mode
   *  mouse pressed event
   *\ param e - mouse click event
   */
  void mousePressEvent_merge(QMouseEvent *e)
  {
    if (e->button() == Qt::LeftButton && is_pure(e->state()))
    {//merge only in case of a left click 
      if ( m_curves_arr->is_empty() )
        return;
		
      setColor(Qt::red);
      Coord_point p(x_real(e->x()) * m_tab_traits.COORD_SCALE ,
                    y_real(e->y()) * m_tab_traits.COORD_SCALE);
      bool       first = true;
      Coord_type min_dist = 0;

      if (first_time_merge)
      {//find the closest mergable half edge to point p 
        first_time_merge = false;
        Halfedge_iterator hei;
        closest_curve = m_curves_arr->halfedges_end();

        for (hei = m_curves_arr->halfedges_begin();
             hei != m_curves_arr->halfedges_end(); ++hei)
        {//find  closest curve to mouse pointer 
          Vertex_iterator   vis = hei->source();
          Vertex_iterator   vit = hei->target();
			 //case the halfedge can't be merged - next iteration         
          if (vis->degree() != 2 && vit->degree() != 2)
            continue;
          X_monotone_curve_2 & xcurve = hei->curve();
          Coord_type dist =
            m_tab_traits.xcurve_point_distance(p, xcurve, this);

          if (first || dist < min_dist)
          {
            min_dist = dist;
            closest_curve = hei;
            first = false;
          }
        }
        if (first)       // we didn't find any "good" curve
        {
          first_time_merge = true;
          return;
        }
        //draw the first half edge to merge with the setColor() chosen above         
        m_tab_traits.draw_xcurve(this , closest_curve->curve() );
        second_curve = m_curves_arr->halfedges_end();
      }   
      else //not first_time_merge
      {
        first_time_merge = true;
        //look for the second halfedge closest to p that is mergable with
        //closest_curve and merge them 
        find_close_curve(closest_curve, second_curve, p, false);
        redraw();
      }
    } else { //not left click event (right click) undo all selections
    		first_time_merge=TRUE;
	 		//repaint all curves to edge_color. 
	 		redraw();    
    	} 
  }

  CGAL::Object locate(const Point_2& pt)
  {
    Walk_point_location* walk_pl;
    if (CGAL::assign(walk_pl, m_point_location))
      return walk_pl->locate(pt);

    Simple_point_location* simple_pl;
    if (CGAL::assign(simple_pl, m_point_location))
      return simple_pl->locate(pt);

    Trap_point_location* trap_pl;
    if (CGAL::assign(trap_pl, m_point_location))
      return trap_pl->locate(pt);

    Lanmarks_point_location* lm_pl;
    if (CGAL::assign(lm_pl, m_point_location))
      return lm_pl->locate(pt);

    // doesnt suppose to reach there
    CGAL_error();
    return CGAL::Object();
  }

  CGAL::Object ray_shoot_up(const Point_2& pt)
  {
    Walk_point_location* walk_pl;
    if (CGAL::assign(walk_pl, m_point_location))
      return walk_pl->ray_shoot_up(pt);

    Simple_point_location* simple_pl;
    if (CGAL::assign(simple_pl, m_point_location))
      return simple_pl->ray_shoot_up(pt);

    Trap_point_location* trap_pl;
    if (CGAL::assign(trap_pl, m_point_location))
      return trap_pl->ray_shoot_up(pt);

    Lanmarks_point_location* lm_pl;
    if (CGAL::assign(lm_pl, m_point_location))
    {
      // QMessageBox::information( this, "Ray shoot down", "Land Marks doesn't
      // support ray shooting");
      return CGAL::Object();
    }

    // doesnt suppose to reach there
    CGAL_error();
    return CGAL::Object();
  }

  CGAL::Object ray_shoot_down(const Point_2& pt)
  {
    Walk_point_location* walk_pl;
    if (CGAL::assign(walk_pl, m_point_location))
      return walk_pl->ray_shoot_down(pt);

    Simple_point_location* simple_pl;
    if (CGAL::assign(simple_pl, m_point_location))
      return simple_pl->ray_shoot_down(pt);

    Trap_point_location* trap_pl;
    if (CGAL::assign(trap_pl, m_point_location))
      return trap_pl->ray_shoot_down(pt);

    Lanmarks_point_location* lm_pl;
    if (CGAL::assign(lm_pl, m_point_location))
    {
      //QMessageBox::information( this, "Ray shoot up", "Land Marks doesn't
      // support ray shooting");
      return CGAL::Object();
    }

    // doesnt suppose to reach there
    CGAL_error();
    return CGAL::Object();
  }

  virtual void change_strategy(Strategy s)
  {
    Walk_point_location* walk_pl = NULL;
    if (CGAL::assign(walk_pl, m_point_location)) delete walk_pl;
    else {
      Simple_point_location* simple_pl = NULL;
      if (CGAL::assign(simple_pl, m_point_location)) delete simple_pl;
      else {
        Trap_point_location* trap_pl = NULL;
        if (CGAL::assign(trap_pl, m_point_location)) delete trap_pl;
        else {
          Lanmarks_point_location* lm_pl = NULL;
          if (CGAL::assign(lm_pl, m_point_location)) delete lm_pl;
        }
      }
    }

    if (s == WALK)
    {
      m_point_location =
        CGAL::make_object(new Walk_point_location(*m_curves_arr));
      return;
    }
    if (s == SIMPLE)
    {
      m_point_location =
        CGAL::make_object(new Simple_point_location(*m_curves_arr));
      return;
    }
    if (s == TRAP)
    {
      QCursor old = cursor();
      setCursor(Qt::WaitCursor);
      m_point_location =
        CGAL::make_object(new Trap_point_location(*m_curves_arr));
      setCursor(old);
      return;
    }
    if (s == LANDMARKS)
    {
     QCursor old = cursor();
      setCursor(Qt::WaitCursor);
      m_point_location =
        CGAL::make_object(new Lanmarks_point_location(*m_curves_arr));
        setCursor(old);
      return;
    }
  }

  virtual bool is_empty()
  {
    return m_curves_arr->is_empty();
  }

  /*Function that is invoked by move or click mouse event functions related to merge. 
   It checks if the parameter closest_curve (first halfedge to merge) is mergable
   with another halfedge (to be stored in second_curve).Second_curve is updated to
   store the mergable halfedge closest to mouse point p.
   If this function is not triggered by a move event the closest_edge and second_edge  
	are merged.	  
   */
  void find_close_curve(Halfedge_iterator &closest_curve,
                        Halfedge_iterator &second_curve,
                        Coord_point &p,
                        bool move_event)
  {
	 //boolean var - if "good" curves were found changed to false   
    bool       first = true;
    Coord_type min_dist = 0;

    for (Halfedge_iterator hei = m_curves_arr->halfedges_begin();
         hei != m_curves_arr->halfedges_end();
         ++hei)
    {
      if (m_curves_arr->are_mergeable(closest_curve, hei))
      {
        X_monotone_curve_2 & xcurve = hei->curve();
        Coord_type dist = m_tab_traits.xcurve_point_distance(p, xcurve , this);
        if (first || dist < min_dist)
        {
          min_dist = dist;
          second_curve = hei;
          first = false;
        }
      }
    }
    if (first)     // didn't find any "good" curve
      return;

    if (!move_event)
    {
      m_curves_arr->merge_edge( closest_curve, second_curve);
    }
  }


  private:

    Face_const_handle get_face(const CGAL::Object& obj)
    {
      Face_const_handle f;
      if (CGAL::assign(f, obj))
        return f;

      Halfedge_const_handle he;
      if (CGAL::assign(he, obj))
        return (he->face());

      Vertex_const_handle v;
      CGAL_assertion(CGAL::assign(v, obj));
      CGAL::assign(v, obj);
      if (v->is_isolated())
        return v->face();
      Halfedge_around_vertex_const_circulator eit = v->incident_halfedges();
      return  (eit->face());
    }

};

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/*! class Segment_tab_traits defines the segment traits
 */
class Segment_tab_traits
{
public:
  typedef Kernel::FT                                FT;
  typedef Arr_seg_list                              Curves_list;
  typedef Seg_arr                                   Arrangement_2;
  typedef Arrangement_2::Face_const_handle          Face_const_handle;
  typedef Arrangement_2::Halfedge_const_handle      Halfedge_const_handle;
  typedef Arrangement_2::Vertex_const_handle        Vertex_const_handle;

  typedef Arrangement_2::Face_handle                Face_handle;
  typedef Arrangement_2::Halfedge_handle            Halfedge_handle;
  typedef Arrangement_2::Vertex_handle              Vertex_handle;
  typedef Seg_traits                                Traits;
  typedef Traits::Point_2                           Point_2;
  typedef Traits::Curve_2                           Curve_2;
  typedef Traits::X_monotone_curve_2                X_monotone_curve_2;

  typedef Arr_seg_iter                              Arr_curve_iter;
  typedef Arr_seg_const_iter                        Arr_curve_const_iter;

  typedef Seg_ccb_halfedge_circulator               Ccb_halfedge_circulator;
  typedef Seg_holes_iterator                        Holes_iterator;
  typedef Arrangement_2::Halfedge_iterator          Halfedge_iterator;
  typedef std::list<Halfedge_iterator>              Hafledge_list;
  typedef Hafledge_list::iterator                   Hafledge_list_iterator;
  typedef Arrangement_2::Vertex_iterator            Vertex_iterator;
  typedef Arrangement_2::Halfedge_around_vertex_circulator
                                          Halfedge_around_vertex_circulator;
  typedef Arrangement_2::Edge_iterator              Edge_iterator;
  typedef Seg_halfedge                              Halfedge;
  typedef Seg_face_iterator                         Face_iterator;

  //point location
  typedef Seg_trap_point_location                   Trap_point_location;
  typedef Seg_simple_point_location                 Simple_point_location;
  typedef Seg_walk_point_location                   Walk_point_location;
  typedef Seg_lanmarks_point_location               Lanmarks_point_location;


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
                                             const X_monotone_curve_2 & cv)
  {
    return (he->source()->point() == cv.source());
  }

  /*! fill_face - fill a face with its color (which is stored at the face)
   * it creates a polyong from the outer boundary of the face and
   * uses CGAL opertaor << of polygons
   */
  void fill_face(Qt_widget_demo_tab<Segment_tab_traits> * w , Face_handle f)
  {
    if (!f->is_unbounded())  // f is not the unbounded face
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
      My_polygon pgn (pts.begin() , pts.end());

      w->setFilled(true);

      // fill the face according to its color (stored at any of her
      // incidents curves)
      if (! f->color().isValid())
        w->setFillColor(def_bg_color);
      else
        w->setFillColor(f->color());

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
      (*w)<<My_polygon(points , points +4 );
      w->setFilled(false);
       w->get_painter().setPen(old_penstyle);
    }
  }


  /*! draw_xcurve - use Qt_Widget operator to draw
   *\ param w - the demo widget
   *\ c - xcurve to be drawen
   */
  void draw_xcurve(Qt_widget_demo_tab<Segment_tab_traits> * w ,
                   X_monotone_curve_2 c )
  {
    (*w) << c;
  }

  /*! Use Qt_Widget operator to draw a portion of an x-monotone curve.
   * \param w The demo widget.
   * \param xcurve The curve to be drawn.
   * \param p_left Defines the left end.
   * \param p_right Defines the right end.
   */
  void draw_xcurve_segment(Qt_widget_demo_tab<Segment_tab_traits> * w,
                           const X_monotone_curve_2 & c,
                           const Point_2 & p_left, const Point_2 & p_right)
  {
    if (m_traits.is_vertical_2_object() (c)) {
      (*w) << c;
      return;
    }

    // Trim the segment, if necessary.
    const Point_2 & p_min = m_traits.construct_min_vertex_2_object()(c);
    const Point_2 & p_max = m_traits.construct_max_vertex_2_object()(c);
    Kernel          ker;
    Kernel::Line_2  l = ker.construct_line_2_object() (p_min, p_max);

    const Point_2 & p1 =
      (ker.compare_x_2_object() (p_left, p_min) == CGAL::LARGER) ?
      Point_2(p_left.x(), ker.compute_y_at_x_2_object()(l, p_left.x())) :
      p_min;

    const Point_2 & p2 =
      (ker.compare_x_2_object() (p_right, p_max) == CGAL::SMALLER) ?
      Point_2(p_right.x(), ker.compute_y_at_x_2_object()(l, p_right.x())) :
      p_max;

    (*w) << ker.construct_segment_2_object() (p1, p2);
  }

  /*! Draw a curve to a Qt widget
   * \param w the demo widget
   * \param c curve to be drawen
   */
  void draw_curve(Qt_widget_demo_tab<Segment_tab_traits> * w , Curve_2 c )
  {
    (*w) << c;
  }

  /*! first_point - a first point of inserted sgment
   */
  void first_point( Coord_point p , Mode  )
  {
    m_p1 = m_p2 = p;
  }

  /*! Obtain the last point of a segment
   */
  void middle_point(Coord_point p,
                    Qt_widget_demo_tab<Segment_tab_traits> * w)
  {
    Coord_kernel      ker;

    if (! ker.equal_2_object() (m_p1, p))
    {
      get_segment( Coord_segment( m_p1 , p ) , w );
      w->active = false;
      //w->redraw();  // not working so I use new_object insted
      w->new_object(make_object(Coord_segment(m_p1 , p)));
    }
  }

  /*! last_point - meaningless for segments
   */
  void last_point( Coord_point  ,
                   Qt_widget_demo_tab<Segment_tab_traits> *  )
  {
    return;
  }

  /*! get_segment - create a new segment, insert him into curves_list
   * and planar map
   */
  void get_segment( Coord_segment coord_seg ,
                    Qt_widget_demo_tab<Segment_tab_traits> * w)
  {
    const Coord_point & coord_source = coord_seg.source();
    const Coord_point & coord_target = coord_seg.target();
    Arr_seg_point_2 source(coord_source.x(), coord_source.y());
    Arr_seg_point_2 target(coord_target.x(), coord_target.y());
    Arr_seg_2 seg (source, target);
    CGAL::insert(*(w->m_curves_arr), seg);
    CGAL::Bbox_2 curve_bbox = seg.bbox();
    w->bbox = w->bbox + curve_bbox;
  }


  /*! xcurve_point_distance - return the distance between a point
   * and a xsegment
   */
  Coord_type xcurve_point_distance(Coord_point p, X_monotone_curve_2 & c ,
                                   Qt_widget_demo_tab<Segment_tab_traits> * )
  {
    const Arr_seg_point_2 & source = c.source();
    const Arr_seg_point_2 & target = c.target();

    Coord_type x1 = CGAL::to_double(source.x());
    Coord_type y1 = CGAL::to_double(source.y());

    Coord_type x2 = CGAL::to_double(target.x());
    Coord_type y2 = CGAL::to_double(target.y());

    Coord_point coord_source(x1 , y1);
    Coord_point coord_target(x2 , y2);
    Coord_segment coord_seg(coord_source, coord_target);
    return CGAL::squared_distance( p, coord_seg);
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
      const Point_2& p = (*vit).point();
      x1 = CGAL::to_double(p.x());
      y1 = CGAL::to_double(p.y());
      dt = w->dist(x1 , y1 , x , y);
      if (first || dt < min_dist)
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
  const X_monotone_curve_2 curve_make_x_monotone(Point_2 p1 , Point_2 p2)
  {
    const Curve_2 cv(p1 , p2);
    X_monotone_curve_2 c;
    CGAL::Object             res;
    CGAL::Oneset_iterator<CGAL::Object> oi(res);
    m_traits.make_x_monotone_2_object()(cv, oi);
    CGAL::assign(c, res);
    return c;
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
  typedef Kernel::FT                                FT;
  typedef Arr_pol_list                              Curves_list;
  typedef Pol_arr                                   Arrangement_2;
  typedef Arrangement_2::Face_const_handle          Face_const_handle;
  typedef Arrangement_2::Halfedge_const_handle      Halfedge_const_handle;
  typedef Arrangement_2::Vertex_const_handle        Vertex_const_handle;
  typedef Pol_traits                                Traits;
  typedef Arr_pol_iter                              Arr_curve_iter;
  typedef Arr_pol_const_iter                        Arr_curve_const_iter;
  typedef Arr_pol_point_2                           Point_2;
  typedef Pol_halfedge_handle                       Halfedge_handle;
  typedef Arr_pol_2                                 Curve_2;
  typedef Arr_xpol_2                                X_monotone_curve_2;
  typedef Pol_face_handle                           Face_handle;
  typedef Pol_ccb_halfedge_circulator               Ccb_halfedge_circulator;
  typedef Pol_holes_iterator                        Holes_iterator;
  typedef Arrangement_2::Halfedge_iterator          Halfedge_iterator;
  typedef std::list<Halfedge_iterator>              Hafledge_list;
  typedef Hafledge_list::iterator                   Hafledge_list_iterator;
  typedef std::vector<Point_2>::iterator            Point_vector_iterator;
  typedef Curve_2::const_iterator                   Curve_const_iterator;
  typedef Arrangement_2::Vertex_iterator            Vertex_iterator;
  typedef Arrangement_2::Halfedge_around_vertex_circulator
    Halfedge_around_vertex_circulator;
  typedef Arrangement_2::Edge_iterator              Edge_iterator;
  typedef Pol_halfedge                              Halfedge;
  typedef  Pol_face_iterator                        Face_iterator;

  //point location
  typedef Pol_trap_point_location                   Trap_point_location;
  typedef Pol_simple_point_location                 Simple_point_location;
  typedef Pol_walk_point_location                   Walk_point_location;
  typedef Pol_lanmarks_point_location               Lanmarks_point_location;


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
                                            const X_monotone_curve_2 & cv)
  {
    return (he->source()->point() == *(cv.begin()));
  }

  /*! fill_face - fill a face with its color (which is stored at the curves)
   * it creates a polyong from the outer boundary of the face and
   * uses CGAL opertaor << of polygons
   */
  void fill_face(Qt_widget_demo_tab<Polyline_tab_traits> * w , Face_handle f)
  {
    if (!f->is_unbounded())  // f is not the unbounded face
    {
      std::list< Coord_point > pts; // holds the points of the polygon
      X_monotone_curve_2::const_iterator           pt_itr;
      X_monotone_curve_2::const_reverse_iterator   pt_rev_itr;
      X_monotone_curve_2 cv;

      /* running with around the outer of the face and generate from it
       * polygon
       */
      Ccb_halfedge_circulator cc=f->outer_ccb();
      do {
        cv = cc->curve();
        if ( curve_has_same_direction (cc) )
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
          for (pt_rev_itr = cv.rbegin() , ++pt_rev_itr; pt_rev_itr != cv.rend();
               ++pt_rev_itr)
          {
            Coord_type x = CGAL::to_double((*pt_rev_itr).x());
            Coord_type y = CGAL::to_double((*pt_rev_itr).y());
            Coord_point coord_source(x , y);
            pts.push_back(coord_source );
          }
        }
        //created from the outer boundary of the face
      } while (++cc != f->outer_ccb());

      // make polygon from the outer ccb of the face 'f'
      My_polygon pgn (pts.begin() , pts.end());

      w->setFilled(true);

      // fill the face according to its color (stored at any of her
      // incidents curves)
      if (! f->color().isValid())
        w->setFillColor(def_bg_color);
      else
        w->setFillColor(f->color());
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
      (*w)<<My_polygon(points , points +4 );
      w->setFilled(false);
       w->get_painter().setPen(old_penstyle);
    }
  }

  /*! draw_xcurve - go over the polyline parts and use Qt_Widget operator
   * to draw
   */
  void draw_xcurve(Qt_widget_demo_tab<Polyline_tab_traits> * w,
                   X_monotone_curve_2 pol )
  {
    Curve_2::const_iterator ps = pol.begin();
    Curve_2::const_iterator pt = ps; ++pt;

    while (pt != pol.end()) {
      const Point_2 & source = *ps;
      const Point_2 & target = *pt;
      Coord_segment coord_seg = convert(source , target);
      *w << coord_seg;
      ++ps; ++pt;
    }
  }

  /*! Use Qt_Widget operator to draw a portion of an x-monotone polyline.
   * \param w The demo widget.
   * \param xcurve The curve to be drawn.
   * \param p_left Defines the left end.
   * \param p_right Defines the right end.
   */
  void draw_xcurve_segment(Qt_widget_demo_tab<Polyline_tab_traits> * w,
                           const X_monotone_curve_2 & pol,
                           const Point_2 & p_left, const Point_2 & p_right)
  {
    if (m_traits.is_vertical_2_object() (pol)) {
      Curve_2::const_iterator pi = pol.begin();
      const Point_2 & source = *pi++;
      const Point_2 & target = *pi;
      Coord_segment coord_seg = convert(source , target);
      (*w) << coord_seg;
      return;
    }

    Curve_2::const_iterator ps = pol.begin();
    Curve_2::const_iterator pt = ps; ++pt;
    Kernel                  ker;
    Kernel::Compare_x_2     comp_x = ker.compare_x_2_object();
    Point_2                 src, trg;

    while (pt != pol.end()) {
      // Skip this segment if it is not in the relevant x-range.
      if (comp_x (p_left, *pt) == CGAL::LARGER)
        continue;

      if (comp_x (p_right, *ps) == CGAL::SMALLER)
        break;

      // Trim the current segment, if necessary.
      Kernel::Line_2      l = ker.construct_line_2_object() (*ps, *pt);

      src = (comp_x (p_left, *ps) == CGAL::LARGER) ?
        Point_2 (p_left.x(), ker.compute_y_at_x_2_object() (l, p_left.x())) :
        *ps;

      trg = (comp_x (p_right, *pt) == CGAL::SMALLER) ?
        Point_2 (p_right.x(), ker.compute_y_at_x_2_object() (l, p_right.x())) :
        *pt;

      Coord_segment coord_seg = convert (src, trg);
      (*w) << coord_seg;
      ++ps; ++pt;
    }
  }


  /*! Draw a curve
   */
  void draw_curve(Qt_widget_demo_tab<Polyline_tab_traits> * w , Curve_2 pol )
  {
    std::list<CGAL::Object> obj_list;
    m_traits.make_x_monotone_2_object()(pol, std::back_inserter(obj_list));
    for(std::list<CGAL::Object>::iterator itr = obj_list.begin();
        itr != obj_list.end();
        ++itr)
    {
      X_monotone_curve_2 cv;

      if (CGAL::assign(cv, *itr))
        draw_xcurve (w, cv);
    }
  }

  /*! first_point - a first point of inserted polyline or a splitter
   */
  void first_point( Coord_point p , Mode m)
  {
    last_of_poly = p;
    if (m == MODE_INSERT)
      points.push_back(Arr_pol_point_2(p.x(),p.y()));
  }

  /*! middle_point - a middle point of a polyline
   */
  void middle_point( Coord_point p ,
                     Qt_widget_demo_tab<Polyline_tab_traits> * w)
  {
    if (last_of_poly == p) return;
    rubber_old = p;

    points.push_back(Arr_pol_point_2(p.x(),p.y()));

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


  /*! xcurve_point_distance - return the distance between a point
   * and a polyline
   */
  Coord_type xcurve_point_distance(Coord_point p, X_monotone_curve_2 & c,
                                   Qt_widget_demo_tab<Polyline_tab_traits> * )
  {
    Curve_const_iterator ps = c.begin();
    Curve_const_iterator pt = ps; pt++;
    bool first = true;
    Coord_type min_dist = 0;

    while (pt != c.end())
    {
      const Point_2 & source = *ps;
      const Point_2 & target = *pt;
      Coord_segment coord_seg = convert(source , target);
      Coord_type dist = CGAL::squared_distance( p, coord_seg);

      if (first || dist < min_dist)
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
      const X_monotone_curve_2& curve = heit->curve();
      Curve_const_iterator cit;
      for (cit = curve.begin(); cit != curve.end(); cit++)
      {
        const Point_2& p = *cit;
        x1 = CGAL::to_double(p.x());
        y1 = CGAL::to_double(p.y());
        dt = w->dist(x1 , y1 , x , y);
        if (first || dt < min_dist)
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
      const Arr_pol_point_2& p = *it;
      x1 = CGAL::to_double(p.x());
      y1 = CGAL::to_double(p.y());
      dt = w->dist(x1 , y1 , x , y);
      if (first || dt < min_dist)
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
  const X_monotone_curve_2 curve_make_x_monotone(Point_2 p1 , Point_2 p2)
  {
    std::vector<Arr_pol_point_2> temp_points;
    temp_points.push_back(p1);
    temp_points.push_back(p2);
    Curve_2 cv(temp_points.begin(), temp_points.end());
    CGAL::Object             res;
    CGAL::Oneset_iterator<CGAL::Object> oi(res);
    m_traits.make_x_monotone_2_object()(cv, oi);
    X_monotone_curve_2 c1 ;
    CGAL::assign(c1, res);
    return c1;
  }


private:

  /*! get_polyline - create a new polyline
   */
  void get_polyline(Qt_widget_demo_tab<Polyline_tab_traits> * w)
  {
    Arr_pol_2  pol (points.begin(), points.end());
    CGAL::insert(*(w->m_curves_arr), pol);
    CGAL::Bbox_2 curve_bbox = pol.bbox();
    w->bbox = w->bbox + curve_bbox;
  }

  /*! convert - convert from Arr_pol_curve to Coord_segment
   */
  Coord_segment convert(Arr_pol_point_2 & source , Arr_pol_point_2 & target)
  {
    Coord_type x1 = CGAL::to_double(source.x());
    Coord_type y1 = CGAL::to_double(source.y());
    Coord_point coord_source(x1, y1);

    Coord_type x2 = CGAL::to_double(target.x());
    Coord_type y2 = CGAL::to_double(target.y());
    Coord_point coord_target(x2, y2);

    return Coord_segment(coord_source, coord_target);
  }

  /*! convert - convert from const Arr_pol_curve to Coord_segment
   */
  Coord_segment convert(const Arr_pol_point_2 & source ,
                        const Arr_pol_point_2 & target)
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
  std::vector<Arr_pol_point_2> points;
};


//////////////////////////////////////////////////////////////////////////////
/*!
 */
#ifdef CGAL_USE_CORE
class Conic_tab_traits
{
public:
  typedef Alg_kernel::FT                            FT;
  typedef Arr_xconic_list                           Curves_list;
  typedef Conic_arr                                 Arrangement_2;
  typedef Arrangement_2::Face_const_handle          Face_const_handle;
  typedef Arrangement_2::Halfedge_const_handle      Halfedge_const_handle;
  typedef Arrangement_2::Vertex_const_handle        Vertex_const_handle;
  typedef Conic_traits                              Traits;
  typedef Arr_xconic_iter                           Arr_curve_iter;
  typedef Arr_xconic_const_iter                     Arr_curve_const_iter;
  typedef Arr_conic_point_2                         Point_2;
  typedef Conic_halfedge_handle                     Halfedge_handle;
  typedef Arrangement_2::Curve_2                    Curve_2;
  typedef Arr_xconic_2                              X_monotone_curve_2;
  typedef Conic_face_handle                         Face_handle;
  typedef Conic_ccb_halfedge_circulator             Ccb_halfedge_circulator;
  typedef Conic_holes_iterator                      Holes_iterator;
  typedef Arrangement_2::Halfedge_iterator          Halfedge_iterator;
  typedef std::list<Halfedge_iterator>              Hafledge_list;
  typedef Hafledge_list::iterator                   Hafledge_list_iterator;
  typedef Arrangement_2::Vertex_iterator            Vertex_iterator;
  typedef Arrangement_2::Halfedge_around_vertex_circulator
    Halfedge_around_vertex_circulator;
  typedef Arrangement_2::Edge_iterator              Edge_iterator;
  typedef Conic_halfedge                            Halfedge;
  typedef Conic_face_iterator                       Face_iterator;

  //point location
  typedef Conic_trap_point_location                 Trap_point_location;
  typedef Conic_simple_point_location               Simple_point_location;
  typedef Conic_walk_point_location                 Walk_point_location;
  typedef Conic_lanmarks_point_location             Lanmarks_point_location;

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
                                            const X_monotone_curve_2 & cv)
  {
    return (he->source()->point() == cv.source());
  }

  /*!
   */
  void fill_face(Qt_widget_demo_tab<Conic_tab_traits> * w , Face_handle f)
  {
    if (! f->is_unbounded())  // f is not the unbounded face
    {
      std::list< Coord_point > pts; // holds the points of the polygon
      /* running with around the outer of the face and generate from it
       * polygon
       */
      Ccb_halfedge_circulator cc=f->outer_ccb();
      do {
        if (w->antenna(cc))
          continue;

        Halfedge_handle he = cc;
        X_monotone_curve_2 c = he->curve();
        // Get the co-ordinates of the curve's source and target.
        double sx = CGAL::to_double(he->source()->point().x()),
               sy = CGAL::to_double(he->source()->point().y()),
               tx = CGAL::to_double(he->target()->point().x()),
               ty = CGAL::to_double(he->target()->point().y());

        Coord_point coord_source(sx / COORD_SCALE, sy / COORD_SCALE);
        Coord_point coord_target(tx / COORD_SCALE, ty / COORD_SCALE);

        if (c.orientation() == CGAL::COLLINEAR)
            pts.push_back(coord_source );
        else
        {
          // If the curve is monotone, than its source and its target has the
          // extreme x co-ordinates on this curve.
          bool is_source_left = (sx < tx);
          int  x_min = is_source_left ? (*w).x_pixel(sx) : (*w).x_pixel(tx);
          int  x_max = is_source_left ? (*w).x_pixel(tx) : (*w).x_pixel(sx);
          double curr_x, curr_y;
          int  x;

          Arr_conic_point_2 px;

          pts.push_back(coord_source );

          if (is_source_left) {
            for (x = x_min + DRAW_FACTOR; x < x_max; x+=DRAW_FACTOR) {
              //= COORD_SCALE)
              curr_x = (*w).x_real(x);
              Alg_kernel   ker;
              Arr_conic_point_2 curr_p(curr_x, 0);
              if (!(ker.compare_x_2_object()(curr_p, c.left()) !=
                      CGAL::SMALLER &&
                    ker.compare_x_2_object()(curr_p, c.right()) !=
                      CGAL::LARGER))
                continue;
              px = c.point_at_x (curr_p);
              curr_y = CGAL::to_double(px.y());
              pts.push_back(Coord_point(curr_x / COORD_SCALE,
                                        curr_y / COORD_SCALE));
            }// for
          }
          else {
            for (x = x_max; x > x_min; x-=DRAW_FACTOR) {
              curr_x = (*w).x_real(x);
              Alg_kernel   ker;
              Arr_conic_point_2 curr_p(curr_x, 0);
              if (!(ker.compare_x_2_object() (curr_p, c.left()) !=
                      CGAL::SMALLER &&
                    ker.compare_x_2_object() (curr_p, c.right()) !=
                      CGAL::LARGER))
                continue;
              px = c.point_at_x (Arr_conic_point_2(curr_x, 0));
              curr_y = CGAL::to_double(px.y());
              pts.push_back(Coord_point(curr_x / COORD_SCALE,
                                        curr_y / COORD_SCALE));
            }// for
          }// else
          pts.push_back(coord_target );
        }
        //created from the outer boundary of the face
      } while (++cc != f->outer_ccb());

      // make polygon from the outer ccb of the face 'f'
      My_polygon pgn (pts.begin() , pts.end());
      QPen old_penstyle = w->get_painter().pen();
      w->get_painter().setPen(Qt::NoPen);
      w->setFilled(true);

      // fill the face according to its color (stored at any of her incidents
      // curves)
      if (! f->color().isValid())
        w->setFillColor(def_bg_color);
      else
        w->setFillColor(f->color());

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
      (*w)<<My_polygon(points , points +4 );
      w->setFilled(false);
       w->get_painter().setPen(old_penstyle);
    }
  }

  /*! draw_xcurve - same as draw_curve
   */
  void draw_xcurve(Qt_widget_demo_tab<Conic_tab_traits> * w,
                   X_monotone_curve_2 c )
  {
    // Get a polyline approximation of the curve.
    const Point_2&  p_min = m_traits.construct_min_vertex_2_object() (c);
    const Point_2&  p_max = m_traits.construct_max_vertex_2_object() (c);

    if (c.orientation() == CGAL::COLLINEAR)
    {
      Coord_point s(CGAL::to_double(p_min.x()), CGAL::to_double(p_min.y()));
      Coord_point t(CGAL::to_double(p_max.x()), CGAL::to_double(p_max.y()));

      Coord_segment seg(s, t);
      *w << seg;
      return;
    }
    const double    x_min = CGAL::to_double (p_min.x());
    const double    x_max = CGAL::to_double (p_max.x());
    const int       ix_min = (*w).x_pixel(x_min);
    const int       ix_max = (*w).x_pixel(x_max);
    unsigned int    n = static_cast<unsigned int> (ix_max - ix_min);

    if (w->x_min() > x_max || w->x_max() < x_min)
      return;

    if (n == 0)
      return;

    CGAL::Bbox_2    c_bbox = c.bbox();

    if (w->y_min() > c_bbox.ymax() || w->y_max() < c_bbox.ymin())
      return;

    std::pair<double, double>  *app_pts = new std::pair<double, double> [n + 1];
    std::pair<double, double>  *end_pts = c.polyline_approximation (n, app_pts);
    std::pair<double, double>  *p_curr = app_pts;
    std::pair<double, double>  *p_next = p_curr + 1;
    Coord_point     ps (p_curr->first, p_curr->second);

    p_curr = app_pts;
    p_next = p_curr + 1;
    do
    {
      Coord_point     pt (p_next->first, p_next->second);

      *w << Coord_segment(ps, pt);
      ps = pt;
      p_curr++;
      p_next++;
    } while (p_next != end_pts);

    delete[] app_pts;
    return;
  }

  /*! Use Qt_Widget operator to draw a portion of an x-monotone conic arc.
   * \param w The demo widget.
   * \param c The curve to be drawn.
   * \param p_left Defines the left end.
   * \param p_right Defines the right end.
   */
  void draw_xcurve_segment(Qt_widget_demo_tab<Conic_tab_traits> * w,
                           const X_monotone_curve_2 & c,
                           const Point_2 & p_left, const Point_2 & p_right)
  {
    // Get a polyline approximation of the curve.
    const Point_2 & p_min = m_traits.construct_min_vertex_2_object() (c);
    const Point_2 & p_max = m_traits.construct_max_vertex_2_object() (c);
    Alg_kernel      ker;

    // Trim the curve, if necessary.
    const Point_2 & p1 =
      (ker.compare_x_2_object() (p_left, p_min) == CGAL::LARGER) ?
      p_left : p_min;

    const Point_2 & p2 =
      (ker.compare_x_2_object() (p_right, p_max) == CGAL::SMALLER) ?
      p_right : p_max;

    const double    x_min = CGAL::to_double (p1.x());
    const double    x_max = CGAL::to_double (p2.x());
    const int       ix_min = (*w).x_pixel(x_min);
    const int       ix_max = (*w).x_pixel(x_max);
    unsigned int    n = static_cast<unsigned int> (ix_max - ix_min);

    if (w->x_min() > x_max || w->x_max() < x_min)
      return;

    if(n == 0)
      return;

    CGAL::Bbox_2    c_bbox = c.bbox();

    if (w->y_min() > c_bbox.ymax() || w->y_max() < c_bbox.ymin())
      return;

    std::pair<double, double>  *app_pts = new std::pair<double, double> [n + 1];
    std::pair<double, double>  *end_pts = c.polyline_approximation (n, app_pts);
    std::pair<double, double>  *p_curr = app_pts;
    std::pair<double, double>  *p_next = p_curr + 1;
    Coord_point     ps (p_curr->first, p_curr->second);

    p_curr = app_pts;
    p_next = p_curr + 1;
    do {
      Coord_point     pt (p_next->first, p_next->second);

      *w << Coord_segment(ps, pt);
      ps = pt;
      p_curr++;
      p_next++;
    } while (p_next != end_pts);

    delete[] app_pts;
    return;
  }

  /*! Draw_a curve
   */
  void draw_curve(Qt_widget_demo_tab<Conic_tab_traits> * w , Curve_2 conic )
  {
    std::list<CGAL::Object> obj_list;
    m_traits.make_x_monotone_2_object()(conic, std::back_inserter(obj_list));
    for(std::list<CGAL::Object>::iterator itr = obj_list.begin();
        itr != obj_list.end();
        ++itr)
    {
      X_monotone_curve_2 cv;

      if (CGAL::assign(cv, *itr))
        draw_xcurve (w, cv);
    }
  }

  ////////////////////////////////////////////////////////////////////////////

  /*! first_point - a first point of inserted sgment
   */
  void first_point( Coord_point p , Mode )
  {
    m_p_old = m_p1 = m_p2 = m_p3 = m_p4 = p;
    num_points = 1;
    first_time = true;
  }

  /*! middle_point - the last point of a segment
   */
  void middle_point( Coord_point p , Qt_widget_demo_tab<Conic_tab_traits> * w)
  {
    Rational r, s, t, u, v, ww;  // The conic coefficients.
    Rational a, b, c, a_sq, b_sq;
    Rational x, y, x1, y1, x0, y0, temp;
    Rational sq_rad;

    x1 = Rational(static_cast<int>(1000 * m_p1.x() + 0.5), 1000);
    y1 = Rational(static_cast<int>(1000 * m_p1.y() + 0.5), 1000);
    x = Rational(static_cast<int>(1000 * p.x() + 0.5), 1000);
    y = Rational(static_cast<int>(1000 * p.y() + 0.5), 1000);

    if (x != x1 || y != y1)
    {
      Arr_conic_2 cv;

      switch (w->conic_type)
      {
       case CIRCLE:
        sq_rad = CGAL::square(x - x1) + CGAL::square(y - y1);
        cv =  Arr_conic_2(Rat_circle_2 (Rat_point_2(x1, y1), sq_rad));
        break;

       case SEGMENT:
        cv =  Arr_conic_2(Rat_segment_2 (Rat_point_2(x,y),
                                                Rat_point_2(x1,y1)));
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

        cv =  Arr_conic_2(r, s, t, u, v, ww);
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
          Rational x2 = Rational(static_cast<int>(1000 * m_p2.x() + 0.5), 1000);
          Rational y2 = Rational(static_cast<int>(1000 * m_p2.y() + 0.5), 1000);
          Rat_kernel ker;
          // the three points of the parabola cannot be collinear (
          if (ker.collinear_2_object()(Rat_point_2(x1,y1),
                                       Rat_point_2(x,y),Rat_point_2(x2,y2)))
          {
            QMessageBox::information( w, "Insert Conic", "Invalid Conic");
            w->active = false;
            *w << m_p1 << m_p2 ;
            return;
          }
          cv =  Arr_conic_2 (Rat_point_2(x1,y1),Rat_point_2(x2,y2),
                                    Rat_point_2(x,y));
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
          Rational x2 = Rational(static_cast<int>(1000 * m_p2.x() + 0.5), 1000);
          Rational y2 = Rational(static_cast<int>(1000 * m_p2.y() + 0.5), 1000);
          Rational x3 = Rational(static_cast<int>(1000 * m_p3.x() + 0.5), 1000);
          Rational y3 = Rational(static_cast<int>(1000 * m_p3.y() + 0.5), 1000);
          Rational x4 = Rational(static_cast<int>(1000 * m_p4.x() + 0.5), 1000);
          Rational y4 = Rational(static_cast<int>(1000 * m_p4.y() + 0.5), 1000);
          cv =  Arr_conic_2 (Rat_point_2(x1,y1),Rat_point_2(x2,y2),
                             Rat_point_2(x3,y3),Rat_point_2(x4,y4),
                             Rat_point_2(x,y));
          if (! cv.is_valid())
          {
            QMessageBox::information( w, "Insert Conic", "Invalid Conic");
            w->active = false;
            *w << m_p1 << m_p2 << m_p3 << m_p4 << p;
            return;
          }
        }
        break;
      }

      CGAL::insert(*(w->m_curves_arr), cv);
      CGAL::Bbox_2 curve_bbox = cv.bbox();
      w->bbox = w->bbox + curve_bbox;
      w->active = false;
      w->new_object(make_object(Coord_segment(m_p1 , p)));
    }
  }

  /*! last_point - meaningless for conics - at least for now
   */
  void last_point( Coord_point  , Qt_widget_demo_tab<Conic_tab_traits> *  )
  {
    return;
  }

  /*! draw_last_segment - call from mouse move event
   */
  void draw_last_segment( Qt_widget_demo_tab<Conic_tab_traits> * w)
  {
    if (w->mode == MODE_SPLIT)
      *w << Coord_segment( m_p1 , m_p_old );
    else
    {
      switch (w->conic_type)
      {
       case CIRCLE:
        *w << Coord_circle(m_p1, pow(m_p1.x() - m_p_old.x(), 2) +
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
    Vertex_iterator vit;
    for (vit = w->m_curves_arr->vertices_begin();
         vit != w->m_curves_arr->vertices_end(); vit++)
    {
      const Point_2& p = (*vit).point();
      x1 = CGAL::to_double(p.x()) / COORD_SCALE;
      y1 = CGAL::to_double(p.y()) / COORD_SCALE;
      dt = w->dist(x1 , y1 , x , y);
      if (first || dt < min_dist)
      {
        min_dist = dt;
        closest = Coord_point(x1 , y1);
        first = false;
      }
    }
    return min_dist;
  }


  /*! xcurve_point_distance - return the distance between
   * a point and a conic
   */
  Coord_type xcurve_point_distance(Coord_point p, X_monotone_curve_2 & c ,
                                   Qt_widget_demo_tab<Conic_tab_traits> * w)
  {
    // Get the co-ordinates of the curve's source and target.
    double sx = CGAL::to_double(c.source().x()),
      sy = CGAL::to_double(c.source().y()),
      tx = CGAL::to_double(c.target().x()),
      ty = CGAL::to_double(c.target().y());

    if (c.orientation() == CGAL::COLLINEAR)
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
        bool is_source_left = (sx < tx);
        int  x_min = is_source_left ? (*w).x_pixel(sx) : (*w).x_pixel(tx);
        int  x_max = is_source_left ? (*w).x_pixel(tx) : (*w).x_pixel(sx);
        double   prev_x = is_source_left ? sx : tx;
        double   prev_y = is_source_left ? sy : ty;
        double   curr_x, curr_y;
        int      x;

        Arr_conic_point_2 px;

        bool         first = true;
        Coord_type   min_dist = 100000000;
        Alg_kernel   ker;


        for (x = x_min + 1; x < x_max; x++)
        {
          curr_x = (*w).x_real(x);
          Arr_conic_point_2 curr_p(curr_x, 0);
          if (!(ker.compare_x_2_object() (curr_p, c.left()) != CGAL::SMALLER &&
                ker.compare_x_2_object() (curr_p, c.right()) != CGAL::LARGER))
            continue;

          px = c.point_at_x(Arr_conic_point_2(curr_x, 0));
          curr_y = CGAL::to_double(px.y());

          Coord_segment coord_seg( Coord_point(prev_x, prev_y) ,
                                   Coord_point(curr_x, curr_y) );
          Coord_type dist = CGAL::squared_distance( p, coord_seg);
          if (first || dist < min_dist)
          {
            first = false;
            min_dist = dist;
          }
          prev_x = curr_x;
          prev_y = curr_y;
        }
        return min_dist;
    }

  }

  /*! curve_make_x_monotone
   */
  const X_monotone_curve_2 curve_make_x_monotone(Point_2 p1 , Point_2 p2)
  {
    Rational x1 (static_cast<int> (1000 * CGAL::to_double(p1.x())+ 0.5), 1000);
    Rational y1 (static_cast<int> (1000 * CGAL::to_double(p1.y())+ 0.5), 1000);

    Rational x2 (static_cast<int> (1000 * CGAL::to_double(p2.x())+ 0.5), 1000);
    Rational y2 (static_cast<int> (1000 * CGAL::to_double(p2.y())+ 0.5), 1000);

    Rat_point_2 my_p1 (x1, y1);
    Rat_point_2 my_p2 (x2, y2);
    Arr_conic_2 cv (Rat_segment_2 (my_p1, my_p2));

    CGAL::Object             res;
    CGAL::Oneset_iterator<CGAL::Object> oi(res);
    m_traits.make_x_monotone_2_object()(cv, oi);
    X_monotone_curve_2 c1;
    CGAL::assign(c1, res);
    return c1;
  }


  Traits m_traits;
  /*! temporary points of the created conic */
  Coord_point m_p_old,m_p1,m_p2,m_p3,m_p4;
  /*! bool flag for hyperbola insertion */
  bool first_time;
  /*! counter for the number of points */
  int num_points;
};
#endif

typedef Qt_widget_demo_tab<Segment_tab_traits>       Qt_widget_segment_tab;
typedef Qt_widget_demo_tab<Polyline_tab_traits>      Qt_widget_polyline_tab;
#ifdef CGAL_USE_CORE
typedef Qt_widget_demo_tab<Conic_tab_traits>         Qt_widget_conic_tab;
#endif

#endif //DEMO_TAB_H
