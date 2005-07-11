// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef OVERLAY_VISITOR_H
#define OVERLAY_VISITOR_H

#include <CGAL/Sweep_line_2/Arr_construction_visitor.h>
#include <CGAL/Arrangement_2/Arr_accessor.h>
#include <CGAL/Unique_hash_map.h> 

CGAL_BEGIN_NAMESPACE



template < class Traits_,
           class Arrangement1_,
           class Arrangement2_,
           class Arrangement_,
           class Event_,
           class Subcurve_,
           class OverlayTraits >
class Overlay_visitor : 
  public Arr_construction_visitor<Traits_, Arrangement_, Event_, Subcurve_>
{
public:

  typedef Traits_                                        Traits; 
  typedef typename Traits::X_monotone_curve_2            X_monotone_curve_2;
  typedef typename Traits::Point_2                       Point_2;
  typedef typename Traits::Curve_info                    Curve_info;

  typedef Arrangement_                                   Arrangement;
  typedef typename Arrangement::Halfedge_handle          Halfedge_handle;
  typedef typename Arrangement::Face_handle              Face_handle;
  typedef typename Arrangement::Vertex_handle            Vertex_handle;

  typedef Arrangement1_                                  Arrangement1;
  typedef typename Arrangement1::Halfedge_const_handle   Halfedge_handle_red;
  typedef typename Arrangement1::Face_const_handle       Face_handle_red;
  typedef typename Arrangement1::Vertex_const_handle     Vertex_handle_red;

  typedef Arrangement2_                                  Arrangement2;
  typedef typename Arrangement2::Halfedge_const_handle   Halfedge_handle_blue;
  typedef typename Arrangement2::Face_const_handle       Face_handle_blue;
  typedef typename Arrangement2::Vertex_const_handle     Vertex_handle_blue;

 
  typedef typename Arrangement::Ccb_halfedge_circulator 
    Ccb_halfedge_circulator;

  typedef Event_                                         Event;
  typedef Subcurve_                                      Subcurve;
  
  typedef Overlay_visitor< Traits,
                           Arrangement1,
                           Arrangement2,
                           Arrangement,
                           Event,
                           Subcurve,
                           OverlayTraits >               Self;
  
  typedef Arr_construction_visitor< Traits,
                                  Arrangement,
                                  Event,
                                  Subcurve >             Base;
  
  typedef typename Base::SubCurveIter                    SubCurveIter;
  typedef typename Base::SubCurveRevIter                 SubCurveRevIter;
  typedef typename Base::SL_iterator                     SL_iterator;

  typedef Unique_hash_map<Halfedge_handle,Curve_info>    Hash_map;
  using Base::m_arr;
  using Base::m_arr_access;


  private:

  //hide default c'tor and assignment operator
  Overlay_visitor (const Self& );
  Self& operator= (const Self& );

public:

  /*! Constructor */
  Overlay_visitor(const Arrangement1& red_arr,
                  const Arrangement2& blue_arr,
                  Arrangement& res_arr,
                  OverlayTraits& overlay_trairs):
    Base(&res_arr),
    m_red_arr(&red_arr),
    m_blue_arr(&blue_arr),
    m_halfedges_map(Curve_info(),
                    // give an initial size for the hash table
                      red_arr.number_of_halfedges() +
                      blue_arr.number_of_halfedges()),
    m_overlay_traits(&overlay_trairs)
  {}


  /*! Destructor */
  virtual ~Overlay_visitor() {}


  bool after_handle_event(Event * event, SL_iterator iter, bool flag)
  {
    bool res = Base::after_handle_event(event, iter, flag);

    SubCurveRevIter rev_iter = event->right_curves_rbegin();
    Subcurve *sc_above = NULL;
    if(iter != this->status_line_end())
      sc_above = (*iter);

    
    if(!sc_above)
    { 
      if(rev_iter != event->right_curves_rend())
      {
        (*rev_iter)->set_above(NULL);
        sc_above = *rev_iter;
        ++rev_iter;     
      }
      else
        return res; // nothing else to do 
    }

    for( ;
         rev_iter != event->right_curves_rend();
         ++rev_iter )
    {
      Subcurve* curr_sc = *rev_iter;

      if(!curr_sc->has_same_color(sc_above))
        curr_sc -> set_above(sc_above);
      else
        curr_sc -> set_above(sc_above->get_above());

      sc_above = curr_sc;
    }
    return res;
  }
 

  void add_subcurve(const X_monotone_curve_2& cv,Subcurve* sc)
  {
    Base::add_subcurve(cv, sc);
  }

    


  void update_event(Event* e,
                    const Point_2& end_point,
                    const X_monotone_curve_2& cv,
                    bool is_left_end)
  {
    Point_2& pt = e->get_point();
    if(pt.is_red_object_null())
    {
      pt.set_red_object(end_point.get_red_object());
    }
    else
      if(pt.is_blue_object_null())
      {
        pt.set_blue_object(end_point.get_blue_object());
      }
  }


  void update_event(Event* e,
                    Subcurve* c1,
                    Subcurve* c2,
                    bool created = false)
  {
    CGAL_assertion(created == true);
  }


  void update_event(Event *e,
                    Subcurve* sc)
  {
    Point_2& pt = e->get_point();

    if(pt.is_red_object_null())
    {
      CGAL_assertion(!pt.is_blue_object_null());
      CGAL_assertion(sc->get_color() == Curve_info::RED);
      Halfedge_handle_red red_he = sc->get_red_halfedge_handle();
      pt.set_red_object(CGAL::make_object(red_he));
    }
    else
      if(pt.is_blue_object_null())
      {
        Halfedge_handle_blue blue_he = sc->get_blue_halfedge_handle();
        pt.set_blue_object(CGAL::make_object(blue_he));
    }
  }


      
                  

  virtual Halfedge_handle insert_in_face_interior
    (const X_monotone_curve_2& cv,
     Subcurve* sc)
  {
    // res is directed from left to right
    Halfedge_handle res = Base::insert_in_face_interior(cv,sc);
    map_halfedge_and_twin(res, true, cv.get_curve_info());
    Subcurve *sc_above = sc->get_above();
    Vertex_handle res_v_left = res->source();
    Vertex_handle res_v_right = res->target();

    //create left vertex
    Event *last_event = 
      reinterpret_cast<Event*>(sc->get_last_event());
    create_vertex(last_event, res_v_left, sc_above);

    //create right vertex
    create_vertex(this ->current_event(), res_v_right, sc_above);

     //update the result edge
    this ->create_edge(sc, res->twin());
    return res;
  }



  virtual Halfedge_handle insert_from_right_vertex
                          (const X_monotone_curve_2& cv,
                           Halfedge_handle he,
                           Subcurve* sc)
  {
    // res is directed from right to left
    Halfedge_handle res = Base::insert_from_right_vertex(cv, he, sc);
    map_halfedge_and_twin(res, false, cv.get_curve_info());
    Subcurve *sc_above = sc->get_above();

    // the new vertex is the left one
    Vertex_handle res_v = res->target();

    Event *last_event = 
      reinterpret_cast<Event*>(sc->get_last_event());
    create_vertex(last_event, res_v, sc_above);

     //update the result edge
    this ->create_edge(sc, res);
    return res;
  }


  virtual Halfedge_handle insert_from_left_vertex
                          (const X_monotone_curve_2& cv,
                           Halfedge_handle he,
                           Subcurve* sc)
  {
    //res is directed from left to right
    Halfedge_handle res = Base::insert_from_left_vertex(cv, he, sc);
    map_halfedge_and_twin(res, true, cv.get_curve_info());
    Subcurve *sc_above = sc->get_above();

     // the new vertex is the right one
    Vertex_handle res_v = res->target();
    create_vertex(this ->current_event(), res_v, sc_above);

    //update the result edge
    this ->create_edge(sc, res->twin());
    return res;
  }


  virtual Halfedge_handle insert_at_vertices(const X_monotone_curve_2& cv,
                                             Halfedge_handle hhandle,
                                             Halfedge_handle prev,
                                             Subcurve* sc,
                                             bool &new_face_created)
  {
    // res is directed from right to left
    Halfedge_handle res = Base::insert_at_vertices(cv, hhandle, prev, sc,
                                                   new_face_created);
    //TODO: add an assertion that res is directed from right to left

    map_halfedge_and_twin(res, false, cv.get_curve_info());

    // new face was created, need to update the new face's data
    if(new_face_created)
    {
      Halfedge_handle_red  red_he  ;
      Halfedge_handle_blue blue_he ;

      // get the new face
      Face_handle new_face = res->face();

      //traverse new_face's boundary
      Ccb_halfedge_circulator ccb_end = new_face->outer_ccb();
      Ccb_halfedge_circulator ccb_circ = ccb_end;
      do
      { 
        //get the current halfedge on the face boundary
        Halfedge_handle he =  ccb_circ->handle();

        CGAL_assertion(m_halfedges_map.is_defined(he));
        const Curve_info& cv_info = m_halfedges_map[he];
        if(cv_info.get_color() == Curve_info::RED)
        {
          red_he = cv_info.get_red_halfedge_handle();
          if(blue_he != Halfedge_handle_blue())
            break;
        }
        else
          if(cv_info.get_color() == Curve_info::BLUE)
          {
            blue_he = cv_info.get_blue_halfedge_handle();
            if(red_he != Halfedge_handle_red())
              break;
          }
          else 
          {
            // overlap
            CGAL_assertion(cv_info.get_color() == Curve_info::PURPLE);
            red_he  = cv_info.get_red_halfedge_handle()  ;
            blue_he = cv_info.get_blue_halfedge_handle() ;
            break;
          }
        ++ccb_circ;
      }
      while(ccb_circ != ccb_end);
      //finished traversing the boundary of the face

      if(red_he != Halfedge_handle_red() && blue_he != Halfedge_handle_blue())
      {
        // red face and blue face intersects (or overlap)
        Face_handle_red red_face = red_he->face();
        Face_handle_blue blue_face = blue_he->face();
        m_overlay_traits->create_face(red_face, blue_face,new_face);
      }
      else
      {
        // red face inside blue face
        if(red_he != Halfedge_handle_red())
        {
          Face_handle_red red_face = red_he->face();
          Face_handle_blue blue_face;
          Subcurve* sc_above = sc->get_above();
          if(!sc_above)
            blue_face = m_blue_arr->unbounded_face();
          else
            blue_face = 
              sc_above->get_blue_halfedge_handle()->face();
          
          m_overlay_traits->create_face(red_face, blue_face,new_face);
        }
        else
        {
          // blue face inside red face
          CGAL_assertion(blue_he != Halfedge_handle_blue() && 
                         red_he == Halfedge_handle_red());
          Face_handle_red red_face;
          Face_handle_blue blue_face = blue_he->face();
          Subcurve* sc_above = sc->get_above();
          if(!sc_above)
            red_face = m_red_arr->unbounded_face();
          else
            red_face = 
              sc_above->get_red_halfedge_handle()->face();
          
          m_overlay_traits->create_face(red_face, blue_face,new_face);
        }
      }
    }

    //update the result edge
    this ->create_edge(sc, res);
    return res;
  }



  virtual Vertex_handle insert_isolated_vertex(const Point_2& pt,
                                               SL_iterator iter)
  {
    Vertex_handle v = Base::insert_isolated_vertex(pt, iter);
    Object red = pt.get_red_object();
    Object blue = pt.get_blue_object();
    Vertex_handle_red     red_v;
    Vertex_handle_blue    blue_v;
    assign(red_v, red);
    assign(blue_v, blue);

    if(!red.is_empty() && !blue.is_empty())
    {
      m_overlay_traits->create_vertex(red_v, blue_v, v);
      return v;
    }
    
    CGAL_assertion(!red.is_empty() || !blue.is_empty());

    Subcurve* sc_above; 
    if(red.is_empty())
    {
      // isolated blue vertex inside red face
      Face_handle_red red_f ;
      if( iter == this ->status_line_end())
      {
        m_overlay_traits->create_vertex(m_red_arr->unbounded_face(),blue_v,v);
        return v;
      }
      sc_above = *iter;
      if(! sc_above)
        red_f = m_red_arr->unbounded_face();
      else
      {
        if(sc_above->get_color() == Curve_info::RED)
          red_f = sc_above->get_red_halfedge_handle()->face();
        else
        {
          sc_above = sc_above->get_above();
          if(!sc_above)
            red_f = m_red_arr->unbounded_face();
          else
            red_f = sc_above->get_red_halfedge_handle()->face();
        }
      }
      m_overlay_traits->create_vertex(red_f, blue_v, v);
      return v;
    }
    
    CGAL_assertion(blue.is_empty());
    // isolated red vertex inside blue face

    Face_handle_blue    blue_f;
    if( iter == this ->status_line_end())
    {
      m_overlay_traits->create_vertex(red_v,m_blue_arr->unbounded_face(),v);
      return v;
    }
    sc_above = *iter;
    if(! sc_above)
      blue_f = m_blue_arr->unbounded_face();
    else
    {
      if(sc_above->get_color() == Curve_info::BLUE)
        blue_f = sc_above->get_blue_halfedge_handle()->face();
      else
      {
        sc_above = sc_above->get_above();
          if(!sc_above)
          blue_f = m_blue_arr->unbounded_face();
        else
          blue_f = sc_above->get_blue_halfedge_handle()->face();
      }
    }
    m_overlay_traits->create_vertex(red_v, blue_f, v);
    return v;
  }
   
      

  // maps halfedge and his twin, right_dir is true iff he is directed from
  // left to right
  void map_halfedge_and_twin(Halfedge_handle he, bool right_dir, 
                             const Curve_info& cv_info)
  {
    // original halfedges that were stored is directed from right to left
    Halfedge_handle_red     red_he_info = cv_info.get_red_halfedge_handle();
    Halfedge_handle_blue    blue_he_info = cv_info.get_blue_halfedge_handle();
    Halfedge_handle_red     red_he_info_twin;
    Halfedge_handle_blue    blue_he_info_twin; 


    if(red_he_info  != Halfedge_handle_red())
      red_he_info_twin = red_he_info->twin();

    if(blue_he_info != Halfedge_handle_blue())
      blue_he_info_twin = blue_he_info->twin();


    //cv_info_twin will have the twin halfedge
    Curve_info cv_info_twin(red_he_info_twin, blue_he_info_twin);

    if(right_dir)
    {
      m_halfedges_map[he] = cv_info_twin;
      m_halfedges_map[he->twin()] = cv_info;
    }
    else
    {
      m_halfedges_map[he] = cv_info;
      m_halfedges_map[he->twin()] = cv_info_twin;
    }
  }


 


  void create_vertex(Event *event, Vertex_handle res_v, Subcurve* sc_above)
  {
    const Point_2& pt = event->get_point();
    CGAL_assertion( !pt.is_red_object_null() || !pt.is_blue_object_null());
    const Object& red_obj  = pt.get_red_object();
    const Object& blue_obj = pt.get_blue_object();
    Vertex_handle_red red_v;
    if(assign(red_v, red_obj))
    {
      Vertex_handle_blue    blue_v;
      if(assign(blue_v, blue_obj))
      {
        // red vertex on blue vertex
        m_overlay_traits ->create_vertex(red_v, blue_v, res_v);
      }
      else
      {
        Halfedge_handle_blue    blue_he;
        if(assign(blue_he, blue_obj))
        {
          //red vertex on blue halfedge
          m_overlay_traits->create_vertex(red_v, blue_he, res_v);
        }
        else
        {
          // red vertex inside blue face
          CGAL_assertion(blue_obj.is_empty());
          Face_handle_blue    blue_f;
          if(!sc_above)
            blue_f = m_blue_arr->unbounded_face();
          else
          {
            blue_f = sc_above ->get_blue_halfedge_handle()->face();
          }
          m_overlay_traits->create_vertex(red_v, blue_f, res_v);
        }
      }
    }
    else
    {
      Halfedge_handle_red    red_he;
      if(assign(red_he, red_obj))
      {
        Halfedge_handle_blue   blue_he;
        if(assign(blue_he, blue_obj))
        {
          //itersection red halfedge and blue halfedge
          m_overlay_traits->create_vertex(red_he, blue_he, res_v);
        }
        else
        {
          Vertex_handle_blue    blue_v;
          CGAL_assertion(assign(blue_v, blue_obj));

          assign(blue_v, blue_obj);
          // blue vertex on red halfedge
          m_overlay_traits->create_vertex(red_he, blue_v, res_v);
        }
      }
      else
      {
        // blue vertex inside red face
        CGAL_assertion(red_obj.is_empty());
        Vertex_handle_blue    blue_v;

        CGAL_assertion(assign(blue_v, blue_obj));
        assign(blue_v, blue_obj);
        Face_handle_red    red_f;
        if(!sc_above)
          red_f = m_red_arr->unbounded_face();
        else
        {
          red_f = sc_above ->get_red_halfedge_handle()->face();
        }
        m_overlay_traits->create_vertex(red_f, blue_v, res_v);
      }
    }
  }




  void create_edge(Subcurve *sc, Halfedge_handle res_he)
  {
    Halfedge_handle_red  red_he;
    Halfedge_handle_blue blue_he;

    // update the result halfedge
    if(sc->get_color() == Curve_info::PURPLE)
    {
       // overlap edge
      red_he = sc->get_red_halfedge_handle();
      blue_he = sc->get_blue_halfedge_handle();
      m_overlay_traits ->create_edge(red_he, blue_he, res_he);
    }
    else
      if(sc->get_color() == Curve_info::RED)
      {
        // red edge on blue face
        red_he = sc->get_red_halfedge_handle();
        Face_handle_blue blue_f;
        Subcurve* sc_above = sc->get_above();
        if(!sc_above)
          blue_f = m_blue_arr->unbounded_face();
        else
          blue_f = sc_above->get_blue_halfedge_handle()->face();
        m_overlay_traits ->create_edge(red_he, blue_f, res_he);
      }
      else
      {
        // blue edge on red face
        CGAL_assertion(sc->get_color() == Curve_info::BLUE);

        blue_he = sc->get_blue_halfedge_handle();
        Face_handle_red red_f;
        Subcurve* sc_above = sc->get_above();
        if(!sc_above)
          red_f = m_red_arr->unbounded_face();
        else
          red_f = sc_above->get_red_halfedge_handle()->face();
        m_overlay_traits ->create_edge(red_f, blue_he, res_he);
      }
  }


  void after_sweep()
  {
    //after sweep finshed, merge the two unbouded_faces from each arrangment
    m_overlay_traits ->create_face(m_red_arr->unbounded_face(),
                                   m_blue_arr->unbounded_face(),
                                   m_arr->unbounded_face());
  }
  


protected:

  const Arrangement1*       m_red_arr;
  const Arrangement2*       m_blue_arr;
  Hash_map                  m_halfedges_map;
  OverlayTraits*            m_overlay_traits;

};



CGAL_END_NAMESPACE

#endif
