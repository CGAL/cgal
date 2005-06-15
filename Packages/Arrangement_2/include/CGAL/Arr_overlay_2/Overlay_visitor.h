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

#include <CGAL/Sweep_line_2/Arr_sweep_line_visitor.h>
#include <CGAL/Arrangement_2/Arr_accessor.h>
#include <CGAL/Unique_hash_map.h> 

CGAL_BEGIN_NAMESPACE

template <class _Arrangement>
class Hash_function
{
private:
  Arr_accessor<_Arrangement>    m_arr_accessor;

public:
  typedef std::size_t                              result_type;
  typedef typename _Arrangement::Halfedge_handle   Halfedge_handle;

  Hash_function(_Arrangement& arr) : m_arr_accessor(arr)
  {}

  std::size_t operator() (Halfedge_handle he) const
  {
    return std::size_t(&*m_arr_accessor.halfedge(he)) / 
      sizeof( *m_arr_accessor.halfedge(he));
  } 
};


template < class _Traits,
           class _Arrangement1,
           class _Arrangement2,
           class _Arrangement,
           class Event,
           class Subcurve,
           class OverlayTraits >
class Overlay_visitor : 
  public Arr_sweep_line_visitor<_Traits, _Arrangement, Event, Subcurve>
{
public:

  typedef _Traits                                        Traits; 
  typedef typename Traits::X_monotone_curve_2            X_monotone_curve_2;
  typedef typename Traits::Point_2                       Point_2;
  typedef typename Traits::Curve_info                    Curve_info;

  typedef _Arrangement                                   Arrangement;
  typedef typename Arrangement::Halfedge_handle          Halfedge_handle;
  typedef typename Arrangement::Face_handle              Face_handle;
  typedef typename Arrangement::Vertex_handle            Vertex_handle;

  typedef _Arrangement1                                  Arrangement1;
  typedef typename Arrangement1::Halfedge_const_handle   Halfedge_handle_red;
  typedef typename Arrangement1::Face_const_handle       Face_handle_red;
  typedef typename Arrangement1::Vertex_const_handle     Vertex_handle_red;

  typedef _Arrangement2                                  Arrangement2;
  typedef typename Arrangement2::Halfedge_const_handle   Halfedge_handle_blue;
  typedef typename Arrangement2::Face_const_handle       Face_handle_blue;
  typedef typename Arrangement2::Vertex_const_handle     Vertex_handle_blue;

 
  typedef typename Arrangement::Ccb_halfedge_circulator 
    Ccb_halfedge_circulator;
  
  typedef Overlay_visitor< Traits,
                           Arrangement1,
                           Arrangement2,
                           Arrangement,
                           Event,
                           Subcurve,
                           OverlayTraits >               Self;
  
  typedef Arr_sweep_line_visitor< Traits,
                                  Arrangement,
                                  Event,
                                  Subcurve >             Base;
  
  typedef typename Event::SubCurveIter                   SubCurveIter;
  typedef typename Event::SubCurveRevIter                SubCurveRevIter;

  typedef Sweep_line_2_impl<Traits,
                            Event,
                            Subcurve,
                            Self,
                            CGAL_ALLOCATOR(int)>         Sweep_line;
  typedef typename Sweep_line::StatusLineIter            StatusLineIter;
  
  typedef Hash_function<Arrangement>                     Hash_function; 
  typedef Unique_hash_map< Halfedge_handle,
                           Curve_info,
                           Hash_function >               Hash_map;
  using Base::m_arr;
  using Base::m_arr_access;
  using Base::m_sweep_line;
  using Base::m_currentEvent;


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
                      blue_arr.number_of_halfedges(),
                    Hash_function(res_arr)),
    m_overlay_traits(&overlay_trairs)
  {}


  /*! Destructor */
  virtual ~Overlay_visitor() {}


  void before_handle_event(Event* event)
  {
    Base::before_handle_event(event);
  }

  bool after_handle_event(Event* event, StatusLineIter iter, bool flag)
  {
    bool res = Base::after_handle_event(event, iter, flag);

    SubCurveRevIter rev_iter = event->right_curves_rbegin();
    Subcurve* sc_above = NULL;
    if( iter != static_cast<Sweep_line*>(m_sweep_line) -> StatusLine_end())
      sc_above = static_cast<Subcurve*>(*iter);
    else
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

  void init_subcurve(Subcurve* sc)
  {
    Base::init_subcurve(sc);
  }

  void add_subcurve(const X_monotone_curve_2& cv,Subcurve* sc)
  {
    Base::add_subcurve(cv, sc);
  }



  virtual Halfedge_handle insert_in_face_interior
    (const X_monotone_curve_2& cv,
     Subcurve* sc)
  {
    // res is directed from left to right
    Halfedge_handle res = Base::insert_in_face_interior(cv,sc);
    map_halfedge_and_twin(res, true, cv.get_curve_info());

    return res;
  }


  virtual Halfedge_handle insert_from_right_vertex
                          (const X_monotone_curve_2& cv,
                           Halfedge_handle he)
  {
    // res is directed from right to left
    Halfedge_handle res = Base::insert_from_right_vertex(cv, he);
    map_halfedge_and_twin(res, false, cv.get_curve_info());

    return res;
  }


  virtual Halfedge_handle insert_from_left_vertex
                          (const X_monotone_curve_2& cv,
                           Halfedge_handle he)
  {
    //res is directed from left to right
    Halfedge_handle res = Base::insert_from_left_vertex(cv, he);
    map_halfedge_and_twin(res, true, cv.get_curve_info());

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
      Face_handle new_face = res.face();

      //traverse new_face's boundary
      Ccb_halfedge_circulator ccb_end = new_face.outer_ccb();
      Ccb_halfedge_circulator ccb_circ = ccb_end;
      do
      { 
        //get the current halfedge on the face boundary
        Halfedge_handle he (*ccb_circ);

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
        Face_handle_red red_face = red_he.face();
        Face_handle_blue blue_face = blue_he.face();
        m_overlay_traits->create_face(red_face, blue_face,new_face);
      }
      else
      {
        // red face inside blue face
        if(red_he != Halfedge_handle_red())
        {
          Face_handle_red red_face = red_he.face();
          Face_handle_blue blue_face;
          Subcurve* sc_above = sc->get_above();
          if(!sc_above)
            blue_face = m_blue_arr->unbounded_face();
          else
            blue_face = 
              sc_above->get_blue_halfedge_handle().face();
          
          m_overlay_traits->create_face(red_face, blue_face,new_face);
        }
        else
        {
          // blue face inside red face
          CGAL_assertion(blue_he != Halfedge_handle_blue() && 
                         red_he == Halfedge_handle_red());
          Face_handle_red red_face;
          Face_handle_blue blue_face = blue_he.face();
          Subcurve* sc_above = sc->get_above();
          if(!sc_above)
            red_face = m_red_arr->unbounded_face();
          else
            red_face = 
              sc_above->get_red_halfedge_handle().face();
          
          m_overlay_traits->create_face(red_face, blue_face,new_face);
        }
      }
    }
    return res;
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
      red_he_info_twin = red_he_info.twin();

    if(blue_he_info != Halfedge_handle_blue())
      blue_he_info_twin = blue_he_info.twin();


    //cv_info_twin will have the twin halfedge
    Curve_info cv_info_twin(red_he_info_twin, blue_he_info_twin);

    if(right_dir)
    {
      m_halfedges_map[he] = cv_info_twin;
      m_halfedges_map[he.twin()] = cv_info;
    }
    else
    {
      m_halfedges_map[he] = cv_info;
      m_halfedges_map[he.twin()] = cv_info_twin;
    }
  }



protected:

  const Arrangement1*       m_red_arr;
  const Arrangement2*       m_blue_arr;
  Hash_map                  m_halfedges_map;
  OverlayTraits*            m_overlay_traits;
};




CGAL_END_NAMESPACE


#endif
