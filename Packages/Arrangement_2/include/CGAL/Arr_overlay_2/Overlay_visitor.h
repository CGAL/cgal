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

#ifndef OVERLAY_VISITOR
#define OVERLAY_VISITOR

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
  typedef std::size_t result_type;
  typedef typename _Arrangement::Halfedge_handle   Halfedge_handle;

  Hash_function(_Arrangement& arr) : m_arr_accessor(arr)
  {}

  std::size_t operator() (Halfedge_handle he) const
  {
    return std::size_t(&*m_arr_accessor.halfedge(he)) / 
      sizeof( *m_arr_accessor.halfedge(he));
  } 
};


template <class _Traits, class _Arrangement, class Event, class Subcurve>
class Ovelay_visitor : 
  public Arr_sweep_line_visitor<_Traits, _Arrangement, Event, Subcurve>
{
public:

  typedef _Traits                                             Traits;
  typedef _Arrangement                                        Arrangement;
  typedef typename Arrangement::Halfedge_handle               Halfedge_handle;
  typedef typename Arrangement::Face_handle                   Face_handle;
  typedef typename Arrangement::Vertex_handle                 Vertex_handle;
  typedef typename Arrangement::Ccb_halfedge_circulator
                                                      Ccb_halfedge_circulator;
  
  typedef Ovelay_visitor<Traits,Arrangement,Event, Subcurve>  Self;
  
  typedef Arr_sweep_line_visitor<Traits,Arrangement,Event, Subcurve>
                                                              Base;
  
  typedef typename Event::SubCurveIter                        SubCurveIter;
  typedef typename Event::SubCurveRevIter                     SubCurveRevIter;

  typedef Sweep_line_2_impl<Traits,
                            Event,
                            Subcurve,
                            Self,
                            CGAL_ALLOCATOR(int)>              Sweep_line;
  typedef typename Sweep_line::StatusLineIter                 StatusLineIter;

   
  typedef typename Traits::X_monotone_curve_2            X_monotone_curve_2;
  typedef typename Traits::Point_2                       Point_2;
  

  typedef typename Traits::Curve_info                    Curve_info;
  typedef Hash_function<Arrangement>                     Hash_function;
  
  typedef Unique_hash_map<Halfedge_handle, Curve_info, Hash_function>
                                                         Hash_map;
  using Base::m_arr;
  using Base::m_arr_access;
  using Base::m_sweep_line;
  using Base::m_currentEvent;


  private:

  Ovelay_visitor (const Self& );
  Self& operator= (const Self& );

public:

  /*! Constructor */
  Ovelay_visitor(Arrangement& res_arr,
                 const Arrangement&  red_arr,
                 const Arrangement&  blue_arr):
    Arr_sweep_line_visitor(res_arr),
    /*m_red_arr(&red_arr),
    m_blue_arr(&blue_arr),*/
    m_halfedges_map(Curve_info(),
                    red_arr.number_of_halfedges() 
                    + 
                    blue_arr.number_of_halfedges(),
                    Hash_function(res_arr))
  {}


  /*! Destructor */
  virtual ~Ovelay_visitor() {}


  void before_handle_event(Event* event)
  {
    Base::before_handle_event(event);
  }

  void after_handle_event(Event* event)
  {
    Base::after_handle_event(event, StatusLineIter iter, bool flag);

    Subcurve* sc_above = NULL;
    if( iter != static_cast<Sweep_line*>(m_sweep_line) -> StatusLine_end())
      sc_above = *iter;

    for(SubCurveRevIter rev_iter = event->right_curves_rbegin();
        rev_iter != event->right_curves_rend();
        ++rev_iter)
    {
      Subcurve* curr_sc = *rev_iter;

      if(!curr_sc->has_same_color(sc_above))
        curr_sc->set_above(above_sc);
      else
        curr_sc->set_above(above_sc->get_above());

      sc_above = curr_sc
    }
  }

  void init_subcurve(Subcurve* sc)
  {
    Base::init_subcurve(sc);
  }

  void add_subcurve(const X_monotone_curve_2& cv,Subcurve* sc)
  {}



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
    Halfedge_handle res = Base::insert_at_vertices(cv, hhandle, prev,
                                                   new_face_created);
    //TODO: add an assertion that res is directed from right to left

    map_halfedge_and_twin(res, false, cv.get_curve_info());

    // new face was created, need to update the new face's data
    if(new_face_created)
    {
      Haldedge null_he(NULL);
      Halfedge_handle red_he  (null_he);
      Halfedge_handle blue_he (null_he);

      // get the new face
      Face_handle new_face = res.face();

      //traverse new_face's boundary
      Ccb_halfedge_circulator ccb_end = new_face.outer_ccb();
      Ccb_halfedge_circulator ccb_circ = ccb_end;
      do
      { 
        //get the current halfedge on the fce boundary
        Halfedge_handle he(ccb_circ.halfedge());

        CGAL_assertion(m_halfedges_map.is_defined(he);
        const Curve_info& cv_info = m_halfedges_map[he];
        if(cv_info.get_color() == Curve_info::RED)
        {
          red_he = cv_info.get_halfedge();
          if(blue_he != null_he)
            break;
        }
        else
          if(cv.info.get_color() == Curve_info::BLUE)
          {
            blue_he = cv_info.get_halfedge();
            if(red_he != null_he)
              break;
          }
          else 
          {
            // overlap
            CGAL_assertion(cv_info.get_color() == Curve_info::PURPLR);
            red_he = cv_info.get_pair_halfedges().first;
            blue_he = cv_info.get_pair_halfedges().second;
            break;
          }
        ++ccb_circ;
      }
      while(ccb_circ != ccb_end);
      //finished traversing the boundary of the face

      if(red_he != null_he && blue_he != null_he)
      {
        Face_handle red_face = red_he.face();
        Face_handle blue_handle = blue.face();
        //TODO : merge the two face's data
      }
      else
      {
        if(red_he != null_he)
        {
          Face_handle red_face = red_he.face();
          Face_handle blue_face = sc->get_above()->get_halfedge_handle().face();
          //TODO : merge the two face's data
        }
        else
        {
          CGAL_assertion(blue_he != null_he && red_he == null_he);
           Face_handle red_face = sc->get_above()->get_halfedge_handle().face();
          Face_handle blue_face = blue.face();
          //TODO : merge the two face's data
        }
      }
    }
  }





  // maps halfedge and his twin, right_dir is true iff he is directed from
  // left to right
  void map_halfedge_and_twin(Halfedge_handle he, bool right_dir, 
                             const Curve_info& cv_info)
  {
    // original halfedge that was stored is directed from right to left
    Halfedge_handle he_info = cv_info.get_halfedge_handle();

    //cv_info_twin will have the twin halfedge
    Curve_info cv_info_twin(he_info.twin(), cv_info.get_color());

    if(right_dir)
    {
      m_halfedges_map[he] = cv_info_twin;
      m_halfedges_map[res->twin()] = cv_info;
    }
    else
    {
      m_halfedges_map[he] = cv_info;
      m_halfedges_map[res->twin()] = cv_info_twin;
    }
  }





protected:

  /*const Arrangement*        m_red_arr;
  const Arrangement*        m_blue_arr;*/
  Hash_map                  m_halfedges_map;
};




CGAL_END_NAMESPACE


#endif
