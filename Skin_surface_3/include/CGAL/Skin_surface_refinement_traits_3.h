// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
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
// $URL$
// $Id$
// 
//
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_SKIN_SURFACE_SUBDIVISION_POLICY_H
#define CGAL_SKIN_SURFACE_SUBDIVISION_POLICY_H

#include <CGAL/intersection_3_1.h>
#include <CGAL/Skin_surface_polyhedral_items_3.h>

CGAL_BEGIN_NAMESPACE

template <class Polyhedron_3, class SkinSurface_3>
class Skin_surface_subdivision_policy_base_3 {
public:
  typedef Polyhedron_3                            Polyhedron;
  typedef SkinSurface_3                           Skin_surface;
  typedef typename Polyhedron::Traits             P_traits;
  typedef typename Skin_surface::TMC_Traits       SS_traits;

  typedef Cartesian_converter<SS_traits, P_traits> T2P_converter;
  typedef Cartesian_converter<P_traits, SS_traits> P2T_converter;

  typedef typename Polyhedron::Vertex_handle      P_vertex_handle;

  typedef typename P_traits::RT          P_rt;
  typedef typename P_traits::Point_3     P_point;
  typedef typename P_traits::Segment_3   P_segment;
  typedef typename P_traits::Line_3      P_line;
  typedef typename P_traits::Vector_3    P_vector;
  typedef typename P_traits::Plane_3     P_plane;

  typedef typename Skin_surface::TMC_Cell_handle   SS_cell_handle;
  typedef typename Skin_surface::TMC_Vertex_handle SS_vertex_handle;

  typedef std::map<SS_vertex_handle,bool>     SS_vertex_map;
  typedef typename SS_vertex_map::iterator    SS_vertex_map_it;

  Skin_surface_subdivision_policy_base_3(Skin_surface const& skin)
    : skin(skin), t2p_converter(), p2t_converter()
  {
    
  }
    
  virtual P_point to_surface(P_vertex_handle vh) = 0;
//  {
//     SS_cell_handle ch = locate(vh->point());
//     return to_surface_along_transversal_segment(vh->point(),ch);
//   }

//   SS_cell_handle locate(const P_point &p) {
//     return skin.locate(p2t_converter(p));
//   }

  virtual P_vector normal(P_vertex_handle vh) = 0;
//   {
//     SS_cell_handle ch = skin.locate(p2t_converter(vh->point()));
//     return ch->surf->gradient(vh->point());
//   }

  virtual ~Skin_surface_subdivision_policy_base_3() {}
protected:
  P_point to_surface_along_transversal_segment(
    P_point const &p, SS_cell_handle ch) {
      
    P_segment transversal_segment = get_transversal_segment(ch,p);
      
    return ch->surf->to_surface(
      transversal_segment.start(), transversal_segment.end());
  }
private:
  P_segment get_transversal_segment(
    SS_cell_handle ch, P_point const &p) {
    // Compute signs on vertices and sort them:
    int nIn = 0;
    int sortedV[4];
    for (int i=0; i<4; i++) {
      if (is_inside(ch,i)) {
        sortedV[nIn] = i; nIn++;
      } else {
        sortedV[3-i+nIn] = i;
      }
    }
    P_point begin_point, end_point;
    Object obj;
    if (nIn==1) {
      begin_point = t2p_converter(ch->vertex(sortedV[0])->point());
      obj = CGAL::intersection(
        P_plane(
          t2p_converter(ch->vertex(sortedV[1])->point()),
          t2p_converter(ch->vertex(sortedV[2])->point()),
          t2p_converter(ch->vertex(sortedV[3])->point())),
        P_line(begin_point, p));
      if ( !assign(end_point, obj) )
        CGAL_assertion_msg(false,"intersection: no intersection.");
    } else if (nIn==2) {
      obj = CGAL::intersection(
        P_plane(
          t2p_converter(ch->vertex(sortedV[2])->point()),
          t2p_converter(ch->vertex(sortedV[3])->point()),
          p),
        P_line(
          t2p_converter(ch->vertex(sortedV[0])->point()),
          t2p_converter(ch->vertex(sortedV[1])->point())));
      if ( !assign(begin_point, obj) )
        CGAL_assertion_msg(false,"intersection: no intersection.");
      obj = CGAL::intersection(
        P_plane(
          t2p_converter(ch->vertex(sortedV[0])->point()),
          t2p_converter(ch->vertex(sortedV[1])->point()),
          p),
        P_line(
          t2p_converter(ch->vertex(sortedV[2])->point()),
          t2p_converter(ch->vertex(sortedV[3])->point())));
      if ( !assign(end_point, obj) )
        CGAL_assertion_msg(false,"intersection: no intersection.");
    } else if (nIn==3) {
      end_point = t2p_converter(ch->vertex(sortedV[3])->point());
      obj = CGAL::intersection(
        P_plane(
          t2p_converter(ch->vertex(sortedV[0])->point()),
          t2p_converter(ch->vertex(sortedV[1])->point()),
          t2p_converter(ch->vertex(sortedV[2])->point())),
        P_line(end_point, p));
      if ( !assign(begin_point, obj) )
        CGAL_assertion_msg(false,"intersection: no intersection.");
    } else {
      CGAL_assertion(false);
    }
    return P_segment(begin_point, end_point);
  }
  
  P_point 
  intersection(P_plane const &plane, P_line const &line) {
    P_point p;

    Object result = CGALi::intersection(plane, line);
    if ( !CGAL::assign(p, result) )
      CGAL_assertion_msg(false,"intersection: no intersection.");
    return p;
  }

  Sign sign(SS_cell_handle ch, int i) {
    return CGAL_NTS sign(
       ch->surf->value(t2p_converter(ch->vertex(i)->point())));
  }  
  bool is_inside(SS_cell_handle ch, int i) {
    //return (sign(ch,i) == POSITIVE);
    SS_vertex_map_it it = triang_vertex_signs.find(ch->vertex(i));
    
    if (it == triang_vertex_signs.end()) {
      bool side = (sign(ch,i) == POSITIVE);
      triang_vertex_signs[ch->vertex(i)] = side;
      CGAL_assertion(triang_vertex_signs[ch->vertex(i)] == side);
      return side;
    } else {
      return it->second;
    }
  }

  
protected:
  Skin_surface const &skin;
  SS_vertex_map triang_vertex_signs;
  T2P_converter const t2p_converter;
  P2T_converter const p2t_converter;
};

template <class Polyhedron_3, class SkinSurface_3>
class Skin_surface_subdivision_policy_with_face_info_3 
  : public Skin_surface_subdivision_policy_base_3<Polyhedron_3, SkinSurface_3> 
{
  typedef Polyhedron_3                           Polyhedron;
  typedef SkinSurface_3                          Skin_surface;
  typedef Skin_surface_subdivision_policy_base_3<Polyhedron, Skin_surface> Base;
public:
  typedef typename Base::P_point             P_point;
  typedef typename Base::P_vector            P_vector;
  typedef typename Base::P_vertex_handle     P_vertex_handle;
  typedef typename Base::SS_cell_handle      SS_cell_handle;

  Skin_surface_subdivision_policy_with_face_info_3(Skin_surface const& skin)
    : Base(skin) {
  }
    
  P_point to_surface(P_vertex_handle vh) {
    SS_cell_handle ch = vh->halfedge()->facet()->triang_ch;
    return Base::to_surface_along_transversal_segment(vh->point(),ch);
  }
  P_vector normal(P_vertex_handle vh) {
    SS_cell_handle ch = vh->halfedge()->facet()->triang_ch;
    return ch->surf->gradient(vh->point());
  }
};

template <class Polyhedron_3, class SkinSurface_3>
class Skin_surface_subdivision_policy_default_3 
  : public Skin_surface_subdivision_policy_base_3<Polyhedron_3, SkinSurface_3> 
{
  typedef Polyhedron_3                           Polyhedron;
  typedef SkinSurface_3                          Skin_surface;
  typedef Skin_surface_subdivision_policy_base_3<Polyhedron, Skin_surface> Base;
public:
  typedef typename Base::P_point             P_point;
  typedef typename Base::P_vector            P_vector;
  typedef typename Base::P_vertex_handle     P_vertex_handle;
  typedef typename Base::SS_cell_handle      SS_cell_handle;

  Skin_surface_subdivision_policy_default_3(Skin_surface &skin) 
    : Base(skin) {}
  
  P_point to_surface(P_vertex_handle vh) 
  {
    SS_cell_handle ch = Base::skin.locate(Base::p2t_converter(vh->point()));
    return to_surface_along_transversal_segment(vh->point(),ch);
  }
  
  P_vector normal(P_vertex_handle vh) 
  {
    SS_cell_handle ch = Base::skin.locate(Base::p2t_converter(vh->point()));
    return ch->surf->gradient(vh->point());
  }
};

template <class Polyhedron_3,
	  class SkinSurface_3>
Skin_surface_subdivision_policy_base_3<Polyhedron_3, SkinSurface_3> *
get_subdivision_policy(Polyhedron_3  &p,
		       SkinSurface_3 &skinsurface) 
{
  typedef Skin_surface_subdivision_policy_default_3<Polyhedron_3, SkinSurface_3>
    Policy;
  return new Policy(skinsurface);
}

template <class P_Traits,
	  class SkinSurface_3>
Skin_surface_subdivision_policy_base_3<
  Polyhedron_3<P_Traits, Skin_surface_polyhedral_items_3<SkinSurface_3> >, 
  SkinSurface_3> *
get_subdivision_policy(Polyhedron_3<
		       P_Traits, 
		       Skin_surface_polyhedral_items_3<SkinSurface_3> > &p,
		       SkinSurface_3 &skinsurface) 
{
  typedef Polyhedron_3<P_Traits, 
    Skin_surface_polyhedral_items_3<SkinSurface_3> >        Polyhedron;
  typedef Skin_surface_subdivision_policy_with_face_info_3<
    Polyhedron, 
    SkinSurface_3>                                          Policy;
  return new Policy(skinsurface);
}

CGAL_END_NAMESPACE

#endif // CGAL_SKIN_SURFACE_SUBDIVISION_TRAITS_H
