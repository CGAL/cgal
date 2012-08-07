// Copyright (c) 2005  Max-Planck-Institute Saarbruecken (Germany).
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
// $URL$ 
// $Id$ 
// 
//
// Author(s)     : Ralf Osbild <osbild@mpi-sb.mpg.de>

#ifndef CGAL_NEF_VERTEX_CYCLE_TO_NEF_3_H
#define CGAL_NEF_VERTEX_CYCLE_TO_NEF_3_H

#include <iostream>
#include <sstream>

// triangulation
#include <CGAL/Nef_3/Exact_triangulation_euclidean_traits_xy_3.h>
#include <CGAL/Nef_3/Exact_triangulation_euclidean_traits_yz_3.h>
#include <CGAL/Nef_3/Exact_triangulation_euclidean_traits_xz_3.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

// Nef polyhedra
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/SNC_structure.h>
#include <CGAL/Nef_3/SNC_constructor.h>
#include <CGAL/Nef_3/SNC_point_locator.h>

namespace CGAL {

/*
template<typename Items>
class Index_assigner {
  Index_assigner() {}
  template<typename Handle>
    void assign_index(Handle& ) const {}
  template<typename Handle>
    void assign_first_index() const {}
  template<typename Handle>
    void assign_new_index() {} 
};

template<> class Index_assigner<CGAL::SNC_indexed_items> {
  int first;
  int index;
 public:
 Index_assigner() : 
  first(Index_generator::get_unique_index()), index(first) {}
  
  template<typename Handle>
    void assign_index(Handle& h) const 
    { h->set_index(index); }
  template<typename Handle>
    void assign_first_index(Handle& h) const 
    { h->set_index(first); }
  template<typename Handle>
    void assign_new_index(Handle& h)
    { h->set_index(); index = h->get_index(); } 
};
*/

//template<typename Items> class Index_matcher;

template<typename Vertex>
class Compare_cpte {

  typedef typename std::pair<Vertex*, Vertex*> VV;

 public:
  Compare_cpte() {}

  bool operator()(const VV p0, const VV p1) const {
    if(p0.first != p1.first)
      return p0.first < p1.first;
    else
      return p0.second < p1.second;
  }  
};

class Compare_face {

 public:
  Compare_face() {}

  template<typename Face>  
    bool operator()(const Face* f0, const Face* f1) const {
    return f0 < f1;
  }
};


template<typename Items,
  typename Edge, typename CompareEdges> class Index_matcher {
 public:
  Index_matcher() {}
  template<typename Handle> 
  void set_index(Handle /*h*/, Edge /*e*/) {}
};

template<typename Edge, typename CompareEdges> 
  class Index_matcher<CGAL::SNC_indexed_items, Edge, CompareEdges> 
{

  bool plusTwin;
  std::map<Edge, int, CompareEdges> edge2int;
  typedef typename std::map<Edge, int, CompareEdges>::iterator e2i_iterator;

 public:
 Index_matcher(bool withTwin) : plusTwin(withTwin) {}

  template<typename Handle> void set_index(Handle h, Edge e) {
    e2i_iterator ei = edge2int.find(e);
    if(ei != edge2int.end()) {
      h->set_index(ei->second);
      if(plusTwin) h->twin()->set_index(h->get_index()+1);
    } else {
      int new_index = Index_generator::get_unique_index();
      edge2int.insert(std::make_pair(e,new_index));
      h->set_index(new_index);
      if(plusTwin)
	h->twin()->set_index();
    }
  }
};

// return value reports success ("true" means nef contains result)
// note: facets are considered to be compact
// CTP - Constrained_triangulation_plus
template <class CTP, class Nef_3, class II>
bool projected_vertex_cycle_to_nef_3 (typename Nef_3::SNC_structure &snc,
      II v_first, II v_last, std::ostringstream &ostr)
{
   typedef typename CTP::Vertex           CTP_vertex;
   typedef typename CTP::Vertex_handle    CTP_vertex_handle;
   typedef typename CTP::Face             CTP_face;
   typedef typename CTP::Face_handle      CTP_face_handle;
   typedef typename CTP::Face_circulator  CTP_face_circulator;
   typedef typename CTP::size_type        CTP_size_type;

   typedef typename Nef_3::SNC_structure            SNC_structure;
   typedef typename SNC_structure::SM_decorator       SM_decorator;
   typedef typename SNC_structure::Vertex_handle      Vertex_handle;
   typedef typename SNC_structure::Sphere_point       Sphere_point;
   typedef typename SNC_structure::Sphere_circle      Sphere_circle;
   typedef typename SNC_structure::SVertex_handle     SVertex_handle;
   typedef typename SNC_structure::SHalfedge_handle   SHalfedge_handle;
   typedef typename SNC_structure::SFace_handle       SFace_handle;

   typedef typename std::pair<CTP_vertex*, CTP_vertex*> Point_pair;
   typedef Compare_cpte<CTP_vertex> Compare_edge;

   typedef typename SNC_structure::Items               Items;

   // declarations and defaults
   II v_it, v_pred_it;
   CTP ctp;
   bool cond;
   CTP_size_type nov;

   snc.clear();

   Index_matcher<Items, Point_pair, Compare_edge> im_edge(false);
   Index_matcher<Items, CTP_face*, Compare_face> im_sedge(true);

   for (nov=0, v_it=v_first; v_it!=v_last; ++v_it)
   {  if ( *v_it != (ctp.insert (*v_it))->point() )
      {  ostr << " -> Different vertices are projected to same location!";
	 return false;
      }
      ++nov;
   }

   if ( ctp.number_of_vertices() != nov )
   {  ostr << " -> Same vertex appears multiple in cycle.";
      return false;
   }
   if ( nov < 3 )
   {  ostr << " -> Not at least 3 vertices.";
      return false;
   }

   // construct constrained triangulation
   v_it = v_first;
   while ( true )
   {  // move iterators
      v_pred_it = v_it;
      ++v_it;

      if ( v_it == v_last ) break; // while-end
      if ( *v_pred_it == *v_it ) continue ; // no self-loops
      ctp.insert_constraint (*v_pred_it, *v_it);
      if ( ctp.number_of_vertices() != nov ) break; // constraints intersect
   }
   if ( v_it == v_last && *v_pred_it != *v_first) // no self-loops
   {  ctp.insert_constraint (*v_pred_it, *v_first);
   }

   // assertion
   if ( ctp.number_of_vertices() != nov )
   {  ostr << " -> Vertex cycle is not simple; edges intersect.";
      return false;
   }
   CGAL_assertion (ctp.is_valid());

   // determine initial positions for the walk-around
   CTP_vertex_handle t_vh, t_vh_0;
   CTP_face_handle t_fh;
   CTP_face_circulator t_fc, t_fc_0, t_fc_1;
   int idx(0);
   {  // search a constrained edge from "outside"
      t_fh = ctp.infinite_face ();
      cond = t_fh->has_vertex(ctp.infinite_vertex(), idx);
      CGAL_assertion ( cond );
      idx = ctp.ccw(idx);
      t_vh = t_fh->vertex(idx);
      // now: *t_vh lies ccw from the infinite vertex

      t_fc = t_fc_0 = ctp.incident_faces(t_vh, t_fh);
      if(t_fc == CTP_face_circulator())
	return false;
      do
      {  if ((cond=t_fc->is_constrained(ctp.cw(t_fc->index(t_vh))))) break;
      } while (--t_fc != t_fc_0)
      ; CGAL_assertion ( cond ); // do-while ends with break
      t_fh = t_fc->neighbor(ctp.cw(t_fc->index(t_vh)));
      // note: if (not sure at this moment!) input is a simple
      // polygon, *t_fh lies inside this polygon, *t_vh is a vertex
      // of *t_fh and the edge of *t_fh starting clockwise at *t_vh
      // is a constraint.
   }

   // walk along constraints and build sphere maps
   t_vh_0 = t_vh;
   do
   {  // circulate faces beginning with *t_fh around *t_vh
      idx = t_fh->index(t_vh);
      CTP_vertex_handle t_target_vh = t_fh->vertex(ctp.cw(idx));

      // create new vertex in SNC
      Vertex_handle p_vh = snc.new_vertex (t_vh->point(), true);
      SM_decorator p_dec (&*p_vh); // &*

      // create new svertex in SM
      Sphere_point dir ( CGAL::ORIGIN+
         (t_target_vh->point() - t_vh->point()) ), dir_0 = dir;
      SVertex_handle p_svh (p_dec.new_svertex (dir)), p_svh_pred=p_svh;
      p_svh->mark() = true;
      if(&*t_target_vh < &*t_vh)
	im_edge.set_index(p_svh, std::make_pair(&*t_target_vh, &*t_vh));
      else
	im_edge.set_index(p_svh, std::make_pair(&*t_vh, &*t_target_vh));

      // consider triangles (implicit edges) that are incident to *t_vh
      t_fc = t_fc_0 = t_fc_1 = ctp.incident_faces(t_vh, t_fh);
      ++t_fc_1;
      do
      {  idx = t_fc->index(t_vh);
	 t_target_vh = t_fc->vertex(ctp.ccw(idx));

	 // assertion
	 if ( ctp.is_infinite(t_target_vh) )
         {  ostr << " -> Vertex cycle is not simple; edges overlap.";
            return false;
	 }

	 // create new svertex in SM
         dir = Sphere_point ( CGAL::ORIGIN+
	    (t_target_vh->point() - t_vh->point()) );

	 // assertion
	 if ( dir == dir_0 )
         {  ostr << " -> Vertex cycle is not simple; edges overlap.";
            return false;
	 }
	 p_svh = SVertex_handle (p_dec.new_svertex(dir));
	 p_svh->mark() = true;
	 if(&*t_target_vh < &*t_vh)
	   im_edge.set_index(p_svh, std::make_pair(&*t_target_vh, &*t_vh));
	 else
	   im_edge.set_index(p_svh, std::make_pair(&*t_vh, &*t_target_vh));
	 
	 // create new sphere edges in SM
         SHalfedge_handle p_seh =
	    p_dec.new_shalfedge_pair (p_svh_pred, p_svh);
         Sphere_circle p_scirc
	    (p_seh->source()->point(), p_seh->target()->point());
         p_seh->circle() = p_scirc;
         p_seh->twin()->circle() = p_scirc.opposite();
         p_seh->mark() = p_seh->twin()->mark() = true;
	 im_sedge.set_index(p_seh, &*t_fc);

	 // constrained edge detected?
	 if ((cond = t_fc->is_constrained(ctp.cw(idx))))
	 {  // create new sface in SM
	    SFace_handle p_fh = p_dec.new_sface ();
	    p_fh->mark() = false;
	    p_dec.link_as_face_cycle(p_seh,p_fh);
	    break;
	 }
	 p_svh_pred = p_svh;
      } while (--t_fc != t_fc_0)
      ; CGAL_assertion ( cond ); // do-while ends with break

      // check
      p_dec.check_integrity_and_topological_planarity();

      // prepare next iteration
      t_fh = t_fc;

      // assertion
      --t_fc;
      while (--t_fc != t_fc_0)
      {  idx = t_fc->index(t_vh);
	 if ( t_fc->is_constrained(ctp.ccw(idx)) )
         {  ostr << " -> Vertex cycle is not simple; edge contains vertex.";
            return false;
         }
      }

      // prepare next iteration
      t_vh = t_fh->vertex(ctp.ccw(t_fh->index(t_vh)));
   } while (t_vh != t_vh_0)
   ; // do-while ends

   return true;
}

// II - input iterator; KN - kernel of normal
// return value reports success ("true" means snc contains result)
template <class Nef_3, class II, class KN>
bool vertex_cycle_to_nef_3 (
   typename Nef_3::SNC_structure &snc, II v_first, II v_last,
   const CGAL::Vector_3<KN> &normal, bool verb = false)
{
   // Constrained_triangulation_plus with Exact_intersections_tag
   typedef typename Nef_3::Kernel          Kernel;
   typedef CGAL::Exact_intersections_tag   Tri_itag;

   // XY projection + triangulation
typedef CGAL::Exact_triangulation_euclidean_traits_xy_3<Kernel>  XY_kernel;
typedef CGAL::Triangulation_vertex_base_2<XY_kernel>             XY_vb;
typedef CGAL::Constrained_triangulation_face_base_2<XY_kernel>   XY_fb;
typedef CGAL::Triangulation_data_structure_2<XY_vb,XY_fb>        XY_ds;
typedef CGAL::Constrained_triangulation_2<XY_kernel, XY_ds, Tri_itag> XY_tri;
typedef CGAL::Constrained_triangulation_plus_2<XY_tri>     XY_tri_plus;

   // XZ projection + triangulation
typedef CGAL::Exact_triangulation_euclidean_traits_xz_3<Kernel>  XZ_kernel;
typedef CGAL::Triangulation_vertex_base_2<XZ_kernel>             XZ_vb;
typedef CGAL::Constrained_triangulation_face_base_2<XZ_kernel>   XZ_fb;
typedef CGAL::Triangulation_data_structure_2<XZ_vb,XZ_fb>        XZ_ds;
typedef CGAL::Constrained_triangulation_2<XZ_kernel, XZ_ds, Tri_itag> XZ_tri;
typedef CGAL::Constrained_triangulation_plus_2<XZ_tri>     XZ_tri_plus;

   // YZ projection + triangulation
typedef CGAL::Exact_triangulation_euclidean_traits_yz_3<Kernel>  YZ_kernel;
typedef CGAL::Triangulation_vertex_base_2<YZ_kernel>             YZ_vb;
typedef CGAL::Constrained_triangulation_face_base_2<YZ_kernel>   YZ_fb;
typedef CGAL::Triangulation_data_structure_2<YZ_vb,YZ_fb>        YZ_ds;
typedef CGAL::Constrained_triangulation_2<YZ_kernel, YZ_ds, Tri_itag> YZ_tri;
typedef CGAL::Constrained_triangulation_plus_2<YZ_tri>     YZ_tri_plus;

   // defaults
   bool is_nef = false;
   char direc='?';
   snc.clear();
   std::ostringstream ostr;

   if ( normal == NULL_VECTOR )
   {  // report it
      ostr << " -> function parameter 'normal' is NULL_VECTOR"
	   << " (this can be a symptom of an error).";
   }

   if ( normal != NULL_VECTOR )
   {  // projection depending on normal vector
      direc='z';
      if ( CGAL::abs(normal.x()) > CGAL::abs(normal.y()) )
      {  if ( CGAL::abs(normal.x()) > CGAL::abs(normal.z()) ) direc='x';
      }
      else
      {  if ( CGAL::abs(normal.y()) > CGAL::abs(normal.z()) ) direc='y';
      }

      // project and triangulate vertices,
      // convert result to Nef polyhedron
      ostr << " Direction of projection is '" << direc << "'.";
      if ( direc == 'x' )
      {  is_nef = projected_vertex_cycle_to_nef_3<YZ_tri_plus,Nef_3> (
                  snc, v_first, v_last, ostr);
      }
      else if ( direc == 'y' )
      {  is_nef = projected_vertex_cycle_to_nef_3<XZ_tri_plus,Nef_3> (
                  snc, v_first, v_last, ostr);
      }
      else
      {  CGAL_assertion ( direc == 'z' );
         is_nef = projected_vertex_cycle_to_nef_3<XY_tri_plus,Nef_3> (
                  snc, v_first, v_last, ostr);
      }
   }

   if ( !is_nef )
   {  // if conversion is unsuccessful so far, try another projection
      if ( !is_nef && direc != 'x' )
      {  ostr << " Now, direction of projection is 'x'.";
	 is_nef = projected_vertex_cycle_to_nef_3<YZ_tri_plus,Nef_3> (
                  snc, v_first, v_last, ostr);
      }
      if ( !is_nef && direc != 'y' )
      {  ostr << " Now, direction of projection is 'y'.";
         is_nef = projected_vertex_cycle_to_nef_3<XZ_tri_plus,Nef_3> (
                  snc, v_first, v_last, ostr);
      }
      if ( !is_nef && direc != 'z' )
      {  ostr << " Now, direction of projection is 'z'.";
         is_nef = projected_vertex_cycle_to_nef_3<XY_tri_plus,Nef_3> (
                  snc, v_first, v_last, ostr);
      }
   }

   // no successful conversion?
   if ( !is_nef && verb )
   {  std::cerr << "\nConversion from vertex cycle to Nef_polyhedron_3"
	 << " was not successful. Error history:"
         << ostr.str().c_str()
	 << " Finally, empty Nef_polyhedron_3 was constructed." << std::endl;
   }

   return is_nef;
}

} //namespace CGAL
#endif // CGAL_NEF_VERTEX_CYCLE_TO_NEF_3_H
