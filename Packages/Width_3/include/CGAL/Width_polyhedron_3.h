// ======================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 2000, September 20
//
// file          : include/CGAL/Width_polyhedron_3.h
// package       : Width_3 (1.6)
// maintainer    : Thomas Herrmann <herrmann@ifor.math.ethz.ch>
// chapter       : Geometric Optimisation
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Thomas Herrmann
// coordinator   : ETH Zuerich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// implementation: 3D Width of a Point Set
// ======================================================================

#ifndef CGAL_WIDTH_POLYHEDRON_3_H
#define CGAL_WIDTH_POLYHEDRON_3_H

#include <CGAL/basic.h>
#include <CGAL/Halfedge_data_structure_using_list.h>
#include <CGAL/Halfedge_data_structure_bases.h>
#include <map>
#include <CGAL/width_assertions.h>

CGAL_BEGIN_NAMESPACE

template<class InputPoint, class InputNormal, class InputPlane, 
  class Width_Traits>
class Width_polyhedron_default_traits_3{
 public:
  typedef InputPoint Point;
  typedef InputNormal Normal;
  typedef InputPlane Plane;
  Width_Traits tco;
  void reverse_normal( Normal& normal) const { 
    tco.inverse_normal(normal); 
  }
  void reverse_plane( Plane& plane) const { 
    tco.opposite_plane(plane); 
  }
  Width_polyhedron_default_traits_3() {}
  ~Width_polyhedron_default_traits_3() {}
};

template<class InputPoint, class Width_Traits>
class Width_vertex_default_base : public CGAL::Vertex_max_base<InputPoint> {
 private:
  typedef Width_Traits WT;
  typedef typename WT::RT RT;
  WT tco;
 public:
  Width_vertex_default_base() {}
  Width_vertex_default_base( const InputPoint& p){
    RT px,py,pz,ph;
    tco.get_point_coordinates(p,px,py,pz,ph);
    pt=tco.make_point(px,py,pz,ph);
  }
};

class Width_halfedge_default_base : public CGAL::Halfedge_min_base {
 protected:
  void* prv;
  void* v;
  void* f;
 public:
  typedef Halfedge_min_base Base;
  typedef CGAL::Tag_true  Supports_halfedge_prev;
  typedef CGAL::Tag_true  Supports_halfedge_vertex;
  typedef CGAL::Tag_true  Supports_halfedge_facet;

  Width_halfedge_default_base() : f(NULL) {}

  void* prev() { return prv;} 
  const void* prev() const { return prv;}
  // the previous halfedge along the facet.

  void* vertex() { return v;} 
  const void* vertex() const { return v;} 
  // the incident vertex.

  void* facet() { return f;} 
  const void* facet() const { return f;} 
  //the facet to the left.

  bool is_border() const { return f == NULL;}
  // is true if `h' is a border halfedge).

  void  set_prev( void* h)        { prv = h;}
  void  set_vertex( void* _v)     { v = _v;}
  void  set_facet( void* _f)      { f = _f;}

  // Avoids unnecessary matchings with base class. (g++ 2.7.2 bug)
  void*       opposite()       { return Base::opposite();}
  const void* opposite() const { return Base::opposite();}
  void*       next()           { return Base::next();}
  const void* next() const     { return Base::next();}

};

template <class InputNormal, class InputPlane, class Width_Traits>
class Width_facet_default_base : public CGAL::Facet_max_base {
 public:
  typedef CGAL::Tag_true     Supports_facet_plane;
  typedef CGAL::Tag_true     Supports_facet_normal;
  typedef InputNormal Normal;
  typedef InputPlane Plane;
  Width_Traits tco;
 protected:
  Plane   pln;
 public:
  Normal normal() const { 
    return tco.orthogonal_vector(pln);
  }
  Plane& plane() { 
    return pln;
  }
  const Plane& plane() const { 
    return pln;
  }
};

template <class InputPolyhedron, class Width_Traits>
class Data_access {
 public:
  typedef InputPolyhedron Polyhedron;
  typedef typename Polyhedron::Vertex Vertex;
  typedef typename Polyhedron::Facet Facet;
  typedef typename Polyhedron::Halfedge Halfedge;
  typedef typename Polyhedron::Facet_handle Facet_handle;
  typedef typename Polyhedron::Vertex_handle Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Polyhedron::Point PolyPoint;
  typedef typename Polyhedron::Plane Plane;
  typedef typename Polyhedron::Vertex_iterator Vertex_iterator;
  typedef typename Polyhedron::Facet_iterator Facet_iterator;
  typedef typename Polyhedron::Halfedge_iterator Halfedge_iterator;
 private:
  typedef typename Width_Traits::RT RT;
  Width_Traits tco;
 public:
  Data_access() {}
  ~Data_access() {}

 private:
  //Precondition: Plane Equation already computed in a deterministic way
  struct Facet_compare {
    bool operator()(const Facet_handle& f, 
		    const Facet_handle& g) const {
      Width_Traits tco;
      Plane fpp=f->plane();
      Plane gpp=g->plane();
      RT fa,fb,fc,fd,ga,gb,gc,gd;
      tco.get_plane_coefficients(fpp,fa,fb,fc,fd);
      tco.get_plane_coefficients(gpp,ga,gb,gc,gd);
      return (fa<ga 
	      || fa==ga && fb<gb 
	      || fa==ga && fb==gb && fc<gc
	      || fa==ga && fb==gb && fc==gc && fd<gd);
    }
  };
  //Precondition: Plane Equation already computed in a deterministic way
  //and a facet is bounded by exactly 3 edges! 
  struct Halfedge_compare {
    bool operator()(const Halfedge_handle& e, 
		    const Halfedge_handle& h) const {
      Width_Traits tco;
      PolyPoint etail=e->opposite()->vertex()->point();
      PolyPoint ehead=e->vertex()->point();
      PolyPoint htail=h->opposite()->vertex()->point();
      PolyPoint hhead=h->vertex()->point();
      RT etx,ety,etz,eth,ehx,ehy,ehz,ehh;
      tco.get_point_coordinates(etail,etx,ety,etz,eth);
      tco.get_point_coordinates(ehead,ehx,ehy,ehz,ehh);
      RT htx,hty,htz,hth,hhx,hhy,hhz,hhh;
      tco.get_point_coordinates(htail,htx,hty,htz,hth);
      tco.get_point_coordinates(hhead,hhx,hhy,hhz,hhh);
      
      return (etx*hth <htx*eth ||
	      etx*hth==htx*eth && ety*hth <hty*eth ||
	      etx*hth==htx*eth && ety*hth==hty*eth && etz*hth <htz*eth ||
	      etx*hth==htx*eth && ety*hth==hty*eth && etz*hth==htz*eth
	      && ehx*hhh <hhx*ehh ||
	      etx*hth==htx*eth && ety*hth==hty*eth && etz*hth==htz*eth 
	      && ehx*hhh==hhx*ehh && ehy*hhh <hhy*ehh ||
	      etx*hth==htx*eth && ety*hth==hty*eth && etz*hth==htz*eth 
	      && ehx*hhh==hhx*ehh && ehy*hhh==hhy*ehh && ehz*hhh<hhz*ehh
	      );
    }
  };

  std::map<const Facet_handle, 
    std::vector<Vertex_handle>,
    Facet_compare> 
    antipodal_vertices;

  std::map<const Halfedge_handle, 
    bool, 
    Halfedge_compare> 
    visited_halfedges;

  std::map<const Halfedge_handle, 
    bool, 
    Halfedge_compare> 
    impassable_halfedges;
  
 public:
  bool is_visited(Halfedge_handle& e) const {
    typename std::map < const Halfedge_handle, bool, 
      Halfedge_compare > ::const_iterator it;
    it=visited_halfedges.find(e);
    CGAL_assertion(it!=visited_halfedges.end());
    DEBUGENDL(VISITED_CHECK,"Visited flag value of edge "
	      <<e->opposite()->vertex()->point()
	      <<" --> ",e->vertex()->point()<<": "<<it->second);
    return it->second;
  }
  void set_visited_flag(Halfedge_handle& e, bool val) {
    DEBUGENDL(VISITED_CHECK,"Set visited flag to: ",val);
    visited_halfedges[e]=val;
  }
  
  bool is_impassable(Halfedge_handle& e) const {
    typename std::map < const Halfedge_handle, 
      bool, 
      Halfedge_compare > ::const_iterator it;
    it=impassable_halfedges.find(e);
    CGAL_assertion(it!=impassable_halfedges.end());
    DEBUGENDL(IMPASSABLE_CHECK,"Impassable flag value of edge ",
	      e->opposite()->vertex()->point()
	      <<" --> "<<e->vertex()->point()<<": "<<it->second);
    return it->second;
  }
  void set_impassable_flag(Halfedge_handle& e, bool val) {
    DEBUGENDL(IMPASSABLE_CHECK,"Set impassable flag to: ",val);
    impassable_halfedges[e]=val;
  }
  
  void set_antipodal_vertices(Facet_handle& f, 
				 std::vector<Vertex_handle>& V) {
    antipodal_vertices[f]=V;
  }
  
  void get_antipodal_vertices(Facet_handle& f, 
			      std::vector<Vertex_handle>& res) const {
    typename std::map<const Facet_handle, 
      std::vector<Vertex_handle>, Facet_compare>::const_iterator it;
    it=antipodal_vertices.find(f);
    CGAL_assertion(it!=antipodal_vertices.end());
    res=it->second;
  }

#if !(defined(CGAL_KERNEL_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
      || defined(NDEBUG))
  int size_of_impassable() {
    return(int(impassable_halfedges.size()));
  }
  int size_of_visited() {
    return(int(visited_halfedges.size()));
  }
  int size_of_antipodal_vertices() {
    return (int(antipodal_vertices.size()));
  }
#endif
};

CGAL_END_NAMESPACE

#endif //WIDTH_POLYHEDRON_3_H
