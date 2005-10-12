// Copyright (c) 2005  Stanford University (USA).
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KDS_KINETIC_DELAUNAY_2_H
#define CGAL_KDS_KINETIC_DELAUNAY_2_H
#include <CGAL/KDS/basic.h>



#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/KDS/Delaunay_triangulation_face_base_2.h>
#include <CGAL/KDS/Delaunay_triangulation_visitor_base_2.h>
#include <CGAL/KDS/Notifying_table_listener_helper.h>
#include <CGAL/KDS/Simulator_kds_listener.h>
#include <CGAL/KDS/internal/tds_2_helpers.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/KDS/Ref_counted.h>
#include <iterator>
#include <map>
#include <set>

CGAL_KDS_BEGIN_NAMESPACE

template <class This>
class Delaunay_edge_failure_event  {
public:
  //! Make sure that the s has been advanced
  Delaunay_edge_failure_event(typename This::Root_stack &s,
			      const typename This::Edge &e,
			      This *kdel): s_(s), e_(e), kdel_(kdel){
    if (!s.empty()) s_.pop();
  }
  void process(const typename This::Time&){
    kdel_->flip(e_);
  }
  void write(std::ostream &out) const {
    out << "Flip " << This::TDS_helper::origin(e_)->point() << "," 
	<< This::TDS_helper::destination(e_)->point();
  }
protected:
  typename This::Root_stack s_;
  const typename This::Edge e_;
  This *kdel_;
};

template <class T>
std::ostream& operator<<(std::ostream &out, const Delaunay_edge_failure_event<T> &e){
  e.write(out);
  return out;
}

//! A 2D kinetic Delaunay triangulation.
/*!  Points are added via the Moving_point_table, so the public
  interface is very limited. See kinetic_Delaunay_2.cc for a useage example.
*/
template <class Simulation_traits, class Visitor= Delaunay_triangulation_visitor_base_2, 
	  class Delaunay
	  = CGAL::Delaunay_triangulation_2<typename Simulation_traits::Instantaneous_kernel,
					   CGAL::Triangulation_data_structure_2<
	    CGAL::Triangulation_vertex_base_2<typename Simulation_traits::Instantaneous_kernel>, 
	    CGAL::KDS::Delaunay_triangulation_face_base_2<Simulation_traits > > > >
class Delaunay_triangulation_2: 
  public Ref_counted<Delaunay_triangulation_2<Simulation_traits, Visitor, Delaunay> >
{
public:
  typedef Delaunay_triangulation_2<Simulation_traits, Visitor, Delaunay> This;

  typedef typename Simulation_traits::Kinetic_kernel Kinetic_kernel;
  typedef typename Simulation_traits::Simulator Simulator;
  typedef typename Simulation_traits::Moving_point_table Moving_point_table;

  typedef typename Moving_point_table::Key Point_key;
  typedef typename Simulator::Event_key Event_key;
  typedef typename Simulator::Root_stack Root_stack;
  typedef typename Simulator::Time Time;

  typedef typename Delaunay::Edge_circulator Edge_circulator;
  typedef typename Delaunay::Face_circulator Face_circulator;
  typedef typename Delaunay::Finite_edges_iterator Finite_edges_iterator;
  //typedef typename Delaunay::Edge_iterator Edge_iterator;

  typedef typename Delaunay::Geom_traits::Point_2 Del_point;
  
  typedef typename Delaunay::Vertex_handle Vertex_handle;
  typedef typename Delaunay::Face_handle Face_handle;
  typedef typename Delaunay::Edge Edge;
  typedef typename Delaunay::All_faces_iterator Face_iterator;
  typedef typename Delaunay::All_edges_iterator Edge_iterator;
  typedef typename Kinetic_kernel::Side_of_oriented_circle_2 SOC;
  typedef typename Kinetic_kernel::Orientation_2 O2;
  typedef internal::Triangulation_data_structure_helper_2<typename Delaunay::Triangulation_data_structure> TDS_helper;

  typedef Delaunay_edge_failure_event<This> Certificate;
  friend class Delaunay_edge_failure_event<This>;

  typedef typename CGAL::KDS::Simulator_kds_listener<typename Simulator::Listener, This> Simulator_listener;
  friend  class CGAL::KDS::Simulator_kds_listener<typename Simulator::Listener, This>;
  typedef typename CGAL::KDS::Notifying_table_listener_helper<typename Moving_point_table::Listener, This> Moving_point_table_listener;
  friend class CGAL::KDS::Notifying_table_listener_helper<typename Moving_point_table::Listener, This>;
public:

  Delaunay_triangulation_2(Simulation_traits st,
			   Visitor w= Visitor()): traits_(st), 
						  siml_(st.simulator_pointer(), this), 
						  motl_(st.moving_point_table_pointer(), this),
						  del_(traits_.instantaneous_kernel_object()),
						  soc_(traits_.kinetic_kernel_object().side_of_oriented_circle_2_object()),
						  o2_(traits_.kinetic_kernel_object().orientation_2_object()),
						  watcher_(w){
  }

  
  //! Just write the objects in order;
  void write(std::ostream &out) const {
    out << del_;
  }
  
  typedef Delaunay Triangulation;
  const Triangulation &triangulation(const typename Simulator::NT &t) const {
    //update_instantaneous_kernel_time();
    del_.geom_traits().set_time(t);
    return del_;
  }

  typedef typename Delaunay::Triangulation_data_structure Triangulation_data_structure;
  const Triangulation_data_structure &triangulation_data_structure() const {
    return del_.tds();
  }

  /*const std::set<Edge>& recent_edges() const {
    return new_edges_;
    }*/


  //! Verify that the current state of the 
  void audit() const {
    audit_structure();
    CGAL_KDS_LOG(CGAL::KDS::LOG_LOTS, *this);
    if (del_.dimension() != 2) return;
    Delaunay sdel(traits_.instantaneous_kernel_object());
    sdel.geom_traits().set_time(traits_.simulator_pointer()->rational_current_time());
    sdel.insert(traits_.moving_point_table_pointer()->keys_begin(), 
		traits_.moving_point_table_pointer()->keys_end());

    CGAL_KDS_LOG(CGAL::KDS::LOG_LOTS, sdel << std::endl);

    CGAL_assertion(del_.dimension() == sdel.dimension());

    for (typename Delaunay::All_vertices_iterator vit = sdel.all_vertices_begin();
	 vit != sdel.all_vertices_end(); ++vit){
      bool found=false;
      //Object_key k= vit->point();
      for (typename Delaunay::All_vertices_iterator vit2= del_.all_vertices_begin();
	   vit2 != del_.all_vertices_end(); ++vit2){
	//Object_key k2= vit2->point();
	if (vit->point() == vit2->point()){
	  found=true;
	  //int d= vit->degree();
	  //int d2= vit2->degree();
	  CGAL_assertion(vit->degree() == vit2->degree());
	}
      }
      CGAL_assertion(found);
    }

  }
  
  void set_has_certificates(bool) {
   // always has certificates
  }
  bool has_certificates() {
    return true;
  }

  void audit_structure() const {
    if (del_.dimension() != 2) return;
    CGAL_KDS_LOG(CGAL::KDS::LOG_LOTS, *this);

    for (typename Delaunay::All_edges_iterator eit= del_.all_edges_begin(); 
	 eit != del_.all_edges_end(); ++eit) {
      Event_key key= TDS_helper::get_undirected_edge_label(*eit);
      CGAL_assertion(key);
      typename Simulator::Root_stack s= compute_failure_time(*eit);
      if (key == traits_.simulator_pointer()->null_event()){
	CGAL_exactness_assertion(s.top() == std::numeric_limits<typename Root_stack::Root>::infinity());
      } else {
	CGAL_exactness_assertion(s.top() != std::numeric_limits<typename Root_stack::Root>::infinity());
      }
    }
  }

  void erase(Point_key k){
    // erase all incident certificates
    Vertex_handle vh= vhs_[k];
    watcher_.delete_vertex(vh);
    Face_circulator fc= vh->incident_faces(), fe=fc;
    if (fc != NULL) {
      do {
	for (unsigned int j=0; j<3; ++j){
	  Edge e(fc, j);
	  Event_key k= TDS_helper::get_undirected_edge_label(e);
	  if (k){
	    traits_.simulator_pointer()->delete_event(k);
	    TDS_helper::set_undirected_edge_label(e, Event_key());
	  }
	}
	++fc;
      } while (fc != fe);
    }
    // remove from triangulation
    del_.geom_traits().set_time(traits_.simulator_pointer()->rational_current_time());
    del_.remove(vh);
    //new_edges_.clear();
    if (del_.dimension()==2){
      std::vector<Face_handle> faces;
      
      del_.get_conflicts(k,std::back_inserter(faces));
      
      for (unsigned int i=0; i< faces.size(); ++i){
	for (unsigned int j=0; j<3; ++j){
	  Edge e(faces[i],j);
	  Event_key k= TDS_helper::get_undirected_edge_label(e);
	  if (!k){
	    // a bit redundant for certificates which don't fail
	    new_certificate(e);
	  } else {
	    //CGAL_assertion(new_edges_.find(TDS_helper::mirror_edge(e)) != new_edges_.end());
	  }
	  //new_edges_.insert(e);
	}
      }
      watcher_.new_faces(faces.begin(), faces.end());
    }
   
  }










  //! The assertion will catch that the object is in the same sorted order
  void set(Point_key k){
    //std::cout << "Object changed " << k << std::endl;

    //new_edges_.clear();

    if (del_.dimension() != 2) {
      CGAL_KDS_LOG(CGAL::KDS::LOG_SOME,"Triangulation is still 1D.\n");
      return;
    }
    Vertex_handle vh= vhs_[k];
    watcher_.change_vertex(vh);
    Edge_circulator ec= vh->incident_edges(), ef=ec;
    if (ec != NULL){
      do {
	rebuild_certificate(*ec);
	++ec;
      } while (ec != ef);
    }
    Face_circulator fc= vh->incident_faces(), fe= fc;
    if (fc != NULL) {
      do {
	int i= fc->index(vh);
	rebuild_certificate(Edge(fc, i));
	++fc;
      } while (fc != fe);
    }
    //write(std::cout);
  }









  
  //! 
  /*!
    Some old certificate edges will be lost, have to find all conflicts first
  */
  void insert(Point_key k){
    bool was_2d= (del_.dimension()==2);

    del_.geom_traits().set_time(traits_.simulator_pointer()->rational_current_time());
    if (was_2d){
      //std::cout << "removing extra certificates.\n";
      std::vector<Face_handle> faces;
      del_.get_conflicts(k, std::back_inserter(faces));
      for (unsigned int i=0; i< faces.size(); ++i){
	Face_handle f= faces[i];
	for (unsigned int j=0; j<3; ++j){
	  Edge e(f, j);
	  Event_key k= TDS_helper::get_undirected_edge_label(e);
	  if (k){
	    traits_.simulator_pointer()->delete_event(k);
	    TDS_helper::set_undirected_edge_label(e, Event_key());
	  }
	}
      }
      watcher_.delete_faces(faces.begin(), faces.end());
    }
   
    vhs_[k]= del_.insert(k);
    watcher_.new_vertex(vhs_[k]);
    // now have to update
    if (!was_2d && del_.dimension()==2){
      //std::cout << "Creating certificates from scratch.\n";
      for (Edge_iterator eit = del_.all_edges_begin(); eit != del_.all_edges_end(); ++eit){
	TDS_helper::set_undirected_edge_label(*eit, Event_key());
	new_certificate(*eit);
      }
      watcher_.new_faces(del_.all_faces_begin(), del_.all_faces_end());
    } else {
      set(k);
    }
    CGAL_expensive_postcondition_code(audit_structure());
    //write(std::cout);
  }









  Edge flip(const Edge &e){
    //CGAL_KDS_LOG(CGAL::KDS::LOG_NONE, "Flipping edge " << TDS_helper::origin(e)->point());
    //CGAL_KDS_LOG(CGAL::KDS::LOG_NONE, TDS_helper::destination(e)->point() << std::endl);
    //CGAL_KDS_LOG(CGAL::KDS::LOG_NONE, " at "  << traits_.simulator()->current_time() << std::endl);
    
    
    
    Face_handle face= e.first;
    int index= e.second;
    int mirror_index = face->mirror_index(index);
    Face_handle mirror_face = face->neighbor(index);

    
    for (unsigned int i=0; i<3; ++i){
      Edge e0(face, i);
      if (e0.second != index){
	traits_.simulator_pointer()->delete_event(TDS_helper::get_undirected_edge_label(e0));
	TDS_helper::set_undirected_edge_label(e0, Event_key());
      } 
      Edge e1(mirror_face, i);
      if (e1.second != mirror_index){
	traits_.simulator_pointer()->delete_event(TDS_helper::get_undirected_edge_label(e1));
	TDS_helper::set_undirected_edge_label(e1, Event_key());
      }
    }

    TDS_helper::set_undirected_edge_label(e, Event_key());
    watcher_.pre_flip(e);
    del_.tds().flip(face,index);
   
   
    // we also know that CGAL preserves the edge index of the flipped edge
    mirror_index = mirror_face->index(face);
    index= face->index(mirror_face);
    
    
    Edge flipped_edge(face,index);
    //CGAL_postcondition(del_.is_face(face));

    mirror_index = face->mirror_index(index);
    mirror_face = face->neighbor(index);
    
    for (unsigned int i=0; i<3; ++i){
      Edge e0(face, i);
      if (!TDS_helper::get_undirected_edge_label(e0)){
	new_certificate(e0);
      } 
      Edge e1(mirror_face, i);
      if (!TDS_helper::get_undirected_edge_label(e1)){
	new_certificate(e1);
      }
    }

    //write(std::cout);
    //new_edges_.clear();
    //new_edges_.insert(flipped_edge);

    CGAL_KDS_LOG(CGAL::KDS::LOG_SOME, "Created " << TDS_helper::origin(flipped_edge)->point());
    CGAL_KDS_LOG(CGAL::KDS::LOG_SOME, TDS_helper::destination(flipped_edge)->point() << std::endl);
    

    CGAL_expensive_postcondition_code(audit_structure());
    watcher_.post_flip(flipped_edge);
    return flipped_edge;
  }
 

  Visitor &visitor() {
    return watcher_;
  }


  const Visitor &visitor() const {
    return watcher_;
  }

protected:
  Simulation_traits traits_;
  Simulator_listener siml_; 
  Moving_point_table_listener motl_; 
  Delaunay del_;
  std::map<Point_key, Vertex_handle> vhs_;
  SOC soc_;
  O2 o2_;
  //std::set<Edge> new_edges_;
  Visitor watcher_;
  
  const typename Moving_point_table::Data& point(Point_key k) const{
    return traits_.moving_point_table_pointer()->at(k);
  }
  

  typename Simulator::Root_stack compute_failure_time(const Edge &e) const {
    Point_key ks[4];
    ks[0]= TDS_helper::origin(e)->point();
    ks[1]= TDS_helper::third_vertex(e)->point();
    ks[2]= TDS_helper::destination(e)->point();
    ks[3]= TDS_helper::mirror_vertex(e)->point();
    
    bool odd_parity=false;
    bool infinity=false;
    for (unsigned int i=0; i<4; ++i){
      if (infinity){
	ks[i-1]=ks[i];
      } else {
	if (!ks[i]){
	  infinity=true;
	  odd_parity= ((i%2)==1);
	}
      }
    }
    Root_stack s;
    if (infinity){
      if (odd_parity){
	std::swap(ks[0], ks[1]);
      }
      typename O2::result_type rt=o2_(point(ks[0]), point(ks[1]), point(ks[2]));
      s= traits_.simulator_pointer()->root_stack_object(rt);
    } else {
      s=traits_.simulator_pointer()->root_stack_object(soc_(point(ks[0]), point(ks[1]),
							    point(ks[2]), point(ks[3])));
    }
    
    return s;
  }

  void new_certificate( const Edge &e){
    CGAL_precondition(!TDS_helper::get_undirected_edge_label(e));
    
    Root_stack s= compute_failure_time(e);
    Event_key k;
    if (!s.empty()) {
      Time t= s.top();
      s.pop();
      k =traits_.simulator_pointer()->new_event(t, Certificate(s, e, this));
    } else {
      k = traits_.simulator_pointer()->null_event();
    }
    TDS_helper::set_undirected_edge_label(e, k);
    //new_edges_.insert(e);
    //std::cout << "Created " << k <<std::endl;
  }

  //! rebuild a certificates
  /*!  I need to check if there is a valid one before since I use
    change_object to initialize the certificates of a new object.
  */
  void rebuild_certificate( const Edge &e){
    if (TDS_helper::get_undirected_edge_label(e)){
      traits_.simulator_pointer()->delete_event(TDS_helper::get_undirected_edge_label(e));
      TDS_helper::set_undirected_edge_label(e,  Event_key());
    }
    new_certificate(e);
  }
};

template <class Sim, class Del, class W> 
std::ostream &operator<<(std::ostream &out, const Delaunay_triangulation_2<Sim, Del, W> &kd){
  kd.write(out);
  return out;
}

CGAL_KDS_END_NAMESPACE
#endif
