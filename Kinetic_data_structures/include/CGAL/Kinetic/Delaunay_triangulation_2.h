// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_KINETIC_DELAUNAY_2_H
#define CGAL_KINETIC_KINETIC_DELAUNAY_2_H
#include <CGAL/Kinetic/basic.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Kinetic/Delaunay_triangulation_face_base_2.h>
#include <CGAL/Kinetic/Delaunay_triangulation_vertex_base_2.h>
#include <CGAL/Kinetic/Delaunay_triangulation_visitor_base_2.h>
#include <CGAL/Kinetic/Active_objects_batch_listener_helper.h>
#include <CGAL/Kinetic/Simulator_kds_listener.h>
#include <CGAL/Kinetic/internal/tds_2_helpers.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Kinetic/Ref_counted.h>
#include <iterator>
#include <CGAL/Kinetic/Event_base.h>
#include <CGAL/Kinetic/Delaunay_triangulation_default_traits_2.h>

CGAL_KINETIC_BEGIN_NAMESPACE

#define CGAL_DELAUNAY_2_DEBUG(x) x

template <class KDel>
struct Delaunay_edge_failure_event: public Event_base<KDel*> {
  typedef Event_base<KDel*>  P;
  Delaunay_edge_failure_event(const typename KDel::Certificate_data &c,
			      const typename KDel::Edge &e,
			      KDel *kdel): P(kdel), c_(c), e_(e) {}
  typename KDel::Certificate_data &certificate() const {
    return c_;
  }
  const typename KDel::Edge edge() const {
    return e_;
  }
  KDel* kdel() {
    return P::kds();
  }

  KDel* kdel() const {
    return P::kds();
  }
  void process() {
    kdel()->flip(e_, c_);
  }

  CGAL::Comparison_result compare_concurrent(typename KDel::Event_key a,
					     typename KDel::Event_key b) const {
    return kdel()->compare_concurrent(a, b);
  }	

  std::ostream& write(std::ostream &out) const
  {
    out << "Flip " << KDel::TDS_helper::origin(edge())->point() << ","
	<< KDel::TDS_helper::destination(edge())->point() 
	<< " to " << KDel::TDS_helper::third_vertex(edge())->point() 
	<< ", " << KDel::TDS_helper::mirror_vertex(edge())->point() ;
    return out;
  }

  typename KDel::Certificate_data c_;
  const typename KDel::Edge e_;
};





//! A 2D kinetic Delaunay triangulation.
/*!  Points are added via the Moving_point_table, so the public
  interface is very limited. See kinetic_Delaunay_2.cc for a useage example.
*/
template <class Simulation_traits_t, 
	  class Visitor= Delaunay_triangulation_visitor_base_2,
	  class Delaunay
	  = CGAL::Delaunay_triangulation_2<typename Simulation_traits_t::Instantaneous_kernel,
					   CGAL::Triangulation_data_structure_2<
	    Delaunay_triangulation_vertex_base_2<typename Simulation_traits_t::Instantaneous_kernel>,
  CGAL::Kinetic::Delaunay_triangulation_face_base_2<Simulation_traits_t > > >,
	  class Delaunay_traits_t= Delaunay_triangulation_default_traits_2<Simulation_traits_t, Delaunay> >
class Delaunay_triangulation_2:
  public Ref_counted<Delaunay_triangulation_2<Simulation_traits_t, Visitor, Delaunay, Delaunay_traits_t> >
{

  typedef CGAL::Delaunay_triangulation_2<typename Simulation_traits_t::Instantaneous_kernel,
					 CGAL::Triangulation_data_structure_2<
    Delaunay_triangulation_vertex_base_2<typename Simulation_traits_t::Instantaneous_kernel>,
    CGAL::Kinetic::Delaunay_triangulation_face_base_2<Simulation_traits_t > > > Basic_Delaunay;

public:
  typedef Delaunay_traits_t Traits;
  typedef Simulation_traits_t Simulation_traits;
  typedef Delaunay_triangulation_2<Simulation_traits, Visitor, Delaunay, Traits> This;

  typedef typename Simulation_traits::Kinetic_kernel Kinetic_kernel;
  typedef typename Simulation_traits::Simulator Simulator;
  typedef typename Simulation_traits::Active_points_2_table Moving_point_table;

  typedef typename Moving_point_table::Key Point_key;
  typedef typename Simulator::Event_key Event_key;
  //typedef typename Simulator::Root_stack Root_stack;

  typedef typename Traits::Triangulation Triangulation;

  typedef typename Triangulation::Edge_circulator Edge_circulator;
  typedef typename Triangulation::Face_circulator Face_circulator;
  typedef typename Triangulation::Finite_edges_iterator Finite_edges_iterator;
  //typedef typename Triangulation::Edge_iterator Edge_iterator;

  typedef typename Triangulation::Geom_traits::Point_2 Del_point;

  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Face_handle Face_handle;
  typedef typename Triangulation::Edge Edge;
  typedef typename Triangulation::All_faces_iterator Face_iterator;
  typedef typename Triangulation::All_edges_iterator Edge_iterator;
  typedef typename Traits::Certificate_data Certificate_data;
  typedef typename Traits::Time Time;
  typedef internal::Triangulation_data_structure_helper_2<typename Triangulation::Triangulation_data_structure> TDS_helper;

  typedef Delaunay_edge_failure_event<This> Event;
  
  //friend class Delaunay_edge_failure_event<This>;
  //friend class Delaunay_hull_edge_failure_event<This>;

  typedef typename CGAL::Kinetic::Simulator_kds_listener<typename Simulator::Listener, This> Simulator_listener;
  friend  class CGAL::Kinetic::Simulator_kds_listener<typename Simulator::Listener, This>;
  typedef typename CGAL::Kinetic::Active_objects_batch_listener_helper<typename Moving_point_table::Listener, This> Moving_point_table_listener;
  friend class CGAL::Kinetic::Active_objects_batch_listener_helper<typename Moving_point_table::Listener, This>;

  struct Compare_edges{
    bool operator()(const Edge &a, const Edge &b) const {
      Point_key a0, a1, b0, b1;
      a0= a.first->vertex((a.second+1)%3)->point();
      a1= a.first->vertex((a.second+2)%3)->point();
      b0= b.first->vertex((b.second+1)%3)->point();
      b1= b.first->vertex((b.second+2)%3)->point();
      if (a0 > a1) std::swap(a0, a1);
      if (b0 > b1) std::swap(b0, b1);
      if (a0 < b0) return true;
      else if (a0 > b0) return false;
      else return a1 < b1;
    }
  };

  void init_data() {
    siml_ = Simulator_listener(traits_.simulator_handle(), this);
    motl_= Moving_point_table_listener(traits_.active_points_2_table_handle(), this);
    has_certificates_=false; 
    clear_stats();
   
    batching_=false;
 
  }

public:

  Delaunay_triangulation_2(Traits st,
			   Triangulation del,
			   Visitor w= Visitor()): 
    traits_(st),
    watcher_(w),
    del_(del) {
    vhs_.resize(del_.number_of_vertices());
    for (typename Triangulation::Vertex_iterator vit = del_.vertices_begin(); vit != del_.vertices_end(); ++vit) {
      CGAL_assertion(vit->point().to_index() < del_.number_of_vertices());
      vhs_[vit->point().to_index()]=vit;
    }
    init_data();
  
    set_has_certificates(true);
  }
  Delaunay_triangulation_2(Simulation_traits st,
			   Visitor w= Visitor()): 
    traits_(st),
    watcher_(w),
    del_(traits_.instantaneous_kernel_object()) {
    init_data();
    set_has_certificates(true);
  }



  //! Just write the objects in order;
  void write(std::ostream &out) const
  {
    out << del_;
  }

  void clear_stats() {
    num_events_=0;
    // num_single_certificates_=0;
  }

  void write_stats(std::ostream &out) const {
    out << "Num events is " << num_events_ << std::endl;
  }

  const Triangulation &triangulation(const typename Simulator::NT &t) const
  {
    //update_instantaneous_kernel_time();
    del_.geom_traits().set_time(t);
    return del_;
  }

  const Triangulation &triangulation() const
  {
    return del_;
  }


  // for Qt_triangulation
  /*Event_key null_event() const {
    return traits_.simulator_handle()->null_event();
    }*/

  /*const Simulation_traits& simulation_traits_object() const {
    return traits_;
    }*/

  typedef typename Triangulation::Triangulation_data_structure Triangulation_data_structure;
  const Triangulation_data_structure &triangulation_data_structure() const
  {
    return del_.tds();
  }

  bool is_batch_editing() const {
    return batching_;
  }

  void set_is_batch_editing(bool tf) {
    if (tf== batching_) return;
    else if (tf) {
      batching_=true;
    } else {
      //unsigned int num_certs= num_certificates_;
      // this is important that it be before update
      batching_=false;

      for (unsigned int i=0; i< batched_certs_.size(); ++i){
	Point_key s=TDS_helper::origin(batched_certs_[i])->point();
	Point_key t=TDS_helper::destination(batched_certs_[i])->point();
	
	std::cout << std::min(s,t) << "--" << std::max(s,t) << std::endl;
	//update_edge(batched_certs_[i]);
      }

      std::sort(batched_certs_.begin(), batched_certs_.end(), Compare_edges());
      /*std::cout << "Batch has " << batched_certs_.size() << " non-unique edges to compute." << std::endl;*/
      batched_certs_.erase(std::unique(batched_certs_.begin(), 
				       batched_certs_.end()), 
			   batched_certs_.end());
      //std::cout << "Batch has " << batched_certs_.size() << " edges to compute." << std::endl;
      for (unsigned int i=0; i< batched_certs_.size(); ++i){
	Point_key s=TDS_helper::origin(batched_certs_[i])->point();
	Point_key t=TDS_helper::destination(batched_certs_[i])->point();
	
	std::cout << std::min(s,t) << "--" << std::max(s,t) << std::endl;
	//update_edge(batched_certs_[i]);
      }
      for (unsigned int i=0; i< batched_certs_.size(); ++i){
	//Point_key s=TDS_helper::origin(batched_certs_[i])->point();
	//Point_key t=TDS_helper::destination(batched_certs_[i])->point();
	
	//std::cout << std::min(s,t) << "--" << std::max(s,t) << std::endl;
	update_edge(batched_certs_[i]);
      }
      
      batched_certs_.clear();
      /*CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME, 
       *traits_.simulator_handle() << std::endl;);*/
      //int dnum= num_certificates_-num_certs;
      //std::cout << "Edit had " << dnum << " certificate computations" << std::endl;

      //audit();
    }
  }
  
  /*const std::set<Edge>& recent_edges() const {
    return new_edges_;
    }*/

  //! Verify that the current state of the
  void audit() const
  {
    if (!has_certificates_) return;
    CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME, "Auditing delaunay" << std::endl);
    
    if (del_.number_of_vertices() < 50) {
      CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_LOTS, *this);
    }
    if (del_.dimension() != 2) return;
    Basic_Delaunay sdel(traits_.instantaneous_kernel_object());
    sdel.geom_traits().set_time(traits_.rational_current_time());
    for (typename Triangulation::Finite_vertices_iterator cit= del_.finite_vertices_begin();
	 cit != del_.finite_vertices_end(); ++cit){
      sdel.insert(cit->point());
    }
    /*    sdel.insert(traits_.active_points_2_table_handle()->keys_begin(),
	  traits_.active_points_2_table_handle()->keys_end());*/

    //CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_LOTS, sdel << std::endl);

    if (del_.dimension() != sdel.dimension()) {
      CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_NONE, "AUDIT FAILURE Dimensions don't match in audit" << std::endl);
      return;
    }
    CGAL_exactness_assertion(del_.dimension() == sdel.dimension());
   
   for (typename Triangulation::All_vertices_iterator vit = del_.all_vertices_begin();
	vit != del_.all_vertices_end(); ++vit) {
     if (vit->point() != Point_key()) {
	if (vit->neighbors_is_changed()) {
	  CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_NONE, "AUDIT FAILURE neighbors is changed " 
			   << vit->point() << std::endl);
	  CGAL_exactness_assertion(!vit->neighbors_is_changed());

	}
	if (vit->neighbors() != vit->degree()) {
	  CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_NONE, "AUDIT FAILURE stored degree is " << vit->neighbors() 
			   << " and actual is " << vit->degree() << " for " << vit->point() << std::endl);
	  CGAL_exactness_assertion(vit->neighbors() == vit->degree());
	}
      }
   }
   for (typename Basic_Delaunay::All_vertices_iterator vit = sdel.all_vertices_begin();
	 vit != sdel.all_vertices_end(); ++vit) {
      bool found=false;
    
      //Object_key k= vit->point();
      for (typename Triangulation::All_vertices_iterator vit2= del_.all_vertices_begin();
	   vit2 != del_.all_vertices_end(); ++vit2) {
	//Object_key k2= vit2->point();
	if (vit->point() == vit2->point()) {
	  found=true;
	  //int d= vit->degree();
	  //int d2= vit2->degree();
	  if (vit->degree() != vit2->degree()) {
	    CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_NONE, "AUDIT FAILURE Degrees don't match in: " 
			     << vit->point() << std::endl);
	  }
	  CGAL_exactness_assertion(vit->degree() == vit2->degree());


	}
      }
      if (!found) {
	CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_NONE, "AUDIT FAILURE Matching vertex not found: " 
			 << vit->point() << std::endl);
      }
      CGAL_exactness_assertion(found);
     
    }

  

    typename Simulation_traits::Instantaneous_kernel ik= traits_.instantaneous_kernel_object();
    ik.set_time(traits_.rational_current_time());

    typename Simulation_traits::Instantaneous_kernel::Side_of_oriented_circle_2 ic2
      = ik.side_of_oriented_circle_2_object();
    for (typename Triangulation::Finite_edges_iterator fit = del_.finite_edges_begin(); 
	 fit != del_.finite_edges_end(); ++fit){
      Point_key k0= fit->first->vertex((fit->second+1)%3)->point();
      Point_key k2= fit->first->vertex((fit->second+2)%3)->point();
      Point_key k3= TDS_helper::mirror_vertex(*fit)->point();
      Point_key k1= TDS_helper::third_vertex(*fit)->point();
      if (k1== Point_key() || k3== Point_key()) continue;
      typename Triangulation::Geom_traits::Current_coordinates cc= del_.geom_traits().current_coordinates_object();
      typedef typename Triangulation::Geom_traits::Current_coordinates::result_type P2;
      P2 p0= cc(k0);
      P2 p1= cc(k1);
      P2 p2= cc(k2);
      P2 p3= cc(k3);
      if (ic2(k0, k1, k2, k3) != CGAL::ON_POSITIVE_SIDE) {
	CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_NONE, "AUDIT FAILURE Failed certificate: " << k0 << " " << k1 << " " 
			 << k2 << " " << k3 << std::endl);
	CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_NONE, "AUDIT FAILURE Points are: " << p0 << ": " << p1 << ": " << p2 
			 << ": " << p3 << std::endl);
      }
      CGAL_exactness_assertion(ic2(k0, k1, k2, k3) == CGAL::ON_POSITIVE_SIDE);
      
    }

    for (typename Triangulation::Edge_iterator fit = del_.edges_begin(); fit != del_.edges_end(); ++fit){
      if (TDS_helper::origin(*fit)->degree()==3 || TDS_helper::destination(*fit)->degree()==3) {
	if (TDS_helper::get_undirected_edge_label(*fit) != Event_key()) {
	  CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_NONE, "AUDIT FAILURE certificate on degree 3 edge: " 
			   << TDS_helper::origin(*fit)->point()
			   << " " <<  TDS_helper::destination(*fit)->point() 
			   << TDS_helper::origin(*fit)->degree() <<  " "
			   << TDS_helper::destination(*fit)->degree() << std::endl);
	} 
	CGAL_exactness_assertion(TDS_helper::get_undirected_edge_label(*fit) == Event_key());
      } else {
	if (TDS_helper::get_undirected_edge_label(*fit) == Event_key()) {
	  CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_NONE, "AUDIT FAILURE no certificate on edge: " 
			   << TDS_helper::origin(*fit)->point()
			   << " " <<  TDS_helper::destination(*fit)->point() << std::endl);
	  CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_NONE, "AUDIT FAILURE degrees are: " 
			   << TDS_helper::origin(*fit)->degree()
			   << " " <<  TDS_helper::destination(*fit)->degree() << std::endl);
	  
	} 
	CGAL_exactness_assertion(TDS_helper::get_undirected_edge_label(*fit) != Event_key());
      }
    }

  }

 

  void set_has_certificates(bool tf) {
    if (tf == has_certificates_){

    } else {
      if (tf==true && del_.dimension()==2) {
	CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME, "DELAUNAY2: Creating certificates."<< std::endl);
	for (typename Triangulation::All_vertices_iterator vit = del_.all_vertices_begin(); 
	     vit != del_.all_vertices_end(); ++vit) {
	  int deg=vit->degree();
	  CGAL_assertion(deg >=3);
	  Vertex_handle vh= vit;
	  vh->set_neighbors(deg);
	  CGAL_DELAUNAY_2_DEBUG(std::cout << "Set degree of " << vit->point() << " to " << deg << std::endl);
	  //vit->set_neighbors_is_changed(false);
	}
	for (Edge_iterator eit = del_.all_edges_begin(); eit != del_.all_edges_end(); ++eit) {
	  TDS_helper::set_undirected_edge_label(*eit, Event_key());
	  update_edge(*eit);
	}
	watcher_.create_faces(del_.all_faces_begin(), del_.all_faces_end());
      } else if (tf==false) { 

	for (typename Triangulation::Edge_iterator it = del_.edges_begin();
	     it != del_.edges_end(); ++it){
	  if (TDS_helper::get_undirected_edge_label(*it) != Event_key()){
	    delete_certificate(*it);
	  }
	}
      } 
      CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME, 
		       *traits_.simulator_handle() << std::endl;);
      has_certificates_=tf;
    }
  }
  bool has_certificates() {
    return has_certificates_;
  }

  void erase(Point_key k) {
    // erase all incident certificates
    Vertex_handle vh= vertex_handle(k);
    if (vh == Vertex_handle()) {
      CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME, "Point " << k << " is not in triangulation on removal."<< std::endl);
      return;
    }
    watcher_.remove_vertex(vh);
    if (has_certificates_) {
      Face_circulator fc= vh->incident_faces(), fe=fc;
      if (fc != NULL) {
	do {
	  for (unsigned int j=0; j<3; ++j) {
	    Edge e(fc, j);
	  Event_key k= TDS_helper::get_undirected_edge_label(e);
	  if (k.is_valid()) {
	    delete_certificate(e);
	  }
	  fc->vertex(j)->set_neighbors_is_changed(true);
	  }
	  ++fc;
	} while (fc != fe);
      }
    }
    // remove from triangulation
    del_.geom_traits().set_time(traits_.rational_current_time());
    del_.remove(vh);
    //new_edges_.clear();
    if (del_.dimension()==2 && has_certificates_) {
      std::vector<Face_handle> faces;

      del_.get_conflicts(k,std::back_inserter(faces));

      for (unsigned int i=0; i< faces.size(); ++i) {
	for (unsigned int j=0; j<3; ++j) {
	  Edge e(faces[i],j);
	  //Event_key k= TDS_helper::get_undirected_edge_label(e);
	  // a bit redundant for certificates which don't fail
	  update_edge(e);
	  //new_edges_.insert(e);
	}
      }
      watcher_.create_faces(faces.begin(), faces.end());
    }

  }

  //! The assertion will catch that the object is in the same sorted order
  void set(Point_key k) {
    //std::cout << "Object changed " << k << std::endl;

    //new_edges_.clear();
    traits_.point_changed(k);
    if (del_.dimension() != 2) {
      CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME,"Triangulation is still 1D.\n");
      return;
    }

    Vertex_handle vh=vertex_handle(k);
    if (vh == Vertex_handle()) {
      CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME, "Point " << k << " is not in triangulation on set."<< std::endl);
      return;
    }
    if (has_certificates_) {
      Edge_circulator ec= vh->incident_edges(), ef=ec;
      if (ec != NULL) {
	do {
	  if (TDS_helper::get_undirected_edge_label(*ec) != Event_key()) {
	    delete_certificate(*ec);
	  }
	  ++ec;
	} while (ec != ef);
      }
      Face_circulator fc= vh->incident_faces(), fe= fc;
      if (fc != NULL) {
	do {
	  int i= fc->index(vh);
	  Edge e(fc, i);
	  if (TDS_helper::get_undirected_edge_label(e) != Event_key()) {
	    delete_certificate(e);
	  }
	  ++fc;
	} while (fc != fe);
      }
    }

   
    watcher_.modify_vertex(vh);

    if (has_certificates_) {
      Edge_circulator ec= vh->incident_edges(), ef=ec;
      if (ec != NULL) {
	do {
	  update_edge(*ec);
	  ++ec;
	} while (ec != ef);
      }
      Face_circulator fc= vh->incident_faces(), fe= fc;
      if (fc != NULL) {
	do {
	  int i= fc->index(vh);
	  update_edge(Edge(fc, i));
	  ++fc;
	} while (fc != fe);
      }
    }
    //write(std::cout);
  }

  //!
  /*!
    Some old certificate edges will be lost, have to find all conflicts first
  */
  void insert(Point_key k) {
    // evil hack
    if (vertex_handle(k) != Vertex_handle()) return;
    CGAL_DELAUNAY_2_DEBUG(std::cout << "Inserting " << k << std::endl);
    bool was_2d= (del_.dimension()==2);

    del_.geom_traits().set_time(traits_.rational_current_time());
    if (was_2d && has_certificates_) {
      //std::cout << "removing extra certificates.\n";
      std::vector<Face_handle> faces;
      del_.get_conflicts(k, std::back_inserter(faces));
      for (unsigned int i=0; i< faces.size(); ++i) {
	Face_handle f= faces[i];
	for (unsigned int j=0; j<3; ++j) {
	  f->vertex(j)->set_neighbors_is_changed(true);
	  Edge e(f, j);
	  Event_key k= TDS_helper::get_undirected_edge_label(e);
	  if (k != Event_key()) {
	    delete_certificate(e);
	  }
	}
      }
      watcher_.remove_faces(faces.begin(), faces.end());

      if (faces.empty()) {
	CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME, "DELAUNAY vertex not successfully inserted " << k << std::endl);
	return;
      }
    }
    set_vertex_handle(k, del_.insert(k));
    CGAL_assertion(vertex_handle(k) != Vertex_handle());
    vertex_handle(k)->set_neighbors(vertex_handle(k)->degree());
    CGAL_DELAUNAY_2_DEBUG(std::cout << "Vertex " << vertex_handle(k)->point() << " has " << vertex_handle(k)->neighbors() << std::endl);
    watcher_.create_vertex(vertex_handle(k));

    // now have to update
    if (!was_2d && del_.dimension()==2) {
      //std::cout << "Creating certificates from scratch.\n";
      if (has_certificates_) {
	has_certificates_=false;
	set_has_certificates(true);
      }
      
    } else if (del_.dimension() == 2 && has_certificates_) {
      Vertex_handle vh= vertex_handle(k);
      Edge_circulator ec= vh->incident_edges(), ef=ec;
      do {
	if (TDS_helper::get_undirected_edge_label(*ec) == Event_key()) {
	  update_edge(*ec);
	}
	++ec;
      } while (ec != ef);
      Face_circulator fc= vh->incident_faces(), fe= fc;
      do {
	int i= fc->index(vh);
	Edge e(fc, i);
	if (TDS_helper::get_undirected_edge_label(e) == Event_key()) {
	  update_edge(e);
	}
	++fc;
      } while (fc != fe);
      
    }
    //write(std::cout);
    //if (del_.dimension()==2) audit();
  }



  Comparison_result compare_concurrent(Event_key a, Event_key b) const {
    Edge ea= traits_.simulator_handle()->template event<Event>(a).edge();
    Edge eb= traits_.simulator_handle()->template event<Event>(b).edge();
    return traits_.compare_concurrent(a, ea, b, eb);
  }




  Edge flip(const Edge &e, Certificate_data cert) {
    ++num_events_;
    CGAL_precondition(!batching_);
    CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME, "DELAUNAY Flipping edge " << TDS_helper::origin(e)->point()
		     << TDS_helper::destination(e)->point() << std::endl);
    //CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_NONE, TDS_helper::destination(e)->point() << std::endl);
    //CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME, " at "  << traits_.simulator()->current_time() << std::endl);

    

    Face_handle face= e.first;
    int index= e.second;
    int mirror_index = face->mirror_index(index);
    Face_handle mirror_face = face->neighbor(index);


    Face_handle bef;
    int bei;
    if (del_.is_edge(TDS_helper::third_vertex(e), TDS_helper::mirror_vertex(e),
		     bef, bei)) {
      // we have a numeric error, lets try to rebuild the neighboring certificates
      CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME,
		       "DELAUNAY ERROR not flipping unflippable edge" << std::endl);
      CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME, 
		       *traits_.simulator_handle() << std::endl;);
      //make this better
      //double ub=to_interval(traits_.simulator_handle()->next_event_time()).second;
      Edge bad_edge(bef, bei);
      Event_key bek= TDS_helper::get_undirected_edge_label(bad_edge);
      
      if (bek == traits_.simulator_handle()->null_event()) {
	CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME,
			 "Dropping the event." << std::endl);
	TDS_helper::set_undirected_edge_label(e, Event_key());
	return e;
      } else {

	double ub = CGAL::to_interval(traits_.simulator_handle()->event_time(bek)).second;
	
	ub= (std::max)(ub+.0000001, 
		       nextafter(ub, (std::numeric_limits<double>::max)()));
	Time t(ub);
	CGAL_precondition(CGAL::compare(t, traits_.simulator_handle()->next_event_time()) == CGAL::LARGER);
	Event_key k =traits_.simulator_handle()->new_event(t, Event(cert, e, this));
	TDS_helper::set_undirected_edge_label(e, k);
	return e;
      }
    }


    TDS_helper::set_undirected_edge_label(e, Event_key());

    for (unsigned int i=0; i<3; ++i) {
      Edge e0(face, i);
      if (TDS_helper::get_undirected_edge_label(e0) != Event_key()) {
	delete_certificate(e0);
      }

      face->vertex(i)->set_neighbors_is_changed(true);
      
      Edge e1(mirror_face, i);
      if (TDS_helper::get_undirected_edge_label(e1) != Event_key()) {
	delete_certificate(e1);
      }

      mirror_face->vertex(i)->set_neighbors_is_changed(true);
    }

   
    watcher_.before_flip(e);
    del_.tds().flip(face,index);

    // we also know that CGAL preserves the edge index of the flipped edge
    mirror_index = mirror_face->index(face);
    index= face->index(mirror_face);

    Edge flipped_edge(face,index);
    Edge mirror_flipped_edge= TDS_helper::mirror_edge(flipped_edge);
    //CGAL_postcondition(del_.is_face(face));
 

    mirror_index = face->mirror_index(index);
    mirror_face = face->neighbor(index);

    /*for (unsigned int i=0; i<3; ++i) {
      Edge e0(face, i);
      update_edge(e0);
      Edge e1(mirror_face, i);
      update_edge(e1);
      }*/
    TDS_helper::set_undirected_edge_label(flipped_edge, traits_.simulator_handle()->null_event());

    update_edge(Edge(flipped_edge.first, (flipped_edge.second+1)%3));
    update_edge(Edge(flipped_edge.first, (flipped_edge.second+2)%3));
    update_edge(Edge(mirror_flipped_edge.first, (mirror_flipped_edge.second+1)%3));
    update_edge(Edge(mirror_flipped_edge.first, (mirror_flipped_edge.second+2)%3));



   {
      std::pair<Time, Certificate_data> sp= traits_.certificate_failure_time(e,cert);
      Event_key k =traits_.simulator_handle()->new_event(sp.first,
							 Event(sp.second, flipped_edge, this));
      TDS_helper::set_undirected_edge_label(flipped_edge, k);
    }
    //write(std::cout);
    //new_edges_.clear();
    //new_edges_.insert(flipped_edge);

    CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME, "Created " << TDS_helper::origin(flipped_edge)->point());
    CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME, TDS_helper::destination(flipped_edge)->point() << std::endl);

    watcher_.after_flip(flipped_edge);
    return flipped_edge;
  }

  Visitor &visitor() {
    return watcher_;
  }

  const Visitor &visitor() const
  {
    return watcher_;
  }

protected:
  Traits traits_;
  Visitor watcher_;
  Triangulation del_;


  Simulator_listener siml_;
  Moving_point_table_listener motl_; 
  std::vector<Vertex_handle> vhs_;

  bool has_certificates_;
  bool batching_;
  std::vector<Edge> batched_certs_;


  mutable unsigned int num_events_;
  //mutable unsigned int num_single_certificates_;

  const typename Traits::Point_2& point(Point_key k) const
  {
    return traits_.point(k);
  }


  Vertex_handle vertex_handle(Point_key k) const {
    if (k.to_index() >= vhs_.size()) return Vertex_handle();
    CGAL_precondition(k.to_index() < vhs_.size());
    return vhs_[k.to_index()];
  }
  void set_vertex_handle(Point_key k, Vertex_handle vh) {
    vhs_.resize(std::max(k.to_index()+1, vhs_.size()));
    vhs_[k.to_index()]=vh;
  }

  void update_vertex(Vertex_handle vh) {
    if (!vh->neighbors_is_changed()) {

    } else {
      int deg= vh->degree();
      if (deg ==3 && vh->neighbors() != 3) {
	CGAL_DELAUNAY_2_DEBUG(std::cout << "Degree 3 for " 
			      << vh->point() << std::endl);
	vh->set_neighbors(deg);
	typename Triangulation::Edge_circulator ec= del_.incident_edges(vh);
	do {
	  if (TDS_helper::get_undirected_edge_label(*ec) != Event_key()) {
	    delete_certificate(*ec);
	  }
	  ++ec;
	} while (ec != del_.incident_edges(vh));
      } else if (vh->neighbors()==3 && deg != 3) {
	CGAL_DELAUNAY_2_DEBUG(std::cout << "Degree 4 for " 
			      << vh->point() << std::endl);
	vh->set_neighbors(deg);
	typename Triangulation::Edge_circulator ec= del_.incident_edges(vh);
	do {
	  if (TDS_helper::get_undirected_edge_label(*ec) == Event_key()) {
	    // check other vertex, it it is not changed and not 3, build
	    Vertex_handle ov= ec->first->vertex((ec->second+1)%3);
	    if (ov== vh) {
	      ov= ec->first->vertex((ec->second+2)%3);
	    }
	    if (!ov->neighbors_is_changed() && ov->neighbors() != 3){
	      new_certificate(*ec);
	      CGAL_DELAUNAY_2_DEBUG(std::cout << "New cert for "
				    << TDS_helper::origin(*ec)->point() << " " 
				    << TDS_helper::destination(*ec)->point() 
				    << std::endl);
	    } else {
	      CGAL_DELAUNAY_2_DEBUG(std::cout << "Not creating cert for "
				    << TDS_helper::origin(*ec)->point() << " " 
				    << TDS_helper::destination(*ec)->point() 
				    << std::endl);
	    }
	  }
	  ++ec;
	} while (ec != del_.incident_edges(vh));
      }  else {
	 vh->set_neighbors(deg);
      }

      CGAL_DELAUNAY_2_DEBUG( std::cout << "Vertex " << vh->point() << " has " << vh->neighbors() << std::endl);
    }
  }


 
  void update_edge(const Edge &e) {
    if (batching_) {
      delete_certificate(e);
      batched_certs_.push_back(canonicalize(e));
    } else {
      CGAL_DELAUNAY_2_DEBUG(std::cout << "Updating edge " 
			    << TDS_helper::origin(e)->point() << " " 
			    << TDS_helper::destination(e)->point() 
			    << std::endl);
      Vertex_handle ov=TDS_helper::origin(e);
      Vertex_handle dv=TDS_helper::destination(e);
      
      update_vertex(ov);
      update_vertex(dv);
      
      CGAL_assertion(!ov->neighbors_is_changed());
      CGAL_assertion(!dv->neighbors_is_changed());
      
      if (TDS_helper::get_undirected_edge_label(e) != Event_key()) {
	CGAL_DELAUNAY_2_DEBUG(std::cout << "Already has event " << std::endl);
	// can't do this since I create all edges around vertex of degree 4 at once
	// CGAL_assertion(0);
      } else if (ov->neighbors() ==3 
		 || dv->neighbors() ==3) {
	CGAL_DELAUNAY_2_DEBUG(std::cout << "One end has 3 " << std::endl);
      } else {
	CGAL_DELAUNAY_2_DEBUG(std::cout << "New certificate" << std::endl);
	new_certificate(e);
      }
    }
  }


  void new_certificate(Edge e) {
   CGAL_precondition(TDS_helper::get_undirected_edge_label(e) == Event_key());
    CGAL_DELAUNAY_2_DEBUG(std::cout << "Making certificate for " << TDS_helper::origin(e)->point() << " " 
			  << TDS_helper::destination(e)->point() 
			  << " which would make " << TDS_helper::mirror_vertex(e)->point() << " " 
			  << TDS_helper::third_vertex(e)->point()
			  << std::endl);
    std::pair<Time, Certificate_data> sp= traits_.certificate_failure_time(e);
    Event_key k= traits_.simulator_handle()->new_event(sp.first, Event(sp.second, e, this));
    TDS_helper::set_undirected_edge_label(e, k);
  }

  void delete_certificate(Edge e) {
    Event_key k=  TDS_helper::get_undirected_edge_label(e);
    traits_.simulator_handle()->delete_event(k);
    TDS_helper::set_undirected_edge_label(e, Event_key());
  }


  Edge canonicalize(Edge e) const {
    if (e.first->neighbor(e.second) < e.first) {
      return TDS_helper::mirror_edge(e);
    } else {
      return e;
    }
  }
};

template <class Sim, class Del, class W, class T>
std::ostream &operator<<(std::ostream &out, const Delaunay_triangulation_2<Sim, Del, W, T> &kd)
{
  kd.write(out);
  return out;
}


CGAL_KINETIC_END_NAMESPACE
#endif
