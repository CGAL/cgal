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
#include <CGAL/Kinetic/Active_objects_listener_helper.h>
#include <CGAL/Kinetic/Simulator_kds_listener.h>
#include <CGAL/Kinetic/internal/tds_2_helpers.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Kinetic/Ref_counted.h>
#include <iterator>
#include <map>
#include <set>

CGAL_KINETIC_BEGIN_NAMESPACE

#define CGAL_DELAUNAY_2_DEBUG(x) 

template <class KDel>
struct Delaunay_2_event_base {
  Delaunay_2_event_base(const typename KDel::Certificate &c,
			const typename KDel::Edge &e,
			KDel *kdel): c_(c), e_(e), d_(kdel) {}
  typename KDel::Certificate &certificate() const {
    return c_;
  }
  const typename KDel::Edge edge() const {
    return e_;
  }
  KDel* kdel() {
    return d_;
  }
  void process() {
    d_->flip(e_, c_);
  }
  typename KDel::Certificate c_;
  const typename KDel::Edge e_;
  KDel *d_;
};

template <class K>
class Delaunay_edge_failure_event: public Delaunay_2_event_base<K>
{
  typedef Delaunay_2_event_base<K> P;
public:
  //! Make sure that the s has been advanced
  Delaunay_edge_failure_event(typename K::SOC_certificate &s,
			      const typename K::Edge &e,
			      K *kdel): P(s, e, kdel) {
  }

  void write(std::ostream &out) const
  {
    out << "Flip " << K::TDS_helper::origin(P::edge())->point() << ","
	<< K::TDS_helper::destination(P::edge())->point();
  }
};

template <class K>
class Delaunay_hull_edge_failure_event: public Delaunay_2_event_base<K> 
{
  typedef Delaunay_2_event_base<K> P;
public:
  //! Make sure that the s has been advanced
  Delaunay_hull_edge_failure_event(typename K::O2_certificate &s,
				   const typename K::Edge &e,
				   K *kdel): P(s,e,kdel) {
    //if (!s.empty()) s_.pop();
  }
  
  void write(std::ostream &out) const
  {
    out << "Flip " << K::TDS_helper::origin(P::edge())->point() << ","
	<< K::TDS_helper::destination(P::edge())->point();
  }
};

template <class T>
std::ostream& operator<<(std::ostream &out, const Delaunay_edge_failure_event<T> &e)
{
  e.write(out);
  return out;
}


template <class T>
std::ostream& operator<<(std::ostream &out, const Delaunay_hull_edge_failure_event<T> &e)
{
  e.write(out);
  return out;
}


//! A 2D kinetic Delaunay triangulation.
/*!  Points are added via the Moving_point_table, so the public
  interface is very limited. See kinetic_Delaunay_2.cc for a useage example.
*/
template <class Simulation_traits_t, class Visitor= Delaunay_triangulation_visitor_base_2,
	  class Delaunay
	  = CGAL::Delaunay_triangulation_2<typename Simulation_traits_t::Instantaneous_kernel,
					   CGAL::Triangulation_data_structure_2<
	    Delaunay_triangulation_vertex_base_2<typename Simulation_traits_t::Instantaneous_kernel>,
	    CGAL::Kinetic::Delaunay_triangulation_face_base_2<Simulation_traits_t > > > >
class Delaunay_triangulation_2:
  public Ref_counted<Delaunay_triangulation_2<Simulation_traits_t, Visitor, Delaunay> >
{
public:
  typedef Simulation_traits_t Simulation_traits;
  typedef Delaunay_triangulation_2<Simulation_traits, Visitor, Delaunay> This;

  typedef typename Simulation_traits::Kinetic_kernel Kinetic_kernel;
  typedef typename Simulation_traits::Simulator Simulator;
  typedef typename Simulation_traits::Active_points_2_table Moving_point_table;

  typedef typename Moving_point_table::Key Point_key;
  typedef typename Simulator::Event_key Event_key;
  //typedef typename Simulator::Root_stack Root_stack;
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
  typedef typename Kinetic_kernel::Positive_side_of_oriented_circle_2 SOC;
  typedef typename Kinetic_kernel::Positive_orientation_2 O2;
  typedef typename SOC::result_type SOC_certificate;
  typedef typename O2::result_type O2_certificate;
  typedef O2_certificate Certificate;
  typedef internal::Triangulation_data_structure_helper_2<typename Delaunay::Triangulation_data_structure> TDS_helper;

  typedef Delaunay_edge_failure_event<This> SOC_event;
  typedef Delaunay_hull_edge_failure_event<This> O2_event;
  
  friend class Delaunay_edge_failure_event<This>;
  friend class Delaunay_hull_edge_failure_event<This>;

  typedef typename CGAL::Kinetic::Simulator_kds_listener<typename Simulator::Listener, This> Simulator_listener;
  friend  class CGAL::Kinetic::Simulator_kds_listener<typename Simulator::Listener, This>;
  typedef typename CGAL::Kinetic::Active_objects_listener_helper<typename Moving_point_table::Listener, This> Moving_point_table_listener;
  friend class CGAL::Kinetic::Active_objects_listener_helper<typename Moving_point_table::Listener, This>;
public:

  Delaunay_triangulation_2(Simulation_traits st,
			   Visitor w= Visitor()): traits_(st),
						  del_(traits_.instantaneous_kernel_object()),
						  soc_(traits_.kinetic_kernel_object().positive_side_of_oriented_circle_2_object()),
						  o2_(traits_.kinetic_kernel_object().positive_orientation_2_object()),
						  watcher_(w),
						  has_certificates_(true){
    siml_ = Simulator_listener(st.simulator_handle(), this);
    motl_= Moving_point_table_listener(st.active_points_2_table_handle(), this);
  }

  //! Just write the objects in order;
  void write(std::ostream &out) const
  {
    out << del_;
  }

  typedef Delaunay Triangulation;
  const Triangulation &triangulation(const typename Simulator::NT &t) const
  {
    //update_instantaneous_kernel_time();
    del_.geom_traits().set_time(t);
    return del_;
  }

  const Simulation_traits& simulation_traits_object() const {
    return traits_;
  }

  typedef typename Delaunay::Triangulation_data_structure Triangulation_data_structure;
  const Triangulation_data_structure &triangulation_data_structure() const
  {
    return del_.tds();
  }

  /*const std::set<Edge>& recent_edges() const {
    return new_edges_;
    }*/

  //! Verify that the current state of the
  void audit() const
  {
    if (!has_certificates_) return;
    CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME, "Auditing delaunay" << std::endl);
    
    CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_LOTS, *this);
    if (del_.dimension() != 2) return;
    Delaunay sdel(traits_.instantaneous_kernel_object());
    sdel.geom_traits().set_time(traits_.simulator_handle()->rational_current_time());
    for (typename Delaunay::Finite_vertices_iterator cit= del_.finite_vertices_begin();
	 cit != del_.finite_vertices_end(); ++cit){
      sdel.insert(cit->point());
    }
    /*    sdel.insert(traits_.active_points_2_table_handle()->keys_begin(),
	  traits_.active_points_2_table_handle()->keys_end());*/

    CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_LOTS, sdel << std::endl);

    if (del_.dimension() != sdel.dimension()) {
      CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_NONE, "AUDIT FAILURE Dimensions don't match in audit" << std::endl);
      return;
    }
    CGAL_exactness_assertion(del_.dimension() == sdel.dimension());
   
   for (typename Delaunay::All_vertices_iterator vit = del_.all_vertices_begin();
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
    for (typename Delaunay::All_vertices_iterator vit = sdel.all_vertices_begin();
	 vit != sdel.all_vertices_end(); ++vit) {
      bool found=false;
    
      //Object_key k= vit->point();
      for (typename Delaunay::All_vertices_iterator vit2= del_.all_vertices_begin();
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

    del_.geom_traits().set_time(traits_.simulator_handle()->rational_current_time());
    typename Simulation_traits::Instantaneous_kernel::Side_of_oriented_circle_2 ic2
      = del_.geom_traits().side_of_oriented_circle_2_object();
    for (typename Delaunay::Finite_edges_iterator fit = del_.finite_edges_begin(); 
	 fit != del_.finite_edges_end(); ++fit){
      Point_key k0= fit->first->vertex((fit->second+1)%3)->point();
      Point_key k2= fit->first->vertex((fit->second+2)%3)->point();
      Point_key k3= TDS_helper::mirror_vertex(*fit)->point();
      Point_key k1= TDS_helper::third_vertex(*fit)->point();
      if (k1== Point_key() || k3== Point_key()) continue;
      typename Simulation_traits::Static_kernel::Point_2 p0= del_.geom_traits().static_object(k0);
      typename Simulation_traits::Static_kernel::Point_2 p1= del_.geom_traits().static_object(k1);
      typename Simulation_traits::Static_kernel::Point_2 p2= del_.geom_traits().static_object(k2);
      typename Simulation_traits::Static_kernel::Point_2 p3= del_.geom_traits().static_object(k3);
      if (ic2(k0, k1, k2, k3) != CGAL::ON_POSITIVE_SIDE) {
	CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_NONE, "AUDIT FAILURE Failed certificate: " << k0 << " " << k1 << " " 
			 << k2 << " " << k3 << std::endl);
	CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_NONE, "AUDIT FAILURE Points are: " << p0 << ": " << p1 << ": " << p2 
			 << ": " << p3 << std::endl);
      }
      CGAL_exactness_assertion(ic2(k0, k1, k2, k3) == CGAL::ON_POSITIVE_SIDE);
      
    }

    for (typename Delaunay::Edge_iterator fit = del_.edges_begin(); fit != del_.edges_end(); ++fit){
      if (TDS_helper::origin(*fit)->degree()==3 || TDS_helper::destination(*fit)->degree()==3
	  /*|| del_.is_edge(TDS_helper::third_vertex(*fit), TDS_helper::mirror_vertex(*fit))*/) {
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
      if (tf==true) {
	CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME, "DELAUNAY2: Creating certificates."<< std::endl);
	for (typename Delaunay::All_vertices_iterator vit = del_.all_vertices_begin(); 
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
      } else {

	for (typename Delaunay::Edge_iterator it = del_.edges_begin();
	     it != del_.edges_end(); ++it){
	  if (TDS_helper::get_undirected_edge_label(*it) != Event_key()){
	    traits_.simulator_handle()->delete_event(TDS_helper::get_undirected_edge_label(*it));
	    TDS_helper::set_undirected_edge_label(*it, Event_key());
	  }
	}
      }
      has_certificates_=tf;
    }
  }
  bool has_certificates() {
    return has_certificates_;
  }

  void erase(Point_key k) {
    // erase all incident certificates
    Vertex_handle vh= vhs_[k];
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
	    traits_.simulator_handle()->delete_event(k);
	    TDS_helper::set_undirected_edge_label(e, Event_key());
	  }
	  fc->vertex(j)->set_neighbors_is_changed(true);
	  }
	  ++fc;
	} while (fc != fe);
      }
    }
    // remove from triangulation
    del_.geom_traits().set_time(traits_.simulator_handle()->rational_current_time());
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

    if (del_.dimension() != 2) {
      CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME,"Triangulation is still 1D.\n");
      return;
    }

    Vertex_handle vh= vhs_[k];
    if (vh == Vertex_handle()) {
      CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME, "Point " << k << " is not in triangulation on set."<< std::endl);
      return;
    }
    if (has_certificates_) {
      Edge_circulator ec= vh->incident_edges(), ef=ec;
      if (ec != NULL) {
	do {
	  if (TDS_helper::get_undirected_edge_label(*ec) != Event_key()) {
	    traits_.simulator_handle()->delete_event(TDS_helper::get_undirected_edge_label(*ec));
	    TDS_helper::set_undirected_edge_label(*ec, Event_key());
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
	    traits_.simulator_handle()->delete_event(TDS_helper::get_undirected_edge_label(e));
	    TDS_helper::set_undirected_edge_label(e, Event_key());
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
    CGAL_DELAUNAY_2_DEBUG(std::cout << "Inserting " << k << std::endl);
    bool was_2d= (del_.dimension()==2);

    del_.geom_traits().set_time(traits_.simulator_handle()->rational_current_time());
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
	    traits_.simulator_handle()->delete_event(k);
	    TDS_helper::set_undirected_edge_label(e, Event_key());
	  }
	}
      }
      watcher_.remove_faces(faces.begin(), faces.end());

      if (faces.empty()) {
	CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME, "DELAUNAY vertex not successfully inserted " << k << std::endl);
	return;
      }
    }

    vhs_[k]= del_.insert(k);
    CGAL_assertion(vhs_[k] != Vertex_handle());
    vhs_[k]->set_neighbors(vhs_[k]->degree());
    CGAL_DELAUNAY_2_DEBUG(std::cout << "Vertex " << vhs_[k]->point() << " has " << vhs_[k]->neighbors() << std::endl);
    watcher_.create_vertex(vhs_[k]);

    // now have to update
    if (!was_2d && del_.dimension()==2) {
      //std::cout << "Creating certificates from scratch.\n";
      if (has_certificates_) {
	has_certificates_=false;
	set_has_certificates(true);
      }
      
    } else if (del_.dimension() == 2 && has_certificates_) {
      Vertex_handle vh= vhs_[k];
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
    if (del_.dimension()==2) audit();
  }

  Edge flip(const Edge &e, Certificate cert) {
    CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME, "DELAUNAY Flipping edge " << TDS_helper::origin(e)->point()
		     << TDS_helper::destination(e)->point() << std::endl);
    //CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_NONE, TDS_helper::destination(e)->point() << std::endl);
    //CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME, " at "  << traits_.simulator()->current_time() << std::endl);

    

    Face_handle face= e.first;
    int index= e.second;
    int mirror_index = face->mirror_index(index);
    Face_handle mirror_face = face->neighbor(index);

    if (del_.is_edge(TDS_helper::third_vertex(e), TDS_helper::mirror_vertex(e))) {
      // we have a numeric error, lets try to rebuild the neighboring certificates
      CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_SOME, "DELAUNAY ERROR not flipping unflippable edge" << std::endl);
      //make this better
      double ub=to_interval(traits_.simulator_handle()->next_event_time()).second;
      ub= (std::max)(ub+.0000001, 
		     nextafter(ub, (std::numeric_limits<double>::max)()));
      Time t=ub;
      //cert.pop_failure_time();
      Event_key k =traits_.simulator_handle()->new_event(t, O2_event(cert, e, this));
      TDS_helper::set_undirected_edge_label(e, k);
      return e;
    }

    for (unsigned int i=0; i<3; ++i) {
      Edge e0(face, i);
      if (e0.second != index && TDS_helper::get_undirected_edge_label(e0) != Event_key()) {
	traits_.simulator_handle()->delete_event(TDS_helper::get_undirected_edge_label(e0));
	TDS_helper::set_undirected_edge_label(e0, Event_key());
      }

      face->vertex(i)->set_neighbors_is_changed(true);
      
      Edge e1(mirror_face, i);
      if (e1.second != mirror_index  && TDS_helper::get_undirected_edge_label(e1) != Event_key()) {
	traits_.simulator_handle()->delete_event(TDS_helper::get_undirected_edge_label(e1));
	TDS_helper::set_undirected_edge_label(e1, Event_key());
      }

      mirror_face->vertex(i)->set_neighbors_is_changed(true);
    }

    TDS_helper::set_undirected_edge_label(e, Event_key());
    watcher_.before_flip(e);
    del_.tds().flip(face,index);

    // we also know that CGAL preserves the edge index of the flipped edge
    mirror_index = mirror_face->index(face);
    index= face->index(mirror_face);

    Edge flipped_edge(face,index);
    //CGAL_postcondition(del_.is_face(face));
    {
      Time t= cert.failure_time();
      cert.pop_failure_time();
      Event_key k =traits_.simulator_handle()->new_event(t, O2_event(cert, flipped_edge, this));
      TDS_helper::set_undirected_edge_label(flipped_edge, k);
    }

    mirror_index = face->mirror_index(index);
    mirror_face = face->neighbor(index);

    for (unsigned int i=0; i<3; ++i) {
      Edge e0(face, i);
      update_edge(e0);
      Edge e1(mirror_face, i);
      update_edge(e1);
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
  Simulation_traits traits_;
  Simulator_listener siml_;
  Moving_point_table_listener motl_;
  Delaunay del_;
  std::map<Point_key, Vertex_handle> vhs_;
  SOC soc_;
  O2 o2_;
  //std::set<Edge> new_edges_;
  Visitor watcher_;
  bool has_certificates_;

  const typename Moving_point_table::Data& point(Point_key k) const
  {
    return traits_.active_points_2_table_handle()->at(k);
  }


  /*void add_certificate_around_vertex(Face_handle f0, Face_handle f1) {
    int i0= (e.second+1)%3;
    Face_handle n0= ce.first->neighbor(i0);
    Vertex_handle v0=ce.first->vertex((e.second+2)%3);
    Edge e0(n0,n0->index(v0));
    
    int i1= (e.second+2)%3;
    Face_handle n1= ce.first->neighbor(i1);
    Vertex_handle v1=ce.first->vertex((e.second+1)%3);
    Edge e1(n1, n1->index(v1));
    CGAL_assertion(v0 != v1);
    CGAL_assertion(n0 != n1);
    
    if (TDS_helper::mirror_edge(e1) == e0 && !TDS_helper::get_undirected_edge_label(e1).is_valid()) {
      CGAL_KINETIC_LOG(CGAL::Kinetic::LOG_LOTS, "DELAUNAY adding certificate to edge " << 
		       TDS_helper::origin(e1)->point() << " " << TDS_helper::destination(e1)->point() << std::endl;);
      new_certificate(e0);
    } else {
      std::cout << "Not adding certificate to edge " 
		<< TDS_helper::origin(e1)->point() << " " << TDS_helper::destination(e1)->point()
		<< " or " << TDS_helper::origin(e0)->point() << " " << TDS_helper::destination(e0)->point();
      if (TDS_helper::get_undirected_edge_label(e1).is_valid()) {
	std::cout << " because it is there";
      }
      std::cout  << std::endl;
    }
    }*/

  bool is_hull_edge(const Edge &e) const {
    return ! TDS_helper::mirror_vertex(e)->point().is_valid()
      || ! TDS_helper::third_vertex(e)->point().is_valid()
      || ! TDS_helper::origin(e)->point().is_valid()
      || ! TDS_helper::destination(e)->point().is_valid();
  }

  SOC_certificate compute_failure_time(const Edge &e) const {
    Point_key ks[4];
    ks[0]= TDS_helper::origin(e)->point();
    ks[1]= TDS_helper::third_vertex(e)->point();
    ks[2]= TDS_helper::destination(e)->point();
    ks[3]= TDS_helper::mirror_vertex(e)->point();
      
      //bool odd_parity=false;
      //bool infinity=false;
      
      SOC_certificate s=soc_(point(ks[0]), point(ks[1]),
			     point(ks[2]), point(ks[3]),
			     traits_.simulator_handle()->current_time(),
			     traits_.simulator_handle()->end_time());
      
      return s;
  }

  O2_certificate compute_hull_failure_time(const Edge &e) const {
    Point_key ks[4];
    ks[0]= TDS_helper::origin(e)->point();
    ks[1]= TDS_helper::third_vertex(e)->point();
    ks[2]= TDS_helper::destination(e)->point();
    ks[3]= TDS_helper::mirror_vertex(e)->point();

      bool odd_parity=false;
      bool infinity=false;
      for (unsigned int i=0; i<4; ++i) {
	if (infinity) {
	  ks[i-1]=ks[i];
	}
	else {
	  if (!ks[i].is_valid()) {
	    infinity=true;
	    odd_parity= ((i%2)==1);
	  }
	}
      }
      if (odd_parity) {
	std::swap(ks[0], ks[1]);
      }
      O2_certificate s=o2_(point(ks[0]), point(ks[1]), point(ks[2]),
			   traits_.simulator_handle()->current_time(),
			   traits_.simulator_handle()->end_time());
      return s;
  }

  void update_vertex(Vertex_handle vh) {
    if (!vh->neighbors_is_changed()) {

    } else {
      int deg= vh->degree();
      if (deg ==3 && vh->neighbors() != 3) {
	CGAL_DELAUNAY_2_DEBUG(std::cout << "Degree 3 for " << vh->point() << std::endl);
	vh->set_neighbors(deg);
	typename Delaunay::Edge_circulator ec= del_.incident_edges(vh);
	do {
	  if (TDS_helper::get_undirected_edge_label(*ec) != Event_key()) {
	    traits_.simulator_handle()->delete_event(TDS_helper::get_undirected_edge_label(*ec));
	  }
	  TDS_helper::set_undirected_edge_label(*ec, Event_key());
	  ++ec;
	} while (ec != del_.incident_edges(vh));
      } else if (vh->neighbors()==3 && deg != 3) {
	CGAL_DELAUNAY_2_DEBUG(std::cout << "Degree 4 for " << vh->point() << std::endl);
	vh->set_neighbors(deg);
	typename Delaunay::Edge_circulator ec= del_.incident_edges(vh);
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


  void new_certificate( const Edge &e) {
    CGAL_precondition(TDS_helper::get_undirected_edge_label(e) == Event_key());
    CGAL_DELAUNAY_2_DEBUG(std::cout << "Making certificate for " << TDS_helper::origin(e)->point() << " " 
			  << TDS_helper::destination(e)->point() 
			  << " which would make " << TDS_helper::mirror_vertex(e)->point() << " " 
			  << TDS_helper::third_vertex(e)->point()
			  << std::endl);
      
    /*if (del_.is_edge(TDS_helper::mirror_vertex(e),
		     TDS_helper::third_vertex(e))) {
      std::cout << "Skipping" << std::endl;
      return;
      }*/
    Event_key k;
    if (static_cast<bool>(is_hull_edge(e))) {
      O2_certificate s= compute_hull_failure_time(e);
      Time t= s.failure_time();
      s.pop_failure_time();
      k =traits_.simulator_handle()->new_event(t, O2_event(s, e, this));
    } else {
      SOC_certificate s= compute_failure_time(e);
      Time t= s.failure_time();
      s.pop_failure_time();
      k =traits_.simulator_handle()->new_event(t, O2_event(s, e, this));
    }
    TDS_helper::set_undirected_edge_label(e, k);
  }

  void update_edge(const Edge &e) {
    CGAL_DELAUNAY_2_DEBUG(std::cout << "Updating edge edge " 
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
    } else if (ov->neighbors() ==3 
	       || dv->neighbors() ==3) {
      CGAL_DELAUNAY_2_DEBUG(std::cout << "One end has 3 " << std::endl);
    } else {
      CGAL_DELAUNAY_2_DEBUG(std::cout << "New certificate" << std::endl);
      new_certificate(e);
    }
  }

  //! rebuild a certificates
  /*!  I need to check if there is a valid one before since I use
    change_object to initialize the certificates of a new object.
  */
  /*void rebuild_certificate( const Edge &e) {
    if (TDS_helper::get_undirected_edge_label(e).is_valid()) {
      traits_.simulator_handle()->delete_event(TDS_helper::get_undirected_edge_label(e));
      TDS_helper::set_undirected_edge_label(e,  Event_key());
      compute_certificate(e);
    } else {
      std::cout << "Not rebuilding for edge " << TDS_helper::origin(e)->point() << " " << TDS_helper::destination(e)->point() 
		<< std::endl;
    }
    }*/
};

template <class Sim, class Del, class W>
std::ostream &operator<<(std::ostream &out, const Delaunay_triangulation_2<Sim, Del, W> &kd)
{
  kd.write(out);
  return out;
}


CGAL_KINETIC_END_NAMESPACE
#endif
