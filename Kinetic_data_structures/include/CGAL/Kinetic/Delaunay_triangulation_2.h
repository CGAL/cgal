// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
#include <CGAL/Kinetic/listeners.h>
#include <CGAL/Kinetic/internal/tds_2_helpers.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Kinetic/Ref_counted.h>
#include <iterator>
#include <CGAL/Kinetic/Event_base.h>
#include <CGAL/Kinetic/Delaunay_triangulation_default_traits_2.h>

namespace CGAL { namespace Kinetic {
//#ifdef NDEBUG
#define CGAL_DELAUNAY_2_DEBUG(x)
/*#else
#define CGAL_DELAUNAY_2_DEBUG(x) x
//#endif*/

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

  void audit(typename KDel::Event_key k) const {
    kdel()->audit_event(k, e_);
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
	<< ", " << kdel()->triangulation_data_structure().mirror_vertex(edge().first, 
                                                                        edge().second)->point() ;
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

  /*typedef typename CGAL::Kinetic::Simulator_kds_listener<typename Simulator::Listener, This> Simulator_listener;
  friend  class CGAL::Kinetic::Simulator_kds_listener<typename Simulator::Listener, This>;
  typedef typename CGAL::Kinetic::Active_objects_batch_listener_helper<typename Moving_point_table::Listener, This> Moving_point_table_listener;
  friend class CGAL::Kinetic::Active_objects_batch_listener_helper<typename Moving_point_table::Listener, This>;*/
  CGAL_KINETIC_DECLARE_BATCH_LISTENERS(typename Simulator,
				       typename Moving_point_table)

  /*struct Compare_edges{
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
    };*/

  void init_data(bool insert) {
    /*siml_ = Simulator_listener(traits_.simulator_handle(), this);
      motl_= Moving_point_table_listener(traits_.active_points_2_table_handle(), this, insert);*/
    CGAL_KINETIC_INITIALIZE_BATCH_LISTENERS(traits_.simulator_handle(),
					    traits_.active_points_2_table_handle(),
					    insert);
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
      CGAL_assertion(vit->point().index() < del_.number_of_vertices());
      vhs_[vit->point().index()]=vit;
    }
    init_data(false);
  
    set_has_certificates(false, 0);
  }

  Delaunay_triangulation_2(Simulation_traits st,
			   Visitor w= Visitor()): 
    traits_(st),
    watcher_(w),
    del_(traits_.instantaneous_kernel_object()) {
    init_data(true);
    set_has_certificates(true, 0);
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
    if (tf) {
      batching_=true;
    } else if (batching_) {

      //unsigned int num_certs= num_certificates_;
      // this is important that it be before update
      batching_=false;

      //std::sort(batched_certs_.begin(), batched_certs_.end(), Compare_edges());
     
      /*batched_certs_.erase(std::unique(batched_certs_.begin(), 
				       batched_certs_.end()), 
				       batched_certs_.end());*/
     
      
      for (unsigned int i=0; i< batched_certs_.size(); ++i){
	//Point_key s=TDS_helper::origin(batched_certs_[i])->point();
	//Point_key t=TDS_helper::destination(batched_certs_[i])->point();
	
	//std::cout << std::min(s,t) << "--" << std::max(s,t) << std::endl;
	update_edge(batched_certs_[i]);
      }
      
      batched_certs_.clear();
      /*CGAL_LOG(CGAL::Kinetic::Log::SOME, 
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
  void audit()const;

  void audit_event(Event_key k, Edge e) const {
    if (get_undirected_edge_label(e) !=k) {
      std::cerr << "AUDIT FAILURE orphan event " << k << std::endl;
    }
    CGAL_assertion(get_undirected_edge_label(e) ==k);
  }

  void set_neighbors_initialized(bool tf) const {
    if (tf) {
      for (typename Triangulation::All_vertices_iterator vit = del_.all_vertices_begin(); 
	   vit != del_.all_vertices_end(); ++vit) {
	vit->set_neighbors(0);
      }
      for (Face_iterator f = del_.all_faces_begin(); f != del_.all_faces_end(); ++f) {
	f->vertex(0)->set_neighbors(f->vertex(0)->neighbors()+1);
	f->vertex(1)->set_neighbors(f->vertex(1)->neighbors()+1);
	f->vertex(2)->set_neighbors(f->vertex(2)->neighbors()+1);
      }
    }
  }

  enum New_certificate_state {HAD_NO_FAILURES=1, HAS_NO_FAILURES = 2, NO_STRUCTURE_CHANGES = 4};

  void set_has_certificates(bool tf, int state=0) {
    if (tf == has_certificates_){

    } else {
      if (tf==true && del_.dimension()==2) {
	CGAL_LOG(CGAL::Log::SOME, "DELAUNAY2: Creating certificates."<< std::endl);
	if (!(state & NO_STRUCTURE_CHANGES)) {
	  set_neighbors_initialized(true);
	}
	if (state & HAS_NO_FAILURES && !(state & HAD_NO_FAILURES)) {
	  for (Face_iterator f = del_.all_faces_begin(); f != del_.all_faces_end(); ++f) {
	    
	    set_directed_edge_label(Edge(f,0), traits_.simulator_handle()->null_event());
	    set_directed_edge_label(Edge(f,1), traits_.simulator_handle()->null_event());
	    set_directed_edge_label(Edge(f,2), traits_.simulator_handle()->null_event());
	    if (f->vertex(0)->neighbors() == 3) {
	      set_directed_edge_label(Edge(f,1), Event_key()); 
	      set_directed_edge_label(Edge(f,2), Event_key()); 
	    }
	    if (f->vertex(1)->neighbors() == 3) {
	      set_directed_edge_label(Edge(f,0), Event_key()); 
	      set_directed_edge_label(Edge(f,2), Event_key()); 
	    }
	    if (f->vertex(2)->neighbors() == 3) {
	      set_directed_edge_label(Edge(f,1), Event_key()); 
	      set_directed_edge_label(Edge(f,0), Event_key()); 
	    }
	    watcher_.create_face(f);
	  }
	} else if (state & HAD_NO_FAILURES  && state & HAS_NO_FAILURES)  {
	  // nothing to do, yeah!!!
	} else {
	  for (Edge_iterator eit = del_.all_edges_begin(); eit != del_.all_edges_end(); ++eit) {
	    set_undirected_edge_label(*eit, Event_key());
	    update_edge(*eit);
	  }
	}
      } else if (tf==false) { 
	if (!(state & HAD_NO_FAILURES)) {
	  for (Face_iterator f = del_.all_faces_begin(); f != del_.all_faces_end(); ++f) {
	    for(unsigned int i=0; i< 3; ++i) {
	      Edge e(f,i);
	      delete_certificate(e);
	    }
	    watcher_.destroy_face(f);
	  }
	}
      } 
      CGAL_LOG(Log::SOME, 
	       *traits_.simulator_handle() << std::endl);
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
      CGAL_LOG(Log::SOME, "Point " << k << " is not in triangulation on removal."<< std::endl);
      return;
    }
    watcher_.pre_remove_vertex(vh);
    if (has_certificates_) {
      Face_circulator fc= del_.incident_faces(vh), fe=fc;
      if (fc != NULL) {
	do {
	  for (unsigned int j=0; j<3; ++j) {
	    Edge e(fc, j);
	    Event_key k= get_undirected_edge_label(e);
	    delete_certificate(e);
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
	  update_neighbors(faces[i]->vertex(j));
	}
      }
      for (unsigned int i=0; i< faces.size(); ++i) {
	for (unsigned int j=0; j<3; ++j) {
	  update_vertex(faces[i]->vertex(j));
	}
      }
      for (unsigned int i=0; i< faces.size(); ++i) {
	for (unsigned int j=0; j<3; ++j) {
	  Edge e(faces[i],j);
	  //Event_key k= get_undirected_edge_label(e);
	  // a bit redundant for certificates which don't fail
	  update_edge(e);
	  //new_edges_.insert(e);
	}
	watcher_.create_face(faces[i]);
      }
      
    }
    watcher_.post_remove_vertex(k);
  }

  //! The assertion will catch that the object is in the same sorted order
  void set(Point_key k) {
    //std::cout << "Object changed " << k << std::endl;

    //new_edges_.clear();
    traits_.point_changed(k);
    if (del_.dimension() != 2) {
      CGAL_LOG(Log::SOME,"Triangulation is still 1D.\n");
      return;
    }

    Vertex_handle vh=vertex_handle(k);
    if (vh == Vertex_handle()) {
      CGAL_LOG(Log::SOME, "Point " << k << " is not in triangulation on set."<< std::endl);
      return;
    }

    if (has_certificates_) {
      Face_handle f= vh->face(), fe= vh->face();
      int i= f->index(vh);
      do {
	//int i= fc->index(vh);
	Edge e0(f, i);
	delete_certificate(e0);
	update_edge(e0);
	Edge e1=Edge(f, (i+1)%3);
	delete_certificate(e1);
	update_edge(e1);
	f= f->neighbor((i+1)%3);
	i= f->index(vh);
      } while (f != fe);
    }

   
    watcher_.change_vertex(vh);

    /*if (has_certificates_) {
      Face_handle f= vh->face(), fe= vh->face();
      int i= f->index(vh);
      do {
	//int i= fc->index(vh);
	Edge e0(f, i);
	update_edge(e0);
	Edge e1=Edge(f, (i+1)%3);
	update_edge(e1);
	f= f->neighbor((i+1)%3);
	i= f->index(vh);
      } while (f != fe);
      }*/
    //write(std::cout);
  }


  void insert(Point_key k) {
    // evil hack
    CGAL_precondition(k.index() >= vhs_.size() || vertex_handle(k) == Vertex_handle());
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
	  Edge e(f, j);
	  Event_key k= get_undirected_edge_label(e);
	  delete_certificate(e);
	}
	watcher_.destroy_face(faces[i]);
      }

      if (faces.empty()) {
	CGAL_LOG(Log::SOME, "DELAUNAY vertex not successfully inserted " << k << std::endl);
	return;
      }
    }
    watcher_.pre_insert_vertex(k);
    Vertex_handle vh= del_.insert(k);
    set_vertex_handle(k, vh);
    
    CGAL_assertion(vertex_handle(k) != Vertex_handle());
    
    CGAL_DELAUNAY_2_DEBUG(std::cout << "Vertex " << vertex_handle(k)->point() 
			  << " has " << vertex_handle(k)->neighbors() << std::endl);
   

    // now have to update
    if (has_certificates_) {
   

      if (!was_2d && del_.dimension()==2) {
	vh->set_neighbors(del_.degree(vh));
	has_certificates_=false;
	set_has_certificates(true, 0); 
      } else if (del_.dimension() == 2) {
	update_neighbors(vh);
	update_vertex(vh);
	// update vertices
	{
	  Face_handle f= vh->face(), fe= vh->face();
	  int i= f->index(vh);
	  do {
	    update_neighbors(f->vertex((i+1)%3));
	    f= f->neighbor((i+1)%3);
	    i= f->index(vh);
	  } while (f != fe);
	} 

	{
	  Face_handle f= vh->face(), fe= vh->face();
	  int i= f->index(vh);
	  do {
	    update_vertex(f->vertex((i+1)%3));
	    f= f->neighbor((i+1)%3);
	    i= f->index(vh);
	  } while (f != fe);
	} 
	// update edges 
	{
	  Face_handle f= vh->face(), fe= vh->face();
	  int i= f->index(vh);
	  do {
	    //int i= fc->index(vh);
	    Edge e0(f, i);
	    update_edge(e0);
	    Edge e1=Edge(f, (i+1)%3);
	    update_edge(e1);
	    f= f->neighbor((i+1)%3);
	    i= f->index(vh);
	  } while (f != fe);
	}
      }
    } else {
      vertex_handle(k)->set_neighbors(del_.degree(vh));
    }
    watcher_.post_insert_vertex(vh);
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
    CGAL_LOG(Log::SOME, "\n\n\n\n\n\nDELAUNAY Flipping edge "
		     << TDS_helper::origin(e)->point()
		     << TDS_helper::destination(e)->point() 
		     << " to get " << TDS_helper::third_vertex(e)->point()
		     << ", " << del_.mirror_vertex(e.first, e.second)->point()<< std::endl);
    //CGAL_LOG(Log::NONE, TDS_helper::destination(e)->point() << std::endl);
    //CGAL_LOG(Log::SOME, " at "  << traits_.simulator()->current_time() << std::endl);

    

    Face_handle face= e.first;
    int index= e.second;
    Face_handle mirror_face = face->neighbor(index);
    int mirror_index =face->neighbor(index)->index(face);
    Edge em(mirror_face,mirror_index);
    CGAL_precondition(mirror_face->neighbor(mirror_index) == face);

    Face_handle bef;
    int bei;

    if (!traits_.is_exact()
	&& del_.is_edge(TDS_helper::third_vertex(e), del_.mirror_vertex(e.first, e.second),
			bef, bei)) {
      // we have a numeric error, lets try to rebuild the neighboring certificates
      CGAL_LOG(Log::SOME,
	       "DELAUNAY ERROR not flipping unflippable edge" << std::endl);
      CGAL_LOG(Log::SOME, 
	       *traits_.simulator_handle() << std::endl);
      //make this better
      //double ub=to_interval(traits_.simulator_handle()->next_event_time()).second;
      Edge bad_edge(bef, bei);
      Event_key bek= get_undirected_edge_label(bad_edge);
      
      if (bek == traits_.simulator_handle()->null_event()) {
	CGAL_LOG(Log::SOME,
			 "Dropping the event." << std::endl);
	set_undirected_edge_label(e, Event_key());
	return e;
      } else {

	double ub = CGAL::to_interval(traits_.simulator_handle()->event_time(bek)).second;
	
	ub= (std::max)(ub+.0000001, 
		       nextafter(ub, (std::numeric_limits<double>::max)()));
	Time t(ub);
	CGAL_precondition(CGAL::compare(t, traits_.simulator_handle()->next_event_time()) == CGAL::LARGER);
	Event_key k =traits_.simulator_handle()->new_event(t, Event(cert, e, this));
	set_undirected_edge_label(e, k);
	return e;
      }
    }
    CGAL_precondition(!del_.is_edge(TDS_helper::third_vertex(e), del_.mirror_vertex(e.first, e.second),
				    bef, bei));

    set_directed_edge_label(e, Event_key());
    set_directed_edge_label(em, Event_key());

    delete_certificate(Edge(face, (index+1)%3));
    delete_certificate(Edge(face, (index+2)%3));
    delete_certificate(Edge(mirror_face, (mirror_index+1)%3));
    delete_certificate(Edge(mirror_face, (mirror_index+2)%3));

   
    watcher_.pre_flip(e);
    del_.tds().flip(face,index);
   
    // we also know that CGAL preserves the edge index of the flipped edge
    mirror_index = mirror_face->index(face);
    index= face->index(mirror_face);

    Edge flipped_edge(face,index);
    Edge mirror_flipped_edge(face->neighbor(index), mirror_index);
    CGAL_assertion(mirror_flipped_edge == mirror_edge(flipped_edge));
    //CGAL_postcondition(del_.is_face(face));

    CGAL_assertion(mirror_index == face->neighbor(index)->index(face));
    CGAL_assertion(mirror_face == face->neighbor(index));

    watcher_.post_flip(flipped_edge);

    decrease_neighbors(flipped_edge.first->vertex(flipped_edge.second));
    increase_neighbors(flipped_edge.first->vertex((flipped_edge.second+1)%3));
    increase_neighbors(flipped_edge.first->vertex((flipped_edge.second+2)%3));
    decrease_neighbors(mirror_flipped_edge.first->vertex(mirror_flipped_edge.second));

    // make sure it is not created by a 3-4 transition
    set_directed_edge_label(flipped_edge, traits_.simulator_handle()->null_event());
    set_directed_edge_label(mirror_flipped_edge, traits_.simulator_handle()->null_event());

    update_decreased_vertex(flipped_edge.first->vertex(flipped_edge.second));
    update_increased_vertex(flipped_edge.first->vertex((flipped_edge.second+1)%3));
    update_increased_vertex(flipped_edge.first->vertex((flipped_edge.second+2)%3));
    update_decreased_vertex(mirror_flipped_edge.first->vertex(mirror_flipped_edge.second));

    update_edge(Edge(flipped_edge.first, (flipped_edge.second+1)%3));
    update_edge(Edge(flipped_edge.first, (flipped_edge.second+2)%3));
    update_edge(Edge(mirror_flipped_edge.first, (mirror_flipped_edge.second+1)%3));
    update_edge(Edge(mirror_flipped_edge.first, (mirror_flipped_edge.second+2)%3));

    {

      Time t; Certificate_data cd;
      if (traits_.certificate_failure_time(flipped_edge,cert, t, cd)){
	Event_key k =traits_.simulator_handle()->new_event(t,
							   Event(cd, flipped_edge, this));
	set_directed_edge_label(flipped_edge, k);
	set_directed_edge_label(mirror_flipped_edge, k);
      } else {
	set_directed_edge_label(flipped_edge, traits_.simulator_handle()->null_event());
	set_directed_edge_label(mirror_flipped_edge, traits_.simulator_handle()->null_event());
      }
    }
    //write(std::cout);
    //new_edges_.clear();
    //new_edges_.insert(flipped_edge);

    CGAL_LOG(Log::SOME, "Created " << TDS_helper::origin(flipped_edge)->point());
    CGAL_LOG(Log::SOME, TDS_helper::destination(flipped_edge)->point() << std::endl);

    return flipped_edge;
  }

  Visitor &visitor() {
    return watcher_;
  }

  const Visitor &visitor() const
  {
    return watcher_;
  }


  bool has_event(const Edge &e) const {
    return get_directed_edge_label(e) != Event_key();
  }

  bool has_finite_event(const Edge &e) const {
    return has_event(e) && get_directed_edge_label(e) != traits_.simulator_handle()->null_event();
  }

protected:
  Traits traits_;
  Visitor watcher_;
  Triangulation del_;


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
    //if (k.index() >= vhs_.size()) return Vertex_handle();
    CGAL_precondition(k.index() < vhs_.size());
    return vhs_[k.index()];
  }
  void set_vertex_handle(Point_key k, Vertex_handle vh) {
    vhs_.resize(std::max BOOST_PREVENT_MACRO_SUBSTITUTION(static_cast<unsigned int>(k.index()+1),
			 static_cast<unsigned int>(vhs_.size())));
    vhs_[k.index()]=vh;
  }

  void update_vertex_to_degree_3(Vertex_handle vh) {
    CGAL_DELAUNAY_2_DEBUG(std::cout << "Degree 3 for " 
			  << vh->point() << std::endl);
    typename Triangulation::Edge_circulator ec= del_.incident_edges(vh);
    do {
      delete_certificate(*ec);
      ++ec;
    } while (ec != del_.incident_edges(vh));
  }

  void update_vertex_to_not_degree_3(Vertex_handle vh) {
    CGAL_DELAUNAY_2_DEBUG(std::cout << "Degree 4 for " 
			      << vh->point() << std::endl);
    typename Triangulation::Edge_circulator ec= del_.incident_edges(vh);
    do {
      if (get_undirected_edge_label(*ec) == Event_key()) {
	// check other vertex, it it is not changed and not 3, build
	Vertex_handle ov= ec->first->vertex((ec->second+1)%3);
	if (ov== vh) {
	  ov= ec->first->vertex((ec->second+2)%3);
	}
	
	if (ov->neighbors() != 3){
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
  }

  void update_neighbors(Vertex_handle vh) {
    unsigned int deg= del_.degree(vh);
    if (deg == vh->neighbors()) return;
    if (deg ==3) {
      vh->set_neighbors(3);
      //update_vertex_to_degree_3(vh);
    } else if (vh->neighbors()==3) {
      vh->set_neighbors(deg);
      //update_vertex_to_not_degree_3(vh);
    }  else {
      vh->set_neighbors(deg);
    }
  }
  
  void decrease_neighbors(Vertex_handle vh) {
    vh->set_neighbors(vh->neighbors()-1);
    CGAL_assertion(neighbors_ok(vh));
    //if (vh->neighbors() == 3) update_vertex_to_degree_3(vh);
  }

  void increase_neighbors(Vertex_handle vh) {
    vh->set_neighbors(vh->neighbors()+1);
    CGAL_assertion(neighbors_ok(vh));
    //if (vh->neighbors() == 4) update_vertex_to_not_degree_3(vh);
  }

  void update_vertex(Vertex_handle vh) {
    CGAL_precondition(neighbors_ok(vh));
    if (vh->neighbors() == 3) {
      update_vertex_to_degree_3(vh);
    } else if (vh->neighbors() == 4) {
      update_vertex_to_not_degree_3(vh);
    }
  }

  void update_decreased_vertex(Vertex_handle vh) {
    if (vh->neighbors() == 3) {
      update_vertex_to_degree_3(vh);
    } 
  }

 void update_increased_vertex(Vertex_handle vh) {
    if (vh->neighbors() == 4) {
      update_vertex_to_not_degree_3(vh);
    } 
  }

  bool neighbors_ok(Vertex_handle vh) const {
    return del_.degree(vh)== vh->neighbors();
    //(vh->neighbors() == vh->degree() && vh->degree() <6 
    //|| vh->degree() >=5 && vh->neighbors() >=5);
  }

  void update_edge_no_batch(const Edge &e) {
    CGAL_DELAUNAY_2_DEBUG(std::cout << "Updating edge " 
			  << TDS_helper::origin(e)->point() << " " 
			  << TDS_helper::destination(e)->point() 
			  << std::endl);
    Vertex_handle ov=TDS_helper::origin(e);
    Vertex_handle dv=TDS_helper::destination(e);
    
    CGAL_precondition(neighbors_ok(ov));
    CGAL_precondition(neighbors_ok(dv));

    if (get_undirected_edge_label(e) != Event_key()) {
      CGAL_DELAUNAY_2_DEBUG(std::cout << "Already has event " << std::endl);
      // can't do this since I create all edges around vertex of degree 4 at once
      // CGAL_error();
    } else if (ov->neighbors() ==3 
	       || dv->neighbors() ==3) {
      CGAL_DELAUNAY_2_DEBUG(std::cout << "One end has 3 " << std::endl);
    } else {
      //CGAL_DELAUNAY_2_DEBUG(std::cout << "New certificate" << std::endl);
      new_certificate(e);
    }
  }



  void update_edge(const Edge &e) {
    if (batching_) {
      //delete_certificate(e);
      batched_certs_.push_back(e);
    } else if (get_undirected_edge_label(e) == Event_key()) {
      update_edge_no_batch(e);
    }
  }

  // return true if hull
  bool points(const Edge &e, Point_key ks[4]) const {
    ks[0]= TDS_helper::origin(e)->point();
    ks[1]= TDS_helper::third_vertex(e)->point();
    ks[2]= TDS_helper::destination(e)->point();
    ks[3]= del_.mirror_vertex(e.first, e.second)->point();

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
    return infinity;
  }


  /*bool is_hull_edge(const Edge &e) const {
    return ! del_.mirror_vertex(e.first, e.second)->point().is_valid()
      || ! TDS_helper::third_vertex(e)->point().is_valid()
      || ! TDS_helper::origin(e)->point().is_valid()
      || ! TDS_helper::destination(e)->point().is_valid();
  }


 
  void edge_points(const Edge &e, Point_key ks[4]) const {
    ks[0]= TDS_helper::origin(e)->point();
    ks[1]= TDS_helper::third_vertex(e)->point();
    ks[2]= TDS_helper::destination(e)->point();
    ks[3]= del_.mirror_vertex(e.first, e.second)->point();
  }

  // very dangerous

  void hull_points(const Edge &e,
		   Point_key ks[4]) const {
    ks[0]= TDS_helper::origin(e)->point();
    ks[1]= TDS_helper::third_vertex(e)->point();
    ks[2]= TDS_helper::destination(e)->point();
    ks[3]= del_.mirror_vertex(e.first, e.second)->point();

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
    }*/


  void new_certificate(Edge e) {
    CGAL_precondition(get_undirected_edge_label(e) == Event_key());
    CGAL_DELAUNAY_2_DEBUG(std::cout << "\nMaking certificate for " << TDS_helper::origin(e)->point() << " " 
			  << TDS_helper::destination(e)->point() 
			  << " which would make " << del_.mirror_vertex(e.first, e.second)->point() << " " 
			  << TDS_helper::third_vertex(e)->point()
			  << std::endl);

    Time t; Certificate_data cd;
    Point_key ks[4];
    if (points(e,ks)) {
      if (traits_.hull_certificate_failure_time(e, ks, t, cd)) {
	Event_key k= traits_.simulator_handle()->new_event(t, Event(cd, e, this));
	set_undirected_edge_label(e, k);
      } else {
	set_undirected_edge_label(e, traits_.simulator_handle()->null_event());
      }
    } else {
      if (traits_.internal_certificate_failure_time(e, ks, t, cd)) {
	Event_key k= traits_.simulator_handle()->new_event(t, Event(cd, e, this));
	set_undirected_edge_label(e, k);
      } else {
	set_undirected_edge_label(e, traits_.simulator_handle()->null_event());
      }
    }
  }

  void delete_certificate(Edge e) {
    CGAL_DELAUNAY_2_DEBUG(std::cout << "Cleaning edge " << TDS_helper::origin(e)->point() << " " << TDS_helper::destination(e)->point() << std::endl);
    Event_key k=  get_undirected_edge_label(e);
    if (k != Event_key()) {
      traits_.simulator_handle()->delete_event(k);
      set_undirected_edge_label(e, Event_key());
    }
  }

  Edge 
  mirror_edge(const Edge e) const 
  {
    return Edge(e.first->neighbor(e.second),
                del_.mirror_index(e.first,  e.second));
  }

  Edge canonicalize(Edge e) const {
    if (e.first->neighbor(e.second) < e.first) {
      return mirror_edge(e);
    } else {
      return e;
    }
  }

  static Event_key get_directed_edge_label(const Edge &e) {
    return e.first->get_edge_label(e.second);
  }

  Event_key get_undirected_edge_label(const Edge &e) const {
#ifndef NDEBUG
    if (get_directed_edge_label(e) != get_directed_edge_label(mirror_edge(e))) {
      std::cerr << "FAILURE Edge from " << TDS_helper::origin(e)->point() << " to " 
		<< TDS_helper::destination(e)->point() << " is screwed." << std::endl;
      std::cerr << get_directed_edge_label(e) << " "
                <<  get_directed_edge_label(mirror_edge(e)) << std::endl;
      CGAL_precondition(get_directed_edge_label(e)
			== get_directed_edge_label(mirror_edge(e)));
    }
#endif
    return e.first->get_edge_label(e.second);
  }

  static void set_directed_edge_label(const Edge &e,
				      Event_key l) {
    e.first->set_edge_label(e.second, l);
  }

  void set_undirected_edge_label(const Edge &e,
					Event_key l) const {
    set_directed_edge_label(e,l);
    set_directed_edge_label(mirror_edge(e),l);
  }

};

template <class Sim, class Del, class W, class T>
std::ostream &operator<<(std::ostream &out, const Delaunay_triangulation_2<Sim, Del, W, T> &kd)
{
  kd.write(out);
  return out;
}

template <class Sim, class Del, class W, class T>
void Delaunay_triangulation_2<Sim, Del, W, T>::audit() const
  {
    if (!has_certificates_) return;
    CGAL_LOG(Log::SOME, "Auditing delaunay" << std::endl);
    
    if (del_.number_of_vertices() < 50) {
      CGAL_LOG(Log::LOTS, *this);
    }
    if (del_.dimension() != 2) return;
    Basic_Delaunay sdel(traits_.instantaneous_kernel_object());
    sdel.geom_traits().set_time(traits_.simulator_handle()->audit_time());
    for (typename Triangulation::Finite_vertices_iterator cit= del_.finite_vertices_begin();
	 cit != del_.finite_vertices_end(); ++cit){
      sdel.insert(cit->point());
    }
    /*    sdel.insert(traits_.active_points_2_table_handle()->keys_begin(),
	  traits_.active_points_2_table_handle()->keys_end());*/

    //CGAL_LOG(Log::LOTS, sdel << std::endl);

    if (del_.dimension() != sdel.dimension()) {
      CGAL_LOG(Log::NONE, "AUDIT FAILURE Dimensions don't match in audit" << std::endl);
      return;
    }
    CGAL_exactness_assertion(del_.dimension() == sdel.dimension());
   
   for (typename Triangulation::All_vertices_iterator vit = del_.all_vertices_begin();
	vit != del_.all_vertices_end(); ++vit) {
     if (vit->point() != Point_key()) {
       
       if (!neighbors_ok(vit)) {
	 CGAL_LOG(Log::NONE, "AUDIT FAILURE stored degree is " << vit->neighbors() 
                  << " and actual is " << del_.degree(vit) << " for " << vit->point() << std::endl);
	 CGAL_exactness_assertion(neighbors_ok(vit));
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
	  if (del_.degree(vit) != del_.degree(vit2)) {
	    CGAL_LOG(Log::NONE, "AUDIT FAILURE Degrees don't match in: " 
			     << vit->point() << std::endl);
	  }
	  CGAL_exactness_assertion(del_.degree(vit) == del_.degree(vit2));


	}
      }
      if (!found) {
	CGAL_LOG(Log::NONE, "AUDIT FAILURE Matching vertex not found: " 
			 << vit->point() << std::endl);
      }
      CGAL_exactness_assertion(found);
     
    }

  

    typename Simulation_traits::Instantaneous_kernel ik= traits_.instantaneous_kernel_object();
    ik.set_time(traits_.simulator_handle()->audit_time());

    typename Simulation_traits::Instantaneous_kernel::Side_of_oriented_circle_2 ic2
      = ik.side_of_oriented_circle_2_object();
    for (typename Triangulation::Finite_edges_iterator fit = del_.finite_edges_begin(); 
	 fit != del_.finite_edges_end(); ++fit){
      Point_key k0= fit->first->vertex((fit->second+1)%3)->point();
      Point_key k2= fit->first->vertex((fit->second+2)%3)->point();
      Point_key k3= del_.mirror_vertex(fit->first, fit->second)->point();
      Point_key k1= TDS_helper::third_vertex(*fit)->point();
      if (k1== Point_key() || k3== Point_key()) continue;
      typename Triangulation::Geom_traits::Current_coordinates cc= del_.geom_traits().current_coordinates_object();
      /*typedef typename Triangulation::Geom_traits::Current_coordinates::result_type P2;
      P2 p0= cc(k0);
      P2 p1= cc(k1);
      P2 p2= cc(k2);
      P2 p3= cc(k3);*/
      if (ic2(k0, k1, k2, k3) != CGAL::ON_POSITIVE_SIDE) {
	CGAL_LOG(Log::NONE, "AUDIT FAILURE Failed certificate: " << k0 << " " << k1 << " " 
			 << k2 << " " << k3 << std::endl);
	CGAL_LOG(Log::NONE, "AUDIT FAILURE Points are: " << cc(k0) << ": " << cc(k1) << ": " << cc(k2) 
			 << ": " << cc(k3) << std::endl);
      }
      CGAL_exactness_assertion(ic2(k0, k1, k2, k3) == CGAL::ON_POSITIVE_SIDE);
      
    }

    for (typename Triangulation::Edge_iterator fit = del_.edges_begin(); fit != del_.edges_end(); ++fit){
      if (get_directed_edge_label(*fit) != 
	  get_directed_edge_label(mirror_edge(*fit))) {
	CGAL_LOG(Log::NONE, "AUDIT FAILURE mismatched labels on " 
			 << TDS_helper::origin(*fit)->point() << " " 
			 << TDS_helper::destination(*fit)->point() 
			 << " front has " << get_directed_edge_label(*fit)
			 << " and back has " << get_directed_edge_label(mirror_edge(*fit))<< std::endl);
      }
      if (del_.degree(TDS_helper::origin(*fit))==3 || del_.degree(TDS_helper::destination(*fit))==3) {
	if (get_undirected_edge_label(*fit) != Event_key()) {
	  CGAL_LOG(Log::NONE, "AUDIT FAILURE certificate on degree 3 edge: " 
			   << TDS_helper::origin(*fit)->point()
			   << " " <<  TDS_helper::destination(*fit)->point() 
                   << del_.degree(TDS_helper::origin(*fit)) <<  " "
                   << del_.degree(TDS_helper::destination(*fit)) << std::endl);
	} 
	CGAL_exactness_assertion(get_undirected_edge_label(*fit) == Event_key());
      } else {
	if (get_undirected_edge_label(*fit) == Event_key()) {
	  CGAL_LOG(Log::NONE, "AUDIT FAILURE no certificate on edge: " 
			   << TDS_helper::origin(*fit)->point()
			   << " " <<  TDS_helper::destination(*fit)->point() << std::endl);
	  CGAL_LOG(Log::NONE, "AUDIT FAILURE degrees are: " 
                   << del_.degree(TDS_helper::origin(*fit))
                   << " " <<  del_.degree(TDS_helper::destination(*fit)) << std::endl);
	  
	} 
	CGAL_exactness_assertion(get_undirected_edge_label(*fit) != Event_key());
      }
    }

  }

} } //namespace CGAL::Kinetic
#endif
