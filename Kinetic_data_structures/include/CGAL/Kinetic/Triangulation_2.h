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

#ifndef CGAL_KINETIC_KINETIC_TRIANGULATION_2_H
#define CGAL_KINETIC_KINETIC_TRIANGULATION_2_H
#include <CGAL/Kinetic/basic.h>

#include <CGAL/Kinetic/Triangulation_face_base_2.h>
#include <CGAL/Kinetic/Triangulation_vertex_base_2.h>
#include <CGAL/Kinetic/Triangulation_visitor_base_2.h>
#include <CGAL/Kinetic/Active_objects_batch_listener_helper.h>
#include <CGAL/Kinetic/Simulator_kds_listener.h>
#include <CGAL/Kinetic/internal/tds_2_helpers.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Kinetic/Ref_counted.h>
#include <iterator>
#include <CGAL/Kinetic/Event_base.h>
#include <CGAL/Kinetic/Triangulation_default_traits_2.h>

namespace CGAL { namespace Kinetic {
#ifdef NDEBUG
#define CGAL_TRIANGULATION_2_DEBUG(x)
#else
#define CGAL_TRIANGULATION_2_DEBUG(x) x
#endif

template <class KDel>
struct Triangulation_event: public Event_base<KDel*> {
  typedef Event_base<KDel*>  P;
  Triangulation_event(const typename KDel::Face_handle &e,
		      KDel *kdel): P(kdel), e_(e) {}
  const typename KDel::Face_handle face() const {
    return e_;
  }
  KDel* kdel() {
    return P::kds();
  }

  KDel* kdel() const {
    return P::kds();
  }
  void process() {
    kdel()->handle_failure(e_);
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
    out << "Face " << e_->vertex(0)->point() << ", "
        << e_->vertex(1)->point() 
        << ", " << e_->vertex(2)->point() ;
    return out;
  }
  const typename KDel::Face_handle e_;
};





//! A 2D kinetic triangulation.
/*!  There is one certificate for each face. For the internal faces
  (inside the convex hull), that certificate is for the face
  collapsing at which point we need to figure out which edge to flip.

  For the faces outside the convex hull, the certificate is for the
  edge that is CCW from the infinite vertex and for whether the other
  vertex on the edge leaves the convex hull.

  The annoying part is how to figure out which edge to flip when an
  internal face collapses.  As far as I can tell this requires
  evaluating which point is in the segment defined by which other two,
  so four 1D orientation predicates.
*/
template <class Simulation_traits_t, 
          class Visitor= Triangulation_visitor_base_2,
          class Tri
          = CGAL::Triangulation_2<typename Simulation_traits_t::Instantaneous_kernel,
                                           CGAL::Triangulation_data_structure_2<
            Triangulation_vertex_base_2<typename Simulation_traits_t::Instantaneous_kernel>,
  CGAL::Kinetic::Triangulation_face_base_2<Simulation_traits_t > > >,
          class Triangulation_traits_t= Triangulation_default_traits_2<Simulation_traits_t, Tri> >
class Triangulation_2:
  public Ref_counted<Triangulation_2<Simulation_traits_t, Visitor, Tri, Triangulation_traits_t> >
{

  typedef CGAL::Triangulation_2<typename Simulation_traits_t::Instantaneous_kernel,
                                         CGAL::Triangulation_data_structure_2<
    Triangulation_vertex_base_2<typename Simulation_traits_t::Instantaneous_kernel>,
    CGAL::Kinetic::Triangulation_face_base_2<Simulation_traits_t > > > Basic_Delaunay;

public:
  typedef Triangulation_traits_t Traits;
  typedef Simulation_traits_t Simulation_traits;
  typedef Triangulation_2<Simulation_traits, Visitor, Tri, Traits> This;

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

  typedef Triangulation_event<This> Event;
  
  //friend class Delaunay_edge_failure_event<This>;
  //friend class Delaunay_hull_edge_failure_event<This>;

  typedef typename CGAL::Kinetic::Simulator_kds_listener<typename Simulator::Listener, This> Simulator_listener;
  friend  class CGAL::Kinetic::Simulator_kds_listener<typename Simulator::Listener, This>;
  typedef typename CGAL::Kinetic::Active_objects_batch_listener_helper<typename Moving_point_table::Listener, This> Moving_point_table_listener;
  friend class CGAL::Kinetic::Active_objects_batch_listener_helper<typename Moving_point_table::Listener, This>;

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
    siml_ = Simulator_listener(traits_.simulator_handle(), this);
    motl_= Moving_point_table_listener(traits_.active_points_2_table_handle(), this, insert);
    has_certificates_=false; 
    clear_stats();
   
    batching_=false;
  }

public:

  Triangulation_2(Traits st,
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
  
    set_has_certificates(true);
  }

  Triangulation_2(Simulation_traits st,
		  Visitor w= Visitor()): 
    traits_(st),
    watcher_(w),
    del_(traits_.instantaneous_kernel_object()) {
    init_data(true);
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
        
        //std::cout << std::min(s,t) << "--" << std::max
		//BOOST_PREVENT_MACRO_SUBSTITUTION(s,t) << std::endl;
        update_face(batched_certs_[i]);
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

  void audit_event(Event_key k, Face_handle e) const {
    if (e->event() !=k) {
      std::cerr << "AUDIT FAILURE orphan event " << k << std::endl;
    }
    CGAL_assertion(e->event() ==k);
  }

  void set_has_certificates(bool tf, bool no_failures=false) {
    if (tf == has_certificates_){

    } else {
      if (tf==true && del_.dimension()==2) {
        CGAL_LOG(CGAL::Kinetic::Log::SOME, "DELAUNAY2: Creating certificates."<< std::endl);
        if (no_failures) {
          for (Face_iterator f = del_.all_faces_begin(); f != del_.all_faces_end(); ++f) {
            f->set_event(traits_.simulator_handle()->null_event());
          }
        }  
        /*for (typename Triangulation::All_vertices_iterator vit = del_.all_vertices_begin(); 
	     vit != del_.all_vertices_end(); ++vit) {
	  int deg=TDS_helper::low_degree(vit, del_.tds());
	  CGAL_assertion(deg >=3);
	  Vertex_handle vh= vit;
	  vh->set_neighbors(deg);
	  CGAL_TRIANGULATION_2_DEBUG(std::cout << "Set degree of " << vit->point() << " to " << deg << std::endl);
          //vit->set_neighbors_is_changed(false);
          }*/
        if (!no_failures) {
          for (Face_iterator eit = del_.all_faces_begin(); eit != del_.all_faces_end(); ++eit) {
            eit->set_event(Event_key());
            update_face(eit);
          }
        }

        watcher_.create_faces(del_.all_faces_begin(), del_.all_faces_end());
      } else if (tf==false) { 
        for (Face_iterator f = del_.all_faces_begin(); f != del_.all_faces_end(); ++f) {
          delete_certificate(f);
        }
      } 
      CGAL_LOG(CGAL::Kinetic::Log::SOME, 
		       *traits_.simulator_handle() << std::endl;);
      has_certificates_=tf;
    }
  }
  bool has_certificates() {
    return has_certificates_;
  }

  void erase(Point_key k) {
    CGAL_error();
#if 0    
    // erase all incident certificates
    Vertex_handle vh= vertex_handle(k);
    if (vh == Vertex_handle()) {
      CGAL_LOG(CGAL::Kinetic::Log::SOME, "Point " << k << " is not in triangulation on removal."<< std::endl);
      return;
    }
    watcher_.remove_vertex(vh);
    Face_handle hull_cert;
    if (has_certificates_) {
      // have to clean up hull events which are not adjacent
      Face_circulator fc= vh->incident_faces(), fe=fc;
      if (fc != NULL) {
	do {
	  delete_certificate(fc);
	  ++fc;
	} while (fc != fe);
      }
      hull_cert= hull_face(k);
      if (hull_cert != Face_handle()) {
	delete_certificate(hull_cert);
      }
    }
    // remove from triangulation
    del_.geom_traits().set_time(traits_.rational_current_time());
    del_.remove(vh);
    //new_edges_.clear();
    if (del_.dimension()==2 && has_certificates_) {
      std::vector<Face_handle> faces;

      del_.get_conflicts(k,std::back_inserter(faces));

      if (hull_cert != Face_handle()) {
	faces.push_back(hull_cert);
      }
      // update non-adjacent hull certificates
      for (unsigned int i=0; i< faces.size(); ++i) {
	update_face(faces[i]);
      }
      
      watcher_.create_faces(faces.begin(), faces.end());
    }
#endif
  }

  //! The assertion will catch that the object is in the same sorted order
  void set(Point_key k) {
    std::cout << "Object changed " << k << std::endl;

    //new_edges_.clear();
    traits_.point_changed(k);
    if (del_.dimension() != 2) {
      CGAL_LOG(CGAL::Kinetic::Log::SOME,"Triangulation is still 1D.\n");
      return;
    }

    Vertex_handle vh=vertex_handle(k);
    if (vh == Vertex_handle()) {
      CGAL_LOG(CGAL::Kinetic::Log::SOME, "Point " << k << " is not in triangulation on set."<< std::endl);
      return;
    }

    if (has_certificates_) {
      Face_circulator fc= vh->incident_faces(), fe=fc;
      if (fc != NULL) {
	do {
	  delete_certificate(fc);
	  update_face(fc);
	  ++fc;
	} while (fc != fe);
      }
      Face_handle hull_cert= hull_face(vh);
      if (hull_cert != Face_handle()) {
	std::cout << "Hull face is " << hull_cert->vertex(0)->point() 
		  << ", " << hull_cert->vertex(1)->point() << ", " 
		  << hull_cert->vertex(2)->point() << std::endl;
	delete_certificate(hull_cert);
	update_face(hull_cert);
      }
    }

   
    watcher_.modify_vertex(vh);

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
    CGAL_TRIANGULATION_2_DEBUG(std::cout << "Inserting " << k << std::endl);
    bool was_2d= (del_.dimension()==2);

    del_.geom_traits().set_time(traits_.rational_current_time());

    Face_handle lf= del_.locate(k);

    if (was_2d && has_certificates_) {
      watcher_.remove_faces(&lf, &lf+1);
      delete_certificate(lf);
    }
    Vertex_handle vh= del_.insert(k, lf);
    set_vertex_handle(k, vh);
    
    CGAL_assertion(vertex_handle(k) != Vertex_handle());
    
    //CGAL_TRIANGULATION_2_DEBUG(std::cout << "Vertex " << vertex_handle(k)->point() << std::endl);
    watcher_.create_vertex(vertex_handle(k));

    // now have to update
    if (has_certificates_) {
      if (!was_2d && del_.dimension()==2) {
	has_certificates_=false;
	set_has_certificates(true); 
      } else if (del_.dimension() == 2) {
	set(k);
      }
    } 
    //write(std::cout);
    //if (del_.dimension()==2) audit();
  }



  Comparison_result compare_concurrent(Event_key a, Event_key b) const {
    //Edge ea= traits_.simulator_handle()->template event<Event>(a).edge();
    //Edge eb= traits_.simulator_handle()->template event<Event>(b).edge();
    return traits_.compare_concurrent(a, b);
  }



  void handle_failure(Face_handle f) {
    int infi=-1;
    for (unsigned int i=0; i< 3; ++i) {
      if (f->vertex(i)->point() == Point_key()) {
	infi=i;
	break;
      }
    }
    if (infi != -1) {
      flip(Edge(f, (infi+2)%3));
    } else {
      // determine which edge to flip
      CGAL::Comparison_result cx01= traits_.compare(f->vertex(0)->point(), f->vertex(1)->point(), 0);
      CGAL::Comparison_result cx12= traits_.compare(f->vertex(1)->point(), f->vertex(2)->point(), 0);
      if (cx01== cx12) {
	CGAL_precondition(cx01 != CGAL::EQUAL);
	flip(Edge(f, 1));
      } else {
	CGAL::Comparison_result cx02= traits_.compare(f->vertex(0)->point(), f->vertex(2)->point(), 0);
	if (cx02 == -cx12) {
	  flip(Edge(f,2));
	} else {
	  flip(Edge(f,0));
	}
      }
    }
  }


  void flip(const Edge &e) {
    ++num_events_;
    CGAL_precondition(!batching_);
    CGAL_LOG(CGAL::Kinetic::Log::SOME, "\n\n\n\n\n\nFlipping edge "
		     << TDS_helper::origin(e)->point()
		     << TDS_helper::destination(e)->point() 
		     << " to get " << TDS_helper::third_vertex(e)->point()
		     << ", " << TDS_helper::mirror_vertex(e)->point()<< std::endl);
    //CGAL_LOG(CGAL::Kinetic::Log::NONE, TDS_helper::destination(e)->point() << std::endl);
    //CGAL_LOG(CGAL::Kinetic::Log::SOME, " at "  << traits_.simulator()->current_time() << std::endl);

    Face_handle face= e.first;
    int index= e.second;
    Face_handle mirror_face = face->neighbor(index);
    int mirror_index =face->mirror_index(index);
    Edge em(mirror_face,mirror_index);
    CGAL_precondition(mirror_face->neighbor(mirror_index) == face);

    Face_handle bef;
    int bei;
 
    CGAL_assertion (!del_.is_edge(TDS_helper::third_vertex(e), TDS_helper::mirror_vertex(e),
				  bef, bei));/* {
      std::cout << "Flipping edge out of the way" << std::endl;
      flip(Edge(bef, bei));
      }*/

    //delete_certificate(face);
    delete_certificate(mirror_face);
    face->set_event(Event_key());

   
    watcher_.before_flip(e);
    del_.tds().flip(face,index);
   
    // we also know that CGAL preserves the edge index of the flipped edge
    //mirror_index = mirror_face->index(face);
    //index= face->index(mirror_face);

    watcher_.after_flip(Edge(face, index));

    update_face(face);
    update_face(mirror_face);
    for (unsigned int i=0; i< 3; ++i) {
      if (face->vertex(i)->point() == Point_key()) {
	if (face->neighbor((i+1)%3) != mirror_face) {
	  delete_certificate(face->neighbor((i+1)%3));
	  update_face(face->neighbor((i+1)%3));
	}
      } 
      if (mirror_face->vertex(i)->point() == Point_key()) {
	if (mirror_face->neighbor((i+1)%3) != face) {
	  delete_certificate(mirror_face->neighbor((i+1)%3));
	  update_face(mirror_face->neighbor((i+1)%3));
	}
      }
    }

    //return flipped_edge;
  }

  Visitor &visitor() {
    return watcher_;
  }

  const Visitor &visitor() const
  {
    return watcher_;
  }

  bool has_event(const Edge &e) const {
    return e.first->event() != Event_key() || e.first->neighbor(e.second)->event() != Event_key();
  }

  bool has_finite_event(const Edge &e) const {
    return e.first->event() != traits_.simulator_handle()->null_event()
    || e.first->neighbor(e.second)->event() != traits_.simulator_handle()->null_event() ;
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
  std::vector<Face_handle> batched_certs_;


  Face_handle hull_face(Vertex_handle v) const {
    Face_circulator fc= v->incident_faces(), fe=fc;
    if (fc != NULL) {
      do {
	int i= fc->index(v);
	if (fc->vertex((i+2)%3)->point() == Point_key()) {
	  return fc->neighbor(i);
	}
	++fc;
      } while (fc != fe);
    }
    return Face_handle();
  }

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
    vhs_.resize(std::max BOOST_PREVENT_MACRO_SUBSTITUTION(k.index()+1, vhs_.size()));
    vhs_[k.index()]=vh;
  }




  void update_face_no_batch(const Face_handle e) {
    CGAL_TRIANGULATION_2_DEBUG(std::cout << "Updating face " 
			  << e->vertex(0)->point() << ", " << e->vertex(1)->point()
			  << ", " << e->vertex(2)->point() << std::endl);

    if (e->event() != Event_key()) {
      CGAL_TRIANGULATION_2_DEBUG(std::cout << "Already has event " << std::endl);
      // can't do this since I create all edges around vertex of degree 4 at once
      // CGAL_error();
    } else {
      //CGAL_TRIANGULATION_2_DEBUG(std::cout << "New certificate" << std::endl);
      new_certificate(e);
    }
  }



  void update_face(const Face_handle e) {
    if (batching_) {
      //delete_certificate(e);
      batched_certs_.push_back(e);
    } else if (e->event() == Event_key()) {
      update_face_no_batch(e);
    }
  }

  // return true if hull
  Vertex_handle points(const Face_handle e, Point_key ks[3]) const {
    ks[0]= e->vertex(0)->point();
    ks[1]= e->vertex(1)->point();
    ks[2]= e->vertex(2)->point();

    
    Vertex_handle cv;
    for (unsigned int i=0; i<3; ++i) {
      if (ks[i]== Point_key()) {
	cv=e->vertex((i+1)%3);
	ks[i]= TDS_helper::mirror_vertex(Edge(e, (i+2)%3))->point();
      }
    }
   
    if (cv != Vertex_handle()) {
      std::swap(ks[0], ks[1]);
    }
    return cv;
  }

  void new_certificate(Face_handle e) {
    CGAL_precondition(e->event() == Event_key());
 
    Time t;
    Point_key ks[3];
    Vertex_handle cv= points(e,ks);
    if (cv != Vertex_handle()) {
      std::cout << "Center is " << cv->point() << std::endl;
      if (cv->degree() ==3) {
	std::cout << "Skipping set " << ks[0] << ", " << ks[1] << ", " << ks[2] << std::endl;
	e->set_event(traits_.simulator_handle()->null_event());
	return;
      }
    }
    if (traits_.certificate_failure_time(e, ks, t)) {
      Event_key k= traits_.simulator_handle()->new_event(t, Event(e, this));
      e->set_event(k);
    } else {
      e->set_event(traits_.simulator_handle()->null_event());
    }
   
  }

  void delete_certificate(Face_handle e) {
    //CGAL_TRIANGULATION_2_DEBUG(std::cout << "Cleaning edge " << TDS_helper::origin(e)->point() << " " << TDS_helper::destination(e)->point() << std::endl);
    Event_key k=  e->event();
    if (k != Event_key()) {
      traits_.simulator_handle()->delete_event(k);
      e->set_event(Event_key());
    }
  }
};

template <class Sim, class Del, class W, class T>
std::ostream &operator<<(std::ostream &out, const Triangulation_2<Sim, Del, W, T> &kd)
{
  kd.write(out);
  return out;
}

template <class Sim, class Del, class W, class T>
void Triangulation_2<Sim, Del, W, T>::audit() const  {
  del_.geom_traits().set_time(traits_.rational_current_time());
  typename Simulation_traits::Instantaneous_kernel::Orientation_2 o2=
    del_.geom_traits().orientation_2_object();
  for (typename Triangulation::All_faces_iterator f = del_.all_faces_begin(); 
       f != del_.all_faces_end(); ++f) {
    Point_key ks[3];
    points(f,ks);
    CGAL_assertion(o2(ks[0],ks[1],ks[2])== CGAL::POSITIVE);
  }
}

} } //namespace CGAL::Kinetic
#endif
