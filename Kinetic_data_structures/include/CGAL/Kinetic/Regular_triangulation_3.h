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

#ifndef CGAL_KINETIC_KINETIC_REGULAR_TRIANGULATION_3_H
#define CGAL_KINETIC_KINETIC_REGULAR_TRIANGULATION_3_H
#include <CGAL/basic.h>

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Kinetic/internal/Delaunay_triangulation_base_3.h>
#include <CGAL/Kinetic/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Kinetic/Regular_triangulation_cell_base_3.h>
#include <CGAL/Kinetic/Regular_triangulation_visitor_base_3.h>
#include <CGAL/Kinetic/Listener.h>
#include <CGAL/Kinetic/Ref_counted.h>

#include <CGAL/Kinetic/Simulator_kds_listener.h>
#include <CGAL/Kinetic/Active_objects_listener_helper.h>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4355) // complaint about using 'this' to
#endif                          // initialize a member


CGAL_KINETIC_BEGIN_INTERNAL_NAMESPACE

template <class KD>
class Regular_3_pop_event: public Delaunay_event_base_3<KD>
{
  typedef Delaunay_event_base_3<KD>  P;
public:
  Regular_3_pop_event(const typename KD::Root_stack &s,
		      const typename KD::Vertex_handle &vh,
		      KD *kd): P(s, kd), vh_(vh) {
  }

  void process() {
    P::kdel()->pop(vh_, P::root_stack());
  }

  typename KD::Vertex_handle vertex() const
  {
    return vh_;
  }

  std::ostream& write(std::ostream &out) const
  {
    out << "Pop " << vh_->point();
    return out;
  }

  void audit(typename KD::Event_key k) const {
    P::kdel()->audit_pop(k, vh_);
  }
  /*virtual bool is_of_type(int tp) const {
    return (tp &type())!=0 || P::is_of_type(tp);
    }
    static int type() {
    return 2;
    }*/

  virtual ~Regular_3_pop_event(){};

protected:
  const typename KD::Vertex_handle vh_;
};

/*template <class K, class S, class VH>
  std::ostream& operator<<(std::ostream &out, const Regular_3_pop_event<K,S,VH> &e)
  {
  e.write(out);
  return out;
  }*/


template <class KD>
class Regular_3_non_vertex_event: public Delaunay_event_base_3<KD>
{
  typedef Delaunay_event_base_3<KD>  P;
public:

  Regular_3_non_vertex_event(const typename KD::Root_stack &s,
			     const typename KD::Point_key &k,
			     const typename KD::Cell_handle &c,
			     KD *kd): P(s,kd), k_(k), cell_(c) {
  }

  Regular_3_non_vertex_event(){}

  std::ostream& write(std::ostream &out) const
  {
    out << "Nothing " << P::vh_->point();
    return out;
  }

  /* virtual bool is_of_type(int tp) const {
     return (tp &type())!=0 || P::is_of_type(tp);
     }
     static int type() {
     return 4;
     }*/

  typename KD::Point_key point() const {return k_;}

  typename KD::Cell_handle cell_handle() const {return cell_;}

  virtual ~Regular_3_non_vertex_event(){};

protected:
  const typename KD::Point_key k_;
  const typename KD::Cell_handle cell_;
};

template <class KD>
class Regular_3_move_event: public Regular_3_non_vertex_event<KD>
{
  typedef Regular_3_non_vertex_event<KD>  P;
public:
  Regular_3_move_event(const typename KD::Root_stack &s,
		       const typename KD::Point_key &k,
		       const typename KD::Cell_handle &c,
		       int dir,
		       KD *kd): P(s,k, c, kd), dir_(dir) {
  }

  void process() {
    P::kdel()->move(P::k_, P::cell_, dir_, P::root_stack());
  }

  std::ostream& write(std::ostream &out) const
  {
    out << "Move " << P::point();
    return out;
  }

  void audit(typename KD::Event_key k) const{
    P::kdel()->audit_move(k, P::point(), P::cell_, dir_);
  }
  /*virtual bool is_of_type(int tp) const {
    return (tp &type())!=0 || P::is_of_type(tp);
    }
    static int type() {
    return 8;
    }*/

  virtual ~Regular_3_move_event(){};

protected:
  int dir_;
};

/*template <class K,  class S, class KK, class C>
  std::ostream& operator<<(std::ostream &out, const Regular_3_move_event<K,S,KK,C> &e)
  {
  e.write(out);
  return out;
  }*/


template <class KD>
class Regular_3_push_event: public Regular_3_non_vertex_event<KD>
{
  typedef Regular_3_non_vertex_event<KD>  P;
public:

  Regular_3_push_event(const typename KD::Root_stack &s,
		       const typename KD::Point_key &k,
		       const typename KD::Cell_handle &c,
		       KD *kd): P(s,k, c, kd) {
  }

  void process() {
    P::kdel()->push(P::k_, P::cell_, P::root_stack());
  }

  std::ostream& write(std::ostream &out) const
  {
    out << "Push " << P::point();
    return out;
  }
  
  void audit(typename KD::Event_key k) const {
    P::kdel()->audit_push(k, P::point(), P::cell_);
  }
  /*virtual bool is_of_type(int tp) const {
    return (tp &type())!=0 || P::is_of_type(tp);
    }
    static int type() {
    return 16;
    }*/

  virtual ~Regular_3_push_event(){};
};

/*template <class K, class S, class KK, class C>
  std::ostream& operator<<(std::ostream &out, const Regular_3_push_event<K,S,KK,C> &e)
  {
  e.write(out);
  return out;
  }*/


template <class Traits>
struct Regular_triangulation_3_types
{
  typedef typename Traits::Active_points_3_table MPT; // here
  typedef typename Traits::Kinetic_kernel KK;
  // typedef CGAL::Triangulation_cell_base_3<typename Traits::Instantaneous_kernel> CFB;

  /*typedef CGAL::Triangulation_cell_base_with_info_3<Delaunay_cache_3<MPT, KK>,
    typename Traits::Instantaneous_kernel, CFB> CFBI;*/

  //typedef Triangulation_labeled_edge_cell_base_3<CFBI, typename Traits::Simulator::Event_key> TFB;
  //typedef Triangulation_labeled_facet_cell_base_3<TFB, typename Traits::Simulator::Event_key> FB;
  /*typedef CGAL::Triangulation_vertex_base_3<typename Traits::Instantaneous_kernel> CVB;
    typedef CGAL::Triangulation_vertex_base_with_info_3<typename Traits::Simulator::Event_key,
    typename Traits::Instantaneous_kernel,
    CVB> LVB;*/
  typedef CGAL::Kinetic::Regular_triangulation_cell_base_3<Traits> FB;
  typedef CGAL::Kinetic::Regular_triangulation_vertex_base_3<Traits> LVB;

  typedef CGAL::Triangulation_data_structure_3<LVB, FB> TDS;

  typedef CGAL::Regular_triangulation_3<typename Traits::Instantaneous_kernel, TDS> Default_triangulation;

  template <class KRT>
  struct Mirroring_visitor {
    template <class Point_key, class Cell_handle>
    void pre_insert_vertex(Point_key v, Cell_handle h) {
      v_.pre_insert_vertex(v, h);
    }

    template <class Vertex_handle>
    void post_insert_vertex(Vertex_handle v) {
      v_.post_insert_vertex(v);
    }

    template <class Vertex_handle>
    void pre_remove_vertex(Vertex_handle v) {
      v_.pre_remove_vertex(v);
    }


    template <class Point_key>
    void post_remove_vertex(Point_key v) {
      v_.post_remove_vertex(v);
    }

    template <class Vertex_handle>
    void change_vertex(Vertex_handle v) {
      v_.change_vertex(v);
    }

    template <class Cell_handle> 
    void create_cell(Cell_handle h) {
      krt_->create_cell(h);
      v_.create_cell(h);
    }

    template <class Cell_handle>
    void destroy_cell(Cell_handle h) {
      krt_->destroy_cell(h);
      v_.destroy_cell(h);
    }

    template <class Edge>
    void pre_edge_flip(const Edge &e) {
      v_.pre_edge_flip(e);
    }
    template <class Edge>
    void post_facet_flip(const Edge& e) {
      v_.post_facet_flip(e);
    }

    template <class Facet>
    void pre_facet_flip(const Facet &f) {
      v_.pre_facet_flip(f);
    }

    template <class Facet>
    void post_edge_flip(const Facet& f) {
      v_.post_edge_flip(f);
    }

    template <class Key, class Cell>
    void pre_move(Key k, Cell h){
      v_.pre_move(k,h);
    }

    template <class Key, class Cell>
    void post_move(Key k, Cell h){
      v_.post_move(k,h);
    }

    Mirroring_visitor(KRT* krt, typename KRT::Visitor v): krt_(krt), v_(v){}

    KRT* krt_;
    typename KRT::Visitor v_;
  };

  //friend class CGAL::Delaunay_triangulation_3<typename P::Instantaneous_kernel, TDS>;

};

CGAL_KINETIC_END_INTERNAL_NAMESPACE

CGAL_KINETIC_BEGIN_NAMESPACE

/*!
  redundant_cells_ maps each cell with redundant points to the ids of the points in that cell

  redundant_points_ maps each redundant point to a certificate
*/
template <class TraitsT,
	  class VisitorT= Regular_triangulation_visitor_base_3,
	  class TriangulationT= typename internal::Regular_triangulation_3_types<TraitsT>::Default_triangulation>
class Regular_triangulation_3:
  public Ref_counted<Regular_triangulation_3<TraitsT, VisitorT, TriangulationT> >
{
private:
  typedef Regular_triangulation_3<TraitsT, VisitorT, TriangulationT> This;

public:
  typedef TraitsT Traits;
  typedef typename Traits::Active_points_3_table::Key Point_key; //here

protected:
  typedef typename Traits::Active_points_3_table MPT; // here
  typedef typename Traits::Simulator Simulator;
  typedef typename Traits::Simulator::Event_key Event_key;
  typedef typename Traits::Simulator::Time Time;

  typedef typename Traits::Kinetic_kernel::Certificate Root_stack;
  typedef TriangulationT Delaunay;

  typedef internal::Regular_triangulation_3_types<TraitsT> Types;

  typedef typename Types::template Mirroring_visitor<This> Delaunay_visitor;
  friend class internal::Regular_triangulation_3_types<TraitsT>::Mirroring_visitor<This>;

  typedef typename Delaunay::Facet Facet;
  typedef typename Delaunay::Edge Edge;
  typedef typename Delaunay::Cell_handle Cell_handle;
  typedef typename Delaunay::Vertex_handle Vertex_handle;
public:
  typedef VisitorT Visitor;


  friend class internal::Delaunay_event_base_3<This>;
  typedef internal::Delaunay_3_edge_flip_event<This> Edge_flip;
  friend class internal::Delaunay_3_edge_flip_event<This>;
  typedef internal::Delaunay_3_facet_flip_event<This> Facet_flip;
  friend class internal::Delaunay_3_facet_flip_event<This>;

  typedef internal::Regular_3_pop_event<This> Pop_event;
  friend class internal::Regular_3_pop_event<This>;

  typedef internal::Regular_3_non_vertex_event<This> Non_vertex_event;
  friend class internal::Regular_3_non_vertex_event<This>;

  typedef internal::Regular_3_move_event<This> Move_event;
  friend class internal::Regular_3_move_event<This> ;

  typedef internal::Regular_3_push_event<This> Push_event;
  friend class internal::Regular_3_push_event<This> ;

protected:
  typedef std::multimap<typename Delaunay::Cell_handle,
			Point_key> RCMap;
  typedef std::map<Point_key, Event_key> RPMap;

  struct Base_traits;
  friend struct Base_traits;

  struct Base_traits: public TraitsT
  {
    typedef This Wrapper;
    typedef TriangulationT Triangulation;
    typedef typename This::Edge_flip Edge_flip;
    typedef typename This::Facet_flip Facet_flip;
    typedef typename TraitsT::Kinetic_kernel::Power_test_3 Side_of_oriented_sphere_3;
    typedef typename TraitsT::Kinetic_kernel::Weighted_orientation_3 Orientation_3;
    typedef typename TraitsT::Active_points_3_table Active_points_3_table; // here

    Side_of_oriented_sphere_3 side_of_oriented_sphere_3_object() const
    {
      //std::cout << "Getting power test" << std::endl;
      return TraitsT::kinetic_kernel_object().power_test_3_object();
    }

    Orientation_3 orientation_3_object() const
    {
      return TraitsT::kinetic_kernel_object().weighted_orientation_3_object();
    }

    Base_traits(This *t, const TraitsT &tr): TraitsT(tr), wr_(t) {}

    Wrapper* wrapper_handle() {
      return wr_;
    }
    const Wrapper* wrapper_handle() const
    {
      return wr_;
    }

    Active_points_3_table* active_points_3_table_handle() const {
      return TraitsT::active_points_3_table_handle(); // here
    }

    Wrapper *wr_;
  };

  typedef internal::Delaunay_triangulation_base_3<Base_traits, Delaunay_visitor> KDel;

  typedef typename CGAL::Kinetic::Simulator_kds_listener<typename TraitsT::Simulator::Listener, This>
  Simulator_listener;
  friend  class CGAL::Kinetic::Simulator_kds_listener<typename TraitsT::Simulator::Listener, This>;
  typedef typename CGAL::Kinetic::Active_objects_listener_helper<typename TraitsT::Active_points_3_table::Listener, This>
  Moving_point_table_listener;// here
  friend class CGAL::Kinetic::Active_objects_listener_helper<typename TraitsT::Active_points_3_table::Listener, This>; // here

public:
  

  Regular_triangulation_3(Traits tr, Visitor v= Visitor()): kdel_(Base_traits(this, tr), Delaunay_visitor(this, v)),
							    listener_(NULL) {
    siml_= Simulator_listener(tr.simulator_handle(), this);
    motl_= Moving_point_table_listener(tr.active_points_3_table_handle(), this); // here
  }


  const Visitor &visitor() const
  {
    return kdel_.visitor().v_;
  }

  typedef TriangulationT Triangulation;
  const Triangulation& triangulation() const
  {
    return kdel_.triangulation();
  }

  struct Listener_core
  {
    typedef typename This::Handle Notifier_handle;
    typedef enum {TRIANGULATION}
      Notification_type;
  };
  friend class Listener<Listener_core>;
  typedef Kinetic::Listener<Listener_core> Listener;


  void audit_move(Event_key k, Point_key pk, Cell_handle , int) const {
    CGAL_assertion(kdel_.vertex_handle(pk) == Vertex_handle());
    CGAL_assertion(redundant_points_.find(pk) != redundant_points_.end());
    CGAL_assertion(redundant_points_.find(pk)->second == k);
  }
  
  void audit_push(Event_key k, Point_key pk, Cell_handle) const {
    CGAL_assertion(kdel_.vertex_handle(pk) == Vertex_handle());
    CGAL_assertion(redundant_points_.find(pk) != redundant_points_.end());
    CGAL_assertion(redundant_points_.find(pk)->second == k);
  }
  void audit_pop(Event_key k, Vertex_handle vh) const {
    CGAL_assertion(vh->info() == k);
  }

  void audit() const
  {
    CGAL_KINETIC_LOG(LOG_LOTS, "Verifying regular.\n");
    if (!has_certificates()) return;
    CGAL_KINETIC_LOG(LOG_LOTS, *this << std::endl);
    //P::instantaneous_kernel().set_time(P::simulator()->audit_time());
    kdel_.audit();
    audit_structure();
    //  RPMap redundant_points_;
    // RCMap redundant_cells_;

    Triangulation tri= triangulation();
    for (typename RPMap::const_iterator it= redundant_points_.begin(); it != redundant_points_.end(); ++it) {
      CGAL_assertion(tri.insert(it->first) == Vertex_handle());
      CGAL_assertion_code(Cell_handle ch= get_cell_handle(it->first));
      CGAL_assertion(triangulation().locate(it->first) == ch);
    }
  
  }

  void write(std::ostream &out) const {
    if (triangulation().dimension() != 3) return;
    kdel_.write(out);
    out << "Redundant points: ";
    for (typename RPMap::const_iterator it= redundant_points_.begin(); it != redundant_points_.end();
	 ++it) {
      out << it->first << " ";
    }
    out << std::endl;
    typename Delaunay::Cell_handle last;
    out << "Redundant cells: ";
    for (typename RCMap::const_iterator it= redundant_cells_.begin(); it != redundant_cells_.end();
	 ++it) {
      if (it->first != last) {
	last= it->first;
	internal::write_cell(last, out);
	out << ": ";
      }
      out << it->second << " ";
    }
    out << std::endl;
  }



  void push(Point_key k, typename Triangulation::Cell_handle h, Root_stack rs) {
    CGAL_KINETIC_LOG(LOG_LOTS, "Pushing " << k << " into cell ");
    CGAL_KINETIC_LOG_WRITE(LOG_LOTS, internal::write_cell(h, LOG_STREAM));
    CGAL_KINETIC_LOG(LOG_LOTS, std::endl);
    
    redundant_points_.erase(k);

    remove_redundant(h,k);
    
    typename Triangulation::Vertex_handle vh= kdel_.insert(k, h);

    handle_vertex(vh, rs);

    on_geometry_changed();
  }

  void pop(typename Triangulation::Vertex_handle vh, const Root_stack &rs) {
    CGAL_KINETIC_LOG(LOG_LOTS, "Popping " << vh->point() << std::endl);
   
    Point_key k= vh->point();
    typename Triangulation::Cell_handle h= kdel_.pop_vertex(vh);

    handle_redundant(k, h, rs);
    /*if (!success) {
      std::cerr << "dropped a vertex when popped.\n";
      redundant_points_[k]=kdel_.simulator()->null_event();
      }*/
    //CGAL_postcondition(success);
    on_geometry_changed();
  }

  void move(Point_key k, typename Triangulation::Cell_handle h, int dir, const Root_stack &) {
    kdel_.visitor().pre_move(k,h);
    CGAL_KINETIC_LOG(LOG_LOTS, "Moving " << k << " from ");
    CGAL_KINETIC_LOG_WRITE(LOG_LOTS, internal::write_cell(h, LOG_STREAM));
    CGAL_KINETIC_LOG(LOG_LOTS, " to ");
    CGAL_KINETIC_LOG_WRITE(LOG_LOTS, internal::write_cell(h->neighbor(dir), LOG_STREAM ));
    CGAL_KINETIC_LOG(LOG_LOTS, std::endl);
    typename Triangulation::Cell_handle neighbor = h->neighbor(dir);
    remove_redundant(h, k);

    bool hinf=false;
    for (unsigned int i=0; i<4; ++i) {
      if (neighbor->vertex(i)== triangulation().infinite_vertex()) {
	hinf=true;
	break;
      }
    }
    if (hinf) {
      insert(k, neighbor);
    } else {
      handle_redundant(k, neighbor);
    }
    kdel_.visitor().post_move(k,neighbor);
  }

  //! remove an object
  /*!
    See if it is redundant. If so, remove it from the list and delete its certificate.
    Otherwise, pass it along.
  */
  void erase(Point_key ) {
    CGAL_assertion(0);
    on_geometry_changed();
  }

  void set(Point_key k) {
    if (!kdel_.has_certificates()) return;
    if (kdel_.vertex_handle(k) != NULL) {
      kdel_.change_vertex(k);
      if (kdel_.is_degree_4(kdel_.vertex_handle(k))) {
	handle_vertex(kdel_.vertex_handle(k));
      }
    } else {
      //kdel_.simulator()->template event<Non_vertex_event>(redundant_points_[k]);
      typename Triangulation::Cell_handle h= get_cell_handle(k);
      kdel_.simulator()->delete_event(redundant_points_[k]);
      redundant_points_.erase(k);
      handle_redundant(k, h);
    }
  }

  void insert(Point_key k, Cell_handle h=Cell_handle()) {
    // almost the same as push
    // if insertion fails, then handle redundant
    CGAL_KINETIC_LOG(LOG_LOTS, "Insert " << k << std::endl);

    typename Triangulation::Vertex_handle vh= kdel_.insert(k, h);
    if (vh == typename Triangulation::Vertex_handle()) {
      if (h==Cell_handle()) {
	h= triangulation().locate(k);
      }
      handle_redundant(k,h);
    } else {
      handle_vertex(vh); 
    }

    on_geometry_changed();
  }


public:
  void set_has_certificates(bool tf) {
    if (tf == has_certificates()) return;
    if (tf==false) {
      for (typename Triangulation::Finite_vertices_iterator vit= triangulation().finite_vertices_begin();
	   vit != triangulation().finite_vertices_end(); ++vit) {
	if (vit->info() != Event_key()) {
	  kdel_.simulator()->delete_event(vit->info());
	  vit->info()=  Event_key();
	}
      }
      for (typename RPMap::iterator it = redundant_points_.begin(); it != redundant_points_.end(); ++it) {
	kdel_.simulator()->delete_event(it->second);
	it->second= Event_key();
      }
      kdel_.set_has_certificates(false);
    }
    else {
      kdel_.set_has_certificates(true);
      if (kdel_.triangulation().dimension()==3) {
	// must be first so the vertex handles are set
	CGAL_KINETIC_LOG(LOG_LOTS, "Setting up certificates.\n");
	for (typename Triangulation::Finite_vertices_iterator vit= triangulation().finite_vertices_begin();
	     vit != triangulation().finite_vertices_end(); ++vit) {
	  if (kdel_.is_degree_4( vit)) {
	    handle_vertex(vit);
	  }
	}
	for (typename RCMap::iterator it= redundant_cells_.begin(); it != redundant_cells_.end(); ++it) {
	  CGAL_KINETIC_LOG(LOG_LOTS, "On init " << it->second << " is redundant" << std::endl);
	  typename Triangulation::Cell_handle h= it->first;
	  handle_redundant(it->second, it->first);
	}
      } else {
	CGAL_KINETIC_LOG(LOG_LOTS, "Triangulation does not have dimension 3.\n");
      }
    }

  }

  bool has_certificates() const
  {
    return kdel_.has_certificates();
  }


protected:

  Cell_handle get_cell_handle(Point_key k) const {
    CGAL_precondition(redundant_points_.find(k) != redundant_points_.end());
    if (redundant_points_.find(k)->second == kdel_.simulator()->null_event()) {
      for (typename RCMap::const_iterator it = redundant_cells_.begin();
	   it != redundant_cells_.end(); ++it){
	if (it->second == k) return it->first;
      }
      CGAL_assertion(0);
      return Cell_handle();
    } else {
      return kdel_.simulator()->template event<Non_vertex_event>(redundant_points_.find(k)->second).cell_handle();
    }
  }

  void set_listener(Listener *l) {
    listener_= l;
  }
  Listener* listener() const
  {
    return listener_;
  }
  void audit_structure() const
  {
    if (!has_certificates()) {
      for (typename RPMap::const_iterator it= redundant_points_.begin(); 
	   it != redundant_points_.end(); ++it) {
	CGAL_assertion(it->second==Event_key());
      }

      for (typename RCMap::const_iterator it= redundant_cells_.begin(); 
	   it != redundant_cells_.end(); ++it) {
	Cell_handle h= it->first;
	Point_key k=it->second;
	CGAL_assertion(triangulation().locate(k)==h);
      }
    } else {
      for (typename Triangulation::Finite_vertices_iterator vit= triangulation().finite_vertices_begin();
	   vit != triangulation().finite_vertices_end(); ++vit) {
	if (triangulation().degree(vit) == 4) {
	  CGAL_assertion_code(Point_key k= vit->point());
	  // it could be infinite
	  // !! for VC
	  CGAL_assertion(vit->info() != Event_key() || !k.is_valid());
	}
	else {
	  CGAL_assertion(vit->info() == Event_key());
	}
	CGAL_assertion(redundant_points_.find(vit->point())== redundant_points_.end());
      }
      CGAL_assertion(triangulation().infinite_vertex()->info() == Event_key());
      for (typename RPMap::const_iterator it= redundant_points_.begin(); it != redundant_points_.end(); ++it) {
	CGAL_assertion(kdel_.vertex_handle(it->first) == Vertex_handle());
	Cell_handle ch= get_cell_handle(it->first);
	typename RCMap::const_iterator beg= redundant_cells_.lower_bound(ch);
	typename RCMap::const_iterator end= redundant_cells_.upper_bound(ch);
	bool found=false;
	for (; beg != end; ++beg) {
	  if (beg->second == it->first) {
	    found=true;
	    break;
	  }
	}
	CGAL_assertion(found);
      } 

      for (typename RCMap::const_iterator it= redundant_cells_.begin(); 
	   it != redundant_cells_.end(); ++it) {
	Point_key pk= it->second;
	Cell_handle ch= it->first;
	CGAL_assertion(redundant_points_.find(pk) != redundant_points_.end());
	Event_key k= redundant_points_.find(pk)->second;

	Cell_handle ech= get_cell_handle(pk);
	CGAL_assertion(ch== ech);
      }
    }
  }
protected:
  //! also much check for vertex_events
  void flip(const typename Triangulation::Edge &edge) {
    typename Triangulation::Facet f= kdel_.flip(edge);

    on_geometry_changed();
  }

  void flip(const typename KDel::Facet &flip_facet) {
    typename Triangulation::Edge edge=  kdel_.flip(flip_facet);

    on_geometry_changed();
  }

  CGAL::Sign orientation(Point_key k, Cell_handle h) const {
    CGAL::Sign ret=CGAL::POSITIVE;
    for (int i=0; i< 4; ++i) {
      typename Triangulation::Facet f(h, i);
      typename Traits::Kinetic_kernel::Weighted_orientation_3 w3;
      w3=kdel_.simulation_traits_object().kinetic_kernel_object().weighted_orientation_3_object();
      CGAL::Sign sn= w3(point(internal::vertex_of_facet(f,0)->point()),
			point(internal::vertex_of_facet(f,1)->point()),
			point(internal::vertex_of_facet(f,2)->point()),
			point(k),
			kdel_.simulation_traits_object().simulator_handle()->current_time());
      if (sn ==CGAL::ZERO) {
	CGAL_KINETIC_LOG(LOG_SOME, "Point " << k << " lies on face ") ;
	CGAL_KINETIC_LOG_WRITE(LOG_SOME, internal::write_facet( f, LOG_STREAM));
	CGAL_KINETIC_LOG(LOG_SOME, "\nPoint trajectory is  " << point(k)  << std::endl) ;
	CGAL_KINETIC_LOG(LOG_SOME, "Triangle 0  " << point(internal::vertex_of_facet(f,0)->point())  << std::endl) ;
	CGAL_KINETIC_LOG(LOG_SOME, "Triangle 1  " << point(internal::vertex_of_facet(f,1)->point())  << std::endl) ;
	CGAL_KINETIC_LOG(LOG_SOME, "Triangle 2  " << point(internal::vertex_of_facet(f,2)->point())  << std::endl) ;
	ret=CGAL::ZERO;
      } else if (sn==CGAL::NEGATIVE) {
	ret=CGAL::NEGATIVE;
	return ret;
      }
    }
    return ret;
  }

  void handle_redundant(Point_key k, Cell_handle h, Root_stack s) {
    CGAL_KINETIC_LOG(LOG_LOTS, "Handle redundant " << k << " ") ;
    CGAL_KINETIC_LOG_WRITE(LOG_LOTS, internal::write_cell( h, LOG_STREAM));
    CGAL_KINETIC_LOG(LOG_LOTS, std::endl);
    CGAL_precondition(orientation(k,h) != CGAL::NEGATIVE);
    CGAL_precondition(redundant_points_[k]==Event_key());


   
    Time pst;
    /*if (!ps.empty()) pst = ps.top();
      else pst= std::numeric_limits<Time>::infinity();*/
    if (s.will_fail()) {
      pst=s.failure_time();
    } else {
      pst= kdel_.simulator()->end_time();
    }

    int first=0;
    for (unsigned int i=0; i< 4; ++i) {
      typename Triangulation::Facet f(h, i);
      // order matters
      Root_stack cs
	= kdel_.orientation_object()(point(internal::vertex_of_facet(f,0)->point()),
				     point(internal::vertex_of_facet(f,1)->point()),
				     point(internal::vertex_of_facet(f,2)->point()),
				     point(k),
				     kdel_.simulator()->current_time(),
				     kdel_.simulator()->end_time());
      if (cs.will_fail() && cs.failure_time() < pst) {
	pst= cs.failure_time();
	s=cs;
	first= i+1;
      }
    }

    if (pst < kdel_.simulator()->end_time()) {
      s.pop_failure_time();
      if (first==0 ) {
	CGAL_KINETIC_LOG(LOG_LOTS, "Making push certificate for " << k << std::endl);
	redundant_points_[k]= kdel_.simulator()->new_event(pst, Push_event(s, k, h, this));
      } else {
	CGAL_KINETIC_LOG(LOG_LOTS, "Making move certificate for " << k << std::endl);
	redundant_points_[k]= kdel_.simulator()->new_event(pst, Move_event(s, k, h, first-1, this));
      }
    } else {
      redundant_points_[k]= kdel_.simulator()->null_event();
    }

    for( typename RCMap::iterator it = redundant_cells_.begin();  it != redundant_cells_.end(); ++it) {
      CGAL_assertion(it->second != k);
    }

    redundant_cells_.insert(typename RCMap::value_type(h, k));
  }

  /*
    if (s.will_fail()) {
    Time t= s.failure_time();
    s.pop_failure_time();
    redundant_points_[k]= kdel_.simulator()->new_event(t, Push_event(s, vh, this));
    } else {
    return kdel_.simulator()->null_event();
    }
  */

  void handle_redundant(Point_key k, typename Triangulation::Cell_handle h) {
    Root_stack ps
      = kdel_.power_test_object()(point(h->vertex(0)->point()),
				  point(h->vertex(1)->point()),
				  point(h->vertex(2)->point()),
				  point(h->vertex(3)->point()),
				  point(k),
				  kdel_.simulator()->current_time(),
				  kdel_.simulator()->end_time());
    handle_redundant(k,h,ps);
  }


  bool try_handle_redundant(Point_key k, typename Triangulation::Cell_handle h) {
    CGAL_KINETIC_LOG(LOG_LOTS, "Trying handle redundant " << k << " ") ;
    CGAL_KINETIC_LOG_WRITE(LOG_LOTS, internal::write_cell( h, LOG_STREAM));
    CGAL_KINETIC_LOG(LOG_LOTS, std::endl);
    if (orientation(k,h) != CGAL::NEGATIVE) {
      handle_redundant(k,h);
      return true;
    } else {
      return false;
    }
  }

  void handle_vertex(typename Triangulation::Vertex_handle vh, Root_stack &s) {
    if (s.will_fail()) {
      Time t= s.failure_time();
      s.pop_failure_time();
      vh->info()= kdel_.simulator()->new_event(t, Pop_event(s, vh, this));
    } else {
      vh->info()= kdel_.simulator()->null_event();
    }
  }

  void handle_vertex(typename Triangulation::Vertex_handle vh) {
    if (vh== triangulation().infinite_vertex()) return;
    CGAL_precondition( internal::has_degree_4(triangulation(), vh));
    CGAL_precondition( vh->info() == Event_key());

    typename Triangulation::Cell_handle ch= vh->cell();
    typename Triangulation::Facet f(ch, ch->index(vh));
    std::vector<typename Triangulation::Vertex_handle> n(4);
   
    for (int i=0; i<3; ++i) {
      n[i]= internal::vertex_of_facet(f,i);
      if (n[i]== triangulation().infinite_vertex()) {
	vh->info() = kdel_.simulator()->null_event();
	return;
      }
    }
    int ind= (f.second+1)%4;
    // some vertex on facet
    n[3] = triangulation().mirror_vertex(ch, ind);
    if (n[3]== triangulation().infinite_vertex()) {
      vh->info() = kdel_.simulator()->null_event();
      return;
    }

    CGAL_KINETIC_LOG(LOG_LOTS, "Making D4 certificate for " << n[0]->point() << n[1]->point()
		     << n[2]->point() << n[3]->point() << " around " << vh->point() << std::endl);

   
    //! The order is switched to invert the predicate since we want it to fail when it goes outside
    Root_stack s
      = kdel_.power_test_object()(point(n[1]->point()),
				  point(n[0]->point()),
				  point(n[2]->point()),
				  point(n[3]->point()),
				  point(vh->point()),
				  kdel_.simulator()->current_time(),
				  kdel_.simulator()->end_time());
    return handle_vertex(vh, s);
  }


  void on_geometry_changed() {
    if (listener_!= NULL) {
      listener_->new_notification(Listener::TRIANGULATION);
    }
    CGAL_KINETIC_LOG(LOG_LOTS, *this);
    audit_structure();
  }

  typename MPT::Data point(Point_key k) const {
    return kdel_.moving_object_table()->at(k);
  }

  // clean vertex events, gather redundant points

  // create vertex events, try to insert redundant points


  void destroy_cell(typename Triangulation::Cell_handle h) {
    for (unsigned int i=0; i<4; ++i) {
      if (h->vertex(i)->info() != Event_key()) {
	kdel_.simulator()->delete_event(h->vertex(i)->info());
	h->vertex(i)->info() == Event_key();
      }
    }
    typename RCMap::iterator beg= redundant_cells_.lower_bound(h);
    typename RCMap::iterator end= redundant_cells_.upper_bound(h);
    for (; beg != end; ++beg) {
      unhandled_keys_.push_back(beg->second);
    }
    redundant_cells_.erase(beg,end);
  }

  void create_cell(typename Triangulation::Cell_handle h) {
    for (unsigned int i=0; i< 4; ++i){
      if (h->vertex(i)->info() == Event_key() && kdel_.is_degree_4(h->vertex(i))){
	handle_vertex(h->vertex(i));
      }
    }
    for (typename std::list<Point_key>::iterator it=unhandled_keys_.begin();
	 it != unhandled_keys_.end(); ++it){
      if (try_handle_redundant(*it, h)) {
	typename std::list<Point_key>::iterator p=it;
	--p;
	unhandled_keys_.erase(it);
	it=p;
      }
    }
  }

  void remove_redundant(typename Triangulation::Cell_handle h, Point_key k) {
    typename RCMap::iterator beg = redundant_cells_.lower_bound(h);
    typename RCMap::iterator end = redundant_cells_.upper_bound(h);
    for (; beg != end; ++beg) {
      if (beg->second == k) {
	redundant_cells_.erase(beg);
	return;
      }
    }
    CGAL_assertion(0);
  }

 

  KDel kdel_;
  Simulator_listener siml_;
  Moving_point_table_listener motl_;
  Listener *listener_;
  RPMap redundant_points_;
  RCMap redundant_cells_;
  std::list<Point_key> unhandled_keys_;
  //typename P::Instantaneous_kernel::Orientation_3 po_;
  // typename P::Kinetic_kernel::Weighted_orientation_3 por_;
};

template <class Traits, class Triang, class Visit>
std::ostream &operator<<(std::ostream &out, const Regular_triangulation_3<Traits, Triang, Visit> &rt)
{
  rt.write(out);
  return out;
}


CGAL_KINETIC_END_NAMESPACE

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif
