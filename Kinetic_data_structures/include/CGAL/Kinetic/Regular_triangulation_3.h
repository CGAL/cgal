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

#include <CGAL/Kinetic/listeners.h>

#include <CGAL/use.h>
#include <CGAL/assertions.h>


#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4355) // complaint about using 'this' to
#endif                          // initialize a member


namespace CGAL { namespace Kinetic { namespace internal {

template <class KD>
class Regular_3_pop_event: public Delaunay_event_base_3<KD, typename KD::Root_stack>
{
  typedef Delaunay_event_base_3<KD, typename KD::Root_stack>  P;
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
class Regular_3_non_vertex_event: public Delaunay_event_base_3<KD, typename KD::Root_stack>
{
  typedef Delaunay_event_base_3<KD, typename KD::Root_stack>  P;
public:

  Regular_3_non_vertex_event(const typename KD::Root_stack &s,
			     const typename KD::Point_key &k,
			     const typename KD::Cell_handle &c,
			     KD *kd): P(s,kd), k_(k), cell_(c) {
  }

  Regular_3_non_vertex_event(){}

  std::ostream& write(std::ostream &out) const
  {
    out << "Nothing ";
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
    out << "Move " << P::point() << " from ";
    internal::write_cell(P::cell_handle() , out);
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
    out << "Push " << P::point() << " into ";
    internal::write_cell(P::cell_handle() , out);
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

 

  //friend class CGAL::Delaunay_triangulation_3<typename P::Instantaneous_kernel, TDS>;

};

} } } //namespace CGAL::Kinetic::internal

namespace CGAL { namespace Kinetic {

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
  typedef typename TraitsT::Instantaneous_kernel Instantaneous_kernel;
public:
  typedef Regular_triangulation_3<TraitsT, VisitorT, TriangulationT> This_RT3;
  typedef TraitsT Traits;
  typedef typename Traits::Active_points_3_table::Key Point_key; //here
  typedef typename Traits::Kinetic_kernel::Certificate Root_stack;
protected:
  typedef typename Traits::Active_points_3_table MPT; // here
  typedef typename Traits::Simulator Simulator;
  typedef typename Traits::Simulator::Event_key Event_key;
  typedef typename Traits::Simulator::Time Time;

  
  typedef TriangulationT Delaunay;

  typedef internal::Regular_triangulation_3_types<TraitsT> Types;

  struct Delaunay_visitor {
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

    Delaunay_visitor(This* krt, VisitorT v): krt_(krt), v_(v){}

    This* krt_;
    VisitorT v_;
  };

  friend struct Delaunay_visitor;

  typedef typename Delaunay::Facet Facet;
  typedef typename Delaunay::Edge Edge;

public:
  typedef typename Delaunay::Cell_handle Cell_handle;
  typedef typename Delaunay::Vertex_handle Vertex_handle;
  typedef VisitorT Visitor;


  friend class internal::Delaunay_event_base_3<This, Root_stack>;
  typedef internal::Delaunay_3_edge_flip_event<This, Root_stack> Edge_flip;
  friend class internal::Delaunay_3_edge_flip_event<This, Root_stack>;
  typedef internal::Delaunay_3_facet_flip_event<This, Root_stack> Facet_flip;
  friend class internal::Delaunay_3_facet_flip_event<This, Root_stack>;

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
    
    typedef TriangulationT Triangulation;
    typedef typename This_RT3::Edge_flip Edge_flip;
    typedef typename This_RT3::Facet_flip Facet_flip;
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

    Base_traits(This_RT3 *t, const TraitsT &tr): TraitsT(tr), wr_(t) {}

    This_RT3* wrapper_handle() {
      return wr_;
    }
    const This_RT3* wrapper_handle() const
    {
      return wr_;
    }

    Active_points_3_table* active_points_3_table_handle() const {
      return TraitsT::active_points_3_table_handle(); // here
    }

    This_RT3 *wr_;
  };

  typedef internal::Delaunay_triangulation_base_3<Base_traits, Delaunay_visitor> KDel;

  CGAL_KINETIC_DECLARE_LISTENERS(typename TraitsT::Simulator,
				 typename Traits::Active_points_3_table)

public:
  

  Regular_triangulation_3(Traits tr, Visitor v= Visitor()): kdel_(Base_traits(this, tr), Delaunay_visitor(this, v)) {
    CGAL_KINETIC_INITIALIZE_LISTENERS(tr.simulator_handle(),
				      tr.active_points_3_table_handle());
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

  CGAL_KINETIC_LISTENER1(TRIANGULATION)
  public:

  void audit_move(Event_key k, Point_key pk, Cell_handle h, int) const {
    CGAL_assertion(kdel_.vertex_handle(pk) == Vertex_handle());
    CGAL_assertion(redundant_points_.find(pk) != redundant_points_.end());
    CGAL_assertion(redundant_points_.find(pk)->second == k);
    audit_redundant(pk, h);
  }
  
  void audit_push(Event_key k, Point_key pk, Cell_handle h) const {
    CGAL_assertion(kdel_.vertex_handle(pk) == Vertex_handle());
    CGAL_assertion(redundant_points_.find(pk) != redundant_points_.end());
    CGAL_assertion(redundant_points_.find(pk)->second == k);
    audit_redundant(pk, h);
  }
  void audit_pop(Event_key k, Vertex_handle vh) const {
    CGAL_assertion_code(if (vh->info() != k) std::cerr << vh->info() << std::endl << k << std::endl);
    CGAL_assertion(vh->info() == k);
  }


  void audit_redundant(Point_key pk, Cell_handle h) const {
    CGAL_LOG(Log::LOTS, "Auditing redundant of " << pk << std::endl);
    CGAL_assertion_code(bool found=false);
    for (typename RCMap::const_iterator cur= redundant_cells_.begin();
	 cur != redundant_cells_.end(); ++cur){
      if (cur->first == h && cur->second == pk) {
	CGAL_assertion_code(found=true);
      } else {
	CGAL_assertion(cur->second != pk);
      }
    }
    CGAL_assertion(found);
  }
  void audit() const
  {
    CGAL_LOG(Log::LOTS, "Verifying regular.\n");
    //if (!has_certificates()) return;
    CGAL_LOG(Log::LOTS, *this << std::endl);
    //P::instantaneous_kernel().set_time(P::simulator()->audit_time());
    kdel_.audit();
    audit_structure();
    //  RPMap redundant_points_;
    // RCMap redundant_cells_;

   
    Triangulation tri= triangulation();
    for (typename RPMap::const_iterator it= redundant_points_.begin(); it != redundant_points_.end(); ++it) {
      //Vertex_handle vh= tri.insert(it->first);
      //if (vh != Vertex_handle()) 
      CGAL_assertion_code(Vertex_handle rvh= tri.insert(it->first));
      CGAL_assertion(rvh == Vertex_handle() || rvh->point() != it->first);
      //if (has_certificates()) {
      CGAL_assertion_code(Cell_handle ch= get_cell_handle(it->first));
      CGAL_assertion_code( typename Instantaneous_kernel::Current_coordinates cco= triangulation().geom_traits().current_coordinates_object());
      CGAL_assertion_code(Cell_handle lh= triangulation().locate(it->first));
      CGAL_assertion(lh == ch
		     || cco(it->first).point() == cco(lh->vertex(0)->point()).point()
		     || cco(it->first).point() == cco(lh->vertex(1)->point()).point()
		     || cco(it->first).point() == cco(lh->vertex(2)->point()).point()
		     || cco(it->first).point() == cco(lh->vertex(3)->point()).point());
	//}
    }
  
  }

  void write(std::ostream &out) const {
    if (triangulation().dimension() != 3) return;
    kdel_.write(out);
    for (typename Triangulation::Finite_vertices_iterator vit= triangulation().finite_vertices_begin();
	 vit != triangulation().finite_vertices_end(); ++vit) {
      if (kdel_.is_degree_4(vit)) {
	out << vit->point() << ": " << vit->info() << std::endl;
      } else if (!kdel_.is_degree_4(vit) && vit->info() != Event_key()) {
	out << vit->point() << "******: " << vit->info() << std::endl;
      }
    }
    out << "Redundant points: ";
    for (typename RPMap::const_iterator it= redundant_points_.begin(); it != redundant_points_.end();
	 ++it) {
      out << it->first << " ";
    }
    out << std::endl;
    typename Delaunay::Cell_handle last;
    out << "Redundant cells: \n";
    for (typename RCMap::const_iterator it= redundant_cells_.begin(); it != redundant_cells_.end();
	 ++it) {
      if (it->first != last) {
	last= it->first;
	internal::write_cell(last, out);
	out << ": ";
      }
      out << it->second << std::endl;
    }
    out << std::endl;
    
  }



  void push(Point_key k, typename Triangulation::Cell_handle h, Root_stack rs) {
    CGAL_LOG(Log::LOTS, "Pushing " << k << " into cell ");
    CGAL_LOG_WRITE(Log::LOTS, internal::write_cell(h, LOG_STREAM));
    CGAL_LOG(Log::LOTS, std::endl);
    
    //redundant_points_.erase(k);

    remove_redundant(h,k);
    
    //typename Triangulation::Vertex_handle vh= kdel_.push_vertex(k, h);
    kdel_.clean_cell(h);
    kdel_.visitor().pre_insert_vertex(k, h);
    // into max dim simplex?
    Vertex_handle vh=kdel_.triangulation().tds().insert_in_cell( h);
    vh->set_point(k);
    kdel_.set_vertex_handle(k, vh);
    handle_vertex(vh, rs);

    std::vector<Cell_handle> ics;
    triangulation().incident_cells(vh, std::back_insert_iterator<std::vector<Cell_handle> >(ics));
    CGAL_postcondition(ics.size() == 4);
    for (unsigned int j=0; j< ics.size(); ++j) {
      Cell_handle cc= ics[j];
      kdel_.handle_new_cell(cc);
    }
    kdel_.visitor().post_insert_vertex(vh);
 
    on_geometry_changed();
  }

  void pop(typename Triangulation::Vertex_handle vh, const Root_stack &rs) {
    CGAL_LOG(Log::LOTS, "Popping " << vh->point() << std::endl);
   
    Point_key k= vh->point();
    vh->info()= Event_key();
    typename Triangulation::Cell_handle h= kdel_.pop_vertex(vh);
    CGAL_precondition(redundant_points_[k]==Event_key());
    redundant_cells_.insert(typename RCMap::value_type(h,k));
    handle_redundant(k, h, rs);
    /*if (!success) {
      std::cerr << "dropped a vertex when popped.\n";
      redundant_points_[k]=kdel_.simulator()->null_event();
      }*/
    //CGAL_postcondition(success);
    on_geometry_changed();
  }

  void move(Point_key k, typename Triangulation::Cell_handle h, int dir, const Root_stack &rs) {
    kdel_.visitor().pre_move(k,h);
    CGAL_LOG(Log::LOTS, "Moving " << k << " from ");
    CGAL_LOG_WRITE(Log::LOTS, internal::write_cell(h, LOG_STREAM));
    CGAL_LOG(Log::LOTS, " to ");
    CGAL_LOG_WRITE(Log::LOTS, internal::write_cell(h->neighbor(dir), LOG_STREAM ));
    CGAL_LOG(Log::LOTS, std::endl);
    typename Triangulation::Cell_handle neighbor = h->neighbor(dir);
    
    bool hinf=false;
    for (unsigned int i=0; i<4; ++i) {
      if (neighbor->vertex(i)== triangulation().infinite_vertex()) {
	hinf=true;
	break;
      }
    }
    if (hinf) {
      //insert(k, neighbor);
      push(k, h, rs);
    } else {
      remove_redundant(h, k);
      CGAL_precondition(redundant_points_[k]==Event_key());
      redundant_cells_.insert(typename RCMap::value_type(neighbor, k));
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
    CGAL_error();
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
      CGAL_precondition(redundant_points_[k]==Event_key());
      handle_redundant(k, h);
    }
  }

  void insert(Point_key k, Cell_handle h) {
    // almost the same as push
    // if insertion fails, then handle redundant
    CGAL_LOG(Log::LOTS, "Inserth " << k << std::endl);
    kdel_.set_instantaneous_time();
    typename Instantaneous_kernel::Current_coordinates cco= triangulation().geom_traits().current_coordinates_object();
    typename Triangulation::Vertex_handle vh;
    if (h!= Cell_handle()) {
      for (unsigned int i=0; i<4; ++i) {
	if (h->vertex(i) != Vertex_handle()
	    && h->vertex(i)->point() != Point_key() 
	    && cco(h->vertex(i)->point()).point() == cco(k).point()) {
	  CGAL_LOG(Log::SOME, "Point " << k << " is on point " 
			   << h->vertex(i)->point() << "\n");
	  vh= h->vertex(i);
	  break;
	}
      }
      if (vh== Vertex_handle()) {
	vh=kdel_.insert(k, h);
      } 
    } else {
      vh= kdel_.insert(k);
    }
   
    
    post_insert(k,vh, h);
  }

  void insert(Point_key k) {
    // almost the same as push
    // if insertion fails, then handle redundant
    CGAL_LOG(Log::LOTS, "Insert " << k << std::endl);

    kdel_.set_instantaneous_time();
    Cell_handle h= triangulation().locate(k);

    return insert(k, h);
  }

  void post_insert(Point_key k, Vertex_handle vh, Cell_handle h) {
    if (vh != Vertex_handle() && vh->point() != k) {
      if (!has_certificates()) {
	unhandled_keys_.push_back(k);
	return;
      } else {
	typename Instantaneous_kernel::Current_coordinates 
	  cco= triangulation().geom_traits().current_coordinates_object();
	if (cco(vh->point()).weight() < cco(k).weight()) {
	  // swap them
	  Point_key rp = kdel_.replace_vertex(vh, k);
	  degen_handle_redundant(rp,vh);
	  if (kdel_.has_certificates() && kdel_.is_degree_4(vh)) {
	    handle_vertex(vh);
	  }
	} else {
	  vh= Vertex_handle();
	  degen_handle_redundant(k, vh);
	  // need to try various cells, not just one
	}
      }
    } else if (vh == Vertex_handle()) {
      CGAL_precondition(triangulation().geom_traits().time() 
			==kdel_.simulator()->current_time());
     
      if (h== Cell_handle()) {
	h = triangulation().locate(k);
      }
      CGAL_precondition(redundant_points_[k]==Event_key());
      redundant_cells_.insert(typename RCMap::value_type(h, k));
      handle_redundant(k,h);
    } else if (kdel_.has_certificates() && kdel_.is_degree_4(vh)){
      handle_vertex(vh); 
    }

    on_geometry_changed();
  }

  void degen_handle_redundant(Point_key k, Vertex_handle vh) {
    std::vector<Cell_handle> ics;
    triangulation().incident_cells(vh, std::back_inserter(ics));
    kdel_.set_instantaneous_time(true);
    for (unsigned int i=0; i< ics.size(); ++i) {
      if (try_handle_redundant(k, ics[i])) return;
    }
    CGAL_error();
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
	CGAL_LOG(Log::LOTS, "Setting up certificates.\n");
	for (typename Triangulation::Finite_vertices_iterator vit= triangulation().finite_vertices_begin();
	     vit != triangulation().finite_vertices_end(); ++vit) {
	  /*if (kdel_.is_degree_4( vit)) {
	    handle_vertex(vit);
	    }*/
	  CGAL_assertion(!kdel_.is_degree_4(vit) || vit->info() != Event_key());
	}
	for (typename RCMap::iterator it= redundant_cells_.begin();
	     it != redundant_cells_.end(); ++it) {
	  CGAL_LOG(Log::LOTS, "On init " << it->second 
			   << " is redundant" << std::endl);
	  CGAL_precondition(redundant_points_[it->second]==Event_key());
	  handle_redundant(it->second, it->first);
	}
	CGAL_assertion(unhandled_keys_.empty());
      } else {
	CGAL_LOG(Log::LOTS, "Triangulation does not have dimension 3.\n");
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
    if (redundant_points_.find(k)->second == kdel_.simulator()->null_event()
	|| !has_certificates()) {
      for (typename RCMap::const_iterator it = redundant_cells_.begin();
	   it != redundant_cells_.end(); ++it){
	if (it->second == k) return it->first;
      }
      CGAL_error();
      return Cell_handle();
    } else {
      return kdel_.simulator()->template event<Non_vertex_event>(redundant_points_.find(k)->second).cell_handle();
    }
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
	Cell_handle lh= triangulation().locate(k);
	if (lh != h) {
	  typename Instantaneous_kernel::Current_coordinates cco= triangulation().geom_traits().current_coordinates_object();
	  bool found=false;
	  for (unsigned int i=0; i<4; ++i) {
	    if (lh->vertex(i)->point() != Point_key() && cco(lh->vertex(i)->point()).point() 
		== cco(k).point()) {
	      found=true;
	    }
	  }
	  CGAL_assertion(found);
          CGAL_USE(found);
	}

	
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
        CGAL_USE(found);
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
    kdel_.flip(edge);

    on_geometry_changed();
  }

  void flip(const typename KDel::Facet &flip_facet) {
    kdel_.flip(flip_facet);

    on_geometry_changed();
  }

  CGAL::Sign orientation(Point_key k, Cell_handle h) const {
    kdel_.set_instantaneous_time();
    int hinf=-1;
    for (int i=0; i< 4; ++i) {
      if (h->vertex(i)->point() == Point_key()) {
	hinf=i;
      }
    }
    CGAL::Sign ret=CGAL::POSITIVE;
    for (int i=0; i< 4; ++i) {
      if (hinf == -1 ||  hinf == i) {
	typename Triangulation::Facet f(h, i);
	typename Traits::Instantaneous_kernel::Orientation_3 o3= triangulation().geom_traits().orientation_3_object();
	
	CGAL::Sign sn= o3((internal::vertex_of_facet(f,0)->point()),
			  (internal::vertex_of_facet(f,1)->point()),
			  (internal::vertex_of_facet(f,2)->point()),
			  (k));
	if (sn ==CGAL::ZERO) {
	  CGAL_LOG(Log::SOME, "Point " << k << " lies on face ") ;
	  CGAL_LOG_WRITE(Log::SOME, internal::write_facet( f, LOG_STREAM));
	  CGAL_LOG(Log::SOME, "\nPoint trajectory is  " << point(k)  << std::endl) ;
	  CGAL_LOG(Log::SOME, "Triangle 0  " << point(internal::vertex_of_facet(f,0)->point())  << std::endl) ;
	  CGAL_LOG(Log::SOME, "Triangle 1  " << point(internal::vertex_of_facet(f,1)->point())  << std::endl) ;
	  CGAL_LOG(Log::SOME, "Triangle 2  " << point(internal::vertex_of_facet(f,2)->point())  << std::endl) ;
	  ret=CGAL::ZERO;
	} else if (sn==CGAL::NEGATIVE) {
	  ret=CGAL::NEGATIVE;
	  return ret;
	}
      }
    }
    return ret;
  }

  void handle_redundant(Point_key k, Cell_handle h, Root_stack s) {
    CGAL_LOG(Log::LOTS, "Handle redundant " << k << " ") ;
    CGAL_LOG_WRITE(Log::LOTS, internal::write_cell( h, LOG_STREAM));
    CGAL_LOG(Log::LOTS, std::endl);
    CGAL_precondition(orientation(k,h) != CGAL::NEGATIVE);
    CGAL_precondition(redundant_points_[k]==Event_key());
    CGAL_assertion_code(bool found=false);
    for( typename RCMap::iterator it = redundant_cells_.begin();  it != redundant_cells_.end(); ++it) {
      CGAL_assertion_code(if(it->second==k) found=true);
      CGAL_assertion(h== it->first || it->second != k);
    }
    CGAL_assertion(found);
    redundant_points_[k]= Event_key();
    if (!kdel_.has_certificates()) return;
 
    Time pst;
    /*if (!ps.empty()) pst = ps.top();
      else pst= std::numeric_limits<Time>::infinity();*/
    if (s.will_fail()) {
      pst=s.failure_time();
    } else {
      pst= kdel_.simulator()->end_time();
    }
    int hinf=-1;
    for (int i=0; i< 4; ++i) {
      if (h->vertex(i)->point() == Point_key()) {
	hinf=i;
      }
    }
    int first=0;
    for (int i=0; i< 4; ++i) {
      if (hinf == -1 || hinf ==i) {
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
    }
    if (pst < kdel_.simulator()->end_time()) {
      s.pop_failure_time();
      if (first==0 ) {
	CGAL_LOG(Log::LOTS, "Making push certificate for " << k << std::endl);
	redundant_points_[k]= kdel_.simulator()->new_event(pst, Push_event(s, k, h, this));
	CGAL_assertion_code(kdel_.simulator()->audit_event(redundant_points_[k]));
	CGAL_assertion_code(kdel_.simulator()->audit_events());
      } else {
	CGAL_LOG(Log::LOTS, "Making move certificate for " << k << std::endl);
	redundant_points_[k]= kdel_.simulator()->new_event(pst, Move_event(s, k, h, first-1, this));
	CGAL_assertion_code(kdel_.simulator()->audit_event(redundant_points_[k]));
	CGAL_assertion_code(kdel_.simulator()->audit_events());
      }
    } else {
      redundant_points_[k]= kdel_.simulator()->null_event();
    }

   
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
    int hinf=-1;
    for (int i=0; i< 4; ++i) {
      if (h->vertex(i)->point() == Point_key()) {
	hinf=i;
      }
    }
    Root_stack ps;
    if (hinf ==-1 ) {
      ps
	= kdel_.power_test_object()(point(h->vertex(0)->point()),
				    point(h->vertex(1)->point()),
				    point(h->vertex(2)->point()),
				    point(h->vertex(3)->point()),
				    point(k),
				    kdel_.simulator()->current_time(),
				    kdel_.simulator()->end_time());
    }
    handle_redundant(k,h,ps);
  }


  bool try_handle_redundant(Point_key k, typename Triangulation::Cell_handle h) {
    CGAL_LOG(Log::LOTS, "Trying handle redundant " << k << " ") ;
    CGAL_LOG_WRITE(Log::LOTS, internal::write_cell( h, LOG_STREAM));
    CGAL_LOG(Log::LOTS, std::endl);
    if (orientation(k,h) != CGAL::NEGATIVE) {
      CGAL_precondition(redundant_points_[k]==Event_key());
      redundant_cells_.insert(typename RCMap::value_type(h, k));
      handle_redundant(k,h);
      return true;
    } else {
      return false;
    }
  }

  void handle_vertex(typename Triangulation::Vertex_handle vh, Root_stack &s) {
    CGAL_LOG(Log::LOTS, "Updating vertex " << vh->point() << std::endl);
    CGAL_precondition(vh->info() == Event_key());
    if (s.will_fail()) {
      Time t= s.failure_time();
      s.pop_failure_time();
      vh->info()= kdel_.simulator()->new_event(t, Pop_event(s, vh, this));
      CGAL_assertion_code(kdel_.simulator()->audit_event(vh->info()));
      CGAL_assertion_code(kdel_.simulator()->audit_events());
      CGAL_assertion_code(kdel_.simulator()->template event<Pop_event>(vh->info()).audit(vh->info()));
    } else {
      vh->info()= kdel_.simulator()->null_event();
    }
  }

  void handle_vertex(typename Triangulation::Vertex_handle vh) {
    CGAL_LOG(Log::LOTS, "Handling vertex " << vh->point() << std::endl);
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

    CGAL_LOG(Log::LOTS, "Making D4 certificate for " << n[0]->point() << n[1]->point()
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
    CGAL_KINETIC_NOTIFY(TRIANGULATION);
    CGAL_LOG(Log::LOTS, *this);
    audit_structure();
  }

  typename MPT::Data point(Point_key k) const {
    return kdel_.moving_object_table()->at(k);
  }

  // clean vertex events, gather redundant points

  // create vertex events, try to insert redundant points


  void destroy_cell(typename Triangulation::Cell_handle h) {
    CGAL_LOG(Log::LOTS, "Cleaning cell " << h->vertex(0)->point()
		     << " " << h->vertex(1)->point() << " " << h->vertex(2)->point()
		     << " " << h->vertex(3)->point() << std::endl);
    for (unsigned int i=0; i<4; ++i) {
      if (h->vertex(i)->info() != Event_key()) {
	CGAL_LOG(Log::LOTS, "Cleaning vertex " << h->vertex(i)->point() << std::endl);
	kdel_.simulator()->delete_event(h->vertex(i)->info());
	h->vertex(i)->info() = Event_key();
      }
    }
    typename RCMap::iterator beg= redundant_cells_.lower_bound(h);
    typename RCMap::iterator end= redundant_cells_.upper_bound(h);
    for (; beg != end; ++beg) {
      unhandled_keys_.push_back(beg->second);
      kdel_.simulator()->delete_event(redundant_points_[beg->second]);
      redundant_points_.erase(beg->second);
    }
    redundant_cells_.erase(redundant_cells_.lower_bound(h),
			   redundant_cells_.upper_bound(h));
  }

  void create_cell(typename Triangulation::Cell_handle h) {
    CGAL_LOG(Log::LOTS, "Creating cell " << h->vertex(0)->point()
		     << " " << h->vertex(1)->point() << " " << h->vertex(2)->point()
		     << " " << h->vertex(3)->point() << std::endl);
    for (unsigned int i=0; i< 4; ++i){
      if (h->vertex(i)->info() == Event_key() && kdel_.is_degree_4(h->vertex(i))){
	handle_vertex(h->vertex(i));
      }
    }
    for ( int i=0; i< static_cast<int>(unhandled_keys_.size()); ++i) {
      kdel_.set_instantaneous_time(true);
      if (try_handle_redundant(unhandled_keys_[i], h)) {
	unhandled_keys_.erase(unhandled_keys_.begin()+i);
	--i;
      }
    }
  }

  void remove_redundant(typename Triangulation::Cell_handle h, Point_key k) {
    redundant_points_.erase(k);
    typename RCMap::iterator beg = redundant_cells_.lower_bound(h);
    typename RCMap::iterator end = redundant_cells_.upper_bound(h);
    for (; beg != end; ++beg) {
      if (beg->second == k) {
	redundant_cells_.erase(beg);
	return;
      }
    }
    
    for (typename RCMap::iterator cur= redundant_cells_.begin();
	 cur != redundant_cells_.end(); ++cur){
      CGAL_ERROR_WRITE( internal::write_cell( cur->first, LOG_STREAM));
      CGAL_ERROR(": " << cur->second);
      if (cur->second == k) {
	CGAL_assertion_code(Cell_handle ch= cur->first);
	CGAL_assertion(ch==h);
	CGAL_error();
      }
    }
    CGAL_error();
  }

 

  KDel kdel_;
  RPMap redundant_points_;
  RCMap redundant_cells_;
  std::vector<Point_key> unhandled_keys_;
  //typename P::Instantaneous_kernel::Orientation_3 po_;
  // typename P::Kinetic_kernel::Weighted_orientation_3 por_;
};

template <class Traits, class Triang, class Visit>
std::ostream &operator<<(std::ostream &out, const Regular_triangulation_3<Traits, Triang, Visit> &rt)
{
  rt.write(out);
  return out;
}


} } //namespace CGAL::Kinetic

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif
