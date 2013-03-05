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
// $Id$ $Date$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_INTERNAL_DELAUNAY_BASE_3_H
#define CGAL_KINETIC_INTERNAL_DELAUNAY_BASE_3_H

#include <CGAL/Kinetic/basic.h>

#include <CGAL/Kinetic/internal/Triangulation_helper_3.h>
#include <CGAL/Kinetic/internal/triangulation_helpers_3.h>
#include <CGAL/Kinetic/Event_base.h>
//#include <CGAL/Kinetic/internal/Delaunay_cache_3.h>

// STL
#include <map>
#include <set>
#include <iterator>
#include <ostream>
#include <iostream>

//extern int too_late__;
//extern int filtered__;

namespace CGAL { namespace Kinetic { namespace internal {

template <class KD, class RS>
class Delaunay_event_base_3: public Event_base<KD*>
{
  typedef Event_base<KD*> P;
public:
  Delaunay_event_base_3(const RS &s,
			KD *kdel):P(kdel), s_(s){
  }
  //! Default constructor
  /*!  This really should not be used, but I need it for the
    Simulator::event() call, due to the apparent gcc compiler bug.
  */
  Delaunay_event_base_3() {

  }
  void process() {
    // for some reason VC insists that this be there
    CGAL_assertion(0 && "Process called in Delaunay_event_base_3");
  }
  const RS& root_stack() const
  {
    return s_;
  }
  std::ostream& write(std::ostream &out) const
  {
    out << "Delaunay_event";
    return out;
  }
  KD* kdel() const {
    return P::kds();
  }
protected:
  const RS s_;
};

/*template <class K, class R>
std::ostream& operator<<(std::ostream &out, const Delaunay_event_base_3<K, R> &e)
{
  e.write(out);
  return out;
  }*/


template <class KD, class RS>
class Delaunay_3_edge_flip_event: public Delaunay_event_base_3<KD, RS>
{
public:
  typedef Delaunay_event_base_3<KD, RS>   P;
  Delaunay_3_edge_flip_event(const RS &s,
			     const typename KD::Edge &e,
			     KD *kdel):P(s, kdel), e_(e) {
#ifndef NDEBUG
    o_= edge_point(e_,0);
    d_= edge_point(e_,1);
#endif
  }
  void process() {
    P::kdel()->flip(e_);
  }

  static typename KD::Point_key edge_point(const typename KD::Edge &e, int i) {
    return vertex_of_edge(e, i)->point();
  }

  std::ostream& write(std::ostream &out) const
  {
    out << "Flip ";
    internal::write_edge(e_, out);
#if 0
    out << "(" << o_ << d_<<")" << std::flush;
    CGAL_postcondition(o_== vertex_of_edge(e_,static_cast<int>(0))->point());
    CGAL_postcondition(d_== vertex_of_edge(e_,1)->point());
#endif
    return out;
  }
 

  void audit(typename KD::Event_key k) const {
    if (e_.first->edge_label(e_.second,e_.third) != k) {
      CGAL_ERROR("Mismatch, for label " << k << " had event " << e_.first->edge_label(e_.second,e_.third));
    }
    CGAL_assertion(e_.first->edge_label(e_.second,e_.third) == k);
  }
protected:
  const typename KD::Edge e_;
#ifndef NDEBUG
  typename KD::Point_key o_, d_;
#endif
};

/*template <class B, class C, class D>
std::ostream& operator<<(std::ostream &out, const Delaunay_3_edge_flip_event<B, C, D> &e)
{
  e.write(out);
  return out;
  }*/


template <class KD, class RS>
class Delaunay_3_facet_flip_event:  public Delaunay_event_base_3<KD, RS>
{
public:
  typedef Delaunay_event_base_3<KD, RS>   P;
  Delaunay_3_facet_flip_event(const RS &s,
			      const typename KD::Facet &e,
			      KD *kdel):  P(s, kdel),e_(e) {
#ifndef NDEBUG
    a_= vertex_of_facet(e_,0)->point();
    b_= vertex_of_facet(e_,1)->point();
    c_= vertex_of_facet(e_,2)->point();
#endif
  }
  void process() {
    P::kdel()->flip(e_);
  }
  std::ostream& write(std::ostream &out) const
  {
    out << "Flip ";
    write_facet(e_, out);
#if 0
    out << "(" << a_ << b_<<c_ << ")";
    CGAL_postcondition(a_== vertex_of_facet(e_,0)->point());
    CGAL_postcondition(b_== vertex_of_facet(e_,1)->point());
    CGAL_postcondition(c_== vertex_of_facet(e_,2)->point());
#endif
    return out;
  }
 
  void audit(typename KD::Event_key k) const {
    if (e_.first->facet_label(e_.second) != k) {
      CGAL_ERROR("Mismatch, for label " << k << " had event " << e_.first->facet_label(e_.second));
    }
    CGAL_assertion(e_.first->facet_label(e_.second) == k);
  }
protected:
  const typename KD::Facet e_;
#ifndef NDEBUG
  typename KD::Point_key a_, b_, c_;
#endif
};

/*template <class T, class C, class D>
std::ostream& operator<<(std::ostream &out, const Delaunay_3_facet_flip_event<T, C, D> &e)
{
  e.write(out);
  return out;
  }*/


template <class TraitsT, class Visitor>
class Delaunay_triangulation_base_3;

template <class TraitsT, class Visitor>
std::ostream &operator<<(std::ostream &out,
			 const Delaunay_triangulation_base_3<TraitsT, Visitor> &dt);

// Common base class for Delaunay and regular triangulations
template <class TraitsT, class Visitor>
class Delaunay_triangulation_base_3
{
public:
  typedef Delaunay_triangulation_base_3<TraitsT, Visitor> This;

  // KDS typedefs
  typedef typename TraitsT::Simulator Simulator;
  typedef typename TraitsT::Active_points_3_table Moving_object_table;
  typedef typename TraitsT::Kinetic_kernel Kinetic_kernel;
  typedef typename TraitsT::Instantaneous_kernel Instantaneous_kernel;
  //typedef typename Simulator::Time Time;
  typedef typename Moving_object_table::Key Point_key;
  typedef typename Moving_object_table::Data Point;
  typedef typename Simulator::Event_key Event_key;
  typedef typename TraitsT::Side_of_oriented_sphere_3::result_type Certificate;
  //typedef typename Simulator::NT NT;

  // Delaunay typedefs
  typedef Triangulation_helper_3<typename TraitsT::Triangulation> Tri;

  // Triangulation members
  typedef typename Tri::Facet Facet;
  typedef typename Tri::Edge Edge;
  typedef typename Tri::Facet_circulator Facet_circulator;
  typedef typename Tri::Cell_circulator Cell_circulator;
  typedef typename Tri::Finite_edges_iterator Finite_edges_iterator;
  typedef typename Tri::Finite_facets_iterator Finite_facets_iterator;
  typedef typename Tri::Edge_iterator Edge_iterator;
  //typedef typename Tri::Facet_iterator Facet_iterator;
  typedef typename Tri::Vertex_handle Vertex_handle;
  typedef typename Tri::Cell_handle Cell_handle;
  typedef typename Tri::All_cells_iterator All_cells_iterator;
  typedef typename Tri::All_edges_iterator All_edges_iterator;
  typedef typename Tri::All_facets_iterator All_facets_iterator;
  typedef typename Tri::Finite_vertices_iterator Finite_vertices_iterator;

  // accessory types
  typedef typename TraitsT::Facet_flip Facet_flip_event;
  typedef typename Facet_flip_event::P Event_base;
  typedef typename TraitsT::Edge_flip Edge_flip_event;
  typedef Tri Triangulation;

  Delaunay_triangulation_base_3(TraitsT tr, Visitor v): tr_(tr),
							triangulation_(tr.instantaneous_kernel_object()),
							soc_(tr.side_of_oriented_sphere_3_object()),
							o3_(tr.orientation_3_object()),
							v_(v) {
    if (0) print();
    has_certificates_=false;
  }

  const TraitsT& simulation_traits_object() const {return tr_;}


  const Visitor& visitor() const
  {
    return v_;
  }
  Visitor& visitor() {
    return v_;
  }

  //! Just write the objects in order;
  template <class Stream>
  void write(Stream &out) const
  {
    //out << triangulation_;
    for (All_cells_iterator cit= triangulation_.all_cells_begin(); cit != triangulation_.all_cells_end();
	 ++cit) {
      Cell_handle h= cit;
      for (unsigned int i=0; i<4; ++i) {
	out << h->vertex(i)->point() << " ";
      }
      out << "--";
      for (unsigned int i=0; i<4; ++i) {
	Facet f(h, i);
	out << triangulation_.label(f) << " ";
      }
      out << std::endl;
    }
    for (All_edges_iterator eit = triangulation_.all_edges_begin(); eit != triangulation_.all_edges_end(); ++eit) {
      triangulation_.write_labeled_edge(*eit, out);
      if (is_degree_3(*eit)) {
	//out << " " << triangulation_.label(*eit);
      }
      else {
	if (has_event(*eit)) {
	  out << " ?? " ;           //<< triangulation_.label(*eit);
	}
	//CGAL_assertion(triangulation_.label(*eit)== Event_key::null());
      }
      out << std::endl;
    }
    for (All_facets_iterator eit = triangulation_.all_facets_begin(); eit != triangulation_.all_facets_end(); ++eit) {
      triangulation_.write_labeled_facet(*eit, out);
      if (!has_degree_3_edge(*eit)) {
	//out << " " << triangulation_.label(*eit);
      }
      else {
	//CGAL_assertion(triangulation_.label(*eit)== Event_key::null());
	if (has_event(*eit)) {
	  out << " ?? ";            //<< triangulation_.label(*eit);
	}
      }
      out << std::endl;
    }
  }
  bool is_degree_3(const Edge &e) const {
    return triangulation_.has_degree_3(e);
  }

  bool is_degree_4(Vertex_handle vh) const {
    return triangulation_.degree(vh) == 4;
  }

  bool has_degree_3_edge(const Facet& f) const {
    return triangulation_.has_degree_3(f);
  }

  bool has_degree_4_vertex(const Facet& f) const {
    return is_degree_4(triangulation_.vertex(f,0))
      || is_degree_4(triangulation_.vertex(f,1))
      || is_degree_4(triangulation_.vertex(f,2));
  }

  bool has_degree_4_vertex(const Edge& f) const {
    return is_degree_4(triangulation_.vertex(f,0))
      || is_degree_4(triangulation_.vertex(f,1));
  }

  bool print() const
  {
    write(std::cout);
    return true;
  }

  /*!
    Delete all incident face and edge events (have to push edge events on to facets).
    Insert point.
    Set edges to never fails. Set outside facets.
  */
  /*Vertex_handle push_vertex(Point_key k, Cell_handle c) {
    clean_cell(c);
    v_.pre_insert_vertex(k, c);
    // into max dim simplex?
    Vertex_handle vh=triangulation_.tds().insert_in_cell( c);
    vh->set_point(k);
    set_vertex_handle(k, vh);
    std::vector<Cell_handle> ics;
    triangulation_.incident_cells(vh, std::back_insert_iterator<std::vector<Cell_handle> >(ics));
    CGAL_postcondition(ics.size() == 4);
    for (unsigned int j=0; j< ics.size(); ++j) {
      Cell_handle cc= ics[j];
      handle_new_cell(cc);
    }
    v_.post_insert_vertex(vh);
    return vh;
    }*/

  Cell_handle pop_vertex(Vertex_handle v) {
    v_.pre_remove_vertex(v);
    CGAL_precondition(is_degree_4(v));
    std::vector<Cell_handle> ics;
    triangulation_.incident_cells(v, std::back_insert_iterator<std::vector<Cell_handle> >(ics));
    Point_key k= v->point();
    for (unsigned int j=0; j< ics.size(); ++j) {
      clean_cell(ics[j]);
    }
   
    set_vertex_handle(v->point(), NULL);

    Cell_handle h= triangulation_.tds().remove_from_maximal_dimension_simplex(v);
    handle_new_cell(h);
    v_.post_remove_vertex(k);

    return h;
  }

  void delete_vertex(Point_key) {
    CGAL_error();
  }


  Vertex_handle change_vertex(Point_key k) {
    if (!has_certificates()) return NULL;
    //triangulation_.geom_traits().set_time(simulator()->current_time_as_nt());
    std::vector<Cell_handle> incident_cells;
    triangulation_.incident_cells(vertex_handle(k), back_inserter(incident_cells));
    for (typename std::vector<Cell_handle>::iterator it= incident_cells.begin(); 
	 it != incident_cells.end();
	 ++it) {
      clean_cell(*it);
    }

    for (typename std::vector<Cell_handle>::iterator it= incident_cells.begin();
	 it != incident_cells.end(); ++it) {
      handle_new_cell(*it);
    }
    v_.change_vertex(vertex_handle(k));
    return vertex_handle(k);
  }

  /*typename Triangulation::Vertex_handle new_vertex_regular(Point_key k, Cell_handle h=NULL) {
    CGAL_precondition(!has_certificates_);
    typename Simulator::NT nt= simulator()->next_time_representable_as_nt();
    CGAL_precondition(simulator()->current_time() == nt);
    //simulator()->set_current_time(nt);
    triangulation_.geom_traits().set_time(nt);
    typename Triangulation::Vertex_handle vh= triangulation_.insert(k, h);
    v_.create_vertex(vh);
    return vh;
    }*/

  typename Triangulation::Vertex_handle insert(Point_key k) {
    /* NT nt= simulator()->next_time_representable_as_nt();
    simulator()->set_current_time(nt);
    triangulation_.geom_traits().set_time(nt);*/
    set_instantaneous_time();
    //typename Simulator::NT nt= simulator()->next_time_representable_as_nt();
    //simulator()->set_current_time(nt);
    //std::cout << "Locating at time " << triangulation_.geom_traits().time() << std::endl;
    Cell_handle h= triangulation_.locate(k);
    /*if (h != Cell_handle() && h->vertex(0) != Vertex_handle()
	       && h->vertex(1) != Vertex_handle()
	       && h->vertex(2) != Vertex_handle()
	       && h->vertex(3) != Vertex_handle()) {
      std::cout << "Kinetic located in " << h->vertex(0)->point() << " "
		<< h->vertex(1)->point() << " "
		<< h->vertex(2)->point() << " "
		<< h->vertex(3)->point() << std::endl;
    } else {
      std::cout << "Kinetic located outside hull" << std::endl;
      }*/
    return insert(k,h);
  }

  //!
  /*!
    Some old certificate edges will be lost, have to find all conflicts first.

    An important problem case is the edges which were previously degree 3 but are no longer
    after insertion. They can now be degree 4 and will have their certificate deleted properly.
    However, this leaves a face with no certificate (the one or more not destroyed by the new vertex).

    This occurs when a degree 3 edges has a facet added but none destroyed--i.e. a boundary edge
    with degree 3.
  */
  typename Triangulation::Vertex_handle insert(Point_key k, Cell_handle h) {
    //CGAL_precondition(h != NULL);

    std::vector<Facet> bfacets;
    std::vector<Cell_handle> cells;
    std::vector<Facet> ifacets;
    //typename Simulator::NT nt= simulator()->next_time_representable_as_nt();
    //CGAL_precondition(simulator()->current_time() == nt);
    //triangulation_.geom_traits().set_time(nt);
    set_instantaneous_time();
    CGAL_precondition(triangulation_.geom_traits().time() == simulator()->current_time());
    Vertex_handle vh;
    v_.pre_insert_vertex(k, h);
    if (triangulation_.dimension() == 3) {
      triangulation_.find_conflicts(k, h, back_inserter(bfacets), 
				    back_inserter(cells),back_inserter(ifacets));
      if (has_certificates_) {
	for (unsigned int i=0; i < cells.size(); ++i) {
	  clean_cell(cells[i]);
	}
      }
      
   
      //! \todo replace by insert_in_hole

      CGAL_LOG(Log::LOTS, "Inserting " << k << " at time " 
		       << triangulation_.geom_traits().time() << "\n");
      vh=triangulation_.insert_in_hole(k, cells.begin(), cells.end(), 
						     bfacets.front().first, bfacets.front().second);
    } else {
      vh= triangulation_.insert(k,h);
    }
    set_vertex_handle(k,vh);
    if (triangulation_.dimension() == 3 && vh != Vertex_handle()) {
      change_vertex(k);
    }
    
    CGAL_expensive_postcondition(audit_structure());
    v_.post_insert_vertex(vh);
    return vh;
  }

  bool has_certificates() const
  {
    return has_certificates_;
  }

  void set_instantaneous_time(bool after=false) const {
    if (!triangulation_.geom_traits().has_time() || triangulation_.geom_traits().time() != simulator()->current_time()) {
      typename Simulator::NT nt= simulator()->next_time_representable_as_nt();
      
      if (simulator()->current_time() == nt) {
	if (after) {
	  triangulation_.geom_traits().set_time_to_after(nt);
	} else {
	  triangulation_.geom_traits().set_time(nt);
	}
      } else {
	CGAL_LOG(Log::SOME, "Warning, insertion of points at non-rational times is slow.\n");
	if (after) {
	  triangulation_.geom_traits().set_time_to_after(simulator()->current_time());
	} else {
	  triangulation_.geom_traits().set_time(simulator()->current_time());
	}
      }
    }
  }

  void set_has_certificates(bool b) {
    if (!has_certificates_ && b) {
      if (triangulation().dimension() == 3) {
	for (All_edges_iterator eit = triangulation_.all_edges_begin();
	     eit != triangulation_.all_edges_end(); ++eit) {
	  CGAL_assertion(!has_event(*eit));
	}
	for (All_facets_iterator eit = triangulation_.all_facets_begin();
	     eit != triangulation_.all_facets_end(); ++eit) {
	  CGAL_assertion(!has_event(*eit));
	}
	create_all_certificates();
	has_certificates_=true;
      }
    } else if (has_certificates_ && !b) {
      destroy_all_certificates();
      has_certificates_=false;
    }
  }

 void create_all_certificates() {
    CGAL_precondition(!has_certificates_);
 
    for (All_edges_iterator eit = triangulation_.all_edges_begin();
	 eit != triangulation_.all_edges_end(); ++eit) {
      if (is_degree_3(*eit) && !has_degree_4_vertex(*eit)) {
	make_certificate(*eit);
      }
    }
    for (All_facets_iterator eit = triangulation_.all_facets_begin();
	 eit != triangulation_.all_facets_end(); ++eit) {
      if (!has_degree_3_edge(*eit)) {
	make_certificate(*eit);
      }
    }
    for (All_cells_iterator cit= triangulation_.all_cells_begin(); 
	 cit != triangulation_.all_cells_end(); ++cit) {
      v_.create_cell(cit);
    }
  }


  void destroy_all_certificates() {
    //vhs_.clear();
    CGAL_precondition(has_certificates_);
    for (All_edges_iterator eit = triangulation_.all_edges_begin();
	 eit != triangulation_.all_edges_end(); ++eit) {
      Event_key k= triangulation_.label(*eit);
      if ( k != Event_key() ) {
	simulator()->delete_event(k);
	triangulation_.set_label(*eit,Event_key());
      }
    }
    for (All_facets_iterator eit = triangulation_.all_facets_begin();
	 eit != triangulation_.all_facets_end(); ++eit) {
      Event_key k= triangulation_.label(*eit);
      if (k != Event_key() ) {
	simulator()->delete_event(k);
	triangulation_.set_label(*eit,Event_key());
      }
      //}
    }
    for (All_cells_iterator cit= triangulation_.all_cells_begin(); 
	 cit != triangulation_.all_cells_end(); ++cit) {
      v_.destroy_cell(cit);
    }
  }
 


  Facet flip(const Edge &edge) {
    v_.pre_edge_flip(edge);
    
    CGAL_LOG(Log::LOTS,"\n\nFlipping edge ");
    CGAL_LOG(Log::LOTS,edge.first->vertex(edge.second)->point() << "--" 
		     << edge.first->vertex(edge.third)->point() << std::endl);
    CGAL_assertion(triangulation_.tds().is_edge(edge.first, edge.second, edge.third) || print());
    if (has_degree_4_vertex(edge)) {
      CGAL_LOG(Log::LOTS,"dropping edge since endpoint is degree 4\n ");
      triangulation_.set_label(edge, simulator()->null_event());
      return Facet();
    }
    CGAL_assertion(!has_degree_4_vertex(edge));

    Vertex_handle poles[2];
    poles[0]= edge.first->vertex(edge.second);
    poles[1]= edge.first->vertex(edge.third);
    int polesi[2];
    polesi[0] = edge.first->index(poles[0]);
    polesi[1]= edge.first->index(poles[1]);

    typename Simulator::Event_key failed_key = triangulation_.label(edge);
    CGAL_precondition(failed_key.is_valid());
    Certificate ore= extract_root_stack(failed_key);
   
    Facet neighboring_facet= triangulation_.opposite(Facet(edge.first, polesi[0]));

    triangulation_.set_label(edge, typename Simulator::Event_key());

    // handle the cross edges to make sure that they are no longer edge flips
    {
      Cell_circulator cc= triangulation_.incident_cells(edge), ce= cc;
      do {
	Cell_handle h=cc;
	clean_cell(h);
      } while (++cc != ce);
    }

    triangulation_.tds().flip_flippable(edge);
    CGAL_expensive_assertion(labeling_is_valid());
    /*if (verbose){
      write_labeled_state();
      }*/

    CGAL_assertion(triangulation_.tds().is_facet(neighboring_facet.first, neighboring_facet.second));
    Facet internal_facet= triangulation_.opposite(neighboring_facet);

    Cell_handle cells[2];
    cells[1]= internal_facet.first;
    polesi[1]= cells[1]->index(poles[1]);
    cells[0]= cells[1]->neighbor(polesi[1]);
    polesi[0]= cells[0]->index(poles[0]);

    Facet middle_facet(cells[0], polesi[0]);
    triangulation_.clear_cell_labels(middle_facet.first);
    triangulation_.clear_cell_labels(triangulation_.opposite(middle_facet).first);

    if (ore.will_fail()) {
      typename Simulator::Time t= ore.failure_time();
      ore.pop_failure_time();
      typename Simulator::Event_key k= simulator()->new_event(t, 
							      Facet_flip_event(ore, middle_facet, tr_.wrapper_handle()));
      
      triangulation_.set_label(middle_facet, k);
    } else {
      triangulation_.set_label(middle_facet, simulator()->null_event());
    }

    for (int c=0; c<2; ++c) {
      handle_new_cell(cells[c]);
    }

    CGAL_postcondition(audit_structure());
    v_.post_edge_flip(middle_facet);
    return middle_facet;
  }











  Edge flip(const Facet &flip_facet) {
    v_.pre_facet_flip(flip_facet);
    Vertex_handle poles[2];
    Facet other_flip_facet= triangulation_.opposite(flip_facet);
    poles[0]= flip_facet.first->vertex(flip_facet.second);
    poles[1]= other_flip_facet.first->vertex(other_flip_facet.second);
    
    CGAL_LOG(Log::LOTS,"\n\nFlipping facet ");
    CGAL_LOG(Log::LOTS,flip_facet.first->vertex((flip_facet.second+1)%4)->point() << "--"
		     << flip_facet.first->vertex((flip_facet.second+2)%4)->point() << "--"
		     << flip_facet.first->vertex((flip_facet.second+3)%4)->point() << std::endl);
    //triangulation_.write_labeled_facet(flip_facet, log_lots() );
    CGAL_LOG(Log::LOTS," with poles " << poles[0]->point() << ", " << poles[1]->point());
    CGAL_LOG(Log::LOTS,std::endl);
    
    CGAL_assertion(triangulation_.tds().is_facet(flip_facet.first, flip_facet.second) || print());

    Event_key failed_key= triangulation_.label(flip_facet);
    Certificate ore= extract_root_stack(failed_key);
    if (ore.will_fail()) {
      CGAL_LOG(Log::LOTS, "Next root of certificate is " << ore.failure_time() << std::endl);
    } else {
      CGAL_LOG(Log::LOTS, "Certificate will never fail" << std::endl);
    }
    triangulation_.set_label(flip_facet, Event_key());

    //typename P::Event_key index= triangulation_.label(flip_facet);

    Facet inside_facet(flip_facet.first, (flip_facet.second+1)%4);
    Facet neighbor_facet= triangulation_.opposite(inside_facet);
    Cell_handle cells[2]= {flip_facet.first, other_flip_facet.first};
    CGAL_assertion(neighbor_facet.first != cells[0] && neighbor_facet.first != cells[1]);

    // go around and change the handler for each edge if it is degree 3. Don't have to look off cell
    for (int c=0; c<2; ++c) {
      clean_cell(cells[c]);
    }
   
   
    triangulation_.tds().flip_flippable(flip_facet);

    CGAL_assertion(triangulation_.tds().is_facet(neighbor_facet.first, neighbor_facet.second));

    Facet a_facet= triangulation_.opposite(neighbor_facet);
    Cell_handle a_cell= a_facet.first;

    Edge central_edge(a_cell, a_cell->index(poles[0]), a_cell->index(poles[1]));

    Cell_circulator cc= triangulation_.incident_cells(central_edge), ce=cc;
    do {
      triangulation_.clear_cell_labels(cc);
    } while (++cc != ce);


    if (ore.will_fail()) {
      typename Simulator::Time t= ore.failure_time();
      ore.pop_failure_time();
      typename Simulator::Event_key k= simulator()->new_event(t, Edge_flip_event(ore, central_edge, tr_.wrapper_handle()));
      triangulation_.set_label(central_edge, k);
    } else {
      triangulation_.set_label(central_edge, simulator()->null_event());
    }
    
    CGAL_expensive_assertion(labeling_is_valid());

    {
      Cell_circulator cc= triangulation_.incident_cells(central_edge), ce=cc;
      do {
	handle_new_cell(cc);
      } while (++cc != ce);
    }

    CGAL_postcondition(audit_structure());
    v_.post_facet_flip(central_edge);
    return central_edge;
  }











  void audit() const
  {
    CGAL_LOG(Log::SOME, "Verifying at time " << simulator()->audit_time() << ".\n");
    set_instantaneous_time();
    //CGAL_precondition(triangulation_.geom_traits().time()== simulator()->audit_time());
    //triangulation_.geom_traits().set_time(simulator().rational_current_time());
    audit_structure();
    triangulation_.is_valid(true);
  }









  const Triangulation& triangulation() const
  {
    return triangulation_;
  }

  Triangulation& triangulation() {
    return triangulation_;
  }

  const Moving_object_table* moving_object_table() const
  {
    return tr_.active_points_3_table_handle();
  }

  Simulator* simulator() {
    return tr_.simulator_handle();
  }

  const Simulator* simulator() const
  {
    return tr_.simulator_handle();
  }

  const Point& point(Point_key k) const
  {
    return moving_object_table()->at(k);
  }

  Vertex_handle vertex_handle(Point_key k) const
  {
    if (k.is_valid() && k.index() < vhs_.size()) {
      return vhs_[k.index()];
    }
    else {
      return NULL;
    }
  }

  void set_vertex_handle(Point_key k, Vertex_handle vh) {
    CGAL_precondition(k != Point_key());
    unsigned int bin= k.index();
    while (vhs_.size() <= bin) {
      vhs_.push_back(Vertex_handle(NULL));
    }
    /*if (vhs_.size() <=bin){
      vhs_.resize(bin+1);
      }*/
    vhs_[k.index()]=vh;
  }

  typename TraitsT::Side_of_oriented_sphere_3 power_test_object() const
  {
    return soc_;
  };
  typename TraitsT::Orientation_3 orientation_object() const
  {
    return o3_;
  }
  Certificate extract_root_stack(Event_key k) const
  {
    //typename Simulator::Event_handle<Event_base> eh(simulator()->event(k, Event_base()));
    //typename Simulator::Root_stack s= eh.pointer()->root_stack();
    return simulator()->template event<Event_base>(k/*, Event_base()*/).root_stack();
  }
  /*
    typename Simulator::Time extract_time(Event_key k) const {
    return simulator()->event<Facet_flip_event::Base>(k)->time();
    }*/



  void make_certificate( const Edge &e,
			 const typename Simulator::Time& st) {
    CGAL_LOG(Log::LOTS, "making certificate for edge ");
    CGAL_LOG_WRITE(Log::LOTS, triangulation_.write_edge(e, LOG_STREAM));
    CGAL_LOG(Log::LOTS, std::endl);
    CGAL_precondition(is_degree_3(e));
    CGAL_precondition_code(Facet_circulator fc= triangulation_.incident_facets(e));
    CGAL_precondition_code(Facet_circulator fe= fc);
    CGAL_precondition_code(do {
      );
			   CGAL_precondition(!has_event(*fc));
			   CGAL_precondition_code(++fc);
			   CGAL_precondition_code(
						  }while(fc != fe));

    CGAL_precondition(is_degree_3(e));
    CGAL_precondition(!has_degree_4_vertex(e));
    CGAL_precondition(!has_event(e));
    
    Certificate s= root_stack(e, st);
   
    if (s.will_fail()) {
      CGAL_LOG(Log::LOTS,"Failure time is " << s.failure_time() << std::endl);
      typename Simulator::Time t= s.failure_time();
      s.pop_failure_time();
      if (s.will_fail()) {
	CGAL_LOG(Log::LOTS, "Next root of this cert is " << s.failure_time() << std::endl);
      }
      typename Simulator::Event_key k=  simulator()->new_event(t, Edge_flip_event(s, e, tr_.wrapper_handle()));
      triangulation_.set_label(e, k);
    } else {
      CGAL_LOG(Log::LOTS,"Certificate will not fail "<< std::endl);
      triangulation_.set_label(e, simulator()->null_event());
    }
  }
  void make_certificate( const Edge &e) {
    make_certificate(e, simulation_traits_object().simulator_handle()->current_time());
  }





  void make_certificate( const Facet &e,
			 const typename Simulator::Time &st) {
    CGAL_precondition(!has_event(e));
    CGAL_LOG(Log::LOTS, "making certificate for facet ");
    CGAL_LOG_WRITE(Log::LOTS, triangulation_.write_facet(e, LOG_STREAM ));
    //triangulation_.write_facet(e, log_lots());
    CGAL_LOG(Log::LOTS,  std::endl);

    CGAL_precondition(!has_degree_3_edge(e));
    for (int i=0; i<3; ++i) {
      CGAL_precondition(triangulation_.label(triangulation_.edge(e, i)) == Event_key());
    }
    Certificate s= root_stack(e, st);
    if (s.will_fail()) {
      typename Simulator::Time t= s.failure_time();
      CGAL_LOG(Log::LOTS, "Failure time is " << t << std::endl);
      s.pop_failure_time();
      if (s.will_fail()) {
	CGAL_LOG(Log::LOTS, "Next root of this cert is " << s.failure_time() << std::endl);
      }
      typename Simulator::Event_key k= simulator()->new_event(t, Facet_flip_event(s, e, tr_.wrapper_handle()));
      triangulation_.set_label(e, k);
    } else {
      CGAL_LOG(Log::LOTS, "Certificate will not fail" << std::endl);
      triangulation_.set_label(e, simulator()->null_event());
    }
  }

   void make_certificate( const Facet &e) {
     make_certificate(e,
		      simulation_traits_object().simulator_handle()->current_time());
   }

 

  template <class Oit>
  void point_keys(const Facet &f, Oit out) const
  {
    int hinf=-1;
    for (unsigned int i=0; i<4; ++i) {
      Point_key k= f.first->vertex(i)->point();
      if (!k.is_valid()) {
	hinf=i;
	break;
      }
    }
    if (hinf==-1) {
      Point_key k= triangulation_.mirror_vertex(f.first, f.second)->point();
      if ( ! k.is_valid() ) {
	hinf=4;
      }
    }
    if (hinf != -1) {
      CGAL_LOG(Log::LOTS, "hinf is " << hinf << std::endl);
      if (hinf ==4) {
	for (unsigned int i=0; i<4; ++i) {
	  Point_key k= f.first->vertex(i)->point();
	  *out= k;
	  ++out;
	}
	return;
      }
      else {
	//Facet ff(f.first, hinf);
	if (hinf%2!=0) {
	  Point_key k= triangulation_.mirror_vertex(f.first, f.second)->point();
	  *out= k;
	  ++out;
	}
	for (int i=0; i<4; ++i) {
	  // CGAL infinite cells seem to be misoriented
	  if (i==hinf) continue;
	  *out= f.first->vertex(i)->point();
	  ++out;
	}
	if (hinf%2==0) {
	  Point_key k= triangulation_.mirror_vertex(f.first, f.second)->point();
	  *out= k;
	  ++out;
	}
      }
    }
    else {
      for (unsigned int i=0; i<4; ++i) {
	Point_key k= f.first->vertex(i)->point();
	*out= k;
	++out;
      }
      Point_key k= triangulation_.mirror_vertex(f.first, f.second)->point();
      *out= k;
      ++out;
    }
  }




private:
  Certificate root_stack(const Edge &e, const typename Simulator::Time &st) const
  {
    return root_stack(*triangulation_.incident_facets(e), st);
  }

  

  void make_no_events(const Edge &e) {
    CGAL_precondition(triangulation_.has_degree_3(e));
    Facet_circulator fc= triangulation_.incident_facets(e);
    Facet_circulator fe= fc;
    do {
      if (has_event(*fc)) {
	Event_key k= triangulation_.label(*fc);
	simulator()->delete_event(k);
	triangulation_.set_label(*fc, Event_key());
      }
    } while(++fc != fe);
  }
public:
  Point_key replace_vertex(Vertex_handle vh, Point_key k) {
    Point_key ok= vh->point();
    vh->point()=k;
    vhs_[ok.index()]= Vertex_handle();
    vhs_[k.index()]= vh;
    return ok;
  }
protected:
  bool has_event(const Edge &e) const
  {
    return triangulation_.label(e) != Event_key();
  }

  bool has_event(const Facet &e) const
  {
    
    return triangulation_.label(e) != Event_key();
  }


  void create_edge_flips(Vertex_handle v) {
    CGAL_precondition(!is_degree_4(v));
    std::vector<Cell_handle> ics;
    triangulation().incident_cells(v, std::back_inserter(ics));
    for (unsigned int i=0; i< ics.size(); ++i) {
      int j=-1;// disable warning
      CGAL_assertion_code(bool ret=)
        ics[i]->has_vertex(v, j); // initializes j
      CGAL_assertion(j != -1);
      CGAL_assertion(ret);
      for (int k=0; k<4 ; ++k) {
	if (k==j) continue;
	Edge e(ics[i], j, k);
	if (is_degree_3(e) && !has_event(e) && !has_degree_4_vertex(e)) {
	  // rather than make_certificate due to ordering dependencies
	  make_edge_flip(e);
	}
      }
    }
  }

  void suppress_edge_flips(Vertex_handle v) {
    CGAL_precondition(is_degree_4(v));
    std::vector<Cell_handle> ics;
    triangulation().incident_cells(v, std::back_inserter(ics));
    for (unsigned int i=0; i< ics.size(); ++i) {
      int j=-1; // keep some dumb compiler happy
      CGAL_assertion_code(bool ret=)
        ics[i]->has_vertex(v, j);
      CGAL_assertion(ret);
      for (int k=0; k<4 ; ++k) {
	if (k==j) continue;
	Edge e(ics[i], j, k);
	if (has_event(e)) {
	  simulator()->delete_event(triangulation_.label(e));
	  triangulation().set_label(e, Event_key());
	}
      }
    }
  }

  void make_edge_flip(Edge &edge) {
    CGAL_LOG(Log::LOTS, "Making edge flip ");
    CGAL_LOG_WRITE(Log::LOTS, triangulation_.write_labeled_edge(edge, LOG_STREAM ));
    CGAL_LOG(Log::LOTS,std::endl);
    CGAL_assertion(triangulation_.has_degree_3(edge));
    typename Simulator::Event_key k= typename Simulator::Event_key();
    Facet_circulator fc= triangulation_.incident_facets(edge), fe=fc;
    do {
      if (has_event(*fc)) {
	CGAL_assertion( k == Event_key() );
	k=triangulation_.label(*fc);
	triangulation_.set_label(*fc, typename Simulator::Event_key());
	Event_key kk=change_to_edge_flip(edge, k);
	triangulation_.set_label(edge, kk);
	simulator()->audit_event(kk);
	return;
      }
    } while (++fc != fe);


    CGAL_LOG(Log::LOTS, "Making up edge event.\n");
    make_certificate(edge);
  }

  void make_not_edge_flip(Edge &edge, Cell_handle h) {
    if (true) {
      CGAL_LOG(Log::LOTS, "Making edge ");
      CGAL_LOG_WRITE(Log::LOTS, triangulation_.write_labeled_edge(edge, LOG_STREAM ));
      CGAL_LOG(Log::LOTS, " not an edge flip.\n");
    }
    CGAL_assertion(is_degree_3(edge) || print());
    typename Simulator::Event_key index = triangulation_.label(edge);
    CGAL_precondition( index != Event_key() );

    triangulation_.set_label(edge, typename Simulator::Event_key());

    Cell_circulator fc=triangulation_.incident_cells(edge), pfc=fc, ef;
    ++fc;
    ef=fc;
    do {
      if( (h != fc) && ( h != pfc) ) {
	for (unsigned int i=0; i< 4; ++i) {
	  if (pfc->neighbor(i) == fc){
	    Facet f(pfc, i);
	    CGAL_precondition(!has_event(f));
	    triangulation_.set_label(f, change_to_facet_flip(f, index));
	    simulator()->audit_event(triangulation_.label(f));
	  }
	}
      }
      ++pfc; ++fc;
    } while (ef != fc);
  }

  template <class C>
  std::back_insert_iterator<C> back_inserter(C &c) const
  {
    return std::back_insert_iterator<C>(c);
  }

  bool contains(typename std::vector<Facet>::iterator beg,
		typename std::vector<Facet>::iterator end, Facet f) {
    Facet of = triangulation_.opposite(f);
    for (; beg != end; ++beg) {
      if (f.first == beg->first) {
	if (f.second == beg->second) return true;
      }
      else if (of.first == beg->first) {
	if (of.second == beg->second) return true;
      }
    }
    return false;
  }

  bool labeling_is_valid() const
  {
    for (All_facets_iterator eit= triangulation_.all_facets_begin(); eit != triangulation_.all_facets_end();
	 ++eit) {
      Facet f= *eit;
      Facet of= triangulation_.opposite(f);
      if (triangulation_.label(f) != triangulation_.label(of)) {
	triangulation_.write_labeled_facet(f, std::cerr);
	std::cerr << " does not match ("<< triangulation_.label(f) << ", "
		  << triangulation_.label(of) << ")\n";
      }
    }
    for (All_edges_iterator eit= triangulation_.all_edges_begin(); eit != triangulation_.all_edges_end();
	 ++eit) {
      Edge e=*eit;
      triangulation_.label(e);
    }
    return true;
  }



  bool audit_structure() const
  {
    if (triangulation_.dimension() != 3) return true;

    std::set<Point_key> pks;
    for (Finite_vertices_iterator eit = triangulation_.finite_vertices_begin();
	  eit != triangulation_.finite_vertices_end(); ++eit) {
      CGAL_assertion_code(Point_key k= eit->point());
      CGAL_assertion(pks.find(k) == pks.end());
      CGAL_assertion_code(pks.insert(k));
    }
    
    if (!has_certificates()) {
      for (All_edges_iterator eit = triangulation_.all_edges_begin();
	   eit != triangulation_.all_edges_end(); ++eit) {
	CGAL_assertion(!has_event(*eit));
      }
      for (All_facets_iterator eit = triangulation_.all_facets_begin();
	   eit != triangulation_.all_facets_end(); ++eit) {
	CGAL_assertion(!has_event(*eit));
      }
    } else {
      CGAL_LOG(Log::SOME, "Auditing structure" << std::endl);
      //print();
      
      for (All_edges_iterator eit = triangulation_.all_edges_begin();
	   eit != triangulation_.all_edges_end(); ++eit) {
	bool isd3= is_degree_3(*eit);
	bool hd4= has_degree_4_vertex(*eit);
	if (!isd3 || hd4) {
	  if (has_event(*eit)) {
	    std::cerr << "Edge should not have certificate ";
	    triangulation_.write_labeled_edge(*eit, std::cerr);
	    std::cerr << std::endl;
	    simulator()->audit_event(triangulation_.label(*eit));
	    CGAL_error();
	  }
	} else if ( isd3) {
	  if (!has_event(*eit)) {
	    std::cerr << "Edge should have certificate ";
	    triangulation_.write_labeled_edge(*eit, std::cerr);
	    std::cerr << std::endl;
	    CGAL_error();
	  } else {
	    simulator()->audit_event(triangulation_.label(*eit));
	  }
	}
      }

      for (All_facets_iterator eit = triangulation_.all_facets_begin();
	   eit != triangulation_.all_facets_end(); ++eit) {
	bool hsd3= has_degree_3_edge(*eit);
	bool hd4= has_degree_4_vertex(*eit);
	if (hsd3 || hd4) {
	  if (has_event(*eit)) {
	    std::cerr << "Facet should not have certificate ";
	    triangulation_.write_labeled_facet(*eit, std::cerr);
	    std::cerr << std::endl;
	    simulator()->audit_event(triangulation_.label(*eit));
	    CGAL_error();
	  }
	} else {
	  if (!has_event(*eit)) {
	    std::cerr << "Facet should have certificate ";
	    triangulation_.write_labeled_facet(*eit, std::cerr);
	    std::cerr << std::endl;
	    CGAL_error();
	  } else {
	    simulator()->audit_event(triangulation_.label(*eit));
	  }
	}
      }
    }
    return true;
  }

  Event_key change_to_edge_flip(const Edge &e, Event_key k) {
    if (k== simulator()->null_event()) return k;
    Certificate s= extract_root_stack(k);
    return simulator()->set_event(k, Edge_flip_event(s, e, tr_.wrapper_handle()));
  }

  Event_key change_to_facet_flip(const Facet &f, Event_key k) {
    if (k== simulator()->null_event()) return k;
    Certificate s= extract_root_stack(k);
    return simulator()->set_event(k, Facet_flip_event(s, f, tr_.wrapper_handle()));
  }

  /*Moving_object_table* moving_object_table() {
    return mpt_.pointer();
    }*/

public:
  void clean_cell(Cell_handle h) {
    CGAL_precondition(has_certificates_);
    for (unsigned int i=0; i< 4; ++i) {
      for (unsigned int j=0; j<i; ++j) {
	Edge e(h, i, j);
	if (has_event(e)) {
	  make_not_edge_flip(e, h);
	  //delete_event(triangulation_.label(e));
	  //triangulation_.set_label(e, Event_key());
	}
      }
    }
    for (unsigned int i=0; i< 4; ++i) {
      Facet f(h, i);
      if (has_event(f)) {
	simulator()->delete_event(triangulation_.label(f));
	triangulation_.set_label(f, Event_key());
      }
    }
    v_.destroy_cell(h);
  }

  void handle_new_cell(Cell_handle h) {
    CGAL_precondition(has_certificates_);
    for (unsigned int i=0; i<4; ++i) {
      for (unsigned int j=0; j< i; ++j) {
	Edge e(h, i, j);
	if (is_degree_3(e) && !has_event(e) && !has_degree_4_vertex(e)) {
	  make_edge_flip(e);
	}
      }
      Facet f(h, i);
      if (!has_event(f) && !has_degree_3_edge(f)) {
	make_certificate(f);
      }
     
    }
    for (unsigned int i=0; i<4; ++i) {
      Vertex_handle vh= h->vertex(i);
      if (is_degree_4(vh)) {
	suppress_edge_flips(vh);
      } else {
	create_edge_flips(vh);
      }
    }

    v_.create_cell(h);
  }
protected:
  void handle_changed_cell(Cell_handle) {
    
  }

  Certificate root_stack(const Facet &f,
			 const typename Simulator::Time &st) const
  {
    std::vector<Point_key> ids;
    point_keys(f, back_inserter(ids));
#ifndef NDEBUG
    std::vector<Point_key> mids;
    Facet of= triangulation_.opposite(f);
    point_keys(of, back_inserter(mids));
#endif
#ifndef NDEBUG
    CGAL_LOG(Log::LOTS, "Creating root_stack for points ");
    for (typename std::vector<Point_key>::const_iterator cit= ids.begin(); cit != ids.end(); ++cit) {
      CGAL_LOG(Log::LOTS, *cit);
    }
    CGAL_LOG(Log::LOTS, std::endl);
#endif
    bool is_const=true;
    for (unsigned int i=0; i<ids.size(); ++i) {
      if (!point(ids[i]).is_constant()) {
	is_const=false;
	break;
      }
    }
    if (is_const) {
      //typename Kinetic_kernel::Certificate_function cf(1.0);
      return Certificate(); //simulator()->root_stack_object(cf);
    }
    else if (ids.size()==4) {
      /*if (point(ids[0]).is_constant()){
      // hack
      std::swap(ids[0], ids[3]);
      std::swap(ids[1], ids[2]);
      }*/
      return o3_(point(ids[0]),
		 point(ids[1]),
		 point(ids[2]),
		 point(ids[3]),
		 st,
		 simulator()->end_time());
    }
    else {
      CGAL_assertion(ids.size()==5);
      /*if (point(ids[0]).is_constant()){
      // hack for linear case
      std::swap(ids[0], ids[4]);
      std::swap(ids[1], ids[2]);
      }*/
      return soc_(point(ids[0]),
		  point(ids[1]),
		  point(ids[2]),
		  point(ids[3]),
		  point(ids[4]),
		  st,
		  simulator()->end_time());
    }
    // Some compilers give warnings without this
    CGAL_postcondition(0);
    return Certificate();
    //return simulator()->root_stack_object(typename TraitsT::Simulator::Function_kernel::Function(0));
  }

  
  TraitsT tr_;
  Triangulation triangulation_;
  std::vector<Vertex_handle> vhs_;
  typename TraitsT::Side_of_oriented_sphere_3 soc_;
  typename TraitsT::Orientation_3 o3_;
  bool has_certificates_;
  Visitor v_;
};

template <class TraitsT, class Visitor>
inline std::ostream &operator<<(std::ostream &out,
				const Delaunay_triangulation_base_3<TraitsT, Visitor> &dt)
{
  dt.write(out);
  return out;
}


} } } //namespace CGAL::Kinetic::internal
#endif
