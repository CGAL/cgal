// Copyright (c) 2016  INRIA Nancy - Grand Est (France).
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
// $URL: $
// $Id:  $
//
//
// Author(s)     : Iordan Iordanov

#ifndef CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_2_H
#define CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_2_H

#include <vector>
#include <map>

#include <CGAL/Periodic_4_hyperbolic_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_triangulation_vertex_base_2.h>
#include <CGAL/Periodic_4_hyperbolic_triangulation_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <iterator>
#include <CGAL/intersections.h>
#include <CGAL/result_of.h>
#include <CGAL/iterator.h>
#include <boost/bind.hpp>

#include <CGAL/Timer.h>

#if defined PROFILING_MODE
	#include <CGAL/Timer.h>
	extern long calls_predicate_identity;
	extern long calls_predicate_non_identity;
	extern double time_predicate_identity;
	extern double time_predicate_non_identity;
	extern double time_remove_dp;
#endif	


namespace CGAL {

  template <  class GT,
  class TDS = Triangulation_data_structure_2<
  Periodic_4_hyperbolic_triangulation_vertex_base_2<GT>,
  Periodic_4_hyperbolic_triangulation_face_base_2<GT> 
  >
  >
  class Periodic_4_hyperbolic_Delaunay_triangulation_2: public Periodic_4_hyperbolic_triangulation_2<GT, TDS> {

	typedef Periodic_4_hyperbolic_Delaunay_triangulation_2<GT, TDS>   Self;
	typedef Periodic_4_hyperbolic_triangulation_2<GT, TDS>            Base;

  public:

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2  
	using Base::cw;
	using Base::ccw;
	using Base::geom_traits;
	using Base::tds;
#endif

	using Base::apply;

	typedef typename Base::Locate_type                         Locate_type;
	typedef typename Base::Geometric_traits                    Geometric_traits;
	typedef typename Base::Triangulation_data_structure        Triangulation_data_structure;
	//typedef typename Base::Int                                 Int;
	typedef typename Base::Hyperbolic_translation              Hyperbolic_translation;
	typedef typename Base::Circle_2                            Circle_2;
	typedef Circle_2                                           Circle;
	typedef typename Base::Point_2                             Point_2;
	typedef Point_2                                            Point;
	typedef typename Base::Segment_2                           Segment_2;
	typedef Segment_2                                          Segment;
	typedef typename Base::Triangle_2                          Triangle_2;
	typedef Triangle_2                                         Triangle;

	typedef typename Base::Periodic_point                      Periodic_point;
	typedef typename Base::Periodic_segment                    Periodic_segment;
	typedef typename Base::Periodic_triangle                   Periodic_triangle;  

	typedef typename Base::Vertex                              Vertex;
	typedef typename Base::Edge                                Edge;
	typedef typename Base::Face                                Face;

	typedef typename Base::Vertex_handle                       Vertex_handle;
	typedef typename Base::Face_handle                         Face_handle;

	typedef typename Base::size_type                           size_type;
	typedef typename Base::difference_type                     difference_type;

	typedef typename Base::Face_iterator                       Face_iterator;
	typedef typename Base::Edge_iterator                       Edge_iterator;
	typedef typename Base::Vertex_iterator                     Vertex_iterator;
	typedef typename Base::Face_circulator                     Face_circulator;
	typedef typename Base::Edge_circulator                     Edge_circulator;
	typedef typename Base::Vertex_circulator                   Vertex_circulator;
	typedef typename Base::Line_face_circulator                Line_face_circulator;
	// typedef typename GT::Construct_inexact_intersection_2 	   Construct_intersection_2;
	typedef typename GT::Construct_hyperbolic_circumcenter_2   Construct_hyperbolic_circumcenter_2;
	typedef typename GT::Construct_segment_2 				   Construct_segment_2;

  private:
	typedef typename GT::FT                                    FT;


	class Dummy_point {
	private:
	  Point _pt;
	  bool  _is_inserted;
	  Vertex_handle _vh;

	public:

	  Dummy_point(FT x, FT y): _pt(x,y), _is_inserted(true) {}
	  Dummy_point(Point p):    _pt(p),   _is_inserted(true) {}

	  Point operator()()      const {  return _pt;          }
	  bool  is_inserted()     const {  return _is_inserted; }
	  Vertex_handle vertex()  const {  return _vh;          }
	  void  set_inserted(bool val)      { _is_inserted = val; }
	  void  set(Point val)              { _pt = val;          }
	  void  set_vertex(Vertex_handle v) { _vh = v;            }
	};

	std::vector<Dummy_point> dummy_points;

  public:

	typedef Point                                            value_type;
	typedef const value_type&                                const_reference;
	typedef Tag_false                                        Weighted_tag;

  private:
	int n_dpt;

	std::vector<Vertex_handle> insert_dummy_points(bool rational = true);
	
	bool is_removable(Vertex_handle v, Delaunay_triangulation_2<GT, TDS>& dt, std::map<Vertex_handle, Vertex_handle>& vmap);

  public:
	Periodic_4_hyperbolic_Delaunay_triangulation_2(Geometric_traits gt) : 
	Periodic_4_hyperbolic_triangulation_2<GT, TDS>(gt) { 
		std::cout << "Inserting dummy points..."; std::cout.flush();
		insert_dummy_points();
		std::cout << "  DONE!" << std::endl;
		n_dpt = 14; 
	}  

	Periodic_4_hyperbolic_Delaunay_triangulation_2(
	  const Circle_2 domain = Circle_2(Point_2(FT(0),FT(0)), FT(1*1)), 
	  const Geometric_traits &gt = Geometric_traits() ) :
	Periodic_4_hyperbolic_triangulation_2<GT, TDS>(domain, gt) { 
		std::cout << "Inserting dummy points..."; std::cout.flush();
		insert_dummy_points();
		std::cout << "  DONE!" << std::endl;
		n_dpt = 14; 
	}

	Periodic_4_hyperbolic_Delaunay_triangulation_2(const Periodic_4_hyperbolic_Delaunay_triangulation_2& tr) :
	Periodic_4_hyperbolic_triangulation_2<GT, TDS>(tr) { 
		std::cout << "Inserting dummy points..."; std::cout.flush();
		insert_dummy_points();
		std::cout << "  DONE!" << std::endl;
		n_dpt = 14;
	}

	Vertex_handle insert(const Point  &p, Face_handle start = Face_handle(), bool batch_insertion = false, bool verified_input = false);

	template < class InputIterator >
	std::ptrdiff_t
	insert(InputIterator first, InputIterator last, bool verified_input = false) {
	  size_type n_initial = this->number_of_vertices();
	  std::vector<Point> points(first, last);
	  spatial_sort(points.begin(), points.end(), geom_traits());
	  Face_handle f;
	  int cnt = 0;
	  int cnt_good = 0;
	  for (typename std::vector<Point>::const_iterator p = points.begin(), end = points.end(); p != end; ++p){
		  cnt++;
		  Vertex_handle v = insert(*p, f, true, verified_input);
	  	  if (v != Vertex_handle()) {
	  	  	cnt_good++;
	  	  	f = v->face();
	  	  }
	  }

	  std::ptrdiff_t ret = this->number_of_vertices() - n_initial;

	  #if defined PROFILING_MODE
	  	CGAL::Timer tmr;
	  	tmr.start();
	  #endif

	  for (int i = 0; i < dummy_points.size(); i++) {
			if (dummy_points[i].is_inserted()) {
				typedef Delaunay_triangulation_2<GT, TDS>           Delaunay;
				Delaunay dt;
				std::map<Vertex_handle, Vertex_handle> vmap;

				if (is_removable(dummy_points[i].vertex(), dt, vmap)) {
					size_type tmp_n = this->number_of_vertices();
					remove(dummy_points[i].vertex());
					CGAL_assertion(this->number_of_vertices() == tmp_n - 1);
					dummy_points[i].set_inserted(false);
				}
			}
		}

		#if defined PROFILING_MODE
			tmr.stop();
			time_remove_dp = tmr.time();
		#endif

		n_dpt = 0;
		for (int i = 0; i < dummy_points.size(); i++) {
			if (dummy_points[i].is_inserted())
				n_dpt++;
		}

	  return ret;
	}

	Face_handle locate(const Point& p, Locate_type& lt, int& li, const Face_handle fh = Face_handle()) const {
	  Hyperbolic_translation lo;
	  return this->hyperbolic_locate(p, lt, li, lo, fh);
	}

	Face_handle locate(const Point& p, const Face_handle fh = Face_handle()) const {
	  Hyperbolic_translation lo;
	  Locate_type lt;
	  int li;
	  return this->hyperbolic_locate(p, lt, li, lo, fh);
	}

 	
 	template<class OutputFaceIterator>
 	void
	find_conflicts(	const Point& 		p, 
					OutputFaceIterator it) const {
		Hyperbolic_translation lo;
		Face_handle hint = this->euclidean_locate(p, lo);
        if (hint != Face_handle()) 
            find_conflicts(hint, p, lo, it);
	}

  	template<class OutputFaceIterator>
  	void 
  	find_conflicts( Face_handle         d, 
				  	const Point&        pt, 
				  	const Hyperbolic_translation&       current_off,
				  	OutputFaceIterator  it ) const;

	Face_handle periodic_euclidean_locate(const Point& p, Hyperbolic_translation& o, const Face_handle fh = Face_handle()) const {
		Locate_type lt;
		int li;
		return this->euclidean_locate(p, lt, li, o, fh);
	}

	Face_handle periodic_locate(const Point& p, Locate_type& lt, int& li, Hyperbolic_translation& lo, const Face_handle fh = Face_handle()) const {
	  return this->hyperbolic_locate(p, lt, li, lo, fh);
	}

	Face_handle periodic_locate(const Point& p, Hyperbolic_translation& lo, const Face_handle fh = Face_handle()) const {
	  return this->hyperbolic_locate(p, lo, fh);
	}

	Point_2 get_dummy_point(int i) const {
		CGAL_triangulation_precondition(0 <= i && i <= dummy_points.size());
	  	return dummy_points[i]();
	}

	void remove(Vertex_handle v);

	bool is_dummy_vertex(Vertex_handle vh) const {
		for (int i = 0; i < dummy_points.size(); i++) {
			if (dummy_points[i]() == vh->point())
				return true;
		}
		return false;
	}

	int number_of_dummy_points() const { return n_dpt; }

protected:

	Bounded_side _side_of_circle( const Face_handle &f, const Point &q,
								  const Hyperbolic_translation &translation ) const {
	  Point p[] = { f->vertex(0)->point(),
					f->vertex(1)->point(),
					f->vertex(2)->point()};
	  
	  Hyperbolic_translation o[]= { translation * f->translation(0),
									translation * f->translation(1),
									translation * f->translation(2)};
	  
	  #if defined PROFILING_MODE
		  if (o[0].is_identity() && o[1].is_identity() && o[2].is_identity()) {
		  	calls_predicate_identity++;
		  } else {
		  	calls_predicate_non_identity++;
		  }
	  #endif

	  	#if defined PROFILING_MODE
		  	CGAL::Timer tmr;
			tmr.start();
		#endif

	  Oriented_side os = this->side_of_oriented_circle( p[0], p[1], p[2], q,
														o[0], o[1], o[2], Hyperbolic_translation());
	  #if defined PROFILING_MODE
	  	tmr.stop();
	  	if (o[0].is_identity() && o[1].is_identity() && o[2].is_identity()) {
	  		time_predicate_identity += tmr.time();
	  	} else {
	  		time_predicate_non_identity += tmr.time();
	  	}
	  #endif

	  if (os == ON_NEGATIVE_SIDE) {
		return ON_UNBOUNDED_SIDE;
	  }
	  else if (os == ON_POSITIVE_SIDE) {
		return ON_BOUNDED_SIDE;
	  }
	  else {
		return ON_BOUNDARY;
	  }
   	}
   	

public:

	Point_2	dual (Face_handle f, Hyperbolic_translation nboff = Hyperbolic_translation()) const {
		Point_2 res = Construct_hyperbolic_circumcenter_2()(	f->vertex(0)->point(), 		f->vertex(1)->point(), 		f->vertex(2)->point(),
				 						 						nboff * f->translation(0), nboff * f->translation(1), nboff * f->translation(2));
		return res;
	}


	Segment_2 dual(const Edge &e) const {
		Segment_2 res = Construct_segment_2()(dual(e.first), dual(e.first->neighbor(e.second), e.first->neighbor_translation(e.second)));
		return res;
	}



	void clear() {
		Base::clear();
		insert_dummy_points();
	}

};  // class Periodic_4_hyperbolic_Delaunay_triangulation_2



template <class Gt, class Tds>
bool
Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt,Tds>::
is_removable(Vertex_handle v, Delaunay_triangulation_2<Gt,Tds>& dt, std::map<Vertex_handle, Vertex_handle>& vmap) {
  
  typedef typename Gt::FT                             FT;
  typedef typename Gt::Circle_2                       Circle;
  typedef Delaunay_triangulation_2<Gt, Tds>           Delaunay;
  typedef typename Delaunay::Finite_faces_iterator    Finite_Delaunay_faces_iterator;

  // This is the exact value of the limit.
  // The systole is 2*acosh(1+sqrt(2)), and we want half of that.
  // The max _diameter_ of the hyperbolic circles must be less than this.
  double factor(0.95); 			// percentage of the limit to consider -- failsafe!
  double lim( factor*acosh(1. + sqrt(2.)) );

  std::vector<Vertex_handle> bdry_verts;
  Face_circulator nbf(tds().incident_faces(v)), done(nbf);
  do {
	int idx = nbf->index(v);
	Hyperbolic_translation off = nbf->translation(idx).inverse();
	off = off*nbf->translation(ccw(idx));
	Vertex_handle thisv = nbf->vertex(ccw(idx));
	bdry_verts.push_back(thisv);
	Point pt = apply(thisv, off);
	Vertex_handle new_v = dt.insert(pt);
	vmap.insert(std::pair<Vertex_handle, Vertex_handle>(new_v, thisv));
  } while (++nbf != done);


  int n_verts = bdry_verts.size();
  double max_diam = 0.;
  for (Finite_Delaunay_faces_iterator fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); fit++) {

	bool is_good = true;
	for (int i = 0; i < 3; i++) {
	  Vertex_handle this_v = vmap[fit->vertex(i)];
	  Vertex_handle prev_v = bdry_verts[n_verts - 1];
	  Vertex_handle curr_v = bdry_verts[0];
	  for (int j = 1; curr_v != this_v; j = (j+1)%n_verts) {
		prev_v = curr_v;
		curr_v = bdry_verts[j];
	  }
	  if (vmap[fit->vertex(ccw(i))] == prev_v) {
		is_good = false;
		break;
	  }
	}
	if (is_good) {
	  Circle c(fit->vertex(0)->point(), 
			   fit->vertex(1)->point(), 
			   fit->vertex(2)->point());
	  typename Gt::Compute_hyperbolic_diameter cdiam;
	  double diam = cdiam(c);
	  if (max_diam < diam) {
		max_diam = diam;
	  }
	}
  }

  if (max_diam < lim) {
	return true;
  } else {
	return false;
  }
}


/*********************************************************************************/


template <class Gt, class Tds>
template <class OutputFaceIterator>
void
Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt,Tds>::
find_conflicts( Face_handle         d, 
				const Point&        pt, 
				const Hyperbolic_translation&       current_off,
				OutputFaceIterator  it ) const {
  if (d->tds_data().is_clear()) {
	if (_side_of_circle(d, pt, current_off) == ON_BOUNDED_SIDE) {
	  d->tds_data().mark_in_conflict();
	  d->store_translations(current_off);
	  it++ = d;
	  for (int jj = 0; jj < 3; jj++) {
		if (d->neighbor(jj)->tds_data().is_clear()) {
		  find_conflicts(d->neighbor(jj), pt, current_off * d->neighbor_translation(jj), it);
		}
	  }
	}
  }
}




template < class Gt, class Tds >
inline
typename Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::
insert(const Point  &p,  Face_handle hint, bool batch_insertion, bool verified_input) {

	Vertex_handle v;

	typedef typename Gt::Side_of_original_octagon Side_of_original_octagon;

	//CGAL::Timer tmr;
	//CGAL::Timer ttmr;
	//ttmr.start();

	CGAL::Bounded_side side = CGAL::ON_BOUNDED_SIDE;

	if (!verified_input) {
		//tmr.reset();
		//tmr.start();
		Side_of_original_octagon check = Side_of_original_octagon();
		side = check(p);
		//tmr.stop();
		//std::cout << "  side check: " << tmr.time() << std::endl;
	}

	if (side != CGAL::ON_UNBOUNDED_SIDE) {

		Hyperbolic_translation loff;
		Locate_type lt;
		int li;
		//tmr.reset();
		//tmr.start();
		Face_handle start = this->euclidean_locate(p, lt, li, loff, hint, true);
		//tmr.stop();
		//std::cout << "  locate: " << tmr.time() << std::endl;
		if (lt == Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::VERTEX) {
			return Vertex_handle();
		}

		std::vector<Face_handle> faces;
		//tmr.reset();
		//tmr.start();
		this->find_conflicts(start, p, loff, std::back_inserter(faces));
		//tmr.stop();
		//std::cout << "  conflicts: " << tmr.time() << std::endl;

		for (int i = 0; i < faces.size(); i++) {
		}

		//tmr.reset();
		//tmr.start();
		v = this->insert_in_hole(p, faces.begin(), faces.end());
		//tmr.stop();
		//std::cout << "  in hole: " << tmr.time() << std::endl;

		//tmr.reset();
		//tmr.start();
		Face_circulator ifc = tds().incident_faces(v), done(ifc);
		do {
			ifc->restore_translations(loff);
			ifc->tds_data().clear();
			ifc->make_canonical();
		} while (++ifc != done);

		Vertex_circulator ivc = tds().incident_vertices(v), done_v(ivc);
		do {
			ivc->remove_translation();
		} while (++ivc != done_v);
		//tmr.stop();
		//std::cout << "  cleanup: " << tmr.time() << std::endl;
		
		if (!batch_insertion) {
			CGAL_triangulation_assertion(this->is_valid(true));

			for (int i = 0; i < dummy_points.size(); i++) {
				if (dummy_points[i].is_inserted()) {
					typedef Delaunay_triangulation_2<Gt, Tds>           Delaunay;
					Delaunay dt;
					std::map<Vertex_handle, Vertex_handle> vmap;

					if (is_removable(dummy_points[i].vertex(), dt, vmap)) {
						remove(dummy_points[i].vertex());
						dummy_points[i].set_inserted(false);
                        n_dpt--;
					}
				}
			}
		}

		//ttmr.stop();
		//std::cout << "  --> total: " << ttmr.time();
		return v;
	}

	return Vertex_handle();
}


//------------------------------------------------------

template < class Gt, class Tds >
void
Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::
remove(Vertex_handle v) {
  typedef Delaunay_triangulation_2<Gt, Tds>           Delaunay;
  Delaunay dt;
  std::map<Vertex_handle, Vertex_handle> vmap;

  if (is_removable(v, dt, vmap)) {
  
	typedef typename Delaunay::Finite_faces_iterator    Finite_Delaunay_faces_iterator;
	typedef std::pair<Face_handle, int>                 Neighbor_pair;
	typedef std::pair<Edge, Neighbor_pair>              Edge_neighbor;


	std::vector<Edge>             bdry_edges;
	std::vector<Vertex_handle>    bdry_verts;
	std::map<Edge, Neighbor_pair> bdry_nbrs;

	Face_circulator nb = tds().incident_faces(v), done(nb);
	std::vector<Face_handle> nbrs;
	do {
	  int idx = nb->index(v);
	  Edge e = Edge(nb, idx);
	  bdry_edges.push_back(e);
	  Face_handle nbf = nb->neighbor(idx);
	  int nidx = 0;
	  if (nbf->neighbor(1) == nb) nidx = 1;
	  if (nbf->neighbor(2) == nb) nidx = 2;
	  CGAL_triangulation_assertion(nbf->neighbor(nidx) == nb);
	  bdry_nbrs.insert(Edge_neighbor(e, Neighbor_pair(nbf, nidx)));
	  bdry_verts.push_back(nb->vertex(ccw(idx)));

	  nb->store_translations(nb->translation(idx).inverse());
	  nbrs.push_back(nb);
	  nb++;
	} while(nb != done);

	for (int i = 0; i < bdry_edges.size(); i++) {
	  Edge e = bdry_edges[i];
	  Face_handle f = e.first;
	  int j = e.second;
	}

	int n_verts = bdry_verts.size();
	std::vector<Face_handle> new_f;
	for (Finite_Delaunay_faces_iterator fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); fit++) {
	  bool is_good = true;
	  for (int i = 0; i < 3; i++) {
		Vertex_handle this_v = vmap[fit->vertex(i)];
		Vertex_handle prev_v = bdry_verts[n_verts - 1];
		Vertex_handle curr_v = bdry_verts[0];
		for (int j = 1; curr_v != this_v; j = (j+1)%n_verts) {
		  prev_v = curr_v;
		  curr_v = bdry_verts[j];
		}
		if (vmap[fit->vertex(ccw(i))] == prev_v) {
		  is_good = false;
		  break;
		}
	  }

	  if (is_good) {
		Face_handle f = tds().create_face();
		for (int j = 0; j < 3; j++) {
		  f->set_vertex(j, vmap[fit->vertex(j)]);
		}
		new_f.push_back(f);
	  }
	}


	int internb = 0;
	int bdrynb = 0;
	for (int i = 0; i < new_f.size(); i++) {
	  for (int k = 0; k < 3; k++) {
		bool found_bdry = false;
		for (int j = 0; j < bdry_verts.size(); j++) { 
		  if (new_f[i]->vertex(ccw(k)) == bdry_verts[j] && 
			  new_f[i]->vertex(cw(k)) == bdry_verts[(j+1)%n_verts]) {
			found_bdry = true;
			Neighbor_pair nb = bdry_nbrs[bdry_edges[j]];
			Face_handle nbf = nb.first;
			int nbidx = nb.second;
			tds().set_adjacency(nbf, nbidx, new_f[i], k);
			bdrynb++;
			break;
		  } 
		}
		if (!found_bdry) {
		  for (int l = 0; l < new_f.size(); l++) {
			if (l == i) continue;
			for (int j = 0; j < 3; j++) {
			  if (new_f[i]->vertex(ccw(k)) == new_f[l]->vertex(cw(j)) &&
				  new_f[i]->vertex(cw(k))  == new_f[l]->vertex(ccw(j)) ) {
				tds().set_adjacency(new_f[i], k, new_f[l], j);
				internb++;
				break;
			  }
			}
		  }
		}
	  }
	}


	for (int j = 0; j < new_f.size(); j++) {
	  for (int i = 0; i < 3; i++) {
		new_f[j]->vertex(i)->set_face(new_f[j]);
	  }
	  new_f[j]->restore_translations();
	  new_f[j]->make_canonical();
	}

	for (int j = 0; j < bdry_edges.size(); j++) {
	  Face_handle f = bdry_edges[j].first;
	  int i = bdry_edges[j].second;
	  f->vertex(ccw(i))->remove_translation();
	}

	for (int i = 0; i < nbrs.size(); i++) {
	  tds().delete_face(nbrs[i]);
	}
	tds().delete_vertex(v);

	CGAL_triangulation_assertion(this->is_valid(true));

  } else { // is not removable
	std::cout << "   -> vertex cannot be removed!" << std::endl;
  }

}


} // namespace CGAL


#include <CGAL/Periodic_4_hyperbolic_triangulation_dummy_14.h>

#endif // CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_2_H
