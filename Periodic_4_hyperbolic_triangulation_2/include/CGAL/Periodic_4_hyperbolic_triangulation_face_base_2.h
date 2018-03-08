// Copyright (c) 1999-2016   INRIA Nancy - Grand Est (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Iordan Iordanov  <Iordan.Iordanov@loria.fr>
//                 


#ifndef CGAL_PERIODIC_4_HYPERBOLIC_TIANGULATION_FACE_BASE_2
#define CGAL_PERIODIC_4_HYPERBOLIC_TIANGULATION_FACE_BASE_2

#include <CGAL/basic.h>
#include <CGAL/Dummy_tds_2.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_data_structure_2.h>

namespace CGAL {

class Face_data {
private:
	bool in_conflict;
	bool on_boundary;
	bool cleared;

public:
	Face_data(): in_conflict(false), on_boundary(false), cleared(true) {}

	void mark_in_conflict() { in_conflict = true; cleared = false;  						}
	void mark_on_boundary() { on_boundary = true; cleared = false;  						}
	void clear() 			{ on_boundary = false; in_conflict = false; cleared = true;  	}
	bool is_in_conflict()	{ return in_conflict; }
	bool is_on_boundary() 	{ return on_boundary; }
	bool is_clear()			{ return cleared;     }
};


template< typename GT, typename TDS = Triangulation_data_structure_2<> >
class Periodic_4_hyperbolic_triangulation_face_base_2 {
public:
	typedef TDS 							Triangulation_data_structure;
	typedef typename TDS::Vertex_handle		Vertex_handle;
	typedef typename TDS::Face_handle		Face_handle;
	typedef typename TDS::Vertex 			Vertex;
	typedef typename TDS::Edge 				Edge;
	typedef typename TDS::Face 				Face;
	typedef Face_data 						TDS_data;
	typedef typename GT::Hyperbolic_translation				Hyperbolic_translation;

	template< typename TDS2 >
	struct Rebind_TDS {
		typedef Periodic_4_hyperbolic_triangulation_face_base_2<GT, TDS2> Other;
	};

private:
	Face_handle 	N[3];
	Vertex_handle	V[3];
	Hyperbolic_translation 			o[3];	// Hyperbolic_translations for vertices
	TDS_data 		_tds_data;
	int 			face_number;

public:

	Periodic_4_hyperbolic_triangulation_face_base_2() 
#ifndef CGAL_CFG_NO_CPP0X_UNIFIED_INITIALIZATION_SYNTAX
	: o{Hyperbolic_translation(), Hyperbolic_translation(), Hyperbolic_translation()}, face_number(-1)
	{}
#else
	{
		set_translations();
		face_number = -1;
	}
#endif

	Periodic_4_hyperbolic_triangulation_face_base_2(
		const Vertex_handle& v0, const Vertex_handle& v1,
		const Vertex_handle& v2) 
#ifndef CGAL_CFG_NO_CPP0X_UNIFIED_INITIALIZATION_SYNTAX
	: V{v0, v1, v2},
	  o{Hyperbolic_translation(), Hyperbolic_translation(), Hyperbolic_translation()}, face_number(-1)
	{
		set_neighbors();
	}
#else
	{
		set_translations();
		set_vertices(v0, v1, v2);
		set_neighbors();
		face_number = -1;
	}
#endif



	Periodic_4_hyperbolic_triangulation_face_base_2(
		const Vertex_handle& v0, const Vertex_handle& v1,
		const Vertex_handle& v2, const Face_handle&   n0,
		const Face_handle&   n1, const Face_handle&   n2) 
#ifndef CGAL_CFG_NO_CPP0X_UNIFIED_INITIALIZATION_SYNTAX
	: V{v0, v1, v2},
	  N{n0, n1, n2},
	  o{Hyperbolic_translation(), Hyperbolic_translation(), Hyperbolic_translation()}, face_number(-1)
	{
		set_neighbors();
	}
#else
	{
		set_translations();
		set_vertices(v0, v1, v2);
		set_neighbors(n0, n1, n2);
		face_number = -1;
	}
#endif



	// ACCESS

	const Vertex_handle& vertex(int i) const {
    	CGAL_triangulation_precondition( i >= 0 && i <= 2 );
    	return V[i];
	}


	const Hyperbolic_translation opposite_translation(const Face_handle& fh) {
		return o[opposite_index(fh)];
	}

	const Vertex_handle& opposite_vertex(const Face_handle& fh) {
		return V[opposite_index(fh)];
	}

	int opposite_index(const Face_handle& fh) {
		if (N[0] == fh) { return 0; }
		if (N[1] == fh) { return 1; }
		CGAL_triangulation_assertion( N[2] == fh );
		return 2;
	}


	bool has_vertex(const Vertex_handle& v) const {
		return (V[0] == v) || (V[1] == v) || (V[2]== v);
	}


	int index(const Vertex_handle& v) const {
        //cout << "Index for face " << face_number << endl;
    	//cout << "  Query is vertex " << v->idx() << endl;
    	if (v == V[0]) { return 0; }
    	if (v == V[1]) { return 1; }
    	CGAL_triangulation_assertion( v == V[2] );
    	return 2;
	}


	TDS_data& tds_data() {
		return _tds_data;
	}

	const Face_handle& neighbor(int i) const {
    	CGAL_triangulation_precondition( i >= 0 && i <= 2);
    	return N[i];
	}


	bool has_neighbor(const Face_handle& n) const {
    	return (N[0] == n) || (N[1] == n) || (N[2] == n);
	}


	bool has_neighbor(const Face_handle& n, int & i) const {
    	if(n == N[0]){ i = 0; return true; }
    	if(n == N[1]){ i = 1; return true; }
   	 	if(n == N[2]){ i = 2; return true; }
    	return false;
	}

	int index(const Face_handle& n) const {
    	if (n == N[0]) return 0;
    	if (n == N[1]) return 1;
    	CGAL_triangulation_assertion( n == N[2] );
    	return 2;
	}

	Hyperbolic_translation translation(int i) const {
		CGAL_triangulation_precondition( i >= 0 && i <= 2 );
		return o[i];
	}


	Hyperbolic_translation neighbor_translation(int i) const {
		CGAL_triangulation_precondition( i >= 0 && i <= 2 );
		int myi = Triangulation_cw_ccw_2::ccw(i);
		Hyperbolic_translation myof = o[myi];

		Hyperbolic_translation nbof;
		bool did_it = false;
		for (int c = 0; c < 3; c++) {
			if (N[i]->vertex(c) == V[myi]) {
				nbof = N[i]->translation(c);
				did_it = true;
				break;
			}
		}

		return (myof - nbof);
	}


	bool has_id_translations() {
		bool b = true;
		for (int i = 0; i < 3 && b; i++) {
			b = (b && o[i].is_identity());
		}
		return b;
	}


	// SETTING

	void set_translations() {
		o[0] = Hyperbolic_translation();
		o[1] = Hyperbolic_translation();
		o[2] = Hyperbolic_translation();
	}

	void set_translations(
		const Hyperbolic_translation& o0, const Hyperbolic_translation& o1, 
		const Hyperbolic_translation& o2) 
	{
		o[0] = o0;
		o[1] = o1;
		o[2] = o2;
	}

	void set_translation(int k, Hyperbolic_translation new_o) {
		CGAL_triangulation_precondition( k >= 0 && k <= 2 );
		o[k] = new_o;
	}

	void set_vertices() {
		V[0] = V[1] = V[2] = Vertex_handle();
	}

	void set_vertex(int k, const Vertex_handle& vh) {
		CGAL_triangulation_precondition( k >= 0 && k <= 2 );
		V[k] = vh;
	}

	void set_vertices(
		const Vertex_handle& v0, const Vertex_handle v1, 
		const Vertex_handle& v2)
	{
		V[0] = v0;
		V[1] = v1;
		V[2] = v2;
	}

	void set_neighbors() {
		N[0] = N[1] = N[2] = Face_handle();
	}

	void set_neighbors(
		const Face_handle& n0, const Face_handle& n1,
		const Face_handle& n2)
	{
		N[0] = n0;
		N[1] = n1;
		N[2] = n2;
	}

	void set_neighbor(int k, const Face_handle& nfh) {
		CGAL_triangulation_precondition( k >= 0 && k <= 2 );
		N[k] = nfh;
	}

	  // CHECKING

  	void make_canonical() {
  		
  		// If all translation are the same, there is an image of the face
  		// inside the original octagon. We store that one. This covers the
  		// simplest case.
  		if (o[0] == o[1] && o[1] == o[2]) {
  			o[0] = Hyperbolic_translation();
  			o[1] = Hyperbolic_translation();
  			o[2] = Hyperbolic_translation();
  			return;
  		}

  		// This covers the cases in which two vertices lie in the same domain.
  		for (int i = 0; i < 3; i++) {
  			int j = (i + 1) % 3;

  			if (o[i] == o[j]) {
  				int k = (i + 2) % 3;
  				if ((o[i].inverse()*o[k]) < (o[k].inverse()*o[i])) {
  					o[k] = o[i].inverse() * o[k];
  					o[i] = Hyperbolic_translation();
  					o[j] = Hyperbolic_translation();
  					return;
  				} else {
  					o[i] = o[k].inverse() * o[i];
  					o[j] = o[k].inverse() * o[j];
  					o[k] = Hyperbolic_translation();
  					return;
  				}
  			} else {
  				continue;
  			}
  		}

  		// Now we know that all vertices lie in different regions.
  		Hyperbolic_translation min(7,2,5);
  		Hyperbolic_translation trans;
  		for (int i = 0; i < 3; i++) {
  			int j = ( i + 1) % 3;  // the index of the 'next' vertex
  			Hyperbolic_translation tmp = o[i].inverse() * o[j];
  			if (tmp < min) {
  				min = tmp;
  				trans = o[i].inverse();
  			}
  		}

  		if (!trans.is_identity()) {
  			o[0] = trans * o[0];
  			o[1] = trans * o[1];
  			o[2] = trans * o[2];
  		}
  	}


  	int dimension() {
  		int d = 2;
  		for (int i = 0; i < 3; i++)
  			if (V[0] == Vertex_handle())
  				d--;

  		if (d < 0)
  			d = 0;

  		return d;
  	} 


  	void apply(Hyperbolic_translation io) {
  		o[0] = io * o[0];
  		o[1] = io * o[1];
  		o[2] = io * o[2];
  	}


  	void store_translations(Hyperbolic_translation noff = Hyperbolic_translation()) {
  		if (noff == Hyperbolic_translation()) {
  			V[0]->store_translation(o[0]);
  			V[1]->store_translation(o[1]);
  			V[2]->store_translation(o[2]);
  		} else {
  			V[0]->store_translation(noff * o[0]);
  			V[1]->store_translation(noff * o[1]);
  			V[2]->store_translation(noff * o[2]);
  		}
  	}

  	void restore_translations(Hyperbolic_translation loff = Hyperbolic_translation()) {
  		for (int i = 0; i < 3; i++) {
  			if (loff.is_identity()) {
  				o[i] = V[i]->get_translation();
  			} else {
  				o[i] = loff.inverse()*V[i]->get_translation();
  			}
  		}  		
  	}

  	void reorient() {
  		// N(eighbors), V(ertices), o(ffsets), no (neighbor translations)

  		int idx0 = 0, idx1 = 1; // the indices to swap

  		swap( N[idx0],  N[idx1]);
  		swap( V[idx0],  V[idx1]);
  		swap( o[idx0],  o[idx1]);

  	}


  	template <class T>
  	void swap(T& a, T& b) {
  		T tmp;
  		tmp = a;
  		a = b;
  		b = tmp;
  	}


  	// For use by Compact_container.
  	void * for_compact_container() const { return N[0].for_compact_container(); }
  	void * & for_compact_container()     { return N[0].for_compact_container(); }

};


template < class TDS >
inline
std::istream&
operator>>(std::istream &is, Periodic_4_hyperbolic_triangulation_face_base_2<TDS> &)
  // non combinatorial information. Default = nothing
{
  return is;
}

template < class TDS >
inline
std::ostream&
operator<<(std::ostream &os,
    const Periodic_4_hyperbolic_triangulation_face_base_2<TDS> &)
  // non combinatorial information. Default = nothing
{
  return os;
}

// Specialization for void.
template<typename GT>
class Periodic_4_hyperbolic_triangulation_face_base_2<GT, void>
{
public:
  typedef Dummy_tds_2                  					Triangulation_data_structure;
  typedef Triangulation_data_structure::Vertex_handle   Vertex_handle;
  typedef Triangulation_data_structure::Face_handle     Face_handle;
  template <typename TDS2>
  struct Rebind_TDS {
    typedef Periodic_4_hyperbolic_triangulation_face_base_2<GT, TDS2> Other;
  };
};


}  // namespace CGAL


#endif  // CGAL_PERIODIC_4_HYPERBOLIC_TIANGULATION_FACE_BASE_2