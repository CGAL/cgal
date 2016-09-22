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
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Iordan Iordanov  <Iordan.Iordanov@loria.fr>

#ifndef CGAL_PERIODIC_4_HYPERBOLIC_TIANGULATION_DS_FACE_BASE_2
#define CGAL_PERIODIC_4_HYPERBOLIC_TIANGULATION_DS_FACE_BASE_2

#include <CGAL/basic.h>
#include <CGAL/Dummy_tds_2.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Hyperbolic_octagon_word_4.h>
//#include <CGAL/Hyperbolic_octagon_word_8.h>
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
class Periodic_4_hyperbolic_triangulation_ds_face_base_2 {
public:
	typedef TDS 							Triangulation_data_structure;
	typedef typename TDS::Vertex_handle		Vertex_handle;
	typedef typename TDS::Face_handle		Face_handle;
	typedef typename TDS::Vertex 			Vertex;
	typedef typename TDS::Edge 				Edge;
	typedef typename TDS::Face 				Face;
	typedef Face_data 						TDS_data;
	typedef unsigned short int 				Int;
	typedef Hyperbolic_octagon_word_4<Int, GT>		Offset;

	template< typename TDS2 >
	struct Rebind_TDS {
		typedef Periodic_4_hyperbolic_triangulation_ds_face_base_2<GT, TDS2> Other;
	};

private:
	Face_handle 	N[3];
	Vertex_handle	V[3];
	Offset 			o[3];	// Offsets for vertices
	TDS_data 		_tds_data;
	int 			face_number;

public:

	Periodic_4_hyperbolic_triangulation_ds_face_base_2() 
#ifndef CGAL_CFG_NO_CPP0X_UNIFIED_INITIALIZATION_SYNTAX
	: o{Offset(), Offset(), Offset()}, face_number(-1)
	{}
#else
	{
		set_offsets();
		face_number = -1;
	}
#endif

	Periodic_4_hyperbolic_triangulation_ds_face_base_2(
		const Vertex_handle& v0, const Vertex_handle& v1,
		const Vertex_handle& v2) 
#ifndef CGAL_CFG_NO_CPP0X_UNIFIED_INITIALIZATION_SYNTAX
	: V{v0, v1, v2},
	  o{Offset(), Offset(), Offset()}, face_number(-1)
	{
		set_neighbors();
	}
#else
	{
		set_offsets();
		set_vertices(v0, v1, v2);
		set_neighbors();
		face_number = -1;
	}
#endif



	Periodic_4_hyperbolic_triangulation_ds_face_base_2(
		const Vertex_handle& v0, const Vertex_handle& v1,
		const Vertex_handle& v2, const Face_handle&   n0,
		const Face_handle&   n1, const Face_handle&   n2) 
#ifndef CGAL_CFG_NO_CPP0X_UNIFIED_INITIALIZATION_SYNTAX
	: V{v0, v1, v2},
	  N{n0, n1, n2},
	  o{Offset(), Offset(), Offset()}, face_number(-1)
	{
		set_neighbors();
	}
#else
	{
		set_offsets();
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


	const Offset opposite_offset(const Face_handle& fh) {
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

	Offset offset(int i) const {
		CGAL_triangulation_precondition( i >= 0 && i <= 2 );
		return o[i];
	}


	Offset neighbor_offset(int i) const {
		int myi = Triangulation_cw_ccw_2::ccw(i);
		Offset myof = o[myi];

		//if (myof.length() > o[Triangulation_cw_ccw_2::cw(i)].length()) {
		//	myi = Triangulation_cw_ccw_2::cw(i);
		//	myof = o[myi];
		//}

		Offset nbof;
		bool did_it = false;
		for (int c = 0; c < 3; c++) {
            //std::cout << "    -- neighbor_offset: N[i]->vertex(c) = " << N[i]->vertex(c)->idx() << ", V[myi] = " << V[myi]->idx() << std::endl;
			if (N[i]->vertex(c) == V[myi]) {
				nbof = N[i]->offset(c);
				did_it = true;
				break;
			}
		}

		if (!did_it) {
			cout << "Do not believe the neighbor offset! It's a lie!!1!" << endl;
		}

		return (myof - nbof);
	}


	bool has_zero_offsets() {
		bool b = true;
		for (int i = 0; i < 3 && b; i++) {
			b = (b && o[i].is_identity());
		}
		return b;
	}


	// SETTING

	void set_number(int n) {
		face_number = n;
	}

	int get_number() {
		return face_number;
	}

	void set_offsets() {
		o[0] = Offset();
		o[1] = Offset();
		o[2] = Offset();
	}

	void set_offsets(
		const Offset& o0, const Offset& o1, 
		const Offset& o2) 
	{
		o[0] = o0;
		o[1] = o1;
		o[2] = o2;
	}

	void set_offset(int k, Offset new_o) {
		o[k] = new_o;
	}

	void set_vertices() {
		V[0] = V[1] = V[2] = Vertex_handle();
	}

	void set_vertex(int k, const Vertex_handle& vh) {
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
		N[k] = nfh;
	}

	  // CHECKING

  	// the following trivial is_valid allows
  	// the user of derived cell base classes
  	// to add their own purpose checking

  	bool is_valid(bool, int) const { 
  		return true; 
  	}

/*
  	bool is_canonical() {

  		for (int i = 0; i < 3; i++) {
  			if (o[i].is_identity()) {
  				Offset lo = o[Triangulation_cw_ccw_2::ccw(i)];
  				//Offset ro = o[Triangulation_cw_ccw_2::cw(i)];
  				//if ( (lo.is_identity() || lo(0) == 1 || lo(0) == 3 || lo(0) == 4 || lo(0) == 6) && (ro.is_identity() || ro(0) == 1 || ro(0) == 3 || ro(0) == 4 || ro(0) == 6 )) {
  					//if (!(ro < lo)) {
  					if ( lo.is_identity() || lo(0) > 3 ) {
  						Offset ro = o[Triangulation_cw_ccw_2::cw(i)];
  						if (!(ro < lo)) {
  							return true;
  						}
  					}
  					//}
  				//}
  			}
  		}

  		return false;
  	}
*/

  	void make_canonical() {
  		//std::cout << "face [" << get_number() << "] original offsets: " << o[0] << ", " << o[1] << ", " << o[2] << endl;

  		int dst = 50;
  		Offset inv;

  		if (o[0] == o[1] && o[1] == o[2]) {
  			o[0] = Offset();
  			o[1] = Offset();
  			o[2] = Offset();
  			//std::cout << "face [" << get_number() << "] offsets: " << o[0] << ", " << o[1] << ", " << o[2] << endl;
  			//std::cout << "----------" << endl << endl;
  			return;
  		}

  		int rd;
  		for (int i = 0; i < 3; i++) {
  			int j = (i+1)%3;
  			Offset tmp;
  			
  			if (o[i].is_identity() && !o[j].is_identity()) {
  				rd = offset_reference_distance(o[j]);
  				tmp = Offset();
  				//cout << "   case 1, rd = " << rd << endl;
  			} else if (!o[i].is_identity()) {
  				tmp = o[i].inverse();
  				Offset cmb = tmp.append(o[j]);
  				if (cmb.is_identity()) {
  					rd = dst;
  				} else {
  					rd = offset_reference_distance(cmb);
  				}
  				//cout << "   case 2, rd = " << rd << endl;
  			}

  			if (rd < dst) {
  				dst = rd;
  				inv = tmp;
  			}
  		}

  		//cout << " chosen: " << inv << " with distance " << dst << endl;
  		o[0] = inv.append(o[0]);
  		o[1] = inv.append(o[1]);
  		o[2] = inv.append(o[2]);

  		//std::cout << "face [" << get_number() << "] offsets: " << o[0] << ", " << o[1] << ", " << o[2] << endl;
  		//std::cout << "----------" << endl << endl;
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


  	void apply(Offset io) {
  		o[0] = io.append(o[0]);
  		o[1] = io.append(o[1]);
  		o[2] = io.append(o[2]);
  	}


  	void store_offsets(Offset noff = Offset()) {
  		V[0]->store_offset(noff.append(o[0]));
  		V[1]->store_offset(noff.append(o[1]));
  		V[2]->store_offset(noff.append(o[2]));
  		//cout << "Now the vertices of face " << get_number() << " store offsets " << V[0]->get_offset() << ", " << V[1]->get_offset() << ", " << V[2]->get_offset() << endl;
  	}

  	void restore_offsets(Offset loff = Offset()) {
  		for (int i = 0; i < 3; i++) {
  			o[i] = loff.inverse().append(V[i]->get_offset());
  		}  		
  	}

  	void reorient() {
  		// N(eighbors), V(ertices), o(ffsets), no (neighbor offsets)

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
operator>>(std::istream &is, Periodic_4_hyperbolic_triangulation_ds_face_base_2<TDS> &)
  // non combinatorial information. Default = nothing
{
  return is;
}

template < class TDS >
inline
std::ostream&
operator<<(std::ostream &os,
    const Periodic_4_hyperbolic_triangulation_ds_face_base_2<TDS> &)
  // non combinatorial information. Default = nothing
{
  return os;
}

// Specialization for void.
template<typename GT>
class Periodic_4_hyperbolic_triangulation_ds_face_base_2<GT, void>
{
public:
  typedef Dummy_tds_2                  					Triangulation_data_structure;
  typedef Triangulation_data_structure::Vertex_handle   Vertex_handle;
  typedef Triangulation_data_structure::Face_handle     Face_handle;
  template <typename TDS2>
  struct Rebind_TDS {
    typedef Periodic_4_hyperbolic_triangulation_ds_face_base_2<GT, TDS2> Other;
  };
};


}  // namespace CGAL


#endif  // CGAL_PERIODIC_4_HYPERBOLIC_TIANGULATION_DS_FACE_BASE_2