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
#include <CGAL/Hyperbolic_word_4.h>
#include <CGAL/Triangulation_data_structure_2.h>

namespace CGAL {

template< typename GT, typename TDS = Triangulation_data_structure_2<> >
class Periodic_4_hyperbolic_triangulation_ds_face_base_2 {
public:
	typedef TDS 							Triangulation_data_structure;
	typedef typename TDS::Vertex_handle		Vertex_handle;
	typedef typename TDS::Face_handle		Face_handle;
	typedef typename TDS::Vertex 			Vertex;
	typedef typename TDS::Face 				Face;
	typedef unsigned short int 				Int;
	typedef Hyperbolic_word_4<Int, GT>		Offset;

	template< typename TDS2 >
	struct Rebind_TDS {
		typedef Periodic_4_hyperbolic_triangulation_ds_face_base_2<GT, TDS2> Other;
	};

private:
	Face_handle 	N[3];
	Vertex_handle	V[3];
	Offset 			o[3];	// Offsets for vertices
	Offset 			no[3]; 	// Offsets for neighboring faces (neighbor face i is opposite to vertex i)
	int 			face_number;

public:

	Periodic_4_hyperbolic_triangulation_ds_face_base_2() 
#ifndef CGAL_CFG_NO_CPP0X_UNIFIED_INITIALIZATION_SYNTAX
	: o{Offset(), Offset(), Offset()}, no{Offset(), Offset(), Offset()}
	{}
#else
	{
		set_offsets();
	}
#endif

	Periodic_4_hyperbolic_triangulation_ds_face_base_2(
		const Vertex_handle& v0, const Vertex_handle& v1,
		const Vertex_handle& v2) 
#ifndef CGAL_CFG_NO_CPP0X_UNIFIED_INITIALIZATION_SYNTAX
	: V{v0, v1, v2},
	  o{Offset(), Offset(), Offset()}, no{Offset(), Offset(), Offset()}
	{
		set_neighbors();
	}
#else
	{
		set_offsets();
		set_vertices(v0, v1, v2);
		set_neighbors();
	}
#endif



	Periodic_4_hyperbolic_triangulation_ds_face_base_2(
		const Vertex_handle& v0, const Vertex_handle& v1,
		const Vertex_handle& v2, const Face_handle&   n0,
		const Face_handle&   n1, const Face_handle&   n2) 
#ifndef CGAL_CFG_NO_CPP0X_UNIFIED_INITIALIZATION_SYNTAX
	: V{v0, v1, v2},
	  N{n0, n1, n2},
	  o{Offset(), Offset(), Offset()}, no{Offset(), Offset(), Offset()}
	{
		set_neighbors();
	}
#else
	{
		set_offsets();
		set_vertices(v0, v1, v2);
		set_neighbors(n0, n1, n2);
	}
#endif



	// ACCESS

	const Vertex_handle& vertex(int i) const {
    	CGAL_triangulation_precondition( i >= 0 && i <= 2 );
    	return V[i];
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
		CGAL_triangulation_precondition( i >= 0 && i <= 2 );
		return no[i];
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
		no[0] = Offset();
		no[1] = Offset();
		no[2] = Offset();
	}

	void set_offsets(
		const Offset& o0, const Offset& o1, 
		const Offset& o2) 
	{
		o[0] = o0;
		o[1] = o1;
		o[2] = o2;
	}

	void set_neighbor_face_offsets(
		const Offset& no0, const Offset& no1, 
		const Offset& no2) 
	{
		no[0] = no0;
		no[1] = no1;
		no[2] = no2;
	}

	void set_offsets(
		const Offset& o0,  const Offset& o1, 
		const Offset& o2,  const Offset& no0,
		const Offset& no1, const Offset& no2) 
	{
		o[0] = o0;
		o[1] = o1;
		o[2] = o2;
		no[0] = no0;
		no[1] = no1;
		no[2] = no2;
	}

	void set_vertices() {
		V[0] = V[1] = V[2] = Vertex_handle();
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

	  // CHECKING

  	// the following trivial is_valid allows
  	// the user of derived cell base classes
  	// to add their own purpose checking

  	bool is_valid(bool, int) const { 
  		return true; 
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