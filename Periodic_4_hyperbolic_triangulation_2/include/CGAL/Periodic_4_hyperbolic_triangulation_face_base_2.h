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
#include <CGAL/Triangulation_ds_face_base_2.h>


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


template< typename GT, typename FB = Triangulation_ds_face_base_2<> >
class Periodic_4_hyperbolic_triangulation_face_base_2 : public FB {
	
public:
	typedef typename FB::Vertex_handle		Vertex_handle;
	typedef typename FB::Face_handle		Face_handle;
	typedef Face_data 						TDS_data;
	typedef typename GT::Hyperbolic_translation				Hyperbolic_translation;

	template< typename TDS2 >
	struct Rebind_TDS {
		typedef typename FB::template Rebind_TDS<TDS2>::Other 				FB2;
		typedef Periodic_4_hyperbolic_triangulation_face_base_2<GT, FB2> 	Other;
	};

private:
	Hyperbolic_translation 			o[3];	// Hyperbolic_translations for vertices
	TDS_data 						_tds_data;

public:

	Periodic_4_hyperbolic_triangulation_face_base_2() 
	: FB()
	{	
		o[0] = Hyperbolic_translation();
		o[1] = Hyperbolic_translation();
		o[2] = Hyperbolic_translation();
	}

	Periodic_4_hyperbolic_triangulation_face_base_2(
		const Vertex_handle& v0, const Vertex_handle& v1,
		const Vertex_handle& v2) 
	: FB(v0, v1, v2)
	{
		o[0] = Hyperbolic_translation();
		o[1] = Hyperbolic_translation();
		o[2] = Hyperbolic_translation();
	}



	Periodic_4_hyperbolic_triangulation_face_base_2(
		const Vertex_handle& v0, const Vertex_handle& v1,
		const Vertex_handle& v2, const Face_handle&   n0,
		const Face_handle&   n1, const Face_handle&   n2) 
	: FB(v0,v1,v2,n0,n1,n2)
	{
		o[0] = Hyperbolic_translation();
		o[1] = Hyperbolic_translation();
		o[2] = Hyperbolic_translation();
	}


	TDS_data& tds_data() {
		return _tds_data;
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
			if (this->neighbor(i)->vertex(c) == this->vertex(myi)) {
				nbof = this->neighbor(i)->translation(c);
				did_it = true;
				break;
			}
		}

		return (myof - nbof);
	}


	void set_translation(int k, Hyperbolic_translation new_o) {
		CGAL_triangulation_precondition( k >= 0 && k <= 2 );
		o[k] = new_o;
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


  	void reorient() {
  		
  		// swaps vertex 0 with vertex 1, and neighbor 0 with neighbor 1
  		FB::reorient();

  		// manually swap translation 0 with translation 1
  		Hyperbolic_translation tmp = o[0];
  		o[0] = o[1];
  		o[1] = tmp;
  	}

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


}  // namespace CGAL


#endif  // CGAL_PERIODIC_4_HYPERBOLIC_TIANGULATION_FACE_BASE_2