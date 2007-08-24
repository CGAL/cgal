// Copyright (c) 2001-2004  ENS of Paris (France).
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
// $URL$
// $Id$
//
// Author(s)     : Pierre Angelier, Michel Pocchiola

#ifndef CGAL_VISIBILITY_COMPLEX_2_EDGE_BASE_H
#define CGAL_VISIBILITY_COMPLEX_2_EDGE_BASE_H

#include <cmath>
#include <map>

CGAL_BEGIN_NAMESPACE
namespace Visibility_complex_2_details {


template< class Vc_ >
class Edge_base;

//------------------------- Visibility_complex_edge_base class ----------------------

template < class Vc_ >
class Edge_base 
    : public Vc_::Gt::Arc_2 
{
public:
    // -------------------------------------------------------------------------
    // Geometry types
    typedef typename Vc_::Gt                    Gt;
    typedef typename Gt::Arc_2                  Arc_2;
    typedef Arc_2                               CA;
    typedef typename Gt::Disk                   Disk;
    typedef typename CA::Disk_handle            Disk_handle;
    // -------------------------------------------------------------------------
    typedef typename Vc_::Vertex                Vertex;
    typedef typename Vc_::Vertex_handle         Vertex_handle;
    typedef typename Vc_::Edge                  Edge;
    typedef typename Vc_::Edge_handle           Edge_handle;
    typedef typename Vc_::Face                  Face;
    typedef typename Vc_::Face_handle           Face_handle;
    typedef typename Vc_::Vertex_const_handle   Vertex_const_handle;
    typedef typename Vc_::Edge_const_handle     Edge_const_handle;
    typedef typename Vc_::Face_const_handle     Face_const_handle;
    typedef Edge_base<Vc_>   Self;
    // -------------------------------------------------------------------------
private:
    bool            sign_;
    Vertex_handle   sup_;
    Vertex_handle   inf_;
// public:    Edge_handle     opposite_;
private:    Face_handle     face[3];

public:
  using CA::object;

    // CONSTRUCTORS ------------------------------------------------------------
    Edge_base() 
	: CA()   , sign_(true) ,
	  sup_(0) , inf_(0)// , opposite_(0) 
    { face[0] = 0; face[1] = 0; face[2] = 0; }
    Edge_base(bool s, Disk_handle p);
    Edge_base(Vertex_handle v0 , Vertex_handle v1 , 
				 Disk_handle p);
    ~Edge_base() {
	if (sup() != 0 && 
            (static_cast<Edge_handle>(this)) == sup()->source_cusp_edge())
	    sup()->set_source_cusp_edge(0);
	if (sup() != 0 && 
            (static_cast<Edge_handle>(this)) == sup()->target_cusp_edge())
	    sup()->set_target_cusp_edge(0);
	if (sup() != 0 && (object() == sup()->source_object() || 
			   object() == sup()->target_object()) 
		       && (static_cast<Edge_handle>(this)) == 
                          sup()->cw_edge(object()))
	    sup()->set_cw_edge(0,object());
	if (inf() != 0 && 
            (object() == inf()->source_object() || 
             object() == inf()->target_object()) && 
            (static_cast<Edge_handle>(this)) == inf()->ccw_edge(object()))
	    inf()->set_ccw_edge(0,object());
	if (face[0] != 0 && 
            (static_cast<Edge_handle>(this)) == face[0]->top_edge())
	    face[0]->set_top_edge(0);
	if (face[1] != 0 && 
            (static_cast<Edge_handle>(this)) ==
            (sign()?face[1]->bottom_edge():face[1]->top_edge()))
          sign()?face[1]->set_bottom_edge(0):face[1]->set_top_edge(0);
	if (face[2] != 0 && 
            (static_cast<Edge_handle>(this)) == face[2]->bottom_edge())
	    face[2]->set_bottom_edge(0);

    }

    Self& operator=(const Self& e) { Arc_2::operator=(e); return *this; }

    Face_handle dl() { return face[0]; }
    Face_handle ul() { return face[2]; }
    Face_handle dr() { return (sign()) ? face[0] : face[1]; }
    Face_handle ur() { return (sign()) ? face[1] : face[2]; }

    Face_const_handle dl() const { return face[0]; }
    Face_const_handle ul() const { return face[2]; }
    Face_const_handle dr() const { return (sign()) ? face[0] : face[1]; }
    Face_const_handle ur() const { return (sign()) ? face[1] : face[2]; }


    void set_adjacent_faces(const Face_handle& f0 , 
			    const Face_handle& f1 , 
			    const Face_handle& f2);

    bool sign()      const { return sign_;     }
    void set_sign(bool b)  { sign_ = b;        }
    void flip_sign()       { sign_ = !sign_;   }

    Vertex_handle       sup()       { return sup_; }
    Vertex_const_handle sup() const { return sup_; }
    void           set_sup(const Vertex_handle& v) {
      typedef typename Vc_::Bitangent_2 Bitangent_2;
      sup_ = v; }
    Vertex_handle       inf()       { return inf_;}
    Vertex_const_handle inf() const { return inf_;}
    void           set_inf(const Vertex_handle& v) { inf_ = v; }

    bool is_regular() const {
      return object();
    }
};



template< class Vc_ >
Edge_base<Vc_>::
Edge_base(bool s, Disk_handle p) : CA(p)
{
    inf_ = 0;
    sup_ = 0;
    face[0] = 0; face[1] = 0; face[2] = 0;
    sign_ = s;
}


template< class Vc_ >
Edge_base<Vc_>::
	    Edge_base(Vertex_handle b0 , Vertex_handle b1 ,
					 Disk_handle p)
    : Arc_2(p,(b0->source_object() == p) ? b0->source() : b0->target(),
		     (b1->source_object() == p) ? b1->source() : b1->target())
{
    inf_ = b0;
    sup_ = b1;
    inf()->set_ccw_edge(static_cast<Edge_handle>(this));
    sup()->set_cw_edge (static_cast<Edge_handle>(this));
    face[0] = 0; face[1] = 0; face[2] = 0;

    sign_ = (inf()->source_object() == p) ? 
	inf()->is_left_xx() : inf()->is_xx_left();

}


template< class Vc_ >
inline void 
Edge_base<Vc_>::set_adjacent_faces(const Face_handle& f0 , 
						      const Face_handle& f1 , 
						      const Face_handle& f2)
{
    face[0] = f0; face[1] = f1; face[2] = f2;
}


// -----------------------------------------------------------------------------
}
CGAL_END_NAMESPACE

#endif // CGAL_VISIBILITY_COMPLEX_2_EDGE_BASE_H
