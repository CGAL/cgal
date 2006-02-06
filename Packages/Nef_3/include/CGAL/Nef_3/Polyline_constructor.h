// Copyright (c) 2005  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_NEF_POLGON_CONSTRUCTOR_H
#define CGAL_NEF_POLGON_CONSTRUCTOR_H

CGAL_BEGIN_NAMESPACE

template<class Nef3, typename forward_iterator>
class Polyline_constructor : public Modifier_base<typename Nef3::SNC_structure> {

  typedef Nef3                                   Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_structure   SNC_structure;
  typedef typename SNC_structure::SM_decorator     SM_decorator;
  typedef typename SNC_structure::Vertex_handle    Vertex_handle;
  typedef typename SNC_structure::SVertex_handle   SVertex_handle;
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;
  typedef typename SNC_structure::SFace_handle     SFace_handle;
  typedef typename SNC_structure::Sphere_point     Sphere_point;
  typedef typename SNC_structure::Sphere_circle    Sphere_circle;
  typedef typename SNC_structure::Kernel           Kernel;
  typedef typename Kernel::Line_3                  Line_3;

  typedef typename std::iterator_traits<forward_iterator>::value_type
    point_iterator_pair;
  typedef typename point_iterator_pair::first_type
    point_iterator;

  forward_iterator begin, end;

 public:
  Polyline_constructor(forward_iterator begin_in, forward_iterator end_in) :
    begin(begin_in), end(end_in) {}

    void create_end_sphere_map(SNC_structure& snc,
			       point_iterator cur,
			       point_iterator prev) {
      Vertex_handle v(snc.new_vertex(*cur, true));
      SM_decorator SM(&*v);
      SVertex_handle sv(v->new_svertex(Sphere_point(ORIGIN+(*prev-*cur)),
				       true));
      SFace_handle sf(v->new_sface());
      SM.link_as_isolated_vertex(sv,sf);
    }

    void create_sphere_map(SNC_structure& snc,
			   point_iterator cur,
			   point_iterator prev,
			   point_iterator next) {
      Vertex_handle v(snc.new_vertex(*cur, true));
      SM_decorator SM(&*v);
      SVertex_handle sv1(v->new_svertex(Sphere_point(ORIGIN+(*prev-*cur)),
					true));
      SVertex_handle sv2(v->new_svertex(Sphere_point(ORIGIN+(*next-*cur)),
					true));      
      SFace_handle sf(v->new_sface());
      SM.link_as_isolated_vertex(sv1,sf);
      SM.link_as_isolated_vertex(sv2,sf);
    }

 public:
    void operator()(SNC_structure& snc) {
      point_iterator pbegin, pend, pnext, pprev;
      for(;begin != end; ++end) {
	pend = begin->second;
	pprev = pnext = pbegin = begin->first;
	++pnext;
	create_end_sphere_map(snc,pbegin,pnext);
	for(++pbegin,++pnext; pnext!=pend; ++pbegin,++pprev,++pnext)
	  create_sphere_map(snc,pbegin,pprev,pnext);
	create_end_sphere_map(snc,pbegin,pprev);
      }
    }
};

CGAL_END_NAMESPACE
#endif // CGAL_NEF_POLGON_CONSTRUCTOR_H
