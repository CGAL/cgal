// Copyright (c) 2005  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_NEF_SHELL_CONSTRUCTOR_H
#define CGAL_NEF_SHELL_CONSTRUCTOR_H

namespace CGAL {

template<class Nef3, typename forward_iterator>
class Shell_constructor : public Modifier_base<typename Nef3::SNC_structure> {

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

  std::vector<Vertex_iterator>   Vertex_of;

  std::istream in;


  bool SNC_io_parser<EW>::check_sep(const char* sep) const {
    char c; 
    do in.get(c); while (isspace(c));
    while (*sep != '\0') { 
      if (*sep != c) {
	in.putback(c);
	return false;
      }
      ++sep; in.get(c);
    }
    in.putback(c);
    return true;  
  }

  void add_sedge(SNC_structure& snc, int prev, int cur, int next) {
    Vertex_handle v = Vertex_of[cur];
    Sphere_point sp1(Vertex_of[prev]->point());
    Sphere_point sp2(Vertex_of[next]->point());
    sp1 = normalized(sp1);
    sp2 = normalized(sp2);

    SVertex_iterator sv;
    SVertex_handle sv1, sv2;
    for(sv = v->svertices_begin(); sv != v->svertices_end(); ++sv)
      if(sv->point() == sp1)
	sv1 = sv;
      else if(sv->point() == sp2)
	sv2 = sv;
	
    SM_decorator(&*v);
    if(sv1 == SVertex_handle())
      sv1 = SM.new_svertex(sp1);
    if(sv2 == SVertex_handle())
      sv2 = SM.new_svertex(sp2);
    SHalfedge_handle se = SM.new_shalfedge_pair(sv1, sv2);
    se->mark() = se->twin()->mark() = true;
  }

  Vertex_handle read_vertex(SNC_structure& snc) {
    RT hx, hy, hz;
    in >> hx >> hy >> hz;;
    return v(snc.new_vertex(Point_3(hx,hy,hz,1),true));  
  }
  
  void read_facet(SNC_structure& snc) {
    int n, prev, cur, next, first, second; 
    in >> n;
    CGAL_assertion(n > 2);
    in >> first;
    in >> second;
    prev = first;
    cur = second;
    n-=2;
    for(;n>0;--n) {
      in >> next;
      add_sedge(snc, prev, cur, next);
      prev = cur;
      cur = next;
    }
    add_sedge(snc, cur, next, first);
    add_sedge(snc, next, first, second);
  }

 public:
  Shell_constructor(std::istream& instr) :
    in(instr) {}
    
    void operator()(SNC_structure& snc) {
      if (!check_sep("OFF"))  
	CGAL_error_msg("OFF header is missing!");
      int nv, nf, x;
      in >> nv;
      in >> nf;
      in >> x;
      Vertex_of.reserve(vn);
      for(x=0; x<nv; ++nv)
	Vertex_of[nv] = read_vertex(snc);
      for(x=0; x<vf; ++nf)
	read_facet(snc);
    }
};

} //namespace CGAL
#endif // CGAL_NEF_SHELL_CONSTRUCTOR_H
