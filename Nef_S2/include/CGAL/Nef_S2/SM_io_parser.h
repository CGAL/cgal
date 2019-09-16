// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Michael Seel  <seel@mpi-sb.mpg.de>
//                 Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_SM_IO_PARSER_H
#define CGAL_SM_IO_PARSER_H

#include <CGAL/license/Nef_S2.h>


#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_2/Object_index.h>
#include <CGAL/Nef_S2/SM_decorator_traits.h>
#include <vector>
#include <iostream>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4355) // complaint about using 'this' to
#endif                          // initialize a member

namespace CGAL {

/*{\Moptions outfile=SM_io_parser.man }*/
/*{\Manpage {SM_io_parser}{Decorator_}{IO of embedded maps}{IO}}*/

/*{\Mdefinition An instance |\Mvar| of the data type |\Mname| is a
decorator to provide input and output of a embedded map.  |\Mtype| is
generic with respect to the |Decorator_| parameter.  |Decorator_| has
to be a decorator model of our |SM_decorator| concept.}*/

/*{\Mgeneralization SM_decorator}*/

template <typename Decorator_>
class SM_io_parser : public Decorator_
{
  typedef Decorator_                         Base;
  typedef typename Decorator_::Sphere_point  Sphere_point;
  typedef typename Decorator_::Sphere_circle Sphere_circle;
  typedef typename Decorator_::Mark          Mark;

  typedef typename Decorator_::Decorator_traits  Decorator_traits;

  typedef typename Decorator_traits::SVertex_iterator     SVertex_iterator;
  typedef typename Decorator_traits::SHalfedge_iterator   SHalfedge_iterator;
  typedef typename Decorator_traits::SFace_iterator       SFace_iterator;
  typedef typename Decorator_traits::SVertex_handle             SVertex_handle;
  typedef typename Decorator_traits::SVertex_const_handle       SVertex_const_handle;
  typedef typename Decorator_traits::SHalfedge_handle           SHalfedge_handle;
  typedef typename Decorator_traits::SHalfedge_const_handle     SHalfedge_const_handle;
  typedef typename Decorator_traits::SFace_const_handle         SFace_const_handle;
  typedef typename Decorator_traits::SFace_handle               SFace_handle;
  typedef typename Decorator_traits::SHalfloop_handle           SHalfloop_handle;
  typedef typename Decorator_traits::SHalfloop_const_handle     SHalfloop_const_handle;
  typedef typename Decorator_traits::SFace_cycle_iterator SFace_cycle_iterator;
  typedef typename Decorator_traits::SHalfedge_around_svertex_circulator
                                     SHalfedge_around_svertex_circulator;


  using Base::is_isolated;
  using Base::first_out_edge;
  using Base::out_edges;

  std::istream& in; std::ostream& out;
  bool verbose;
  // a reference to the IO object
  CGAL::Object_index<SVertex_const_handle>   VI;
  CGAL::Object_index<SHalfedge_const_handle> EI;
  CGAL::Object_index<SFace_const_handle>     FI;
  std::vector<SVertex_handle>        SVertex_of;
  std::vector<SHalfedge_handle>      Edge_of;
  std::vector<SFace_handle>          SFace_of;
  SHalfloop_handle                   Loop_of[2];
  // object mapping for input
  std::size_t vn,en,ln,fn,i;
  // the number of objects

  bool check_sep(const char* sep);
  void print_vertex(SVertex_handle) const;
  void print_edge(SHalfedge_handle) const;
  void print_loop(SHalfloop_const_handle) const;
  void print_face(SFace_handle) const;

  bool read_vertex(SVertex_handle);
  bool read_edge(SHalfedge_handle);
  bool read_loop(SHalfloop_handle);
  bool read_face(SFace_handle);

  void debug_vertex(SVertex_handle) const;
  void debug_edge(SHalfedge_handle) const;
  void debug_loop(SHalfloop_const_handle) const;

public:
/*{\Mcreation 3}*/
SM_io_parser(std::istream& is, const Base& D);
/*{\Mcreate creates an instance |\Mvar| of type |\Mname|
to input |H| from |is|.}*/

SM_io_parser(std::ostream& os, const Base& D);
/*{\Mcreate creates an instance |\Mvar| of type |\Mname|
to output |H| to |os|.}*/

/*{\Moperations 2 3}*/
void print() const;
/*{\Mop prints |H| to |os|.}*/
void read();
/*{\Mop reads |H| from |is|.}*/
void debug() const;
void print_faces() const;

std::string index(SVertex_const_handle v) const 
{ return VI(v,verbose); }
std::string index(SHalfedge_const_handle e) const 
{ return EI(e,verbose); }
std::string index(SHalfloop_const_handle l) const 
{ if (verbose)  return (l==this->shalfloop()? "l0" : "l1");
  else return (l==this->shalfloop()? "0" : "1");
}
std::string index(SFace_const_handle f) const 
{ return FI(f,verbose); }

static void dump(const Decorator_& D, std::ostream& os = std::cerr);
/*{\Mstatic prints the plane map decorated by |D| to |os|.}*/

}; // SM_io_parser<Decorator_>


template <typename Decorator_>
SM_io_parser<Decorator_>::
SM_io_parser(std::istream& iin, const Base& H) :
  Base(H), in(iin), out(std::cout), verbose(0), 
  vn(0), en(0), ln(0), fn(0)
{ this->clear(); }

template <typename Decorator_>
SM_io_parser<Decorator_>::
SM_io_parser(std::ostream& iout, const Base& D) 
  : Base(D), in(std::cin), out(iout), 
  VI(this->svertices_begin(),this->svertices_end(),'v'),
  EI(this->shalfedges_begin(),this->shalfedges_end(),'e'),
  FI(this->sfaces_begin(),this->sfaces_end(),'f'),
  vn(this->number_of_svertices()), 
  en(this->number_of_shalfedges()), 
  ln(this->number_of_shalfloops()),
  fn(this->number_of_sfaces())
{ verbose = (get_mode(out) != CGAL::IO::ASCII &&
             get_mode(out) != CGAL::IO::BINARY);
}


//-----------------------------------------------------------------------------
// OUTPUT AND INPUT:
//-----------------------------------------------------------------------------

template <typename Decorator_>
bool SM_io_parser<Decorator_>::check_sep(const char* sep)
{
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

template <typename Decorator_>
void SM_io_parser<Decorator_>::print_vertex(SVertex_handle v) const
{
  // syntax: index { isolated incident_object, mark, point }
  out << index(v) << " { ";
  if ( is_isolated(v) ) out << "1 " << index(v->incident_sface());
  else                  out << "0 " << index(first_out_edge(v));
  out  << ", " << v->mark() << ", " << v->point() <<  "}\n";
}

template <typename Decorator_>
bool SM_io_parser<Decorator_>::read_vertex(SVertex_handle v)
{ 
  // precondition: nodes exist
  // syntax: index { isolated incident_object, mark, point}
  int n; bool iso; int f; Mark m; Sphere_point p; 
  if ( !(in >> n) ||
       !check_sep("{") ||
       !(in >> iso) ||
       !(in >> f) ||
       !check_sep(",") ||
       !(in >> m) ||
       !check_sep(",") ||
       !(in >> p) ||
       !check_sep("}") ) return false;
 
  if (iso) set_face(v,SFace_of[f]);
  else     set_first_out_edge(v,Edge_of[f]);
  v->mark() = m; v->point() = p;
  return true; 
}

template <typename Decorator_>
void SM_io_parser<Decorator_>::print_edge(SHalfedge_handle e) const
{ // syntax: index { twin, prev, next, source, face, mark, circle }

  Decorator_ D;
  out << index(e) << " { "
      << index(e->twin()) << ", " 
      << index(e->sprev()) << ", " << index(e->snext()) << ", "
      << index(e->source()) << ", " << index(e->incident_sface()) << ", "
      << e->mark() << ", " << e->circle() << " }\n";
}

template <typename Decorator_>
bool SM_io_parser<Decorator_>::read_edge(SHalfedge_handle e)
{ // syntax: index { twin, prev, next, source, face, mark, circle }
  int n, eo, epr, ene, v, f; bool m; Sphere_circle k;
  if ( !(in >> n) ||
       !check_sep("{") ||
       !(in >> eo) || !check_sep(",") ||
       !(in >> epr) || !check_sep(",") ||
       !(in >> ene) || !check_sep(",") ||
       !(in >> v) || !check_sep(",") ||
       !(in >> f) || !check_sep(",") ||
       !(in >> m) || !check_sep(",") ||
       !(in >> k) || !check_sep("}") )
    return false;
  CGAL_assertion_msg 
     (eo >= 0 && eo < en && epr >= 0 && epr < en && ene >= 0 && ene < en &&
      v >= 0 && v < vn && f >= 0 && f < fn ,
      "wrong index in read_edge");
  
  // precond: features exist!
  CGAL_assertion(EI[e->twin()]);
  set_prev(e,Edge_of[epr]);
  set_next(e,Edge_of[ene]);
  set_source(e,SVertex_of[v]);
  set_face(e,SFace_of[f]);
  e->mark() = m;
  e->circle() = k;
  return true;
}

template <typename Decorator_>
void SM_io_parser<Decorator_>::print_loop(SHalfloop_const_handle l) const
{ // syntax: index { twin, face, mark, circle }
  out << index(l) << " { "
      << index(l->twin()) << ", " 
      << index(l->incident_sface()) << ", "
      << l->mark() << ", " << l->circle() << " }\n";
}

template <typename Decorator_>
bool SM_io_parser<Decorator_>::read_loop(SHalfloop_handle l)
{ // syntax: index { twin, face, mark, circle }
  int n, lo, f; bool m; Sphere_circle k;
  if ( !(in >> n) ||
       !check_sep("{") ||
       !(in >> lo) || !check_sep(",") ||
       !(in >> f) || !check_sep(",") ||
       !(in >> m) || !check_sep(",") ||
       !(in >> k) || !check_sep("}") )
    return false;
  CGAL_assertion_msg(
    (lo >= 0 && lo < 2 && f >= 0 && f < fn),"wrong index in read_edge");
  
  set_face(l,SFace_of[f]);
  l->mark() = m;
  l->circle() = k;
  return true;
}


template <typename Decorator_>
void SM_io_parser<Decorator_>::print_face(SFace_handle f) const
{ // syntax: index { fclist, ivlist, loop, mark }
  out << index(f) << " { "; 
  SFace_cycle_iterator it;
  CGAL_forall_sface_cycles_of(it,f)
    if ( it.is_shalfedge() ) out << index(SHalfedge_handle(it)) << ' ';
  out << ", ";
  CGAL_forall_sface_cycles_of(it,f)
    if ( it.is_svertex() ) out << index(SVertex_handle(it)) << ' ';
  out << ", ";
  CGAL_forall_sface_cycles_of(it,f)
    if ( it.is_shalfloop() ) out << index(SHalfloop_handle(it));
  out << ", " << f->mark() << " }\n";
}

template <typename Decorator_>
bool SM_io_parser<Decorator_>::read_face(SFace_handle f)
{ // syntax: index { fclist, ivlist, loop, mark }
  int n, ei, vi, li; Mark m;
  if ( !(in >> n) || !check_sep("{") ) return false;
  while (in >> ei) { 
    CGAL_assertion_msg(ei >= 0 && ei < en, 
                           "wrong index in face cycle list.");
    store_sm_boundary_object(Edge_of[ei],f);
  } in.clear();
  if (!check_sep(",")) { return false; }
  while (in >> vi) { 
    CGAL_assertion_msg(vi >= 0 && vi < vn, 
                           "wrong index in iso vertex list.");
    store_sm_boundary_object(SVertex_of[vi],f);
  } in.clear();
  if (!check_sep(",")) { return false; }
  while (in >> li) { 
    CGAL_assertion_msg(li >= 0 && li < 2, 
                           "wrong index in iso vertex list.");
    store_sm_boundary_object(Loop_of[li],f);
  } in.clear();
  if (!check_sep(",") || !(in >> m) || !check_sep("}") ) 
    return false;
  f->mark() = m;
  return true;
}

template <typename Decorator_>
void SM_io_parser<Decorator_>::print() const
{
  out << "Sphere_map_2" << std::endl;
  out << "vertices "  << vn << std::endl;
  out << "edges "     << en << std::endl;
  out << "loops "     << ln << std::endl;
  out << "faces "     << fn << std::endl;
  if (verbose) 
    out << "/* index { isolated ? face : edge, mark, point } */" << std::endl;
  SVertex_iterator vit;
  CGAL_forall_svertices(vit,*this) print_vertex(vit);
  if (verbose) 
    out << "/* index { twin, prev, next, source, face, mark, circle } */" 
	<< std::endl;
  SHalfedge_iterator eit;
  CGAL_forall_shalfedges(eit,*this) print_edge(eit);
  if (verbose) 
    out << "/* index { twin, face, mark, circle } */" << std::endl;
  if ( this->has_shalfloop() ) 
    { print_loop(this->shalfloop()); print_loop(this->shalfloop()->twin()); }
  if (verbose) 
    out << "/* index { fclist, ivlist, loop, mark } */" << std::endl;
  SFace_iterator fit;
  CGAL_forall_sfaces(fit,*this) print_face(fit);
  out.flush();
  if (verbose) debug();
}

template <typename Decorator_>
void SM_io_parser<Decorator_>::read() 
{
  if ( !check_sep("Sphere_map_2") )  
    CGAL_error_msg("SM_io_parser::read: no embedded_PM header.");
  if ( !(check_sep("vertices") && (in >> vn)) ) 
    CGAL_error_msg("SM_io_parser::read: wrong vertex line.");
  if ( !(check_sep("edges") && (in >> en) && (en%2==0)) )
    CGAL_error_msg("SM_io_parser::read: wrong edge line.");
  if ( !(check_sep("loops") && (in >> ln)) )
    CGAL_error_msg("SM_io_parser::read: wrong loop line.");
  if ( !(check_sep("faces") && (in >> fn)) )
    CGAL_error_msg("SM_io_parser::read: wrong face line.");

  SVertex_of.resize(vn);
  Edge_of.resize(en);
  SFace_of.resize(fn);
  for(i=0; i<vn; i++)  SVertex_of[i] =   this->new_svertex();
  for(i=0; i<en; i++) 
    if (i%2==0) Edge_of[i] = this->new_shalfedge_pair();
    else Edge_of[i] = Edge_of[i-1]->twin();
  for(i=0; i<fn; i++)  SFace_of[i] =     this->new_sface();
  if ( ln == 2 ) { 
    Loop_of[0] = this->new_shalfloop_pair(); 
    Loop_of[1] = this->shalfloop()->twin(); 
  }

  for(i=0; i<vn; i++) {
    if (!read_vertex(SVertex_of[i]))
      CGAL_error_msg("SM_io_parser::read: error in node line");
  }
  for(i=0; i<en; i++) {
    if (!read_edge(Edge_of[i]))
      CGAL_error_msg("SM_io_parser::read: error in edge line");
  }
  if ( ln == 2 ) {
    read_loop(Loop_of[0]); read_loop(Loop_of[1]);
  }
  for(i=0; i<fn; i++) {
    if (!read_face(SFace_of[i]))
      CGAL_error_msg("SM_io_parser::read: error in face line");
  }
}

//-----------------------------------------------------------------------------
// VERBOSE OUTPUT:
// note that we output the index of the objects which is stored in them
// this is NOT the member index as produced by the forall loops
//-----------------------------------------------------------------------------

template <typename Decorator_>
void SM_io_parser<Decorator_>::debug_vertex(SVertex_handle v) const
{ 
  out << index(v) << "[" << v->mark() << "," << v->point() << "]" << std::endl; 
}

template <typename Decorator_>
void SM_io_parser<Decorator_>::debug_edge(SHalfedge_handle e) const
{ 
  out << index(e)
      << "(" << index(e->source()) << "," << index(e->target()) << ") "
      << index(e->twin()) << " " << index(e->incident_sface())
      << " ["<< e->mark() << "," << e->circle() << "] " << std::endl;
}

template <typename Decorator_>
void SM_io_parser<Decorator_>::debug_loop(SHalfloop_const_handle l) const
{ 
  out << index(l) << " "
      << index(l->twin()) << " " << index(l->incident_sface())
      << " ["<< l->mark() << "] " << l->circle() << std::endl;
}


template <typename Decorator_>
void SM_io_parser<Decorator_>::debug() const
{ 
  out << "\nDEBUG Plane_map\n";
  out << "Vertices:  " << this->number_of_svertices() << "\n";
  out << "SHalfedges: " << this->number_of_shalfedges() << "\n";
  out << "Loop:      " << this->number_of_shalfloops() << "\n";
  SVertex_iterator vit; 
  CGAL_forall_svertices(vit,*this) {
    if ( is_isolated(vit) ) continue;
    SHalfedge_around_svertex_circulator hcirc(out_edges(vit)), hend(hcirc);
    debug_vertex(vit);
    CGAL_For_all(hcirc,hend) { out << "  "; debug_edge(hcirc); }
  }
  if ( this->has_shalfloop() ) 
    { debug_loop(this->shalfloop()); debug_loop(this->shalfloop()->twin()); }
  out << std::endl;
}

template <typename Decorator_>
void SM_io_parser<Decorator_>::print_faces() const
{ 
  out << "\nFACES\n";
  out << "Vertices:  " << this->number_of_svertices() << "\n";
  out << "SHalfedges: " << this->number_of_shalfedges() << "\n";
  out << "Loop:      " << this->number_of_shalfloops() << "\n";
  SHalfedge_iterator e;
  Unique_hash_map<SHalfedge_iterator,bool> Done(false);
  CGAL_forall_shalfedges(e,*this) {
    if ( Done[e] ) continue;
    typename Base::SHalfedge_around_sface_circulator c(e), ce = c;
    out << "face cycle\n";
    CGAL_For_all(c,ce) 
    { Done[c]=true; out << "  "; debug_vertex(c->source()); }
  }
  if ( this->has_shalfloop() ) 
    { debug_loop(this->shalfloop()); debug_loop(this->shalfloop()->twin()); }
  out << std::endl;
}

template <typename Decorator_>
void SM_io_parser<Decorator_>::dump(const Decorator_& D, std::ostream& os)
{ SM_io_parser<Decorator_> Out(os,D);
  Out.print();
  Out.print_faces();
}



} //namespace CGAL


#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif //CGAL_SM_IO_PARSER_H
