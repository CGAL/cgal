// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_PM_IO_PARSER_H
#define CGAL_PM_IO_PARSER_H

#include <CGAL/license/Nef_2.h>


#include <CGAL/Nef_2/PM_decorator.h>
#include <CGAL/Nef_2/Object_index.h>
#include <vector>

namespace CGAL {

/*{\Moptions outfile=PM_io_parser.man }*/
/*{\Manpage {PM_io_parser}{PMDEC}{IO of plane maps}{IO}}*/

/*{\Mdefinition An instance |\Mvar| of the data type |\Mname| is a
decorator to provide input and output of a plane map.  |\Mtype| is
generic with respect to the |PMDEC| parameter.  |PMDEC| has to be a
decorator model of our |PM_decorator| concept.}*/

/*{\Mgeneralization PM_decorator}*/

template <typename PMDEC>
class PM_io_parser : public PMDEC
{
  typedef PMDEC                     Base;
  typedef typename PMDEC::Plane_map Plane_map;
  typedef typename PMDEC::Point     Point;
  typedef typename PMDEC::Mark      Mark;

  typedef typename PMDEC::Vertex_iterator   Vertex_iterator;
  typedef typename PMDEC::Halfedge_iterator Halfedge_iterator;
  typedef typename PMDEC::Face_iterator     Face_iterator;
  typedef typename PMDEC::Vertex_handle     Vertex_handle;
  typedef typename PMDEC::Halfedge_handle   Halfedge_handle;
  typedef typename PMDEC::Face_handle       Face_handle;

  using Base::clear;
  using Base::vertices_begin;
  using Base::vertices_end;
  using Base::halfedges_begin;
  using Base::halfedges_end;
  using Base::faces_begin;
  using Base::faces_end;
  using Base::number_of_vertices;
  using Base::number_of_halfedges;
  using Base::number_of_faces;
  using Base::new_vertex;
  using Base::new_face;
  using Base::new_halfedge_pair_without_vertices;
  using Base::is_isolated;
  using Base::face;
  using Base::mark;
  using Base::point;
  using Base::twin;
  using Base::previous;
  using Base::next;
  using Base::source;
  using Base::target;
  using Base::out_edges;
  using Base::halfedge;

  std::istream& in;
  std::ostream& out;
  bool verbose;
  // a reference to the IO object
  CGAL::Object_index<Vertex_handle>   VI;
  CGAL::Object_index<Halfedge_handle> EI;
  CGAL::Object_index<Face_handle>     FI;
  std::vector<Vertex_handle>    Vertex_of;
  std::vector<Halfedge_handle>  Halfedge_of;
  std::vector<Face_handle>      Face_of;
  // object mapping for input
  std::size_t vn,en,fn,i;
  // the number of objects


  void print_vertex(Vertex_handle) const;
  void print_hedge(Halfedge_handle) const;
  void print_face(Face_handle) const;

  bool read_vertex(Vertex_handle);
  bool read_hedge(Halfedge_handle);
  bool read_face(Face_handle);

  void debug_vertex(Vertex_handle) const;
  void debug_hedge(Halfedge_handle) const;

public:
/*{\Mcreation 3}*/
PM_io_parser(std::istream& is, Plane_map& H)
/*{\Mcreate creates an instance |\Mvar| of type |\Mname|
   to input |H| from |is|.}*/
    : Base(H), in(is), out(std::cout), verbose(0), vn(0), en(0), fn(0)
        { this->clear(); }


PM_io_parser(std::ostream& os, const Plane_map& H)
/*{\Mcreate creates an instance |\Mvar| of type |\Mname|
to output |H| to |os|.}*/
: Base(const_cast<Plane_map&>(H)), in(std::cin), out(os), 
  VI(Base::vertices_begin(), Base::vertices_end(),'v'),
  EI(Base::halfedges_begin(),Base::halfedges_end(),'e'),
  FI(Base::faces_begin(),Base::faces_end(),'f'),
  vn(Base::number_of_vertices()), 
  en(Base::number_of_halfedges()), 
  fn(Base::number_of_faces())
{ verbose = (get_mode(out) != CGAL::IO::ASCII &&
             get_mode(out) != CGAL::IO::BINARY);
}


PM_io_parser(std::ostream& os, const PMDEC& D)
: Base(D), in(std::cin), out(os), 
  VI(Base::vertices_begin(),Base::vertices_end(),'v'),
  EI(Base::halfedges_begin(),Base::halfedges_end(),'e'),
  FI(Base::faces_begin(),Base::faces_end(),'f'),
  vn(Base::number_of_vertices()), 
  en(Base::number_of_halfedges()), 
  fn(Base::number_of_faces())
{ verbose = (get_mode(out) != CGAL::IO::ASCII &&
             get_mode(out) != CGAL::IO::BINARY);
}


/*{\Moperations 2 3}*/

bool check_sep(const char* sep);

void print() const;
/*{\Mop prints |H| to |os|.}*/
void read();
/*{\Mop reads |H| from |is|.}*/
void debug() const;
std::string index(Vertex_handle v) const { return VI(v,verbose); }
std::string index(Halfedge_handle e) const { return EI(e,verbose); }
std::string index(Face_handle f) const { return FI(f,verbose); }

static void dump(const PMDEC& D, std::ostream& os = std::cerr);
/*{\Mstatic prints the plane map decorated by |D| to |os|.}*/

}; // PM_io_parser<PMDEC>


//-----------------------------------------------------------------------------
// OUTPUT AND INPUT:
//-----------------------------------------------------------------------------

template <typename PMDEC>
bool PM_io_parser<PMDEC>::check_sep(const char* sep)
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

template <typename PMDEC>
void PM_io_parser<PMDEC>::print_vertex(Vertex_handle v) const
{
  // syntax: index { isolated incident_object, mark, point }
  out << index(v) << " { ";
  if ( is_isolated(v) ) out << "1 " << index(face(v));
  else                  out << "0 " << index(v->halfedge());
  out  << ", " << mark(v) << ", " << point(v) <<  "}\n";
}

template <typename PMDEC>
bool PM_io_parser<PMDEC>::read_vertex(Vertex_handle v)
{ 
  // precondition: nodes exist
  // syntax: index { mark, point, isolated object }
  int n; bool iso; int f; Mark m; Point p; 
  if ( !(in >> n) ||
       !check_sep("{") ||
       !(in >> iso) ||
       !(in >> f) ||
       !check_sep(",") ||
       !(in >> m) ||
       !check_sep(",") ||
       !(in >> p) ||
       !check_sep("}") ) return false;
 
  if (iso) v->set_face(Face_of[f]);
  else     v->set_halfedge(Halfedge_of[f]);
  mark(v) = m; point(v) = p;
  return true; 
}

template <typename PMDEC>
void PM_io_parser<PMDEC>::print_hedge(Halfedge_handle e) const
{ // syntax: index { opposite, prev, next, vertex, face, mark }
  out << index(e) << " { "
      << index(twin(e)) << ", " 
      << index(previous(e)) << ", " << index(next(e)) << ", "
      << index(target(e)) << ", " <<   index(face(e)) << ", "
      << mark(e) << " }\n";
}

template <typename PMDEC>
bool PM_io_parser<PMDEC>::read_hedge(Halfedge_handle e)
{ // syntax: index { opposite, prev, next, vertex, face, mark }
  int n, eo, epr, ene, v, f; bool m; 
  if ( !(in >> n) ||
       !check_sep("{") ||
       !(in >> eo) || !check_sep(",") ||
       !(in >> epr) || !check_sep(",") ||
       !(in >> ene) || !check_sep(",") ||
       !(in >> v) || !check_sep(",") ||
       !(in >> f) || !check_sep(",") ||
       !(in >> m) || !check_sep("}") )
    return false;
  CGAL_assertion_msg 
     (eo >= 0 || (std::size_t) eo < en || epr >= 0 || (std::size_t) epr < en || ene >= 0 || (std::size_t) ene < en ||
      v >= 0 || (std::size_t) v < vn || f >= 0 || (std::size_t) f < fn ,
      "wrong index in read_hedge");
  
  // precond: objects exist!
  CGAL_assertion(EI[e->opposite()]);
  e->set_prev(Halfedge_of[epr]);
  e->set_next(Halfedge_of[ene]);
  e->set_vertex(Vertex_of[v]);
  e->set_face(Face_of[f]);
  mark(e) = m;
  return true;
}


template <typename PMDEC>
void PM_io_parser<PMDEC>::print_face(Face_handle f) const
{ // syntax: index { halfedge, fclist, ivlist, mark }
  out << index(f) << " { " << index(halfedge(f)) << ", ";
  typedef typename std::list<Halfedge_handle>::iterator lhiterator;
  lhiterator hit, hend = f->fc_end();
  for(hit = f->fc_begin(); hit!=hend; ++hit) 
    out << index(*hit) << ' ';
  out << ", ";
  typedef typename std::list<Vertex_handle>::iterator lviterator;
  lviterator vit, vend = f->iv_end();
  for(vit = f->iv_begin(); vit!=vend; ++vit) 
    out << index(*vit) << ' ';
  out << ", " << mark(f) << " }\n";
}

template <typename PMDEC>
bool PM_io_parser<PMDEC>::read_face(Face_handle f)
{ // syntax: index { halfedge, fclist, ivlist, mark }
  int n, ei, vi; Mark m;
  if ( !(in >> n) || !check_sep("{") ) return false;
  if ( !(in >> ei) || !check_sep(",") ) return false;
  if (ei >= 0) f->set_halfedge(Halfedge_of[ei]);
  while (in >> ei) { 
    CGAL_assertion_msg(ei >= 0 && (std::size_t) ei < en, "wrong index in face cycle list.");
    f->store_fc(Halfedge_of[ei]);
  } in.clear();
  if (!check_sep(",")) { return false; }
  while (in >> vi) { 
    CGAL_assertion_msg(vi >= 0 && (std::size_t) vi < vn, "wrong index in iso vertex list.");
    f->store_iv(Vertex_of[vi]);
  } in.clear();
  if (!check_sep(",") || !(in >> m) || !check_sep("}") ) 
    return false;
  mark(f) = m;
  return true;
}

template <typename PMDEC>
void PM_io_parser<PMDEC>::print() const
{
  out << "Plane_map_2" << std::endl;
  out << "vertices "  << vn << std::endl;
  out << "halfedges " << en << std::endl;
  out << "faces "     << fn << std::endl;
  if (verbose) 
    out << "/* index { isolated ? face : halfedge, mark, point } */" 
        << std::endl;
  Vertex_iterator vit, vend = this->vertices_end();
  for(vit = this->vertices_begin(); vit!=vend; ++vit) print_vertex(vit);
  if (verbose) 
    out << "/* index { opposite, prev, next, vertex, face, mark } */" 
        << std::endl;
  Halfedge_iterator eit, eend = this->halfedges_end();
  for(eit = this->halfedges_begin(); eit!=eend; ++eit) print_hedge(eit);
  if (verbose) 
    out << "/* index { halfedge, fclist, ivlist, mark } */" 
        << std::endl;
  Face_iterator fit, fend = this->faces_end();
  for(fit = this->faces_begin(); fit!=fend; ++fit) print_face(fit);
  out.flush();
  if (verbose) debug();
}

template <typename PMDEC>
void PM_io_parser<PMDEC>::read() 
{
  if ( !check_sep("Plane_map_2") )  
    CGAL_error_msg("PM_io_parser::read: no embedded_PM header.");
  if ( !(check_sep("vertices") && (in >> vn)) ) 
    CGAL_error_msg("PM_io_parser::read: wrong node line.");
  if ( !(check_sep("halfedges") && (in >> en) && (en%2==0)) )
    CGAL_error_msg("PM_io_parser::read: wrong edge line.");
  if ( !(check_sep("faces") && (in >> fn)) )
    CGAL_error_msg("PM_io_parser::read: wrong face line.");

  Vertex_of.resize(vn);
  Halfedge_of.resize(en);
  Face_of.resize(fn);

  for(i=0; i<vn; i++)  Vertex_of[i] =   this->new_vertex();
  for(i=0; i<en; i++) 
    if (i%2==0) Halfedge_of[i] = this->new_halfedge_pair_without_vertices();
    else Halfedge_of[i] = twin(Halfedge_of[i-1]);
  for(i=0; i<fn; i++)  Face_of[i] =     this->new_face();

  for(i=0; i<vn; i++) {
    if (!read_vertex(Vertex_of[i]))
      CGAL_error_msg("PM_io_parser::read: error in node line");
  }
  for(i=0; i<en; i++) {
    if (!read_hedge(Halfedge_of[i]))
      CGAL_error_msg("PM_io_parser::read: error in halfedge\
      line");
  }
  for(i=0; i<fn; i++) {
    if (!read_face(Face_of[i]))
      CGAL_error_msg("PM_io_parser::read: error in face line");
  }
}

//-----------------------------------------------------------------------------
// VERBOSE OUTPUT:
// note that we output the index of the objects which is stored in them
// this is NOT the member index as produced by the forall loops
//-----------------------------------------------------------------------------

template <typename PMDEC>
void PM_io_parser<PMDEC>::debug_vertex(Vertex_handle v) const
{ 
  out << index(v) << "[" << mark(v) << "," << point(v) << "]" 
      << std::endl; 
}

template <typename PMDEC>
void PM_io_parser<PMDEC>::debug_hedge(Halfedge_handle e) const
{ 
  out << index(e)
      << "(" << index(source(e)) << "," << index(target(e)) << ") "
      << index(face(e)) << " " << index(twin(e))
      << " ["<< mark(e) << "]" <<" "<<&*(face(e)) << std::endl;
}


template <typename PMDEC>
void PM_io_parser<PMDEC>::debug() const
{ 
  out << "\nDEBUG Plane_map\n"
      << "Vertices:  " << this->number_of_vertices() << "\n"
      << "Halfedges: " << this->number_of_halfedges() << "\n";
  Vertex_iterator vit,vend = this->vertices_end(); 
  for( vit = this->vertices_begin(); vit != vend; ++vit ) {
    if ( is_isolated(vit) ) continue;
    typename Base::Halfedge_around_vertex_circulator
      hcirc = out_edges(vit), hend = hcirc;
    debug_vertex(vit);
    CGAL_For_all(hcirc,hend) { out << "  "; debug_hedge(hcirc); }
  }
  out << std::endl;
}

template <typename PMDEC>
void PM_io_parser<PMDEC>::dump(const PMDEC& D, std::ostream& os)
{ PM_io_parser<PMDEC> Out(os,D);
  Out.print();
}



} //namespace CGAL
#endif //CGAL_PM_IO_PARSER_H
