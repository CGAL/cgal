// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : include/CGAL/Nef_3/SNC_SM_io_parser.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// maintainer    : Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// SNC_SM_io_parser.h              input and output of sphere maps
// ============================================================================
#ifndef CGAL_SNC_SM_IO_PARSER_H
#define CGAL_SNC_SM_IO_PARSER_H

#include <CGAL/Nef_3/SNC_SM_decorator.h>
#include <CGAL/Nef_2/Object_index.h>
#include <vector>

using std::endl;

CGAL_BEGIN_NAMESPACE

/*{\Moptions outfile=SNC_SM_io_parser.man }*/
/*{\Manpage {SNC_SM_io_parser}{Refs_}{IO of local graphs}{IO}}*/

/*{\Mdefinition An instance |\Mvar| of the data type |\Mname| is a
decorator to provide input and output of a sphere map associated to a
vertex of an SNC.  |\Mtype| is generic with respect to the |Refs_|
parameter. |Refs_| has to be a model of our |SNC| concept.}*/

/*{\Mgeneralization SNC_SM_decorator<Refs_>}*/

template <typename Refs_>
class SNC_SM_io_parser : public SNC_SM_decorator<Refs_>
{
  typedef SNC_SM_decorator<Refs_> Base;
  typedef typename Refs_::Vertex_handle Vertex_handle;
  typedef typename Refs_::Sphere_point  Sphere_point;
  typedef typename Refs_::Sphere_circle Sphere_circle;
  typedef typename Refs_::Mark          Mark;

  typedef typename Refs_::SVertex_iterator   SVertex_iterator;
  typedef typename Refs_::SHalfedge_iterator SHalfedge_iterator;
  typedef typename Refs_::SFace_iterator     SFace_iterator;
  typedef typename Refs_::SVertex_handle     SVertex_handle;
  typedef typename Refs_::SHalfedge_handle   SHalfedge_handle;
  typedef typename Refs_::SFace_handle       SFace_handle;
  typedef typename Refs_::SHalfloop_handle       SHalfloop_handle;
  typedef typename Refs_::SFace_cycle_iterator SFace_cycle_iterator;
  std::istream& in; std::ostream& out;
  // a reference to the IO object
  bool verbose;
  CGAL::Object_index<SVertex_handle>   VI;
  CGAL::Object_index<SHalfedge_handle> EI;
  CGAL::Object_index<SFace_handle>     FI;
  std::vector<SVertex_handle>    Vertex_of;
  std::vector<SHalfedge_handle>  Edge_of;
  std::vector<SFace_handle>      Face_of;
  SHalfloop_handle               Loop_of[2];
  // object mapping for input
  int vn, en, ln, fn, i;
  // the number of objects


  bool check_sep(char* sep);
  void print_vertex(SVertex_handle) const;
  void print_edge(SHalfedge_handle) const;
  void print_loop(SHalfloop_handle) const;
  void print_face(SFace_handle) const;

  bool read_vertex(SVertex_handle);
  bool read_edge(SHalfedge_handle);
  bool read_loop(SHalfloop_handle);
  bool read_face(SFace_handle);

  void debug_vertex(SVertex_handle) const;
  void debug_edge(SHalfedge_handle) const;
  void debug_loop(SHalfloop_handle l) const;

public:
/*{\Mcreation 3}*/
SNC_SM_io_parser(std::istream& is, Vertex_handle v);
/*{\Mcreate creates an instance |\Mvar| of type |\Mname|
to input the local graph |M| of |v| from |is|.}*/

SNC_SM_io_parser(std::ostream& os, Vertex_handle v);
/*{\Mcreate creates an instance |\Mvar| of type |\Mname|
to output the local graph |M| of |v| to |os|.}*/

/*{\Moperations 2 3}*/
void print() const;
/*{\Mop prints |M| to |os|.}*/
void read();
/*{\Mop reads |M| from |is|.}*/
void debug() const;
std::string index(SVertex_handle v) const 
{ return VI(v,verbose); }
std::string index(SHalfedge_handle e) const 
{ return EI(e,verbose); }
std::string index(SHalfloop_handle l) const 
{ if (verbose) return (l==shalfloop()? "l0" : "l1");
  else return (l==shalfloop()? "0" : "1");
}
std::string index(SFace_handle f) const 
{ return FI(f,verbose); }

static void dump(Vertex_handle v, std::ostream& os = std::cerr);
/*{\Mstatic prints the plane map decorated by |D| to |os|.}*/

}; // SNC_SM_io_parser<Refs_>


template <typename Refs_>
SNC_SM_io_parser<Refs_>::SNC_SM_io_parser(std::istream& iin, Vertex_handle v) :
  Base(v), in(iin), out(std::cout), verbose(0), 
  vn(0), en(0), ln(0), fn(0) 
{ clear(); }

template <typename Refs_>
SNC_SM_io_parser<Refs_>::SNC_SM_io_parser(std::ostream& iout, Vertex_handle v) 
  : Base(v), in(std::cin), out(iout), 
  VI(svertices_begin(),svertices_end(),'v'),
  EI(shalfedges_begin(),shalfedges_end(),'e'),
  FI(sfaces_begin(),sfaces_end(),'f'),
  vn(number_of_svertices()), 
  en(number_of_shalfedges()), 
  ln(number_of_shalfloops()),
  fn(number_of_sfaces())
{ verbose = (out.iword(CGAL::IO::mode) != CGAL::IO::ASCII &&
             out.iword(CGAL::IO::mode) != CGAL::IO::BINARY);
}


//-----------------------------------------------------------------------------
// OUTPUT AND INPUT:
//-----------------------------------------------------------------------------

template <typename Refs_>
bool SNC_SM_io_parser<Refs_>::check_sep(char* sep)
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

template <typename Refs_>
void SNC_SM_io_parser<Refs_>::print_vertex(SVertex_handle v) const
{
  // syntax: index { isolated incident_object, mark, point }
  out << index(v) << " { ";
  if ( is_isolated(v) ) out << "1 " << index(face(v));
  else                  out << "0 " << index(first_out_edge(v));
  out  << ", " << mark(v) << ", " << point(v) <<  "}\n";
}

template <typename Refs_>
bool SNC_SM_io_parser<Refs_>::read_vertex(SVertex_handle v)
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
 
  if (iso) set_face(v,Face_of[f]);
  else     set_first_out_edge(v,Edge_of[f]);
  mark(v) = m; point(v) = p;
  return true; 
}

template <typename Refs_>
void SNC_SM_io_parser<Refs_>::print_edge(SHalfedge_handle e) const
{ // syntax: index { twin, prev, next, source, face, mark, circle }
  out << index(e) << " { "
      << index(twin(e)) << ", " 
      << index(previous(e)) << ", " << index(next(e)) << ", "
      << index(source(e)) << ", " << index(face(e)) << ", "
      << mark(e) << ", " << circle(e) << " }\n"; 
}

template <typename Refs_>
bool SNC_SM_io_parser<Refs_>::read_edge(SHalfedge_handle e)
{ // syntax: index { twin, prev, next, source, face, mark, circle }
  int n, eo, epr, ene, v, f; bool m; Sphere_circle k;  
  if ( !(in >> n) ||
       !check_sep("{") ||
       !(in >> eo) || !check_sep(",") ||
       !(in >> epr) || !check_sep(",") ||
       !(in >> ene) || !check_sep(",") ||
       !(in >> v) || !check_sep(",") ||
       !(in >> f) || !check_sep(",") ||
       !(in >> m) || !check_sep("}") ||
       !(in >> k) || !check_sep("}") )
    return false;
  CGAL_nef3_assertion_msg 
     (eo >= 0 && eo < en && epr >= 0 && epr < en && ene >= 0 && ene < en &&
      v >= 0 && v < vn && f >= 0 && f < fn ,
      "wrong index in read_edge");
  
  // precond: features exist!
  CGAL_nef3_assertion(EI[twin(e)]);
  set_prev(e,Edge_of[epr]);
  set_next(e,Edge_of[ene]);
  set_source(e,Vertex_of[v]);
  set_face(e,Face_of[f]);
  mark(e) = m;
  circle(e) = k;
  return true;
}

template <typename Refs_>
void SNC_SM_io_parser<Refs_>::print_loop(SHalfloop_handle l) const
{ // syntax: index { twin, face, mark, circle }
  out << index(l) << " { "
      << index(twin(l)) << ", " 
      << index(face(l)) << ", "
      << mark(l) << ", " << circle(l) << " }\n";
}

template <typename Refs_>
bool SNC_SM_io_parser<Refs_>::read_loop(SHalfloop_handle l)
{ // syntax: index { twin, face, mark, circle }
  int n, lo, f; bool m; Sphere_circle k;  
  if ( !(in >> n) ||
       !check_sep("{") ||
       !(in >> lo) || !check_sep(",") ||
       !(in >> f) || !check_sep(",") ||
       !(in >> m) || !check_sep("}") ||
       !(in >> k) || !check_sep("}") )
    return false;
  CGAL_nef3_assertion_msg(
    (lo >= 0 && lo < 2 && f >= 0 && f < fn),"wrong index in read_edge");
  
  set_face(l,Face_of[f]);
  mark(l) = m;
  circle(l) = k;
  return true;
}


template <typename Refs_>
void SNC_SM_io_parser<Refs_>::print_face(SFace_handle f) const
{ // syntax: index { fclist, ivlist, loop, mark }
  out << index(f) << " { "; 
  SFace_cycle_iterator it;
  CGAL_nef3_forall_sface_cycles_of(it,f)
    if ( it.is_shalfedge() ) out << index(SHalfedge_handle(it)) << ' ';
  out << ", ";
  CGAL_nef3_forall_sface_cycles_of(it,f)
    if ( it.is_svertex() ) out << index(SVertex_handle(it)) << ' ';
  out << ", ";
  CGAL_nef3_forall_sface_cycles_of(it,f)
    if ( it.is_shalfloop() ) out << index(SHalfloop_handle(it));
  out << ", " << mark(f) << " }\n";
}

template <typename Refs_>
bool SNC_SM_io_parser<Refs_>::read_face(SFace_handle f)
{ // syntax: index { fclist, ivlist, loop, mark }
  int n, ei, vi, li; Mark m;
  if ( !(in >> n) || !check_sep("{") ) return false;
  while (in >> ei) { 
    CGAL_nef3_assertion_msg(ei >= 0 && ei < en, "wrong index in face cycle list.");
    store_boundary_object(Edge_of[ei],f);
  } in.clear();
  if (!check_sep(",")) { return false; }
  while (in >> vi) { 
    CGAL_nef3_assertion_msg(vi >= 0 && vi < vn, "wrong index in iso vertex list.");
    store_boundary_object(Vertex_of[vi],f);
  } in.clear();
  if (!check_sep(",")) { return false; }
  while (in >> li) { 
    CGAL_nef3_assertion_msg(li >= 0 && li < 2, "wrong index in iso vertex list.");
    store_boundary_object(Loop_of[li],f);
  } in.clear();
  if (!check_sep(",") || !(in >> m) || !check_sep("}") ) 
    return false;
  mark(f) = m;
  return true;
}

template <typename Refs_>
void SNC_SM_io_parser<Refs_>::print() const
{
  out << "Sphere Map"  << endl;
  out << "svertices "  << vn << endl;
  out << "shalfedges " << en << endl;
  out << "shalfloops " << ln << endl;
  out << "sfaces "     << fn << endl;
  if (verbose) 
    out << "/* index { isolated ? face : edge, mark, point } */" << endl;
  SVertex_iterator vit;
  CGAL_nef3_forall_svertices(vit,*this) print_vertex(vit);
  if (verbose) 
    out << "/* index { twin, prev, next, source, face, mark, circle } */" 
	<< endl;
  SHalfedge_iterator eit;
  CGAL_nef3_forall_shalfedges(eit,*this) print_edge(eit);
  if (verbose) 
    out << "/* index { twin, face, mark, circle } */" << endl;
  if ( has_loop() ) 
  { print_loop(shalfloop()); print_loop(twin(shalfloop())); }
  if (verbose) 
    out << "/* index { fclist, ivlist, loop, mark } */" << endl;
  SFace_iterator fit;
  CGAL_nef3_forall_sfaces(fit,*this) print_face(fit);
  out.flush();
  if (verbose) debug();
}

template <typename Refs_>
void SNC_SM_io_parser<Refs_>::read() 
{
  if ( !check_sep("Sphere Map") )  
    CGAL_nef3_assertion_msg(0,"SNC_SM_io_parser::read: no embedded_PM header.");
  if ( !(check_sep("svertices") && (in >> vn)) ) 
    CGAL_nef3_assertion_msg(0,"SNC_SM_io_parser::read: wrong vertex line.");
  if ( !(check_sep("shalfedges") && (in >> en) && (en%2==0)) )
    CGAL_nef3_assertion_msg(0,"SNC_SM_io_parser::read: wrong edge line.");
  if ( !(check_sep("shalfloops") && (in >> ln)) )
    CGAL_nef3_assertion_msg(0,"SNC_SM_io_parser::read: wrong loop line.");
  if ( !(check_sep("sfaces") && (in >> fn)) )
    CGAL_nef3_assertion_msg(0,"SNC_SM_io_parser::read: wrong face line.");

  Vertex_of.reserve(vn);
  Edge_of.reserve(en);
  Face_of.reserve(fn);
  for(i=0; i<vn; i++)  Vertex_of[i] = new_vertex();
  for(i=0; i<en; i++) 
    if (i%2==0) Edge_of[i] = new_edge_pair();
    else Edge_of[i] = twin(Edge_of[i-1]);
  for(i=0; i<fn; i++)  Face_of[i] = new_face();
  if ( ln == 2 ) 
  { Loop_of[0] = new_loop_pair(); Loop_of[1] = twin(shalfloop()); }

  for(i=0; i<vn; i++) {
    if (!read_vertex(Vertex_of[i]))
      CGAL_nef3_assertion_msg(0,"SNC_SM_io_parser::read: error in node line");
  }
  for(i=0; i<en; i++) {
    if (!read_edge(Edge_of[i]))
      CGAL_nef3_assertion_msg(0,"SNC_SM_io_parser::read: error in edge line");
  }
  if ( ln == 2 ) {
    read_loop(Loop_of[0]); read_loop(Loop_of[1]);
  }
  for(i=0; i<fn; i++) {
    if (!read_face(Face_of[i]))
      CGAL_nef3_assertion_msg(0,"SNC_SM_io_parser::read: error in face line");
  }
}

//-----------------------------------------------------------------------------
// VERBOSE OUTPUT:
// note that we output the index of the objects which is stored in them
// this is NOT the member index as produced by the forall loops
//-----------------------------------------------------------------------------

template <typename Refs_>
void SNC_SM_io_parser<Refs_>::debug_vertex(SVertex_handle v) const
{ 
  out << index(v) << "[" << mark(v) << "," << point(v) << "]" << endl; 
}

template <typename Refs_>
void SNC_SM_io_parser<Refs_>::debug_edge(SHalfedge_handle e) const
{ 
  out << index(e)
      << "(" << index(source(e)) << "," << index(target(e)) << ") "
      << index(twin(e)) << " " << index(face(e))
      << " ["<< mark(e) << "] " << circle(e) << endl;
}

template <typename Refs_>
void SNC_SM_io_parser<Refs_>::debug_loop(SHalfloop_handle l) const
{ 
  out << index(l) << " "
      << index(twin(l)) << " " << index(face(l))
      << " ["<< mark(l) << "] " << circle(l) << endl;
}


template <typename Refs_>
void SNC_SM_io_parser<Refs_>::debug() const
{ 
  out << "\nDEBUG Sphere Map\n";
  out << "svertices:  " << number_of_svertices()  << "\n";
  out << "shalfedges: " << number_of_shalfedges() << "\n";
  out << "shalfloops: " << number_of_shalfloops() << "\n";
  SVertex_iterator vit; 
  CGAL_nef3_forall_svertices(vit,*this) {
    if ( is_isolated(vit) ) continue;
    typename Base::SHalfedge_around_svertex_circulator
      hcirc = out_edges(vit), hend = hcirc;
    debug_vertex(vit);
    CGAL_For_all(hcirc,hend) { out << "  "; debug_edge(hcirc); }
  }
  if ( has_loop() ) 
  { debug_loop(shalfloop()); debug_loop(twin(shalfloop())); }
  out << endl;
}

template <typename Refs_>
void SNC_SM_io_parser<Refs_>::dump(Vertex_handle v, std::ostream& os)
{ SNC_SM_io_parser<Refs_> Out(os,v); Out.print();
  os << std::endl; }



CGAL_END_NAMESPACE
#endif //CGAL_SNC_SM_IO_PARSER_H

