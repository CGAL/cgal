// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Nef_2/PM_io_parser.h
// package       : Nef_2 
// chapter       : Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: Input and Output of extended plane map
// ============================================================================

#ifndef CGAL_PM_IO_PARSER_H
#define CGAL_PM_IO_PARSER_H

#include <CGAL/Nef_2/PM_decorator.h>
#include <CGAL/Nef_2/Object_index.h>
#include <vector>

CGAL_BEGIN_NAMESPACE

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
  int vn,en,fn,i;
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
PM_io_parser(std::istream& is, Plane_map& H);
/*{\Mcreate creates an instance |\Mvar| of type |\Mname|
   to input |H| from |is|.}*/

PM_io_parser(std::ostream& os, const Plane_map& H);
/*{\Mcreate creates an instance |\Mvar| of type |\Mname|
to output |H| to |os|.}*/

PM_io_parser(std::ostream&, const PMDEC& D);

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


template <typename PMDEC>
PM_io_parser<PMDEC>::PM_io_parser
  (std::istream& iin, Plane_map& H) :
  Base(H), in(iin), out(std::cout), verbose(0), vn(0), en(0), fn(0)
{ clear(); }

template <typename PMDEC>
PM_io_parser<PMDEC>::PM_io_parser
  (std::ostream& iout, const Plane_map& H) :
  Base(const_cast<Plane_map&>(H)), in(std::cin), out(iout), 
  VI(vertices_begin(),vertices_end(),'v'),
  EI(halfedges_begin(),halfedges_end(),'e'),
  FI(faces_begin(),faces_end(),'f'),
  vn(number_of_vertices()), 
  en(number_of_halfedges()), 
  fn(number_of_faces())
{ verbose = (out.iword(CGAL::IO::mode) != CGAL::IO::ASCII &&
             out.iword(CGAL::IO::mode) != CGAL::IO::BINARY);
}

template <typename PMDEC>
PM_io_parser<PMDEC>::PM_io_parser
  (std::ostream& iout, const PMDEC& D) :
  Base(D), in(std::cin), out(iout), 
  VI(vertices_begin(),vertices_end(),'v'),
  EI(halfedges_begin(),halfedges_end(),'e'),
  FI(faces_begin(),faces_end(),'f'),
  vn(number_of_vertices()), 
  en(number_of_halfedges()), 
  fn(number_of_faces())
{ verbose = (out.iword(CGAL::IO::mode) != CGAL::IO::ASCII &&
             out.iword(CGAL::IO::mode) != CGAL::IO::BINARY);
}


//-----------------------------------------------------------------------------
// OUTPUT AND INPUT:
//-----------------------------------------------------------------------------
#ifdef __BORLANDC__
#define ISSPACENS std::
#else
#define ISSPACENS 
#endif

template <typename PMDEC>
bool PM_io_parser<PMDEC>::check_sep(const char* sep)
{
  char c; 
  do in.get(c); while (ISSPACENS isspace(c));
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
     (eo >= 0 || eo < en || epr >= 0 || epr < en || ene >= 0 || ene < en ||
      v >= 0 || v < vn || f >= 0 || f < fn ,
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
    CGAL_assertion_msg(ei >= 0 && ei < en, "wrong index in face cycle list.");
    f->store_fc(Halfedge_of[ei]);
  } in.clear();
  if (!check_sep(",")) { return false; }
  while (in >> vi) { 
    CGAL_assertion_msg(vi >= 0 && vi < vn, "wrong index in iso vertex list.");
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
  Vertex_iterator vit, vend = vertices_end();
  for(vit = vertices_begin(); vit!=vend; ++vit) print_vertex(vit);
  if (verbose) 
    out << "/* index { opposite, prev, next, vertex, face, mark } */" 
        << std::endl;
  Halfedge_iterator eit, eend = halfedges_end();
  for(eit = halfedges_begin(); eit!=eend; ++eit) print_hedge(eit);
  if (verbose) 
    out << "/* index { halfedge, fclist, ivlist, mark } */" 
        << std::endl;
  Face_iterator fit, fend = faces_end();
  for(fit = faces_begin(); fit!=fend; ++fit) print_face(fit);
  out.flush();
  if (verbose) debug();
}

template <typename PMDEC>
void PM_io_parser<PMDEC>::read() 
{
  if ( !check_sep("Plane_map_2") )  
    CGAL_assertion_msg(0,"PM_io_parser::read: no embedded_PM header.");
  if ( !(check_sep("vertices") && (in >> vn)) ) 
    CGAL_assertion_msg(0,"PM_io_parser::read: wrong node line.");
  if ( !(check_sep("halfedges") && (in >> en) && (en%2==0)) )
    CGAL_assertion_msg(0,"PM_io_parser::read: wrong edge line.");
  if ( !(check_sep("faces") && (in >> fn)) )
    CGAL_assertion_msg(0,"PM_io_parser::read: wrong face line.");

  Vertex_of.reserve(vn);
  Halfedge_of.reserve(en);
  Face_of.reserve(fn);
  for(i=0; i<vn; i++)  Vertex_of[i] =   new_vertex();
  for(i=0; i<en; i++) 
    if (i%2==0) Halfedge_of[i] = new_halfedge_pair_without_vertices();
    else Halfedge_of[i] = twin(Halfedge_of[i-1]);
  for(i=0; i<fn; i++)  Face_of[i] =     new_face();

  for(i=0; i<vn; i++) {
    if (!read_vertex(Vertex_of[i]))
      CGAL_assertion_msg(0,"PM_io_parser::read: error in node line");
  }
  for(i=0; i<en; i++) {
    if (!read_hedge(Halfedge_of[i]))
      CGAL_assertion_msg(0,"PM_io_parser::read: error in halfedge\
      line");
  }
  for(i=0; i<fn; i++) {
    if (!read_face(Face_of[i]))
      CGAL_assertion_msg(0,"PM_io_parser::read: error in face line");
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
      << "Vertices:  " << number_of_vertices() << "\n"
      << "Halfedges: " << number_of_halfedges() << "\n";
  Vertex_iterator vit,vend = vertices_end(); 
  for( vit = vertices_begin(); vit != vend; ++vit ) {
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



CGAL_END_NAMESPACE
#endif //CGAL_PM_IO_PARSER_H

