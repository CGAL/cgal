#ifndef CGAL_SM_IO_PARSER_H
#define CGAL_SM_IO_PARSER_H

#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_2/Object_index.h>
#include <vector>
#include <iostream>

CGAL_BEGIN_NAMESPACE

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

  typedef typename Decorator_::Vertex_iterator     Vertex_iterator;
  typedef typename Decorator_::Halfedge_iterator   Halfedge_iterator;
  typedef typename Decorator_::Face_iterator       Face_iterator;
  typedef typename Decorator_::Vertex_handle       Vertex_handle;
  typedef typename Decorator_::Halfedge_handle     Halfedge_handle;
  typedef typename Decorator_::Face_handle         Face_handle;
  typedef typename Decorator_::Halfloop_handle     Halfloop_handle;
  typedef typename Decorator_::Face_cycle_iterator Face_cycle_iterator;
  std::istream& in; std::ostream& out;
  bool verbose;
  // a reference to the IO object
  CGAL::Object_index<Vertex_handle>   VI;
  CGAL::Object_index<Halfedge_handle> EI;
  CGAL::Object_index<Face_handle>     FI;
  std::vector<Vertex_handle>        Vertex_of;
  std::vector<Halfedge_handle>      Edge_of;
  std::vector<Face_handle>          Face_of;
  Halfloop_handle                   Loop_of[2];
  // object mapping for input
  int vn,en,ln,fn,i;
  // the number of objects

  bool check_sep(char* sep);
  void print_vertex(Vertex_handle) const;
  void print_edge(Halfedge_handle) const;
  void print_loop(Halfloop_handle) const;
  void print_face(Face_handle) const;
  void print_init_points() const;

  bool read_vertex(Vertex_handle);
  bool read_edge(Halfedge_handle);
  bool read_loop(Halfloop_handle);
  bool read_face(Face_handle);
  bool read_init_points() const;

  void debug_vertex(Vertex_handle) const;
  void debug_edge(Halfedge_handle) const;
  void debug_loop(Halfloop_handle) const;

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

std::string index(Vertex_handle v) const 
{ return VI(v,verbose); }
std::string index(Halfedge_handle e) const 
{ return EI(e,verbose); }
std::string index(Halfloop_handle l) const 
{ if (verbose)  return (l==halfloop()? "l0" : "l1");
  else return (l==halfloop()? "0" : "1");
}
std::string index(Face_handle f) const 
{ return FI(f,verbose); }

static void dump(const Decorator_& D, std::ostream& os = cerr);
/*{\Mstatic prints the plane map decorated by |D| to |os|.}*/

}; // SM_io_parser<Decorator_>


template <typename Decorator_>
SM_io_parser<Decorator_>::
SM_io_parser(std::istream& iin, const Decorator_& H) :
  Base(H), in(iin), out(std::cout), verbose(0), 
  vn(0), en(0), fn(0), ln(0)
{ clear(); }

template <typename Decorator_>
SM_io_parser<Decorator_>::
SM_io_parser(std::ostream& iout, const Decorator_& D) 
  : Base(D), in(std::cin), out(iout), 
  VI(vertices_begin(),vertices_end(),'v'),
  EI(halfedges_begin(),halfedges_end(),'e'),
  FI(faces_begin(),faces_end(),'f'),
  vn(number_of_vertices()), 
  en(number_of_halfedges()), 
  ln(number_of_halfloops()),
  fn(number_of_faces())
{ verbose = (out.iword(CGAL::IO::mode) != CGAL::IO::ASCII &&
             out.iword(CGAL::IO::mode) != CGAL::IO::BINARY);
}


//-----------------------------------------------------------------------------
// OUTPUT AND INPUT:
//-----------------------------------------------------------------------------

template <typename Decorator_>
bool SM_io_parser<Decorator_>::check_sep(char* sep)
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
void SM_io_parser<Decorator_>::print_vertex(Vertex_handle v) const
{
  // syntax: index { isolated incident_object, mark, point }
  out << index(v) << " { ";
  if ( is_isolated(v) ) out << "1 " << index(face(v));
  else                  out << "0 " << index(first_out_edge(v));
  out  << ", " << mark(v) << ", " << point(v) <<  "}\n";
}

template <typename Decorator_>
bool SM_io_parser<Decorator_>::read_vertex(Vertex_handle v)
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

template <typename Decorator_>
void SM_io_parser<Decorator_>::print_edge(Halfedge_handle e) const
{ // syntax: index { twin, prev, next, source, face, mark, circle }
  out << index(e) << " { "
      << index(twin(e)) << ", " 
      << index(previous(e)) << ", " << index(next(e)) << ", "
      << index(source(e)) << ", " << index(face(e)) << ", "
      << mark(e) << ", " << circle(e) << " }\n";
}

template <typename Decorator_>
bool SM_io_parser<Decorator_>::read_edge(Halfedge_handle e)
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
  CGAL_nef_assertion_msg 
     (eo >= 0 && eo < en && epr >= 0 && epr < en && ene >= 0 && ene < en &&
      v >= 0 && v < vn && f >= 0 && f < fn ,
      "wrong index in read_edge");
  
  // precond: features exist!
  CGAL_nef_assertion(EI[twin(e)]);
  set_prev(e,Edge_of[epr]);
  set_next(e,Edge_of[ene]);
  set_source(e,Vertex_of[v]);
  set_face(e,Face_of[f]);
  mark(e) = m;
  circle(e) = k;
  return true;
}

template <typename Decorator_>
void SM_io_parser<Decorator_>::print_loop(Halfloop_handle l) const
{ // syntax: index { twin, face, mark, circle }
  out << index(l) << " { "
      << index(twin(l)) << ", " 
      << index(face(l)) << ", "
      << mark(l) << ", " << circle(l) << " }\n";
}

template <typename Decorator_>
bool SM_io_parser<Decorator_>::read_loop(Halfloop_handle l)
{ // syntax: index { twin, face, mark, circle }
  int n, lo, f; bool m; Sphere_circle k;
  if ( !(in >> n) ||
       !check_sep("{") ||
       !(in >> lo) || !check_sep(",") ||
       !(in >> f) || !check_sep(",") ||
       !(in >> m) || !check_sep(",") ||
       !(in >> k) || !check_sep("}") )
    return false;
  CGAL_nef_assertion_msg(
    (lo >= 0 && lo < 2 && f >= 0 && f < fn),"wrong index in read_edge");
  
  set_face(l,Face_of[f]);
  mark(l) = m;
  circle(l) = k;
  return true;
}


template <typename Decorator_>
void SM_io_parser<Decorator_>::print_face(Face_handle f) const
{ // syntax: index { fclist, ivlist, loop, mark }
  out << index(f) << " { "; 
  Face_cycle_iterator it;
  CGAL_forall_face_cycles_of(it,f)
    if ( it.is_halfedge() ) out << index(Halfedge_handle(it)) << ' ';
  out << ", ";
  CGAL_forall_face_cycles_of(it,f)
    if ( it.is_vertex() ) out << index(Vertex_handle(it)) << ' ';
  out << ", ";
  CGAL_forall_face_cycles_of(it,f)
    if ( it.is_halfloop() ) out << index(Halfloop_handle(it));
  out << ", " << mark(f) << " }\n";
}

template <typename Decorator_>
void SM_io_parser<Decorator_>::print_init_points() const
{ 
  out << mark_of_halfsphere(-1) << " "
      << mark_of_halfsphere(+1) << "\n"; 
}



template <typename Decorator_>
bool SM_io_parser<Decorator_>::read_face(Face_handle f)
{ // syntax: index { fclist, ivlist, loop, mark }
  int n, ei, vi, li; Mark m;
  if ( !(in >> n) || !check_sep("{") ) return false;
  while (in >> ei) { 
    CGAL_nef_assertion_msg(ei >= 0 && ei < en, 
                           "wrong index in face cycle list.");
    store_boundary_object(Edge_of[ei],f);
  } in.clear();
  if (!check_sep(",")) { return false; }
  while (in >> vi) { 
    CGAL_nef_assertion_msg(vi >= 0 && vi < vn, 
                           "wrong index in iso vertex list.");
    store_boundary_object(Vertex_of[vi],f);
  } in.clear();
  if (!check_sep(",")) { return false; }
  while (in >> li) { 
    CGAL_nef_assertion_msg(li >= 0 && li < 2, 
                           "wrong index in iso vertex list.");
    store_boundary_object(Loop_of[li],f);
  } in.clear();
  if (!check_sep(",") || !(in >> m) || !check_sep("}") ) 
    return false;
  mark(f) = m;
  return true;
}

template <typename Decorator_>
bool SM_io_parser<Decorator_>::read_init_points() const
{ Mark m_neg, m_pos;
  if ( ! (in >> m_neg >> m_pos) ) return false;
  mark_of_halfsphere(-1) = m_neg;
  mark_of_halfsphere(+1) = m_pos;
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
  Vertex_iterator vit;
  CGAL_forall_vertices(vit,*this) print_vertex(vit);
  if (verbose) 
    out << "/* index { twin, prev, next, source, face, mark, circle } */" 
	<< std::endl;
  Halfedge_iterator eit;
  CGAL_forall_halfedges(eit,*this) print_edge(eit);
  if (verbose) 
    out << "/* index { twin, face, mark, circle } */" << std::endl;
  if ( has_loop() ) 
  { print_loop(halfloop()); print_loop(twin(halfloop())); }
  if (verbose) 
    out << "/* index { fclist, ivlist, loop, mark } */" << std::endl;
  Face_iterator fit;
  CGAL_forall_faces(fit,*this) print_face(fit);
  if (verbose) 
    out << "/* mark at y-/y+ */" << std::endl;
  print_init_points();  
  out.flush();
  if (verbose) debug();
}

template <typename Decorator_>
void SM_io_parser<Decorator_>::read() 
{
  if ( !check_sep("Plane_map_2") )  
    CGAL_nef_assertion_msg(0,"SM_io_parser::read: no embedded_PM header.");
  if ( !(check_sep("vertices") && (in >> vn)) ) 
    CGAL_nef_assertion_msg(0,"SM_io_parser::read: wrong vertex line.");
  if ( !(check_sep("edges") && (in >> en) && (en%2==0)) )
    CGAL_nef_assertion_msg(0,"SM_io_parser::read: wrong edge line.");
  if ( !(check_sep("loops") && (in >> ln)) )
    CGAL_nef_assertion_msg(0,"SM_io_parser::read: wrong loop line.");
  if ( !(check_sep("faces") && (in >> fn)) )
    CGAL_nef_assertion_msg(0,"SM_io_parser::read: wrong face line.");

  Vertex_of.reserve(vn);
  Edge_of.reserve(en);
  Face_of.reserve(fn);
  for(i=0; i<vn; i++)  Vertex_of[i] =   new_vertex();
  for(i=0; i<en; i++) 
    if (i%2==0) Edge_of[i] = new_edge_pair_without_vertices();
    else Edge_of[i] = twin(Edge_of[i-1]);
  for(i=0; i<fn; i++)  Face_of[i] =     new_face();
  if ( ln == 2 ) 
  { Loop_of[0] = new_loop(); Loop_of[1] = twin(loop()); }

  for(i=0; i<vn; i++) {
    if (!read_vertex(Vertex_of[i]))
      CGAL_nef_assertion_msg(0,"SM_io_parser::read: error in node line");
  }
  for(i=0; i<en; i++) {
    if (!read_edge(Edge_of[i]))
      CGAL_nef_assertion_msg(0,"SM_io_parser::read: error in edge line");
  }
  if ( ln == 2 ) {
    read_loop(Loop_of[0]); read_loop(Loop_of[1]);
  }
  for(i=0; i<fn; i++) {
    if (!read_face(Face_of[i]))
      CGAL_nef_assertion_msg(0,"SM_io_parser::read: error in face line");
  }
  if (!read_init_points())
    CGAL_nef_assertion_msg(0,"SM_io_parser::read: error in init point line");
}

//-----------------------------------------------------------------------------
// VERBOSE OUTPUT:
// note that we output the index of the objects which is stored in them
// this is NOT the member index as produced by the forall loops
//-----------------------------------------------------------------------------

template <typename Decorator_>
void SM_io_parser<Decorator_>::debug_vertex(Vertex_handle v) const
{ 
  out << index(v) << "[" << mark(v) << "," << point(v) << "]" << std::endl; 
}

template <typename Decorator_>
void SM_io_parser<Decorator_>::debug_edge(Halfedge_handle e) const
{ 
  out << index(e)
      << "(" << index(source(e)) << "," << index(target(e)) << ") "
      << index(twin(e)) << " " << index(face(e))
      << " ["<< mark(e) << "," << circle(e) << "] " << std::endl;
}

template <typename Decorator_>
void SM_io_parser<Decorator_>::debug_loop(Halfloop_handle l) const
{ 
  out << index(l) << " "
      << index(twin(l)) << " " << index(face(l))
      << " ["<< mark(l) << "] " << circle(l) << std::endl;
}


template <typename Decorator_>
void SM_io_parser<Decorator_>::debug() const
{ 
  out << "\nDEBUG Plane_map\n";
  out << "Vertices:  " << number_of_vertices() << "\n";
  out << "Halfedges: " << number_of_halfedges() << "\n";
  out << "Loop:      " << number_of_halfloops() << "\n";
  Vertex_iterator vit; 
  CGAL_forall_vertices(vit,*this) {
    if ( is_isolated(vit) ) continue;
    typename Base::Halfedge_around_vertex_circulator
      hcirc = out_edges(vit), hend = hcirc;
    debug_vertex(vit);
    CGAL_For_all(hcirc,hend) { out << "  "; debug_edge(hcirc); }
  }
  if ( has_loop() ) 
  { debug_loop(halfloop()); debug_loop(twin(halfloop())); }
  out << std::endl;
}

template <typename Decorator_>
void SM_io_parser<Decorator_>::print_faces() const
{ 
  out << "\nFACES\n";
  out << "Vertices:  " << number_of_vertices() << "\n";
  out << "Halfedges: " << number_of_halfedges() << "\n";
  out << "Loop:      " << number_of_halfloops() << "\n";
  Halfedge_iterator e;
  Unique_hash_map<Halfedge_iterator,bool> Done(false);
  CGAL_forall_halfedges(e,*this) {
    if ( Done[e] ) continue;
    typename Base::Halfedge_around_face_circulator c(e), ce = c;
    out << "face cycle\n";
    CGAL_For_all(c,ce) 
    { Done[c]=true; out << "  "; debug_vertex(source(c)); }
  }
  if ( has_loop() ) 
  { debug_loop(halfloop()); debug_loop(twin(halfloop())); }
  out << std::endl;
}

template <typename Decorator_>
void SM_io_parser<Decorator_>::dump(const Decorator_& D, std::ostream& os)
{ SM_io_parser<Decorator_> Out(os,D);
  Out.print();
  Out.print_faces();
}



CGAL_END_NAMESPACE
#endif //CGAL_SM_IO_PARSER_H

