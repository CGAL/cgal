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
// file          : include/CGAL/Nef_3/SNC_io_parser.h
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
// SNC_io_parser.h         input and output of SNCs
// ============================================================================
#ifndef CGAL_SNC_IO_PARSER_H
#define CGAL_SNC_IO_PARSER_H

#include <CGAL/basic.h>
#include <CGAL/Nef_3/SNC_SM_decorator.h>
#include <CGAL/Nef_3/SNC_structure.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_2/Object_index.h>
#include <vector>

CGAL_BEGIN_NAMESPACE

template <typename SNC_structure_>
class SNC_io_parser : public SNC_decorator<SNC_structure_>
{ typedef SNC_structure_ SNC_structure;
  typedef CGAL::SNC_io_parser<SNC_structure_> Self;
  typedef CGAL::SNC_decorator<SNC_structure_> Base;
  typedef CGAL::SNC_SM_decorator<SNC_structure_> SM_decorator;
  std::istream& in; std::ostream& out;
  bool verbose;
  CGAL::Object_index<Vertex_iterator>   VI;
  CGAL::Object_index<Halfedge_iterator> EI;
  CGAL::Object_index<Halffacet_iterator>    FI;
  CGAL::Object_index<Volume_iterator>   CI;
  CGAL::Object_index<SHalfedge_iterator> SEI;
  CGAL::Object_index<SHalfloop_handle>   SLI;
  CGAL::Object_index<SFace_iterator>     SFI;
  std::vector<Vertex_iterator>   Vertex_of;
  std::vector<Halfedge_iterator> Edge_of;
  std::vector<Halffacet_iterator>    Halffacet_of;
  std::vector<Volume_iterator>   Volume_of;
  std::vector<SVertex_iterator>  SVertex_of; 
  std::vector<SHalfedge_iterator> SEdge_of;
  std::vector<SHalfloop_iterator> SLoop_of;
  std::vector<SFace_iterator>     SFace_of;
  long i,vn,en,fn,cn,sen,sln,sfn;

public:
  #define USING(t) typedef typename SNC_structure_::t t
  USING(Vertex_iterator); USING(Vertex_handle);
  USING(Halfedge_iterator); USING(Halfedge_handle);
  USING(Halffacet_iterator); USING(Halffacet_handle);
  USING(Volume_iterator); USING(Volume_handle);
  USING(SVertex_iterator); USING(SVertex_handle);
  USING(SHalfedge_iterator); USING(SHalfedge_handle);
  USING(SFace_iterator); USING(SFace_handle);
  USING(SHalfloop_iterator); USING(SHalfloop_handle);
  USING(Object_iterator); USING(Object_handle);
  USING(SObject_handle);
  USING(SFace_cycle_iterator);
  USING(Halffacet_cycle_iterator);
  USING(Shell_entry_iterator);
  USING(Point_3);
  USING(Plane_3);
  USING(Vector_3);
  USING(Sphere_point);
  USING(Sphere_segment);
  USING(Mark);
  #undef USING
  typedef void* GenPtr;

public:
  SNC_io_parser(std::istream& is, SNC_structure& W);
  SNC_io_parser(std::ostream& os, SNC_structure& W);

  std::string index(Vertex_iterator v) const 
  { return VI(v,verbose); }
  std::string index(Halfedge_iterator e) const 
  { return EI(e,verbose); }
  std::string index(Halffacet_iterator f) const 
  { return FI(f,verbose); }
  std::string index(Volume_iterator c) const 
  { return CI(c,verbose); }
  std::string index(SHalfedge_iterator e) const 
  { return SEI(e,verbose); }
  std::string index(SHalfloop_handle l) const 
  { return SLI(l,verbose); }
  std::string index(SFace_iterator f) const 
  { return SFI(f,verbose); }
  std::string index(Object_iterator o) const
  { if( o == 0 )
      return string("undef");
    Vertex_iterator v;
    Halfedge_iterator e;
    Halffacet_iterator f;
    Volume_iterator c;
    SHalfedge_iterator se;
    SHalfloop_handle sl;
    SFace_iterator sf;
    if( assign( v, *o))
      return index(v);
    else if( assign( e, *o))
      return index(e);
    else if( assign( f, *o))
      return index(f);
    else if( assign( c, *o))
      return index(c);
    else if( assign( se, *o))
      return index(se);
    else if( assign( sl, *o))
      return index(sl);
    else if( assign( sf, *o))
      return index(sf);
    return string("unknown object");
  }

  bool check_sep(char* sep) const;
  void print_vertex(Vertex_handle) const;
  void print_edge(Halfedge_handle) const;
  void print_facet(Halffacet_handle) const;
  void print_volume(Volume_handle) const;
  void print_sedge(SHalfedge_handle) const;
  void print_sloop(SHalfloop_handle) const;
  void print_sface(SFace_handle) const;
  void print() const;
  void print_local_graph(Vertex_handle) const;

  bool read_vertex(Vertex_handle) const;
  bool read_edge(Halfedge_handle) const;
  bool read_facet(Halffacet_handle) const;
  bool read_volume(Volume_handle) const;
  bool read_svertex(SVertex_handle) const;
  bool read_sedge(SHalfedge_handle) const;
  bool read_sloop(SHalfloop_handle) const;
  bool read_sface(SFace_handle) const;
  void read();

  static void dump(SNC_structure& W, std::ostream& os = std::cerr)
  { Self O(os,W); O.print(); }

};

template <typename EW>
SNC_io_parser<EW>::SNC_io_parser(std::istream& is, SNC_structure& W) : 
  Base(W), in(is), out(std::cout) 
{ CGAL_assertion(W.empty());
  verbose = false; }

template <typename EW>
SNC_io_parser<EW>::SNC_io_parser(std::ostream& os, SNC_structure& W) : 
  Base(W), in(std::cin), out(os),
  VI(W.vertices_begin(),W.vertices_end(),'V'),
  EI(W.halfedges_begin(),W.halfedges_end(),'E'),
  FI(W.halffacets_begin(),W.halffacets_end(),'F'),
  CI(W.volumes_begin(),W.volumes_end(),'C'),
  SEI(W.shalfedges_begin(),W.shalfedges_end(),'e'),
  SLI(W.shalfloops_begin(),W.shalfloops_end(),'l'),
  SFI(W.sfaces_begin(),W.sfaces_end(),'f'),
  vn(W.number_of_vertices()), 
  en(W.number_of_halfedges()), 
  fn(W.number_of_halffacets()),
  cn(W.number_of_volumes()),
  sen(W.number_of_shalfedges()),
  sln(W.number_of_shalfloops()),
  sfn(W.number_of_sfaces())
{ verbose = (out.iword(CGAL::IO::mode) != CGAL::IO::ASCII &&
             out.iword(CGAL::IO::mode) != CGAL::IO::BINARY); 
  VI[W.vertices_end()]=-2;
  EI[W.halfedges_end()]=-2;
  FI[W.halffacets_end()]=-2;
  CI[W.volumes_end()]=-2;
  SEI[W.shalfedges_end()]=-2;
  SLI[W.shalfloops_end()]=-2;
  SFI[W.sfaces_end()]=-2;
}


template <typename EW>
bool SNC_io_parser<EW>::check_sep(char* sep) const
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


template <typename EW>
void SNC_io_parser<EW>::print() const
{ 
  out << "Selective Nef Complex" << std::endl;
  out << "vertices   " << sncp()->number_of_vertices() << std::endl;
  out << "halfedges  " << sncp()->number_of_halfedges() << std::endl;
  out << "facets     " << sncp()->number_of_halffacets() << std::endl;
  out << "volumes    " << sncp()->number_of_volumes() << std::endl;
  out << "shalfedges " << sncp()->number_of_shalfedges() << std::endl;
  out << "shalfloops " << sncp()->number_of_shalfloops() << std::endl;
  out << "sfaces     " << sncp()->number_of_sfaces() << std::endl;

  if (verbose) 
    out << "/* Vertex: index { svs sve ses see sfs sfe sl,"
        << " mark, point } */\n";
  Vertex_iterator v;
  CGAL_nef3_forall_vertices(v,*sncp()) print_vertex(v);
  if (verbose) 
  out << "/* Edge: index { twin, source, isolated incident_object,"
      << " mark } */\n";
  Halfedge_iterator e;
  CGAL_nef3_forall_halfedges(e,*sncp()) print_edge(e);
  if (verbose) 
  out << "/* Facet: index { fclist, ivlist, mark, plane } */\n";
  Halffacet_iterator f;
  CGAL_nef3_forall_halffacets(f,*sncp()) print_facet(f);
  if (verbose) 
  out << "/* Volume: index { shlist, mark } */\n";
  Volume_iterator c;
  CGAL_nef3_forall_volumes(c,*sncp()) print_volume(c);
  if (verbose) 
  out << "/* SEdge: index { twin, sprev, snext, source, sface,"
      << " prev, next, facet } */\n";
  SHalfedge_iterator se;
  CGAL_nef3_forall_shalfedges(se,*sncp()) print_sedge(se);
  if (verbose) 
  out << "/* SLoop: index { twin, sface, facet } */" << std::endl;
  SHalfloop_iterator sl;
  CGAL_nef3_forall_shalfloops(sl,*sncp()) print_sloop(sl);
  if (verbose) 
  out << "/* SFace: index { fclist, ivlist, sloop, volume } */" << std::endl;
  SFace_iterator sf;
  CGAL_nef3_forall_sfaces(sf,*sncp()) print_sface(sf);
}

template <typename EW>
void SNC_io_parser<EW>::read()
{ 
  if ( !check_sep("Selective Nef Complex") )  
    CGAL_assertion_msg(0,"SNC_io_parser::read: no embedded_PM header.");
  if ( !(check_sep("vertices") && (in >> vn)) ) 
    CGAL_assertion_msg(0,"SNC_io_parser::read: wrong vertex line.");
  if ( !(check_sep("halfedges") && (in >> en) && (en%2==0)) )
    CGAL_assertion_msg(0,"SNC_io_parser::read: wrong edge line.");
  if ( !(check_sep("facets") && (in >> fn) && (fn%2==0)) )
    CGAL_assertion_msg(0,"SNC_io_parser::read: wrong facet line.");
  if ( !(check_sep("volumes") && (in >> cn)) )
    CGAL_assertion_msg(0,"SNC_io_parser::read: wrong volume line.");
  if ( !(check_sep("shalfedges") && (in >> sen)) )
    CGAL_assertion_msg(0,"SNC_io_parser::read: wrong sedge line.");
  if ( !(check_sep("shalfloops") && (in >> sln)) )
    CGAL_assertion_msg(0,"SNC_io_parser::read: wrong sloop line.");
  if ( !(check_sep("sfaces") && (in >> sfn)) )
    CGAL_assertion_msg(0,"SNC_io_parser::read: wrong sface line.");

  Vertex_of.reserve(vn);
  Edge_of.reserve(en);
  Halffacet_of.reserve(fn);
  Volume_of.reserve(cn);
  SEdge_of.reserve(sen);
  SLoop_of.reserve(sln);
  SFace_of.reserve(sfn);

  for(i=0; i<vn; i++)  Vertex_of[i] = sncp()->new_vertex_only();
  for(i=0; i<en; i++)  Edge_of[i] = sncp()->new_halfedge_only();
  for(i=0; i<fn; i++)  Halffacet_of[i] = sncp()->new_halffacet_only();
  for(i=0; i<cn; i++)  Volume_of[i] = sncp()->new_volume_only();
  for(i=0; i<sen; i++) SEdge_of[i] = sncp()->new_shalfedge_only();
  for(i=0; i<sln; i++) SLoop_of[i] = sncp()->new_shalfloop_only();
  for(i=0; i<sfn; i++) SFace_of[i] = sncp()->new_sface_only();

  for(i=0; i<vn; i++) {
    if (!read_vertex(Vertex_of[i]))
      CGAL_assertion_msg(0,"SNC_io_parser::read: error in node line");
  }
  for(i=0; i<en; i++) {
    if (!read_edge(Edge_of[i]))
      CGAL_assertion_msg(0,"SNC_io_parser::read: error in edge line");
  }
  for(i=0; i<fn; i++) {
    if (!read_facet(Halffacet_of[i]))
      CGAL_assertion_msg(0,"SNC_io_parser::read: error in facet line");
  }
  for(i=0; i<cn; i++) {
    if (!read_volume(Volume_of[i]))
      CGAL_assertion_msg(0,"SNC_io_parser::read: error in volume line");
  }
  for(i=0; i<sen; i++) {
    if (!read_sedge(SEdge_of[i]))
      CGAL_assertion_msg(0,"SNC_io_parser::read: error in sedge line");
  }
  for(i=0; i<sln; i++) {
    if (!read_sloop(SLoop_of[i]))
      CGAL_assertion_msg(0,"SNC_io_parser::read: error in sloop line");
  }
  for(i=0; i<sfn; i++) {
    if (!read_sface(SFace_of[i]))
      CGAL_assertion_msg(0,"SNC_io_parser::read: error in sface line");
  }
}


template <typename EW>
void SNC_io_parser<EW>::print_vertex(Vertex_handle v) const
{ // syntax: index { svs sve ses see sfs sfe sl, mark, point }
  out << index(v) << " { " 
      << index(v->svertices_begin()) << " "
      << index(v->svertices_last()) << " "
      << index(v->shalfedges_begin()) << " "
      << index(v->shalfedges_last()) << " "
      << index(v->sfaces_begin()) << " "
      << index(v->sfaces_last()) << " "
      << index(v->shalfloop()) << ", "
      << mark(v) << ", " << point(v) << " }" << std::endl;
}

template <typename EW>
bool SNC_io_parser<EW>::
read_vertex(Vertex_handle v) const
{ // syntax: index { svs sve ses see sfs sfe sl, mark, point }
  int n, svs, sve, ses, see, sfs, sfe, sl;  Mark m;  Point_3 p;
  if ( !(in >> n) ||
       !check_sep("{") ||
       !(in >> svs) || !(in >> sve) || 
       !(in >> ses) || !(in >> see) || 
       !(in >> sfs) || !(in >> sfe) || 
       !(in >> sl) ||
       !check_sep(",") ||
       !(in >> m) || 
       !check_sep(",") ||
       !(in >> p) || 
       !check_sep("}") ) return false;
  CGAL_assertion_msg(
      Vertex_of[n] == v &&
      svs >= 0 && svs < en && sve >= 0 && sve < en && 
      ses >= 0 && ses < sen && see >= 0 && see < sen && 
      sfs >= 0 && sfs < sfn && sfe >= 0 && sfe < sfn && 
      sl >= 0 && sl < sln , "wrong index in read_vertex");

  v->svertices_begin_ = Edge_of[svs]; v->svertices_last_ = Edge_of[sve];
  v->shalfedges_begin_ = SEdge_of[ses]; v->shalfedges_last_ = SEdge_of[see];
  v->sfaces_begin_ = SFace_of[sfs]; v->sfaces_last_ = SFace_of[sfe];
  v->shalfloop_ = SLoop_of[sl];
  mark(v) = m; point(v) = p;
  return true; 
}


template <typename EW>
void SNC_io_parser<EW>::print_edge(Halfedge_handle e) const
{ // syntax: index { twin, source, isolated incident_object, mark }
  SM_decorator D(vertex(e));
  out << index(e) << " { " << index(twin(e)) << ", "
      << index(source(e)) << ", ";
  if ( D.is_isolated(e) ) out << "1 " << index(D.face(e));
  else out << "0 " << index(D.first_out_edge(e));
  out << ", " << mark(e) << " }"<< std::endl;
}

template <typename EW>
bool SNC_io_parser<EW>::
read_edge(Halfedge_handle e) const
{ // syntax: index { twin, source, isolated incident_object, mark }
  int n, et, vs, ef, efm; bool iso; Mark m;
  if ( !(in >> n) ||
       !check_sep("{") ||
       !(in >> et) || !check_sep(",") ||
       !(in >> vs) || !check_sep(",") ||
       !(in >> iso) || !(in >> ef) || !check_sep(",") ||
       !(in >> m) || !check_sep("}") )
    return false;
  
  if (iso) efm=sfn; else efm=sen;
  CGAL_assertion_msg (
      Edge_of[n] == e &&
      et >= 0 && et < en && vs >= 0 && vs < vn && 
      ef >= 0 && et < efm , "wrong index in read_edge");
  
  e->twin_ = Edge_of[et];
  e->center_vertex_ = Vertex_of[vs];
  if ( iso ) e->incident_sface_ = SFace_of[ef];
  else       e->out_sedge_ = SEdge_of[ef];
  mark(e) = m;
  return true;
}

template <typename EW>
void SNC_io_parser<EW>::print_facet(Halffacet_handle f) const
{ // syntax: index { fclist, ivlist, mark, plane }
  out << index(f) << " { "; 
  Halffacet_cycle_iterator it; 
  CGAL_nef3_forall_facet_cycles_of(it,f)
    if ( it.is_shalfedge() ) out << index(SHalfedge_handle(it)) << ' ';
  out << ", ";
  CGAL_nef3_forall_facet_cycles_of(it,f)
    if ( it.is_shalfloop() ) out << index(SHalfloop_handle(it));
  out << ", " << mark(f) << ", " << plane(f) << " }\n";
}

template <typename EW>
bool SNC_io_parser<EW>::
read_facet(Halffacet_handle f) const
{ // syntax: index { fclist, ivlist, loop, mark }
  int n, ei, li; Mark m;
  if ( !(in >> n) || !check_sep("{") ) return false;
  while (in >> ei) { 
    CGAL_assertion_msg(ei >= 0 && ei < sen, 
      "wrong index in facet cycle list.");
    store_boundary_object(SEdge_of[ei],f);
  } in.clear();
  if (!check_sep(",")) { return false; }
  while (in >> li) { 
    CGAL_assertion_msg(li >= 0 && li < sln, 
      "wrong index in facet cycle list.");
    store_boundary_object(SLoop_of[li],f);
  } in.clear();
  if (!check_sep(",") || !(in >> m) || !check_sep("}") ) 
    return false;
  mark(f) = m;
  return true;
}

template <typename EW>
void SNC_io_parser<EW>::print_volume(Volume_handle c) const
{ // syntax: index { shlist, mark }
  out << index(c) << " { "; 
  Shell_entry_iterator it;
  CGAL_nef3_forall_shells_of(it,c)
    out << index(SFace_handle(it)) << ' ';
  out << ", " << mark(c) << " }\n";
}

template <typename EW>
bool SNC_io_parser<EW>::
read_volume(Volume_handle c) const
{ // syntax: index { shlist, mark }
  int n, fi; Mark m;
  if ( !(in >> n) || !check_sep("{") ) return false;
  while (in >> fi) { 
    CGAL_assertion_msg(fi >= 0 && fi < sfn, 
      "wrong index in shell list.");
    store_boundary_object(SFace_of[fi],c);
  } in.clear();
  if (!check_sep(",") || !(in >> m) || !check_sep("}") ) 
    return false;
  mark(c) = m;
  return true;
}


template <typename EW>
void SNC_io_parser<EW>::
print_sedge(SHalfedge_handle e) const
{ // syntax: index { twin, sprev, snext, source, sface, prev, next, facet }
  SM_decorator D(vertex(e));
  out << index(e) << " { "
      << index(D.twin(e)) << ", " 
      << index(D.previous(e)) << ", " << index(D.next(e)) << ", "
      << index(D.source(e)) << ", " << index(D.face(e)) << ", "
      << index(previous(e)) << ", " << index(next(e)) << ", "
      << index(facet(e)) << " }\n";
}

template <typename EW>
bool SNC_io_parser<EW>::
read_sedge(SHalfedge_handle e) const
{ // syntax: index { twin, sprev, snext, source, sface, prev, next, facet }
  int n, et, sp, sn, vs, sf, ep, en, ft; 
  if ( !(in >> n) ||
       !check_sep("{") ||
       !(in >> et) || !check_sep(",") ||
       !(in >> sp) || !check_sep(",") ||
       !(in >> sn) || !check_sep(",") ||
       !(in >> vs) || !check_sep(",") ||
       !(in >> sf) || !check_sep(",") ||
       !(in >> ep) || !check_sep(",") ||
       !(in >> en) || !check_sep(",") ||
       !(in >> ft) || !check_sep("}") )
    return false;
  CGAL_assertion_msg 
     (et >= 0 && et < sen && sp >= 0 && sp < sen && 
      sn >= 0 && sn < sen && vs >= 0 && vs < en && 
      sf >= 0 && sf < sfn && ep >= 0 && ep < sen &&
      en >= 0 && en < sen && ft >= 0 && ft < fn ,
      "wrong index in read_sedge");
  
  // precond: features exist!
  CGAL_assertion(SEdge_of[n]==e);
  e->sprev_ = SEdge_of[sp];
  e->snext_ = SEdge_of[sn];
  e->source_ = Edge_of[vs];
  e->incident_sface_ = SFace_of[sf];
  e->prev_ = SEdge_of[ep];
  e->next_ = SEdge_of[en];
  e->incident_facet_ = Halffacet_of[ft];
  return true;
}

template <typename EW>
void SNC_io_parser<EW>::
print_sloop(SHalfloop_handle l) const
{ // syntax: index { twin, sface, facet }
  SM_decorator D(vertex(l));
  out << index(l) << " { "
      << index(D.twin(l)) << ", " << index(D.face(l)) << ", "
      << index(facet(l)) << " }\n";
}

template <typename EW>
bool SNC_io_parser<EW>::
read_sloop(SHalfloop_handle l) const
{ // syntax: index { twin, sface, facet }
  int n, lt, sf, ft; 
  if ( !(in >> n) ||
       !check_sep("{") ||
       !(in >> lt) || !check_sep(",") ||
       !(in >> sf) || !check_sep(",") ||
       !(in >> ft) || !check_sep("}") )
    return false;
  CGAL_assertion_msg 
     (lt >= 0 && lt < sen && sf >= 0 && sf < sfn && 
      ft >= 0 && ft < fn ,
      "wrong index in read_sedge");
  
  CGAL_assertion(SLoop_of[n]==l);
  l->twin_ = SLoop_of[lt];
  l->incident_sface_ = SFace_of[sf];
  l->incident_facet_ = Halffacet_of[ft];
  return true;
}


template <typename EW>
void SNC_io_parser<EW>::
print_sface(SFace_handle f) const
{ // syntax: index { vertex, fclist, ivlist, sloop, volume }
  out << index(f) << " { " << index(f->center_vertex_) << ", "; 
  SFace_cycle_iterator it;
  CGAL_nef3_forall_sface_cycles_of(it,f)
    if ( it.is_shalfedge() ) out << index(SHalfedge_handle(it)) << ' ';
  out << ", ";
  CGAL_nef3_forall_sface_cycles_of(it,f)
    if ( it.is_svertex() ) out << index(SVertex_handle(it)) << ' ';
  out << ", ";
  CGAL_nef3_forall_sface_cycles_of(it,f)
    if ( it.is_shalfloop() ) out << index(SHalfloop_handle(it));
  out << ", " << index(volume(f)) << " }\n";
}


template <typename EW>
bool SNC_io_parser<EW>::
read_sface(SFace_handle f) const
{ // syntax: index { vertex, fclist, ivlist, sloop, volume }
  int n, vc, ei, vi, li, c;
  if ( !(in >> n) || !check_sep("{") ||
       !(in >> vc) || !check_sep(",") ) 
    return false;
  CGAL_assertion(vc >= 0 && vc < vn);
  f->center_vertex_ = Vertex_of[vc];
  SM_decorator D(Vertex_of[vc]);
  while (in >> ei) { 
    CGAL_assertion_msg(ei >= 0 && ei < sen, 
      "wrong index in sface cycle list.");
    D.store_boundary_object(SEdge_of[ei],f);
  } in.clear();
  while (in >> vi) { 
    CGAL_assertion_msg(vi >= 0 && vi < en, 
      "wrong index in sface cycle list.");
    D.store_boundary_object(Edge_of[vi],f);
  } in.clear();
  if (!check_sep(",")) { return false; }
  while (in >> li) { 
    CGAL_assertion_msg(li >= 0 && li < sln, 
      "wrong index in sface cycle list.");
    D.store_boundary_object(SLoop_of[li],f);
  } in.clear();
  if (!check_sep(",") || !(in >> c) || !check_sep("}") ) 
    return false;
  f->incident_volume_ = Volume_of[c];
  return true;
}


template <typename EW>
void SNC_io_parser<EW>::print_local_graph(Vertex_handle v) const
{ SM_decorator D(v);
  out << "Local Graph " 
      << D.number_of_vertices() << " " << D.number_of_edges() << " "
      << D.number_of_loops() << " " << D.number_of_faces() << " "
      << std::endl;
  if (verbose) 
    out << "/* index { twin, source, isolated incident_object, mark } */\n";
  SVertex_iterator vit;
  CGAL_nef3_forall_svertices_of(vit,v) print_edge(vit);
  if (verbose) 
    out << "/* index { twin, sprev, snext, source, sface,"
        << " prev, next, facet } */\n";
  SHalfedge_iterator eit;
  CGAL_nef3_forall_shalfedges_of(eit,v) print_sedge(eit);
  if (verbose) 
    out << "/* index { twin, sface, facet } */" << std::endl;
  if ( D.has_loop() ) 
  { print_sloop(D.loop()); print_sloop(D.twin(D.loop())); }
  if (verbose) 
    out << "/* index { fclist, ivlist, sloop, volume } */" << std::endl;
  SFace_iterator fit;
  CGAL_nef3_forall_sfaces_of(fit,v) print_sface(fit);
  out.flush();
}



CGAL_END_NAMESPACE
#endif //CGAL_SNC_IO_PARSER_H

