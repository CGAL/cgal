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
// file          : include/CGAL/Nef_3/SNC_SM_checker.h
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
// SNC_SM_checker.h                checking functions
// ============================================================================
#ifndef CGAL_SNC_SM_CHECKER_H
#define CGAL_SNC_SM_CHECKER_H

#include <CGAL/basic.h>
#include <CGAL/SNC_SM_const_decorator.h>
CGAL_BEGIN_NAMESPACE

/*{\Moptions outfile=SNC_SM_checker.man }*/
/*{\Manpage {SNC_SM_checker}{Decorator}{Plane map checking}{}}*/

/*{\Mdefinition An instance |\Mvar| of the data type |\Mname| is a
decorator to check the structure of a sphere map. It is generic with
respect to two template concepts.  |Decorator_| has to be a decorator
model of our |SNC_SM_const_decorator| concept.}*/

/*{\Mgeneralization Decorator}*/

template <typename Decorator_> 
class SNC_SM_checker : public Decorator_
{ typedef Decorator_ Base;
public:
/*{\Mtypes 3}*/
typedef Decorator_ Const_decorator;
/*{\Mtypemember equals |Decorator_|.}*/

#define USING(t) typedef typename Base::t t
USING(Vertex_handle);
USING(SVertex_handle);
USING(SHalfedge_handle);
USING(SVertex_const_iterator);
USING(SHalfedge_const_iterator);
USING(SHalfedge_around_svertex_const_circulator);
USING(SHalfedge_around_sface_const_circulator);
USING(Sphere_point);
USING(Sphere_segment);
USING(Sphere_direction);
#undef USING

/*{\Mcreation 3}*/
SNC_SM_checker(Vertex_handle v) : Base(v) {}
/*{\Mcreate constructs a plane map checker working on |v|.}*/

SNC_SM_checker(const Base& D) : Base(D) {}

/*{\Moperations 2 }*/
Sphere_direction direction(SHalfedge_const_handle e) const
{ return Sphere_direction(circle(e)); }

void check_order_preserving_embedding(Vertex_const_handle v) const;
/*{\Mop checks if the embedding of the targets of the edges in
the adjacency list |A(v)| is counter-clockwise order-preserving with 
respect to the order of the edges in |A(v)|.}*/

void check_order_preserving_embedding() const;
/*{\Mop checks if the embedding of all vertices of |P| is 
counter-clockwise order-preserving with respect to the adjacency
list ordering of all vertices.}*/

void check_is_triangulation() const;
/*{\Mop checks if |P| is a triangulation.}*/

}; // SNC_SM_checker<Decorator_>


template <typename Decorator_>
void SNC_SM_checker<Decorator_>::
check_order_preserving_embedding(SVertex_const_handle v) const
{
  std::ostrstream error_status;
  CGAL::set_pretty_mode ( error_status );
  SHalfedge_const_handle ef = first_out_edge(v) ,e=ef,en,enn;
  error_status << "check_order_preserving_embedding\n";
  error_status << "vertex " << PH(v) << endl;
  Sphere_point p = point(v);
  if ( e != SHalfedge_const_handle() ) {
    while ( true ) {
      en = cyclic_adj_succ(e);
      enn = cyclic_adj_succ(en);
      if (en == ef) break;
      error_status << "  -> " << point(target(e));
      error_status << " " << point(target(en)) << " ";
      error_status << " " << point(target(enn)) << endl;
      if ( !strictly_ordered_ccw_at(p,direction(e),direction(en),
				    direction(enn)) ||
           !strictly_ordered_ccw_at(p,direction(e),direction(en),
				    direction(ef)) ) {
        error_status << "ccw order violate!" << endl << '\0';
        CGAL_nef3_assertion_msg(0,error_status.str());
      }
      e = en;
    }
  }
  error_status.freeze(0);
}

template <typename Decorator_>
void SNC_SM_checker<Decorator_>::
check_order_preserving_embedding() const
{
  SVertex_const_iterator v;
  CGAL_nef3_forall_svertices_of(v,center_vertex())
    check_order_preserving_embedding(v);
}

template <typename Decorator_>
void SNC_SM_checker<Decorator_>::
check_is_triangulation() const
{
  check_integrity_and_topological_planarity(false);
  CGAL_nef3_assertion(number_of_connected_components() == 1);
  check_order_preserving_embedding();

  std::ostrstream error_status;
  CGAL::set_pretty_mode ( error_status );
  error_status << "check_is_triangulation\n";
  SHalfedge_const_iterator e;
  CGAL_nef3_forall_shalfedges_of(e,center_vertex()) {
    SHalfedge_around_sface_const_circulator hit(e), hend(hit); 
    int edges_in_face_cycle=0;
    CGAL_For_all(hit,hend) {
      error_status << PH(hit);
      ++edges_in_face_cycle;
    }
    CGAL_nef3_assertion_msg(edges_in_face_cycle==3,error_status.str());
  }
  error_status.freeze(0);
}


CGAL_END_NAMESPACE
#endif // CGAL_SNC_SM_CHECKER_H

