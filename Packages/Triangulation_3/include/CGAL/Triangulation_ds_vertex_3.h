// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Triangulation_ds_vertex_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

// vertex of a combinatorial triangulation of any dimension <=3

#ifndef CGAL_TRIANGULATION_DS_VERTEX_3_H
#define CGAL_TRIANGULATION_DS_VERTEX_3_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_ds_iterators_3.h>
#include <CGAL/Triangulation_short_names_3.h>

CGAL_BEGIN_NAMESPACE

template <class Vb, class Cb >
class  Triangulation_ds_vertex_3 
  : public Vb
{
  typedef Triangulation_ds_vertex_3<Vb,Cb> Vertex;
  typedef Triangulation_ds_cell_3<Vb,Cb> Cell;

public:

  Triangulation_ds_vertex_3()
    : Vb()
  { set_order_of_creation(); }
    
  Cell* cell() const
  {
    return (Cell *) (Vb::cell());
  }
    
  void set_cell(Cell* c)
  {
    Vb::set_cell(c);
  }

  bool is_valid(bool verbose = false, int level = 0) const;

  // used for symbolic perturbation in remove_vertex for Delaunay
  // undocumented
  int get_order_of_creation() const
  {
      return _order_of_creation;
  }

private:
  void set_order_of_creation()
  {
    static int nb=-1; 
    _order_of_creation = ++nb;
  }

  int _order_of_creation;
};

template < class VH>
class Vertex_tds_compare_order_of_creation {
public:
  bool operator()(VH u, VH v) const {
    return u->get_order_of_creation() < v->get_order_of_creation();
  }
};

template <class Vb, class Cb >
bool
Triangulation_ds_vertex_3<Vb,Cb>::is_valid
(bool verbose, int level) const
{
  bool result = Vb::is_valid(verbose,level) && cell()->has_vertex(this);
  if ( ! result ) {
    if ( verbose )
      std::cerr << "invalid vertex" << std::endl;
    CGAL_triangulation_assertion(false);
  }
  return result;
}

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_DS_VERTEX_3_H
