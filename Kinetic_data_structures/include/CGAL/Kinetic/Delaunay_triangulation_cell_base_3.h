// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_KINETIC_DELAUNAY_CELL_BASE_3_H
#define CGAL_KINETIC_KINETIC_DELAUNAY_CELL_BASE_3_H

#include <CGAL/Kinetic/basic.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

namespace CGAL { namespace Kinetic {
//! A class to track labels of edges of faces in a triangulation
template <class SimulationTraits, class Cell_base= CGAL::Triangulation_cell_base_3<typename SimulationTraits::Instantaneous_kernel> >
class Delaunay_triangulation_cell_base_3: public Cell_base
{
private:
  typedef typename Cell_base::Triangulation_data_structure   TDS;
public:
  typedef TDS                            Triangulation_data_structure;
  typedef typename TDS::Cell_handle     Cell_handle;
  typedef typename TDS::Vertex_handle   Vertex_handle;
  typedef typename Cell_base::Geom_traits Traits;

  typedef typename SimulationTraits::Simulator::Event_key Edge_label;
  typedef Edge_label Facet_label;
  Delaunay_triangulation_cell_base_3(): Cell_base() {
  }

  Delaunay_triangulation_cell_base_3(Vertex_handle v0, Vertex_handle v1,
				     Vertex_handle v2, Vertex_handle v3): Cell_base(v0, v1, v2, v3) {
  }

  Delaunay_triangulation_cell_base_3(Vertex_handle v0, Vertex_handle v1,
				     Vertex_handle v2, Vertex_handle v3,
				     Cell_handle f0, Cell_handle f1,
				     Cell_handle f2, Cell_handle f3): Cell_base(v0,v1,v2, v3, f0,f1,f2, f3) {
  }

  template < typename TDS3 >
  struct Rebind_TDS
  {
    typedef typename Cell_base::template Rebind_TDS<TDS3>::Other Cb3;
    typedef Delaunay_triangulation_cell_base_3<SimulationTraits, Cb3>  Other;
  };

  //! Set the label for edge i
  void set_edge_label(unsigned int i, unsigned int j, const Edge_label l) {
    elabels_[edge_index(i,j)]=l;
  }

  //! Get the label
  Edge_label edge_label(unsigned int i, unsigned int j) const
  {
    return elabels_[edge_index(i,j)];
  }

  //! Set the label for edge i
  void set_facet_label(unsigned int i, const Facet_label l) {
    CGAL_assertion(i<4);
    flabels_[i]=l;
  }

  //! Get the label
  Facet_label facet_label(unsigned int i) const
  {
    CGAL_assertion(i<4);
    return flabels_[i];
  }

protected:

  int edge_index(unsigned int i, unsigned int j) const
  {
    CGAL_precondition(i<4);
    CGAL_precondition(j<4);
    CGAL_precondition(i!= j);
    //if (i <j) std::swap(i,j);
    unsigned int sum= i+j;
    if ((std::max)(i,j) != 3) --sum;
    CGAL_assertion( sum <6);
    return sum;
  }

  Edge_label elabels_[6];
  Facet_label flabels_[4];
};

template <class Tr, class Fb>
std::ostream &operator<<(std::ostream &out,
			 const Delaunay_triangulation_cell_base_3<Tr, Fb> &f)
{
  out << static_cast<const Fb&>(f);
  //out << " [" << f.edge_label(1,0) << ", " << f.edge_label(2,0) << ", " << f.edge_label(2,1) << ", "
  //   << f.edge_label(3,0) << ", " << f.edge_label(3,1) << ", " << f.edge_label(3,2)<< "]";
  return out;
}


} } //namespace CGAL::Kinetic
#endif
