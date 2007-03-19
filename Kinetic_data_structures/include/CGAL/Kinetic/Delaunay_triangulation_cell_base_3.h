// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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

CGAL_KINETIC_BEGIN_NAMESPACE
//! A class to track labels of edges of faces in a triangulation
template <class SimulationTraits, 
  class Cell_base= CGAL::Triangulation_cell_base_3<typename SimulationTraits::Instantaneous_kernel> >
class Delaunay_triangulation_cell_base_3: public Cell_base
{
private:
  typedef typename Cell_base::Triangulation_data_structure   TDS;
public:
  typedef TDS                            Triangulation_data_structure;
  typedef typename TDS::Cell_handle     Cell_handle;
  typedef typename TDS::Vertex_handle   Vertex_handle;
  typedef typename Cell_base::Geom_traits Traits;

  void clear_elabels() {
    for (unsigned int i=0; i< 4; ++i) {
      elabels_[i]=false;
    }
  }

  typedef typename SimulationTraits::Simulator::Event_key Edge_label;
  typedef Edge_label Facet_label;
  Delaunay_triangulation_cell_base_3(): Cell_base() {
    clear_elabels();
  }

  Delaunay_triangulation_cell_base_3(Vertex_handle v0, Vertex_handle v1,
				     Vertex_handle v2, Vertex_handle v3): Cell_base(v0, v1, v2, v3) {
    clear_elabels();
  }

  Delaunay_triangulation_cell_base_3(Vertex_handle v0, Vertex_handle v1,
				     Vertex_handle v2, Vertex_handle v3,
				     Cell_handle f0, Cell_handle f1,
				     Cell_handle f2, Cell_handle f3): Cell_base(v0,v1,v2, v3, f0,f1,f2, f3) {
    clear_elabels();
  }

  template < typename TDS3 >
  struct Rebind_TDS
  {
    typedef typename Cell_base::template Rebind_TDS<TDS3>::Other Cb3;
    typedef Delaunay_triangulation_cell_base_3<SimulationTraits, Cb3>  Other;
  };

  //! Set the label for edge i
  void set_edge_certificate(unsigned int i, unsigned int j, const Edge_label l) {
    std::pair<int,int> lp= both_edge_indices(i,j);
    CGAL_precondition(flabels_[lp.first]== Facet_label() || elabels_[lp.first]);
    elabels_[lp.first]=true;
    flabels_[lp.first]=l;
    CGAL_precondition(flabels_[lp.second]== Facet_label() || elabels_[lp.second]);
    elabels_[lp.second]=true;
    flabels_[lp.second]=l;
  }

  bool has_edge_certificate(unsigned int i, unsigned int j) const {
    std::pair<int,int> l= both_edge_indices(i,j);
    return elabels_[l.first] && elabels_.[l.second];
  }


  //! Get the label
  Edge_label edge_label(unsigned int i, unsigned int j) const
  {
    CGAL_precondition(has_edge_label(i,j));
    CGAL_precondition(elabels_[first_edge_index(i,j)]);
    CGAL_precondition(flabels_[first_edge_index(i,j) ]
		      == flabels_[second_edge_index(i,j)]);
    return flabels_[first_edge_index(i,j)];
  }

  //! Set the label for edge i
  void set_facet_label(unsigned int i, const Facet_label l) {
    CGAL_assertion(elabels_[i]== false || flabels_[i] == Edge_label());
    CGAL_assertion(i<4);
    flabels_[i]=l;
    elabels_[i]=false;
  }

  //! Get the label
  Facet_label facet_label(unsigned int i) const
  {
    CGAL_assertion(elabels_[i]== false || flabels_[i] == Edge_label());
    CGAL_assertion(i<4);
    return flabels_[i];
  }

protected:

  /*
    need a hash which takes any tripple of faces to a 

    1: nope
    2: 0,1
    3: 0,2
    4: 0,3
    5: prime
    6: 1,2
    7: prime
    8: 1,3
    9: square
    10: nope
    11: prime
    12: 2,3
   */

  int first_edge_index(unsigned int i, unsigned int j) const {
    const static int lu0[]={-1,-1,0,0,0,-1,1,-1,1,-1,-1,-1,2};
    int r= lu0[(i+1)*(j+1)];
    CGAL_postcondition(r != -1);
    return r;
  }

  int second_edge_index(unsigned int i, unsigned int j) const {
    const static int lu1[]={-1,-1,1,2,3,-1,2,-1,3,-1,-1,-1,3};
    int r= lu1[(i+1)*(j+1)];
    CGAL_postcondition(r != -1);
    return r;
  }
  std::pair<int,int> both_edge_indices(unsigned int i, unsigned int j) const {
    CGAL_precondition(i != j);
    int a= first_edge_index(i,j);
    int b= second_edge_index(i,j);
    CGAL_postcondition(a != -1);
    CGAL_postcondition(b != -1);
    CGAL_postcondition(a != b);
    return std::make_pair(a,b);
  }

  //Edge_label elabels_[6];
  Facet_label flabels_[4];
  bool elabels_[4];
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


CGAL_KINETIC_END_NAMESPACE
#endif
