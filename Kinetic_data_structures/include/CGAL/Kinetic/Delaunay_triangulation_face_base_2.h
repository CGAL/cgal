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

#ifndef CGAL_KINETIC_KINETIC_DELAUNAY_FACE_BASE_2_H
#define CGAL_KINETIC_KINETIC_DELAUNAY_FACE_BASE_2_H

#include <CGAL/Kinetic/basic.h>
#include <CGAL/Triangulation_face_base_2.h>

namespace CGAL { namespace Kinetic {
//! A class to track labels of edges of faces in a triangulation
template <class SimulationTraits, class Face_base= CGAL::Triangulation_face_base_2<typename SimulationTraits::Instantaneous_kernel> >
class Delaunay_triangulation_face_base_2: public Face_base
{
private:
  typedef typename Face_base::Triangulation_data_structure   TDS;
public:
  typedef TDS                            Triangulation_data_structure;
  typedef typename TDS::Face_handle     Face_handle;
  typedef typename TDS::Vertex_handle   Vertex_handle;
  typedef typename Face_base::Geom_traits Traits;

  typedef typename SimulationTraits::Simulator::Event_key Edge_label;
  Delaunay_triangulation_face_base_2(): Face_base() {
  }

  Delaunay_triangulation_face_base_2(Vertex_handle v0, Vertex_handle v1,
				     Vertex_handle v2): Face_base(v0, v1, v2) {
  }

  Delaunay_triangulation_face_base_2(Vertex_handle v0, Vertex_handle v1,
				     Vertex_handle v2,
				     Face_handle f0, Face_handle f1,
				     Face_handle f2): Face_base(v0,v1,v2, f0,f1,f2) {
  }

  //! Set the label for edge i
  void set_edge_label(int i, const Edge_label l) {
    CGAL_assertion(i>=0 && i<3);
    _labels[i]=l;
  }

  //! Get the label
  Edge_label get_edge_label(int i) const
  {
    CGAL_assertion(i>=0 && i<3);
    return _labels[i];
  }

  /*typedef Edge_label* Edge_label_iterator;
  typedef const Edge_label* Const_edge_label_iterator;

  Edge_label_iterator begin_edge_labels() {
    return &_labels[0];
  }
  Edge_label_iterator end_edge_labels() {
    return &_labels[0]+3;
  }
  Const_edge_label_iterator begin_edge_labels() const
  {
    return &_labels[0];
  }
  Const_edge_label_iterator end_edge_labels() const
  {
    return &_labels[0]+3;
    }*/

  template < typename TDS2 >
  struct Rebind_TDS
  {
    typedef typename Face_base::template Rebind_TDS<TDS2>::Other Fb2;
    typedef Delaunay_triangulation_face_base_2<SimulationTraits, Fb2>  Other;
  };
protected:
  Edge_label _labels[3];
};

template <class T, class Fb>
std::ostream &operator<<(std::ostream &out,
			 const Delaunay_triangulation_face_base_2<T, Fb> &f)
{
  out << static_cast<const Fb&>(f);
  out << " (" << f.get_edge_label(0) << ", " << f.get_edge_label(1) << ", " << f.get_edge_label(2) << ")";
  return out;
}


} } //namespace CGAL::Kinetic
#endif
