// Copyright (c) 2003-2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Steve Oudot, David Rey, Mariette Yvinec, Laurent Rineau, Andreas Fabri

#ifndef CGAL_SURFACE_MESH_CELL_BASE_3_H
#define CGAL_SURFACE_MESH_CELL_BASE_3_H

#include <CGAL/Complex_2_in_triangulation_cell_base_3.h>

#include <string>

namespace CGAL {
  
  template < class GT, class Cb=Triangulation_cell_base_3<GT> > 
  class Surface_mesh_cell_base_3 : 
    public Complex_2_in_triangulation_cell_base_3<GT, Cb> 
  { 
    typedef Complex_2_in_triangulation_cell_base_3<GT, Cb> Base;
  public:
    typedef Surface_mesh_cell_base_3 <GT, Cb> Self;

    template < class TDS3 >
    struct Rebind_TDS {
      typedef typename Cb::template Rebind_TDS<TDS3>::Other  Cb3;
      typedef Surface_mesh_cell_base_3 <GT, Cb3> Other;
    };
    
    
    typedef typename GT::Point_3 Point;
    
    typedef typename Base::Triangulation_data_structure Tds;
    typedef typename Tds::Vertex_handle Vertex_handle;
    typedef typename Tds::Cell_handle Cell_handle;
    
    
  private:
    
    // -- fields added to this class --
    Point tab_surface_center_facets [4];
    char c_visited; // the lower 4 bits for 'facet_visited', 
                    // the higher 4 bits for 'visits'
    
  public:
    
    // Constructors

    Surface_mesh_cell_base_3() 
      : Base(), c_visited(0) 
    {}
    
    Surface_mesh_cell_base_3 (Vertex_handle v0,
			      Vertex_handle v1,
			      Vertex_handle v2,
			      Vertex_handle v3) 
      : Base (v0, v1, v2, v3), c_visited(0)  
    {}
    
    Surface_mesh_cell_base_3 (Vertex_handle v0,
			      Vertex_handle v1,
			      Vertex_handle v2,
			      Vertex_handle v3,
			      Cell_handle n0,
			      Cell_handle n1,
			      Cell_handle n2,
			      Cell_handle n3) 
      : Base (v0, v1, v2, v3, n0, n1, n2, n3), c_visited(0)  
    {}
    
    
    
    // Access functions
    
    // Facets
    
    bool is_facet_visited (const int facet) const {
      CGAL_assertion (facet>=0 && facet <4);
      return c_visited & (1 << facet);;
    }

    bool visited (const int facet) const {
      CGAL_assertion (facet>=0 && facet <4);
      return c_visited & (1 << (facet+4));;
    }
    
    const Point& get_facet_surface_center(const int facet) const {
      CGAL_assertion (facet>=0 && facet <4);
      return(tab_surface_center_facets[facet]);
  }
    
    // Setting functions
    
    // Facets
    
    void set_facet_visited (const int facet) {
      CGAL_assertion (facet>=0 && facet <4);
      c_visited |= (1 << facet);
    }
    
    void set_visited (const int facet) {
      CGAL_assertion (facet>=0 && facet <4);
      c_visited |= (1 << (4+facet));
    }
    
    void reset_visited (const int facet) {
      CGAL_assertion (facet>=0 && facet <4);
      c_visited &= (15 & ~(1 << facet));
      c_visited &= ((15 | (15<<4)) & ~(1 << (4+facet)));
    }
    
    void set_facet_surface_center(const int facet, 
			  const Point& p) {
      CGAL_assertion (facet>=0 && facet <4);
      tab_surface_center_facets[facet]=p;
    }

#ifdef CGAL_MESH_3_IO_H
    static
    std::string io_signature()
    {
      return 
        Get_io_signature<Complex_2_in_triangulation_cell_base_3<GT, Cb> >()();
    }
#endif
  };  // end Surface_mesh_cell_base_3

template < class GT, class Cb >
std::istream&
operator>>(std::istream &is,
           Surface_mesh_cell_base_3<GT, Cb> &v)
{
  return is >> 
    static_cast<Complex_2_in_triangulation_cell_base_3<GT, Cb>&>(v);
}

template < class GT, class Cb >
std::ostream&
operator<<(std::ostream &os,
           const Surface_mesh_cell_base_3<GT, Cb> &v)
{
  return os <<
    static_cast<const Complex_2_in_triangulation_cell_base_3<GT, Cb>&>(v);
}

}  // namespace CGAL


#endif  // CGAL_SURFACE_MESH_CELL_BASE_3_H

