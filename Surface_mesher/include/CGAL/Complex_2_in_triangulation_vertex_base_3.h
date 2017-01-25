// Copyright (c) 2003-2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Steve Oudot, David Rey, Mariette Yvinec, Laurent Rineau, Andreas Fabri




#ifndef CGAL_COMPLEX_2_IN_TRIANGULATION_VERTEX_BASE_3_H
#define CGAL_COMPLEX_2_IN_TRIANGULATION_VERTEX_BASE_3_H

#include <CGAL/license/Surface_mesher.h>



#include <CGAL/Triangulation_vertex_base_3.h>

namespace CGAL {

  template < class GT, class Vb = Triangulation_vertex_base_3 <GT> > 
  class Complex_2_in_triangulation_vertex_base_3 : public Vb {    
    
  public:
    typedef Complex_2_in_triangulation_vertex_base_3 <GT, Vb> Self;

    template < class TDS3 >
    struct Rebind_TDS {
      typedef typename Vb::template Rebind_TDS<TDS3>::Other  Vb3;
      typedef Complex_2_in_triangulation_vertex_base_3 <GT, Vb3> Other;
    };

    typedef typename Vb::Triangulation_data_structure Tds;
    typedef typename Tds::Vertex_handle Vertex_handle;
    typedef typename Tds::Cell_handle Cell_handle;
    typedef typename Tds::Facet Facet;

  private:
    int number_of_incident_facets_;
    int number_of_components_; // number of components in the adjacency
                               // graph of incident facets (in complex)
    bool cache_validity;
  public:
    // Constructors

    Complex_2_in_triangulation_vertex_base_3()
      : Vb(), 
        number_of_incident_facets_(0),
        number_of_components_(0),
        cache_validity(false)
    {}

    bool is_c2t3_cache_valid() const {
      return cache_validity;
    }

    void invalidate_c2t3_cache()
    {
      cache_validity = false;
    }

    void set_c2t3_cache(const int i, const int j)
    {
      number_of_incident_facets_ = i;
      number_of_components_ = j;
      cache_validity = true;
    }

    int cached_number_of_incident_facets() const
    {
      return number_of_incident_facets_;
    }
    
    int cached_number_of_components() const
    {
      return number_of_components_;
    }

#ifdef CGAL_MESH_3_IO_H
    static
    std::string io_signature()
    {
      return Get_io_signature<Vb>()();
    }
#endif

  };  // end Complex_2_in_triangulation_vertex_base_3


}  // namespace CGAL


#endif  // CGAL_COMPLEX_2_IN_TRIANGULATION_CELL_BASE_3_H
