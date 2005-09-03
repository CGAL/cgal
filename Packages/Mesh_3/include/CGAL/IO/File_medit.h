// Copyright (c) 2004-2005  INRIA Sophia-Antipolis (France).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_IO_FILE_MEDIT_H
#define CGAL_IO_FILE_MEDIT_H

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/Point_traits.h>

namespace CGAL {

template <class C2T3>
void
output_pslg_to_medit (std::ostream& os, const C2T3& c2t3)
{
  typedef typename C2T3::Triangulation_3 Tr;
  typedef typename Tr::Finite_cells_iterator Finite_cells_iterator;
  typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Point Point;
  typedef Point_traits<Point> P_traits;
  typedef typename P_traits::Bare_point Bare_point;

  const Tr& tr = c2t3.triangulation();

  // Header.

  os << "MeshVersionFormatted 1" << std::endl
     << "Dimension 3" << std::endl;

  // Vertices
  
  os << "Vertices" << std::endl
     << tr.number_of_vertices() << std::endl;

  os << std::setprecision(20);
 
  std::map<Vertex_handle, int> V;
  int inum = 1;
  for( Finite_vertices_iterator vit = tr.finite_vertices_begin();
       vit != tr.finite_vertices_end();
       ++vit) 
  {
    V[vit] = inum++;
    Point p = static_cast<Point>(vit->point());
    os << p.x() << " " << p.y() << " " << p.z() << " " 
       << p.surface_index() << std::endl; // reference
  }

  // Facets
  os << "Triangles" << std::endl
     << number_of_facets_on_surface(tr) << std::endl;
  for( Finite_facets_iterator fit = tr.finite_facets_begin(); 
       fit != tr.finite_facets_end(); ++fit)
  {
    int surface_index = 0;
    if (c2t3.complex_subface_type(fit->first,fit->second)
        != C2T3::NOT_IN_COMPLEX)
    {
      for (int i=0; i<4; i++)
        if (i != (*fit).second)
        {
          const Vertex_handle& vh = (*fit).first->vertex(i);
          surface_index = vh->point().surface_index();
          os << V[vh] << " ";
        }
      
      os << surface_index << std::endl; //ref
    }
  }

  // Tetrahedra
  os << "Tetrahedra" << std::endl
     << number_of_cells_in_domain(tr) << std::endl;
  for( Finite_cells_iterator cit = tr.finite_cells_begin(); 
       cit != tr.finite_cells_end(); ++cit)
    if( cit->is_in_domain() )
    {
      for (int i=0; i<4; i++)
        os << V[cit->vertex(i)] << " ";
      os << "1" << std::endl; //ref
    }
  
  // End
  os << "End" << std::endl;

} // end output_pslg_to_medit

namespace Mesh_3 { namespace details {
  
  class Debug {
    bool active;
    std::ostream* out;
    const std::string header;
  public:
    Debug(bool debug,
          std::ostream* debug_str = &std::cout,
          const std::string header_string = "")
      : active(debug), out(debug_str), header(header_string)
    {
    }

    template <class T>
    Debug& apply(T& t)
    {
      if(active)
        *out << header << t;
      return *this;
    }

    operator bool() const
    {
      return false;
    }
  }; // end internal class Debug
    
  template <typename T>
  Debug& operator<<(Debug& d, const T& t)
  {
    d.apply(t);
    return d;
  }
} // end namespace detailes
} // end namespace Mesh_3

template < class C2T3>
bool
input_pslg_from_medit (std::istream& is, 
                       C2T3 & c2t3,
                       bool debug = false, 
                       std::ostream* debug_str = &std::cout) {
  typedef typename C2T3::Triangulation_3 Tr;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Point Point;
  typedef Point_traits<Point> P_traits;
  typedef typename P_traits::Bare_point Bare_point;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Vertex_handle Vertex_handle;

  
  Mesh_3::details::Debug debug_stream(debug,
                                      debug_str,
                                      "CGAL::input_pslg_from_medit()"
                                      "input error:");

  Tr& tr = c2t3.triangulation();

  tr.clear();

  // Header.
  std::string temp_string;
  int temp_int;

  is >> temp_string;
  if( !is || temp_string != "MeshVersionFormatted")
    return debug_stream << "  not a medit file!\n";

  is >> temp_int; // version
  if( !is || temp_int!=1 ) 
    return debug_stream << "  wrong medit version number: " << temp_int
                        << "\n"
                        << "  should be 1!\n";

  is >> temp_string;
  if( !is || temp_string != "Dimension" )
    return debug_stream << "  \"Dimension\" expected!\n";

  is >> temp_int; // dimension
  if( !is || temp_int != 3 )
    return debug_stream << "  dimension 3 expected!\n";

  // Vertices
  is >> temp_string;
  if( !is || temp_string != "Vertices" )
    return debug_stream << "  \"Vertices\" expected!\n";
  
  int number_of_vertices;
  if( ! (is >> number_of_vertices) )
    return debug_stream << "  number of vertices expected!\n";
  
  Cell_handle c_hint = Cell_handle();
  std::vector<Vertex_handle> V(number_of_vertices);
  for(int i = 0; i < number_of_vertices; ++i)
  {
    P_traits point_convert;

    Bare_point bp;
    is >> bp;
    if( !is )
      return debug_stream << "  error in Bare_point.operator>>()\n";

    int index;
    is >> index;
    if( !is )
      return debug_stream << "  index expected!\n";

    Point p = point_convert.point(bp);
    p.set_surface_index(index);

    Vertex_handle v = tr.insert(p, c_hint);
    V[i] = v;
    c_hint = v->cell();
  }

  // Facets
  is >> temp_string;
  if( !is || temp_string != "Triangles" )
    return debug_stream << "  \"Triangles\" expected!\n";

  int number_of_facets_on_surface;
  is >> number_of_facets_on_surface;
  if( !is )
    return debug_stream << "  number of facets expected!\n";

  for(int i = 0; i < number_of_facets_on_surface; ++i)
  {
    int i1, i2, i3;
    is >> i1 >> i2 >> i3 >> temp_int;

    Cell_handle ch;
    int i, j, k;
    if( !is || !tr.is_facet(V[i1-1], V[i2-1], V[i3-1], ch, i, j, k) )
      return debug_stream << "  cannot read facet!\n";

    const int index = 6-i-j-k;

    c2t3.set_in_complex(ch, index);
  }
  
  // Tetrahedra
  is >> temp_string;
  
  if( !is || temp_string != "Tetrahedra" )
    return debug_stream << "  \"Tetrahedra\" expected!\n";

  int number_of_cells_in_domain;
  is >> number_of_cells_in_domain;
  if( !is )
    return debug_stream << "  number of cells expected!\n";
  
  for(int i = 0; i < number_of_cells_in_domain; ++i)
  {
    int i1, i2, i3, i4;
    
    is >> i1 >> i2 >> i3 >> i4 >> temp_int;

    Cell_handle ch;
    if( !is || !tr.is_cell(V[i1-1], V[i2-1], V[i3-1],V[i4-1], ch) )
      return debug_stream << "  cannot read cell!\n";
    ch->set_in_domain(true);
  }
  
  // End
  is >> temp_string;

  if( is && temp_string == "End" )
    return true;
  else
    return debug_stream << "  \"End\" expected!\n"; 
} // end input_pslg_from_medit

template < class Tr>
int number_of_cells_in_domain(const Tr& T) {
  int result=0;
  for (typename Tr::Finite_cells_iterator cit = T.finite_cells_begin(); 
       cit != T.finite_cells_end(); ++cit)
    if (cit->is_in_domain ())
      ++result;
  return result;
}

} // end namespace CGAL

#endif // CGAL_IO_FILE_MEDIT_H
