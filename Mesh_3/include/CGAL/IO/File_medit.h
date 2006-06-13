// Copyright (c) 2004-2006  INRIA Sophia-Antipolis (France).
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
#include <CGAL/utility.h>

namespace CGAL {

template <class C2T3>
void
output_to_medit (std::ostream& os, const C2T3& c2t3)
{
  typedef typename C2T3::Triangulation Tr;
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
    os << CGAL::to_double(p.x()) << " " 
       << CGAL::to_double(p.y()) << " " 
       << CGAL::to_double(p.z()) << " " 
       << p.surface_index() << std::endl; // reference
  }

  // Facets
  os << "Triangles" << std::endl
     << c2t3.number_of_facets() << std::endl;
  for( Finite_facets_iterator fit = tr.finite_facets_begin(); 
       fit != tr.finite_facets_end(); ++fit)
  {
    int surface_index = 0;
    if (c2t3.face_status(fit->first,fit->second)
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

} // end output_to_medit

namespace Mesh_3 { namespace details {
  
  class Debug {
    bool active;
    std::ostream* out;
    const std::string header;

    class Debug_aux {
      Debug* debug;
    public:
      Debug_aux(Debug* d) : debug(d) {}

      template <typename T>
      Debug_aux& operator<<(const T& t)
      {
        debug->apply(t, false);
        return *this;
      }

      operator bool() const
      {
        return false;
      }
    };

    Debug_aux aux;

  public:
    Debug(bool debug,
          std::ostream* debug_str = &std::cout,
          const std::string header_string = "")
      : active(debug), out(debug_str), header(header_string),
        aux(this)
    {
    }

    template <class T>
    void apply(T& t, bool with_header = true)
    {
      if(active)
      {
        if(with_header)
          *out << header;
        *out << t;
      }
    }

    template <typename T>
    Debug_aux& operator<<(const T& t)
    {
      apply(t);
      return aux;
    }

    operator bool() const
    {
      return false;
    }
  }; // end internal class Debug
    
//   template <typename T>
//   Debug_aux& operator<<(Debug& d, const T& t)
//   {
//     d.apply(t);
//     return aux;
//   }
} // end namespace details
} // end namespace Mesh_3

  /**************************************************************/
  /*********NOT POSSIBLE !!!!!***********************************/
  /**************************************************************/

// template < class C2T3>
// bool
// input_from_medit (std::istream& is, 
//                   C2T3 & c2t3,
//                   bool debug = false, 
//                   std::ostream* debug_str = &std::cerr) {
//   typedef typename C2T3::Triangulation Tr;
//   typedef typename Tr::Triangulation_data_structure Tds;
//   typedef typename Tr::Vertex_handle Vertex_handle;
//   typedef typename Tr::Point Point;
//   typedef Point_traits<Point> P_traits;
//   typedef typename P_traits::Bare_point Bare_point;
//   typedef typename Tr::Cell_handle Cell_handle;
//   typedef typename Tr::Vertex_handle Vertex_handle;

  
//   Mesh_3::details::Debug debug_stream(debug,
//                                       debug_str,
//                                       "CGAL::input_from_medit()\n"
//                                       "input error:");

//   Tr& tr = c2t3.triangulation();
//   Tds& tds = tr.tds();

//   tr.clear();

//   tds.set_dimension(3);

//   // Header
//   std::string temp_string;
//   int temp_int;

//   is >> temp_string;
//   if( !is || temp_string != "MeshVersionFormatted")
//     return debug_stream << "  not a medit file\n";

//   is >> temp_int; // version
//   if( !is || temp_int!=1 ) 
//     return debug_stream << "  wrong medit version number: " << temp_int
//                         << "\n"
//                         << "  should be 1\n";

//   is >> temp_string;
//   if( !is || temp_string != "Dimension" )
//     return debug_stream << "  \"Dimension\" expected\n";

//   is >> temp_int; // dimension
//   if( !is || temp_int != 3 )
//     return debug_stream << "  dimension 3 expected\n";

//   // Vertices
//   is >> temp_string;
//   if( !is || temp_string != "Vertices" )
//     return debug_stream << "  \"Vertices\" expected\n";
  
//   int number_of_vertices;
//   if( ! (is >> number_of_vertices) )
//     return debug_stream << "  number of vertices expected\n";
  
//   std::vector<Vertex_handle> V(number_of_vertices+1);

//   V[0] = tr.infinite_vertex();

//   for(int i = 0; i < number_of_vertices; ++i)
//   {
//     P_traits point_convert;

//     Bare_point bp;
//     is >> bp;
//     if( !is )
//       return debug_stream << "  point " << i 
//                           << ", error in Bare_point.operator>>()\n";

//     int index;
//     is >> index;
//     if( !is )
//       return debug_stream << "  index expected\n";

//     Point p = point_convert.point(bp);
//     p.set_surface_index(index);

//     Vertex_handle v = tds.create_vertex();
//     v->set_point(p);
//     V[i+1] = v;
//   }

//   // Facets
//   is >> temp_string;
//   if( !is || temp_string != "Triangles" )
//     return debug_stream << "  \"Triangles\" expected\n";

//   int number_of_facets_on_surface;
//   is >> number_of_facets_on_surface;
//   if( !is )
//     return debug_stream << "  number of facets expected\n";

//   std::vector< CGAL::Triple<int, int, int> > 
//     facets_vector(number_of_facets_on_surface);

//   for(int i = 0; i < number_of_facets_on_surface; ++i)
//   {
//     int i1, i2, i3;
//     is >> i1 >> i2 >> i3 >> temp_int;

//     if( !is )
//       return debug_stream << "  cannot read facet " << i << "\n";

//     facets_vector[i] = CGAL::make_triple(i1, i2, i3);
//   }
  
//   // Tetrahedra
//   is >> temp_string;
  
//   if( !is || temp_string != "Tetrahedra" )
//     return debug_stream << "  \"Tetrahedra\" expected\n";

//   int number_of_cells_in_domain;
//   is >> number_of_cells_in_domain;
//   if( !is )
//     return debug_stream << "  number of cells expected\n";
  
//   std::vector<Cell_handle> C(number_of_cells_in_domain);

//   for(int i = 0; i < number_of_cells_in_domain; ++i)
//   {
//     Cell_handle c = create_cell();
//     for (int k=0; k<=3; ++k) {
//       int ik;
//       is >> ik;
//       c->set_vertex(k, V[ik]);
//       V[ik]->set_cell(c);
//     }
//     C[i] = c;
//   }
//   for(int i = 0; i < number_of_cells_in_domain; ++i)
//   {
//     is >> i1 >> i2 >> i3 >> i4 >> temp_int;

//     Cell_handle ch;
//     if( !is )
//       return debug_stream << "  cannot read cell " << i << "\n";
//     if( !tr.is_cell(V[i1], V[i2], V[i3],V[i4], ch) )
//       debug_stream << "  cell " << i << " is not Delaunay:\n"
//                    << i1 << " " << i2 << " " << i3 << " " << i4 << "\n";
//     ch->set_in_domain(true);
//   }

//   // Inserts surface facets, using facets_vector

//   for(int i = 0; i < number_of_facets_on_surface; ++i)
//   {
//     Cell_handle ch;
//     int i, j, k;
//     if( !tr.is_facet(V[facets_vector[i].first],
//                      V[facets_vector[i].second],
//                      V[facets_vector[i].third],
//                      ch, i, j, k) )
//       return debug_stream << "  facet " << i << "is not Delaunay\n";

//     const int index = 6-i-j-k;

//     c2t3.set_in_complex(ch, index);
//   }
  
//   // End
//   is >> temp_string;

//   if( is && temp_string == "End" )
//     return true;
//   else
//     return debug_stream << "  \"End\" expected\n"; 
// } // end input_from_medit

template < class Tr>
int number_of_cells_in_domain(const Tr& T) {
  int result=0;
  for (typename Tr::Finite_cells_iterator cit = T.finite_cells_begin(); 
       cit != T.finite_cells_end(); ++cit)
    if (cit->is_in_domain ())
      ++result;
  return result;
}

template <class C2T3>
int number_of_facets_on_surface_with_index(const C2T3& c2t3,
                                           const unsigned int surface_index)
{
  typedef typename C2T3::Triangulation Tr;

  const Tr& tr = c2t3.triangulation();

  int count = 0;

  for(typename Tr::Finite_facets_iterator fit = tr.finite_facets_begin();
      fit != tr.finite_facets_end();
      ++fit)
  {
    const typename Tr::Cell_handle& cell = fit->first;
    const int index = fit->second;

    if(c2t3.face_status(cell, index) != C2T3::NOT_IN_COMPLEX)
    {
      const typename Tr::Vertex_handle& va = cell->vertex((index+1)&3);
      const typename Tr::Vertex_handle& vb = cell->vertex((index+1)&3);
      const typename Tr::Vertex_handle& vc = cell->vertex((index+1)&3);

      const unsigned int va_index = va->point().surface_index();
      const unsigned int vb_index = vb->point().surface_index();
      const unsigned int vc_index = vc->point().surface_index();

      if(va_index == surface_index &&
         vb_index == surface_index &&
         vc_index == surface_index)
      {
        ++count;
      }
    }
  }
  return count;
} // end number_of_facets_on_surface_with_index(Tr, int) 


} // end namespace CGAL

#endif // CGAL_IO_FILE_MEDIT_H
