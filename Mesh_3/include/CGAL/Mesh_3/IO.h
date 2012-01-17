// Copyright (c) 2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESH_3_IO_H
#define CGAL_MESH_3_IO_H

#include <CGAL/IO/File_medit.h> // for Debug
#include <iostream>
#include <map>

#include <string>

namespace CGAL { 

template <class T>
struct Get_io_signature
{
  std::string operator()() 
  {
    return T::io_signature();
  }
};

template <class Kernel>
struct Get_io_signature<Point_3<Kernel> >
{
  std::string operator()() {
    return "Point_3";
  }
};

template <class Point, typename FT>
struct Get_io_signature<Weighted_point<Point, FT> >
{
  std::string operator()() {
    return std::string("Weighted_point<") + Get_io_signature<Point>()() + ">";
  }
};

#ifdef CGAL_TRIANGULATION_3_H
template <class Gt, class Vb, class Cb>
struct
Get_io_signature<Triangulation_3<Gt, Triangulation_data_structure_3<Vb, Cb> > >
{
  std::string operator()() {
    return std::string("Triangulation_3(") +
      Get_io_signature<typename Gt::Point_3>()() +
      ",Vb(" + Get_io_signature<Vb>()() +
      "),Cb(" + Get_io_signature<Cb>()() +
      "))";
  }
};
#endif

#ifdef CGAL_DELAUNAY_TRIANGULATION_3_H
template <class Gt, class Tds>
struct 
Get_io_signature<Delaunay_triangulation_3<Gt, Tds> >
{
  std::string operator()() {
    return Get_io_signature<Triangulation_3<Gt, Tds> >()();
  }
};
#endif

#ifdef CGAL_REGULAR_TRIANGULATION_3_H
template <class Gt, class Tds>
struct
Get_io_signature<Regular_triangulation_3<Gt, Tds> >
{
  std::string operator()() {
    return Get_io_signature<Triangulation_3<Gt, Tds> >()();
  }
};
#endif

#ifdef CGAL_TRIANGULATION_VERTEX_BASE_3_H
template <class Gt>
struct Get_io_signature<Triangulation_vertex_base_3<Gt> >
{
  std::string operator()() {
    return "Tcb_3";
  }
};
#endif

#ifdef CGAL_TRIANGULATION_CELL_BASE_3_H
template <class Gt>
struct
Get_io_signature<Triangulation_cell_base_3<Gt> >
{
  std::string operator()() {
    return "Tcb_3";
  }
};
#endif

#ifdef CGAL_REGULAR_TRIANGULATION_CELL_BASE_3_H
template <class Gt, class Cb, class Container>
struct
Get_io_signature<Regular_triangulation_cell_base_3<Gt, Cb, Container> >
{
  std::string operator()() {
    return "Tcb_3";
  }
};
#endif

namespace Mesh_3 {

/** Ouput a C2T3 corresponding to a mesh to a std::ostream.
    This function assumes that the dimension of the triangulation is 3.
*/
template <class C2T3>
void output_mesh(std::ostream& os, const C2T3& c2t3)
{
//   typedef typename C2T3::Triangulation Tr;
//   typedef typename Tr::Cell_iterator Cell_iterator;
//   typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
//   typedef typename Tr::Vertex_handle Vertex_handle;
//   typedef typename Tr::Point Point;

//   const Tr& tr = c2t3.triangulation();

  bool ascii = is_ascii(os);
  set_ascii_mode(os);

  if(ascii)
    os << "format: ";
  else
    os << "binaryformat: ";

  os << Get_io_signature<C2T3>()() << "\n";

  if(ascii)
    set_ascii_mode(os);
  else
    set_binary_mode(os);

  os << c2t3;

//   int number_of_vertices = tr.number_of_vertices();

//   os << number_of_vertices;
//   if (is_ascii(os))
//     os << std::endl;

//   std::map<Vertex_handle, int> indices_of_vertices;

//   indices_of_vertices[tr.infinite_vertex()] = 0;

//   int i = 0;
//   for (Finite_vertices_iterator vit=tr.finite_vertices_begin();
//        vit!=tr.finite_vertices_end();
//        ++vit)
//   {
//     indices_of_vertices[vit] = ++i;
//     os << *vit;
//     if (is_ascii(os))
// 	os << std::endl;
//   }

//   tr.tds().print_cells(os, indices_of_vertices);

//   for(Cell_iterator cit=tr.cells_begin(); cit != tr.cells_end(); ++cit)
//   { 
//     os << *cit;
//     if (is_ascii(os))
//       os << std::endl;
//   }
} // end output_mesh

template <class C2T3>
bool
input_mesh(std::istream& is, 
           C2T3 & c2t3,
           bool debug = false, 
           std::ostream* debug_str = &std::cout)
{
//   typedef typename C2T3::Triangulation Tr;
//   typedef typename Tr::Triangulation_data_structure Tds;
//   typedef typename Tr::Vertex_handle Vertex_handle;  
//   typedef typename Tr::Cell_handle Cell_handle;  

  details::Debug debug_stream(debug,
                              debug_str,
                              "CGAL::Mesh_3::input_mesh()"
                              " input error: ");

  bool ascii = is_ascii(is);
  bool will_be_ascii = ascii;
  set_ascii_mode(is);

  std::string format;
  is >> format;
  if( format == "binaryformat:" )
  {
    if(ascii)
      *debug_str << "Warning: switching stream to binary mode.\n";
    will_be_ascii = false;
  }
  else
    if( format == "format:" )
    {
      if(!ascii)
        *debug_str << "Warning: switching stream to ascii mode.\n";
      will_be_ascii = true;
    }
    else
    {
      debug_stream << "Bad file format!\n";
      debug_stream << "expected \"format: \" or \"binaryformat\", found \"" << format << "\" instead.\n";
      return false;
    }
  
  {
    char ret;
    is.get(ret);
    if( ret != ' ' )
      return debug_stream << "Expected ' ', found '" << ret << "'\n";
  }

  is >> format;
  if( !is || format != Get_io_signature<C2T3>()() )
  {
    debug_stream << "bad format \"" << format << "\"\n";
    debug_stream << "expected format \"" << Get_io_signature<C2T3>()() << "\"\n";
    return false;
  }

   
  {
    char ret;
    is.get(ret);
    if( ret != '\n' )
      return debug_stream << "Expected '\n', found '" << ret << "'\n";
  }

  if(will_be_ascii)
    set_ascii_mode(is);
  else
    set_binary_mode(is);

  if( ! (is >> c2t3) )
    return debug_stream << "Cannot read the triangulation.\n";
  else
    return true;

//   Tds& tds = tr.tds();

//   tr.clear();

//   tds.set_dimension(3);

//   int number_of_vertices;
//   is >> number_of_vertices;

//   if( !is)
//     return debug_stream << "Cannot read number of vertices!\n";

//   std::map<int, Vertex_handle> vertex_handles;

//   vertex_handles[0] = tr.infinite_vertex();

//   for(int i = 1; i <= number_of_vertices; ++i)
//   {
//     vertex_handles[i] = tds.create_vertex();
//     is >> *vertex_handles[i];
//     if( !is )
//       return debug_stream << "Cannot read vertex:" << i << "!\n";
//   }

//   std::map<int, Cell_handle > C;

//   int m;
//   tds.read_cells(is, vertex_handles, m, C);

//   if( !is )
//     return debug_stream << "Cannot read cells!\n";

//   std::cerr << "Reading infos\n";
//   for (int j=0 ; j < m; j++)
//   {
//     is >> *(C[j]);
//     if( !is )
//       return debug_stream << "Cannot read informations for cell:"
//                           << j << "!\n";
//   }

//    CGAL_triangulation_assertion( tr.is_valid(true, 3) );

//   return is;
} // end input_mesh

} // end namespace Mesh_3
} // end namespace CGAL


#endif // CGAL_MESH_3_IO_H
