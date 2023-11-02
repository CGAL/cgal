// Copyright (c) 1999-2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion
//                 Clement Jamin

#ifndef CGAL_IO_TRIANGULATION_RAW_IOSTREAM_3_H
#define CGAL_IO_TRIANGULATION_RAW_IOSTREAM_3_H

#include <CGAL/license/Triangulation_3.h>

#include <iostream>
#include <vector>

#include <CGAL/Unique_hash_map.h>
#include <CGAL/IO/io.h>
#include <CGAL/Mesh_3/Triangulation_3_fwd.h>

namespace CGAL {

namespace IO {

namespace internal {

template < class GT, class Tds_, class Lds_ >
bool construct_vertices_map(const Triangulation_3<GT, Tds_, Lds_>& tr,
                            Unique_hash_map<typename Triangulation_3<GT, Tds_>::Vertex_handle, std::size_t>& vertices_map)
{
  typedef Triangulation_3<GT, Tds_>               Triangulation;
  typedef typename Triangulation::size_type       size_type;
  typedef typename Triangulation::Vertex_iterator Vertex_iterator;

  size_type i = 0;
  vertices_map[tr.infinite_vertex()] = 0;
  for(Vertex_iterator it = tr.vertices_begin(), end = tr.vertices_end(); it != end; ++it)
  {
    vertices_map[it] = i++;
  }

  CGAL_assertion(i == tr.number_of_vertices()+1);

  return true;
}

template < class GT, class Tds_, class Lds_ >
bool write_header(std::ostream& os, const Triangulation_3<GT, Tds_, Lds_>& tr)
{
  // Writes:
  // - the dimension
  if(IO::is_ascii(os))
  {
    os << tr.dimension() << "\n";
  }
  else
  {
    write(os, tr.dimension());
  }
  return (bool)os;
}

template < class GT, class Tds_, class Lds_ >
bool read_header(std::istream& is, Triangulation_3<GT, Tds_, Lds_>& tr)
{
  // Reads:
  // - the dimension
  int dimension;
  if(IO::is_ascii(is))
  {
    is >> dimension;
  }
  else
  {
    read(is, dimension);
  }
  if(dimension > 3 || dimension < -2)
    return false;
  tr.tds().set_dimension(dimension);
  return (bool)is;
}

template < class GT, class Tds_, class Lds_ >
bool write_vertices(std::ostream& os, const Triangulation_3<GT, Tds_, Lds_>& tr)
{
  // Writes:
  // - the number of finite vertices
  // - the non combinatorial information on vertices (point, etc)
  typedef Triangulation_3<GT, Tds_>               Triangulation;
  typedef typename Triangulation::size_type       size_type;
  typedef typename Triangulation::Vertex_iterator Vertex_iterator;

  size_type n = tr.number_of_vertices();

  if(IO::is_ascii(os))
  {
    os << n << '\n';
  }
  else
  {
    write(os, n);
  }

  if (n == 0)
    return false;

  for(Vertex_iterator it = ++tr.vertices_begin(), end = tr.vertices_end(); it != end; ++it)
  {
    os << *it;
    if(IO::is_ascii(os))
      os << '\n';
  }

  return (bool)os;
}

template < class GT, class Tds_, class Lds_,
           typename Tr_src = Triangulation_3<GT, Tds_, Lds_>,
           typename ConvertVertex = Nullptr_t >
struct read_vertices_functor
{
  bool operator()(std::istream& is, Triangulation_3<GT, Tds_, Lds_>& tr,
                  std::vector< typename Triangulation_3<GT, Tds_>::Vertex_handle >& vertices_handles,
                  typename Triangulation_3<GT, Tds_>::size_type n,
                  ConvertVertex convert_vertex = ConvertVertex())
  {
    // Reads:
    // - the non combinatorial information on vertices (point, etc)
    // The vertices are read with the type Tr_src::Vertex
    // and converted to type Triangulation_3<GT, Tds_, Lds_>::Vertex
    // using the ConvertVertex convertor.
    typedef typename Tr_src::Vertex Vertex1;

    for(std::size_t i=1; i <= n; i++)
    {
      Vertex1 v;
      if(!(is >> v))
        return false;
      vertices_handles[i] = tr.tds().create_vertex();
      convert_vertex(v, *vertices_handles[i]);
    }

    return (bool)is;
  }
};

template < class GT, class Tds_, class Lds_, class Tr_src >
struct read_vertices_functor < GT, Tds_, Lds_, Tr_src, Nullptr_t >
{
  bool operator()(std::istream& is, Triangulation_3<GT, Tds_, Lds_>& tr,
                  std::vector< typename Triangulation_3<GT, Tds_>::Vertex_handle >& vertices_handles,
                  typename Triangulation_3<GT, Tds_>::size_type n,
                  Nullptr_t _ = Nullptr_t())
  {
    // Reads:
    // - the non combinatorial information on vertices (point, etc)
    for(std::size_t i=1; i <= n; i++)
    {
      vertices_handles[i] = tr.tds().create_vertex();
      if(!(is >> *vertices_handles[i]))
        return false;
    }
    return (bool)is;
  }
};


template < class GT, class Tds_, class Lds_,
           typename Tr_src = Triangulation_3<GT, Tds_, Lds_>,
           typename ConvertVertex = Nullptr_t >
bool read_vertices(std::istream& is, Triangulation_3<GT, Tds_, Lds_>& tr,
                   std::vector< typename Triangulation_3<GT, Tds_>::Vertex_handle >& vertices_handles,
                   ConvertVertex convert_vertex = ConvertVertex())
{
  // Reads:
  // - the number of finite vertices
  // - the non combinatorial information on vertices (point, etc)
  typedef Triangulation_3<GT, Tds_>         Triangulation;
  typedef typename Triangulation::size_type size_type;

  size_type n;
  if(IO::is_ascii(is))
  {
    is >> n;
  }
  else
  {
    read(is, n);
  }
  if((n+1) > vertices_handles.max_size())
    return false;

  vertices_handles.resize(n+1);
  vertices_handles[0] = tr.infinite_vertex(); // the infinite vertex is numbered 0

  if (!read_vertices_functor<GT, Tds_, Lds_, Tr_src, ConvertVertex>()(is, tr, vertices_handles, n, convert_vertex))
    return false;

  return (bool)is;
}

struct Cell_accessor
{
  template <class Iterator>
  auto operator()(const Iterator& it) { return it; }
};

struct Cell_accessor_first
{
  template <class Iterator>
  auto operator()(const Iterator& it) { return it->first; }
};

template <class Iterator, class Cell_accessor>
bool write_cell_info(std::ostream& os, Iterator start, const Iterator& end)
{
  Cell_accessor accessor;
  if(IO::is_ascii(os))
  {
    for(; start != end; ++start)
    {
      os << *(accessor(start)) << '\n';
    }
  }
  else
  {
    for(; start != end; ++start)
    {
      os << *(accessor(start));
    }
  }
  return (bool)os;
}

template < class GT, class Tds_, class Lds_ >
bool write_cells(std::ostream& os, const Triangulation_3<GT, Tds_, Lds_>& tr,
                 const Unique_hash_map<typename Triangulation_3<GT, Tds_>::Vertex_handle, std::size_t>& vertices_map)
{
  // Writes:
  // [Call to Tds_::print_cells]
  // - the number of cells
  // - the cells by the indices of their vertices in the preceding list
  //   of vertices, plus the non combinatorial information on each cell
  // - the neighbors of each cell by their index in the preceding list of cells
  // [Cells other info]
  // - when dimension < 3 : the same with faces of maximal dimension

  typedef Triangulation_3<GT, Tds_>              Triangulation;
  typedef typename Triangulation::Cell_iterator  Cell_iterator;
  typedef typename Triangulation::Edge_iterator  Edge_iterator;
  typedef typename Triangulation::Facet_iterator Facet_iterator;

  tr.tds().print_cells(os, vertices_map);
  if (!os)
    return false;

  // Write the non combinatorial information on the cells
  // using the << operator of Cell.
  // Works because the iterator of the tds traverses the cells in the
  // same order as the iterator of the triangulation
  switch(tr.dimension())
  {
    case 3:
    {
      return write_cell_info<Cell_iterator, Cell_accessor>(os, tr.cells_begin(), tr.cells_end());
      break;
    }
    case 2:
    {
      return write_cell_info<Facet_iterator, Cell_accessor_first>(os, tr.facets_begin(), tr.facets_end());
      break;
    }
    case 1:
    {
      return write_cell_info<Edge_iterator, Cell_accessor_first>(os, tr.edges_begin(), tr.edges_end());
      break;
    }
  }
  return false;
}


template < class GT, class Tds_, class Lds_,
           typename Tr_src = Triangulation_3<GT, Tds_, Lds_>,
           typename ConvertCell = Nullptr_t >
struct read_cells_functor
{
  bool operator()(std::istream& is, Triangulation_3<GT, Tds_, Lds_>& tr,
                  const std::vector< typename Triangulation_3<GT, Tds_>::Vertex_handle >& vertices_handles,
                  typename Triangulation_3<GT, Tds_>::size_type m, std::vector< typename Triangulation_3<GT, Tds_>::Cell_handle >& C,
                  ConvertCell convert_cell = ConvertCell())
  {
    // Reads:
    // [Cells other info]
    // - when dimension < 3 : the same with faces of maximal dimension
    // The cell are read with the type Tr_src::Cell
    // and converted to type Triangulation_3<GT, Tds_, Lds_>::Cell
    // using the ConvertCell convertor.

    typedef typename Tr_src::Cell Cell1;

    for(std::size_t i=0 ; i < m; i++)
    {
      Cell1 c;
      if(!(is >> c))
        return false;
      convert_cell(c, *C[i]);
    }

    return (bool)is;
  }
};

template < class GT, class Tds_, class Lds_, class Tr_src >
struct read_cells_functor < GT, Tds_, Lds_, Tr_src, Nullptr_t >
{
  bool operator()(std::istream& is, Triangulation_3<GT, Tds_, Lds_>& tr,
                  const std::vector< typename Triangulation_3<GT, Tds_>::Vertex_handle >& vertices_handles,
                  typename Triangulation_3<GT, Tds_>::size_type m, std::vector< typename Triangulation_3<GT, Tds_>::Cell_handle >& C,
                  Nullptr_t _ = Nullptr_t())
  {
    // Reads:
    // - the non combinatorial information on vertices (point, etc)
    for(std::size_t i=0 ; i < m; i++)
    {
      if(!(is >> *C[i]))
        return false;
    }

    return (bool)is;
  }
};

template < class GT, class Tds_, class Lds_,
           typename Tr_src = Triangulation_3<GT, Tds_, Lds_>,
           typename ConvertCell = Nullptr_t >
bool read_cells(std::istream& is, Triangulation_3<GT, Tds_, Lds_>& tr,
                const std::vector< typename Triangulation_3<GT, Tds_>::Vertex_handle >& vertices_handles,
                ConvertCell convert_cell = ConvertCell())
{
  // Reads:
  // [Call to Tds_::read_cells]
  // - the number of cells
  // - the cells by the indices of their vertices in the preceding list
  //   of vertices, plus the non combinatorial information on each cell
  // - the neighbors of each cell by their index in the preceding list of cells
  // [Cells other info]
  // - when dimension < 3 : the same with faces of maximal dimension

  typedef Triangulation_3<GT, Tds_>                 Triangulation;
  typedef typename Triangulation::Cell_handle      Cell_handle;

  std::vector< Cell_handle > C;

  std::size_t m;
  tr.tds().read_cells(is, vertices_handles, m, C);
  if(!is)
    return false;

  if (!read_cells_functor<GT, Tds_, Lds_, Tr_src, ConvertCell>()(is, tr, vertices_handles, m, C, convert_cell))
    return false;

  return (bool)is;
}

} // end namespace internal


template < class GT, class Tds_, class Lds_ >
bool export_triangulation_3(std::ostream& os, const Triangulation_3<GT, Tds_, Lds_>& tr,
                            Unique_hash_map<typename Triangulation_3<GT, Tds_>::Vertex_handle, std::size_t>& vertices_map)
{
  // Writes:
  // [Header]
  //   - the dimension
  // [Vertices]
  //   - the number of finite vertices
  //   - the non combinatorial information on vertices (point, etc)
  // [Cells]
  //     - the number of cells
  //     - the cells by the indices of their vertices in the preceding list
  //       of vertices, plus the non combinatorial information on each cell
  //   [Cells combinatorial information]
  //     - the neighbors of each cell by their index in the preceding list of cells
  //   [Cells other info]
  //     - when dimension < 3 : the same with faces of maximal dimension

  return internal::write_header(os, tr)
      && internal::write_vertices(os, tr)
      && internal::construct_vertices_map(tr, vertices_map)
      && internal::write_cells(os, tr, vertices_map);
}

template < class GT, class Tds_, class Lds_ >
bool export_triangulation_3(std::ostream& os, const Triangulation_3<GT, Tds_, Lds_>& tr)
{
  Unique_hash_map<typename Triangulation_3<GT, Tds_>::Vertex_handle, std::size_t> vertices_map;
  return export_triangulation_3(os, tr, vertices_map);
}


template < class GT, class Tds_, class Lds_,
           typename Tr_src = Triangulation_3<GT, Tds_, Lds_>,
           typename ConvertVertex = Nullptr_t,
           typename ConvertCell = Nullptr_t >
bool import_triangulation_3(std::istream& is, Triangulation_3<GT, Tds_, Lds_>& tr,
                            std::vector<typename Triangulation_3<GT, Tds_>::Vertex_handle >& vertices_handles,
                            ConvertVertex convert_vertex = ConvertVertex(),
                            ConvertCell convert_cell = ConvertCell())
{
  // Reads:
  // [Header]
  //   - the dimension
  // [Vertices]
  //   - the number of finite vertices
  //   - the non combinatorial information on vertices (point, etc)
  // [Cells]
  //     - the number of cells
  //     - the cells by the indices of their vertices in the preceding list
  //       of vertices, plus the non combinatorial information on each cell
  //   [Cells combinatorial information]
  //     - the neighbors of each cell by their index in the preceding list of cells
  //   [Cells other info]
  //     - when dimension < 3 : the same with faces of maximal dimension

  tr.tds().clear(); // infinite vertex deleted
  tr.set_infinite_vertex(tr.tds().create_vertex());

  bool result = internal::read_header(is, tr)
             && internal::read_vertices<GT, Tds_, Lds_, Tr_src, ConvertVertex>
                                       (is, tr, vertices_handles, convert_vertex)
             && internal::read_cells<GT, Tds_, Lds_, Tr_src, ConvertCell>
                                    (is, tr, vertices_handles, convert_cell);

  CGAL_assertion(result && tr.is_valid(false));
  return result;
}

template < class GT, class Tds_, class Lds_,
           typename Tr_src = Triangulation_3<GT, Tds_, Lds_>,
           typename ConvertVertex = Nullptr_t,
           typename ConvertCell = Nullptr_t >
bool import_triangulation_3(std::istream& is, Triangulation_3<GT, Tds_, Lds_>& tr,
                            ConvertVertex convert_vertex = ConvertVertex(),
                            ConvertCell convert_cell = ConvertCell())
{
  std::vector<typename Triangulation_3<GT, Tds_>::Vertex_handle > vertices_handles;
  return import_triangulation_3<GT, Tds_, Lds_, Tr_src, ConvertVertex, ConvertCell>
                               (is, tr, vertices_handles, convert_vertex, convert_cell);
}


} // end namespace IO

} // end namespace CGAL

#endif // CGAL_IO_TRIANGULATION_RAW_IOSTREAM_3_H
