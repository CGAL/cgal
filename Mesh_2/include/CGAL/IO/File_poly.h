// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_FILE_POLY_H
#define CGAL_FILE_POLY_H

#include <CGAL/license/Mesh_2.h>
#include <CGAL/IO/OFF/File_scanner_OFF.h>

namespace CGAL {

namespace IO {

//the function that reads a Shewchuk Triangle .poly file
template <class CDT, class OutputIterator>
void
read_triangle_poly_file(CDT& t, std::istream &f,
                        OutputIterator seeds)
{
  typedef typename CDT::Vertex_handle Vertex_handle;
  typedef typename CDT::Point Point;

  t.clear();

  unsigned int number_of_points;
  skip_comment_OFF(f);
  f >> number_of_points;
  skip_until_EOL(f);
  skip_comment_OFF(f);

  // read vertices
  std::vector<Vertex_handle> vertices(number_of_points);
  for(unsigned int i = 0; i < number_of_points; ++i)
    {
      unsigned int j;
      double x, y;
      f >> j >> iformat(x) >> iformat(y);
      Point p(x, y);
      skip_until_EOL(f); skip_comment_OFF(f);
      vertices[--j] = t.insert(p);
    }

  // read segments
  unsigned int number_of_segments;
  f >> number_of_segments;
  skip_until_EOL(f); skip_comment_OFF(f);
  for(unsigned int k = 0; k < number_of_segments; ++k)
    {
      unsigned int l, v1, v2;
      f >> l >> v1 >> v2;
      skip_until_EOL(f); skip_comment_OFF(f);
      t.insert_constraint(vertices[--v1], vertices[--v2]);
    }

  // read holes
  unsigned int number_of_holes;
  f >> number_of_holes;
  for(unsigned int m = 0; m < number_of_holes; ++m)
    {
      unsigned int n;
      Point p;
      f >> n >> p;
      skip_until_EOL(f); skip_comment_OFF(f);
      *seeds++ = p;
    }
}

template <class CDT>
inline
void
read_triangle_poly_file(CDT& t, std::istream &f)
{
  read_triangle_poly_file(t, f, Emptyset_iterator());
}

//the function that write a Shewchuk Triangle .poly file
template <class CDT, typename InputIterator>
void
write_triangle_poly_file(const CDT& t, std::ostream &f,
                         InputIterator begin, InputIterator end)
{
  typedef typename CDT::Vertex_handle Vertex_handle;
  typedef typename CDT::Finite_vertices_iterator
    Finite_vertices_iterator;
  typedef typename CDT::Finite_edges_iterator
    Finite_edges_iterator;

  std::map<Vertex_handle, unsigned int> index;

  // write vertices
  f << "# Shewchuk Triangle .poly file, produced by the CGAL::Mesh_2 package"
    << std::endl
    << "# Neither attributes nor boundary markers are used." << std::endl
    << t.number_of_vertices() << " " << 2 << " "
    << 0 << " " << 0 << std::endl;

  f << std::endl;

  unsigned int vertices_counter = 0;
  for(Finite_vertices_iterator vit = t.finite_vertices_begin();
      vit != t.finite_vertices_end();
      ++vit)
    {
      f << ++vertices_counter << " " << vit->point() << std::endl;
      index[vit] = vertices_counter;
    }

  f << std::endl;

  // write constrained edges

  int number_of_constrained_edges = 0;
  for(Finite_edges_iterator it = t.finite_edges_begin();
      it != t.finite_edges_end();
      ++it)
    if(it->first->is_constrained(it->second))
      ++number_of_constrained_edges;

  f << number_of_constrained_edges << " " << 0 << std::endl;
  unsigned int edges_counter = 0;

  for(Finite_edges_iterator eit = t.finite_edges_begin();
      eit != t.finite_edges_end();
      ++eit)
    if(eit->first->is_constrained(eit->second))
      f << ++edges_counter << " "
        << index[eit->first->vertex(t.cw(eit->second))] << " "
        << index[eit->first->vertex(t.ccw(eit->second))]
        << std::endl;

  f << std::endl;


  // write seeds, assuming that the seeds unmarks faces
  f << std::distance(begin, end) << std::endl;
  unsigned int seeds_counter = 0;
  for(InputIterator sit = begin;
      sit!=end; ++sit)
    f << ++seeds_counter << " " << *sit << std::endl;
}

//the same without holes.
template <class CDT>
void
write_triangle_poly_file(const CDT& t, std::ostream &f)
{
  std::list<int> l;

  write_triangle_poly_file(t, f, l.begin(), l.end());
}

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE
using IO::read_triangle_poly_file;
using IO::write_triangle_poly_file;
#endif

} // end namespace CGAL

#endif // CGAL_FILE_POLY_H
