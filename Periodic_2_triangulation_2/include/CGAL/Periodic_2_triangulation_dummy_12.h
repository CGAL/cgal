// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Nico Kruithof <Nico@nghk.nl>

#ifdef CGAL_INCLUDE_FROM_PERIODIC_2_TRIANGULATION_2_H

#include <vector>

template < class GT, class Tds >
inline std::vector<typename Periodic_2_triangulation_2<GT, Tds>::Vertex_handle >
Periodic_2_triangulation_2<GT, Tds>::insert_dummy_points()
{
  clear();

  Vertex_handle vertices[12];
  // 6 faces per row, 4 rows
  Face_handle faces[24];

  // Initialise vertices:
  for (int i = 0; i < 4; i++)
    {
      for (int j = 0; j < 3; j++)
        {
          // Initialise virtual vertices out of the domain for debugging
          vertices[3 * i + j] = _tds.create_vertex();
          Point p(j * (1.0 / 3.0) + i * (1.0 / 6.0), i * (1.0 / 4.0) );
          p = Point((p.x() > FT(0.9375) ? (std::max)( p.x() - 1, FT(0) ) : p.x()),
                    p.y());
          p = Point((_domain.xmax() - _domain.xmin()) * p.x(),
                    (_domain.xmax() - _domain.xmin()) * p.y());
          p = Point(p.x() + _domain.xmin(),
                    p.y() + _domain.ymin());
          vertices[3 * i + j]->set_point(p);
        }
    }

  // Create faces:
  for (int i = 0; i < 24; i++)
    {
      faces[i] = _tds.create_face();
    }

  // bottom row
  faces[0]->set_vertices(vertices[0], vertices[1], vertices[3]);
  faces[1]->set_vertices(vertices[1], vertices[2], vertices[4]);
  faces[2]->set_vertices(vertices[2], vertices[0], vertices[5]);
  faces[3]->set_vertices(vertices[0], vertices[3], vertices[5]);
  faces[4]->set_vertices(vertices[1], vertices[4], vertices[3]);
  faces[5]->set_vertices(vertices[2], vertices[5], vertices[4]);
  // second row
  faces[6]->set_vertices(vertices[3], vertices[4], vertices[6]);
  faces[7]->set_vertices(vertices[4], vertices[5], vertices[7]);
  faces[8]->set_vertices(vertices[5], vertices[3], vertices[8]);
  faces[9]->set_vertices(vertices[3], vertices[6], vertices[8]);
  faces[10]->set_vertices(vertices[4], vertices[7], vertices[6]);
  faces[11]->set_vertices(vertices[5], vertices[8], vertices[7]);
  // third row
  faces[12]->set_vertices(vertices[6], vertices[7], vertices[9]);
  faces[13]->set_vertices(vertices[7], vertices[8], vertices[10]);
  faces[14]->set_vertices(vertices[8], vertices[6], vertices[11]);
  faces[15]->set_vertices(vertices[6], vertices[9], vertices[11]);
  faces[16]->set_vertices(vertices[7], vertices[10], vertices[9]);
  faces[17]->set_vertices(vertices[8], vertices[11], vertices[10]);
  // fourth row
  faces[18]->set_vertices(vertices[9], vertices[10], vertices[2]);
  faces[19]->set_vertices(vertices[10], vertices[11], vertices[0]);
  faces[20]->set_vertices(vertices[11], vertices[9], vertices[1]);
  faces[21]->set_vertices(vertices[9], vertices[2], vertices[1]);
  faces[22]->set_vertices(vertices[10], vertices[0], vertices[2]);
  faces[23]->set_vertices(vertices[11], vertices[1], vertices[0]);

  faces[0]->set_neighbors(faces[4], faces[3], faces[23]);
  faces[1]->set_neighbors(faces[5], faces[4], faces[21]);
  faces[2]->set_neighbors(faces[3], faces[5], faces[22]);
  faces[3]->set_neighbors(faces[8], faces[2], faces[0]);
  faces[4]->set_neighbors(faces[6], faces[0], faces[1]);
  faces[5]->set_neighbors(faces[7], faces[1], faces[2]);

  faces[6]->set_neighbors(faces[10], faces[9], faces[4]);
  faces[7]->set_neighbors(faces[11], faces[10], faces[5]);
  faces[8]->set_neighbors(faces[9], faces[11], faces[3]);
  faces[9]->set_neighbors(faces[14], faces[8], faces[6]);
  faces[10]->set_neighbors(faces[12], faces[6], faces[7]);
  faces[11]->set_neighbors(faces[13], faces[7], faces[8]);

  faces[12]->set_neighbors(faces[16], faces[15], faces[10]);
  faces[13]->set_neighbors(faces[17], faces[16], faces[11]);
  faces[14]->set_neighbors(faces[15], faces[17], faces[9]);
  faces[15]->set_neighbors(faces[20], faces[14], faces[12]);
  faces[16]->set_neighbors(faces[18], faces[12], faces[13]);
  faces[17]->set_neighbors(faces[19], faces[13], faces[14]);

  faces[18]->set_neighbors(faces[22], faces[21], faces[16]);
  faces[19]->set_neighbors(faces[23], faces[22], faces[17]);
  faces[20]->set_neighbors(faces[21], faces[23], faces[15]);
  faces[21]->set_neighbors(faces[1], faces[20], faces[18]);
  faces[22]->set_neighbors(faces[2], faces[18], faces[19]);
  faces[23]->set_neighbors(faces[0], faces[19], faces[20]);

  set_offsets(faces[0], 0, 0, 0);
  set_offsets(faces[1], 0, 0, 0);
  set_offsets(faces[2], 0, 2, 0);
  set_offsets(faces[3], 2, 2, 0);
  set_offsets(faces[4], 0, 0, 0);
  set_offsets(faces[5], 0, 0, 0);
  set_offsets(faces[6], 0, 0, 0);
  set_offsets(faces[7], 0, 0, 0);
  set_offsets(faces[8], 0, 2, 2);
  set_offsets(faces[9], 0, 0, 0);
  set_offsets(faces[10], 0, 0, 0);
  set_offsets(faces[11], 0, 2, 0);
  set_offsets(faces[12], 0, 0, 0);
  set_offsets(faces[13], 0, 2, 0);
  set_offsets(faces[14], 0, 0, 0);
  set_offsets(faces[15], 0, 0, 0);
  set_offsets(faces[16], 0, 0, 0);
  set_offsets(faces[17], 2, 2, 0);
  set_offsets(faces[18], 0, 0, 1);
  set_offsets(faces[19], 0, 2, 3);
  set_offsets(faces[20], 0, 0, 1);
  set_offsets(faces[21], 0, 1, 1);
  set_offsets(faces[22], 0, 3, 1);
  set_offsets(faces[23], 0, 1, 1);

  vertices[0]->set_face(faces[0]);
  vertices[1]->set_face(faces[1]);
  vertices[2]->set_face(faces[2]);
  vertices[3]->set_face(faces[3]);
  vertices[4]->set_face(faces[4]);
  vertices[5]->set_face(faces[5]);
  vertices[6]->set_face(faces[6]);
  vertices[7]->set_face(faces[7]);
  vertices[8]->set_face(faces[8]);
  vertices[9]->set_face(faces[12]);
  vertices[10]->set_face(faces[13]);
  vertices[11]->set_face(faces[14]);

  _tds.set_dimension(2);
  _cover = make_array(1, 1);

  std::vector<Vertex_handle> ret_vector(12);
  for (int i = 0; i < 12; i++)
    {
      ret_vector[i] = vertices[i];
    }

  CGAL_assertion(is_valid());

  return ret_vector;
}

#endif // CGAL_INCLUDE_FROM_PERIODIC_2_TRIANGULATION_2_H
