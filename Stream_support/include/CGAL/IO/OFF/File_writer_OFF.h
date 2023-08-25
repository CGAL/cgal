// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_IO_OFF_FILE_WRITER_OFF_H
#define CGAL_IO_OFF_FILE_WRITER_OFF_H

#include <CGAL/IO/binary_file_io.h>
#include <CGAL/IO/OFF/File_header_OFF.h>

#include <iostream>
#include <cstddef>

namespace CGAL {

class File_writer_OFF
{
  std::ostream* m_os;
  File_header_OFF m_header;

public:
  File_writer_OFF(bool verbose = false) : m_os(nullptr), m_header(verbose) {}
  File_writer_OFF(const File_header_OFF& h) : m_os(nullptr), m_header(h) {}

  std::ostream& out() { return *m_os; }
  File_header_OFF& header() { return m_header; }
  const File_header_OFF& header() const { return m_header; }

  void write_header(std::ostream& os,
                    std::size_t vertices,
                    std::size_t /*halfedges*/,
                    std::size_t facets,
                    bool colors = false,
                    bool normals = false,
                    bool textures = false)
  {
    m_os = &os;

    m_header.set_vertices(vertices);
    m_header.set_facets(facets);
    m_header.set_normals(normals);
    m_header.set_colors(colors);
    m_header.set_textures(textures);

    // Print header.
    out() << m_header;
  }

  void write_footer()
  {
    out() << std::endl;
  }

  void write_vertex(const double x, const double y, const double z)
  {
    if(m_header.binary())
    {
      I_Binary_write_big_endian_float32(out(), float(x));
      I_Binary_write_big_endian_float32(out(), float(y));
      I_Binary_write_big_endian_float32(out(), float(z));
    }
    else
    {
      out() << '\n' << IO::oformat(x) << ' ' << IO::oformat(y) << ' ' << IO::oformat(z);
    }
  }

  void write_vertex_normal(const double x, const double y, const double z)
  {
    if(m_header.binary())
    {
      I_Binary_write_big_endian_float32(out(), float(x));
      I_Binary_write_big_endian_float32(out(), float(y));
      I_Binary_write_big_endian_float32(out(), float(z));
    }
    else
    {
      out() << ' ' << ' ' << IO::oformat(x) << ' ' << IO::oformat(y) << ' ' << IO::oformat(z);
    }
  }

  void write_vertex_color(const double r, const double g, const double b)
  {
    if(m_header.binary())
    {
      I_Binary_write_big_endian_float32(out(), float(r));
      I_Binary_write_big_endian_float32(out(), float(g));
      I_Binary_write_big_endian_float32(out(), float(b));
    }
    else
    {
      out() << ' ' << ' ' << IO::oformat(r) << ' ' << IO::oformat(g) << ' ' << IO::oformat(b);
    }
  }

  void write_vertex_texture(const double tx, const double ty)
  {
    if(m_header.binary())
    {
      I_Binary_write_big_endian_float32(out(), float(tx));
      I_Binary_write_big_endian_float32(out(), float(ty));
    }
    else
    {
      out() << ' ' << ' ' << IO::oformat(tx) << ' ' << IO::oformat(ty);
    }
  }

  void write_facet_header()
  {
    if(m_header.ascii())
    {
      out() << '\n';
    }
  }

  void write_facet_begin(std::size_t no)
  {
    if(m_header.binary())
      I_Binary_write_big_endian_integer32(out(), static_cast<std::int32_t>(no));
    else
      out() << no << ' ';
  }

  void write_facet_vertex_index(std::size_t index)
  {
    if(m_header.binary())
      I_Binary_write_big_endian_integer32(out(), static_cast<std::int32_t>(index));
    else
      out() << ' ' << index;
  }

  void write_face_color(const double r, const double g, const double b)
  {
    write_vertex_color(r, g, b);
  }

  void write_facet_end()
  {
    if(m_header.binary())
      I_Binary_write_big_endian_integer32(out(), 0);
    else
      out() << '\n';
  }
};

} //namespace CGAL

#endif // CGAL_IO_OFF_FILE_WRITER_OFF_H
