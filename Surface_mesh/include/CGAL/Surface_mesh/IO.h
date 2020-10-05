// Copyright (C) 2001-2005 by Computer Graphics Group, RWTH Aachen
// Copyright (C) 2011 by Graphics & Geometry Group, Bielefeld University
// Copyright (C) 2014 GeometryFactory
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//

#ifndef CGAL_SURFACE_MESH_IO_H
#define CGAL_SURFACE_MESH_IO_H

#include <CGAL/license/Surface_mesh.h>


#include <CGAL/disable_warnings.h>

//== INCLUDES =================================================================

#include <string>
#include <fstream>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <stdexcept>


#include <boost/array.hpp>

#include <CGAL/assertions.h>
#include <CGAL/use.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Surface_mesh/Properties.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

//=============================================================================

namespace CGAL {

namespace internal {
// helper function
  template <typename T> void read(std::istream& in, T& t)
{
    in.read(reinterpret_cast<char*>(&t), sizeof(t));
}

template <typename Point_3>
bool read_off_binary(Surface_mesh<Point_3>& mesh,
                     std::istream& in,
                     const bool has_normals,
                     const bool has_texcoords)
{
    typedef Surface_mesh<Point_3> Mesh;
    typedef typename Kernel_traits<Point_3>::Kernel K;
    typedef typename K::Vector_3 Vector_3;
    typedef typename K::Vector_2 Vector_2;
    typedef typename K::Vector_3 Normal;
    typedef typename K::Vector_3 Texture_coordinate;

    unsigned int       i, j, idx;
    unsigned int       nV, nF, nE;
    Point_3            p;
    Vector_3           n, c;
    Vector_2           t;
    typename Mesh::Vertex_index  v;

    // properties
    typename Mesh::template Property_map<typename Mesh::Vertex_index, Normal>              normals;
    typename Mesh::template Property_map<typename Mesh::Vertex_index, Texture_coordinate>  texcoords;
    if (has_normals)   normals   = mesh.template add_property_map<typename Mesh::Vertex_index, Normal>("v:normal").first;
    if (has_texcoords) texcoords = mesh.template add_property_map<typename Mesh::Vertex_index, Texture_coordinate>("v:texcoord").first;


    // #Vertice, #Faces, #Edges
    internal::read(in, nV);
    internal::read(in, nF);
    internal::read(in, nE);
    mesh.clear();
    mesh.reserve(nV, (std::max)(3*nV, nE), nF);


    // read vertices: pos [normal] [color] [texcoord]
    for (i=0; i<nV && in.good(); ++i)
    {
        // position
        internal::read(in, p);
        v = mesh.add_vertex(p);

        // normal
        if (has_normals)
        {
            internal::read(in, n);
            normals[v] = n;
        }

        // tex coord
        if (has_texcoords)
        {
            internal::read(in, t);
            texcoords[v] = Vector_3(t[0], t[1], 0.0);
        }
    }


    // read faces: #N v[1] v[2] ... v[n-1]
    std::vector<typename Mesh::Vertex_index> vertices;
    for (i=0; i<nF; ++i)
    {
        internal::read(in, nV);
        vertices.resize(nV);
        for (j=0; j<nV; ++j)
        {
            internal::read(in, idx);
            vertices[j] = typename Mesh::Vertex_index(idx);
        }
        if(!mesh.add_face(vertices).is_valid()) {
          // adding a face did not succeed, stop reading the rest
          return false;
        }
    }

    return true;
}

template <typename Point_3>
bool read_off_ascii(Surface_mesh<Point_3>& mesh,
                    std::istream& in,
                    const bool has_normals,
                    const bool has_texcoords)
{
    typedef Surface_mesh<Point_3> Mesh;
    typedef typename Kernel_traits<Point_3>::Kernel K;
    typedef typename K::Vector_3 Vector_3;
    typedef typename K::Vector_3 Normal;
    typedef typename K::Vector_3 Texture_coordinate;

    boost::array<double, 3> buffer;
    std::string             line;
    unsigned int            i, j, idx;
    unsigned int            nV, nF, nE;
    typename Mesh::Vertex_index   v;

    // properties
    typename Mesh::template Property_map<typename Mesh::Vertex_index, Normal>                 normals;
    typename Mesh::template Property_map<typename Mesh::Vertex_index, Texture_coordinate>     texcoords;

    if (has_normals)   normals   = mesh.template add_property_map<typename Mesh::Vertex_index, Normal>("v:normal").first;
    if (has_texcoords) texcoords = mesh.template add_property_map<typename Mesh::Vertex_index, Texture_coordinate>("v:texcoord").first;

    char c;
    do {
      c = in.get();
      if(c == '#'){
        getline(in,line);
      } else {
        in.putback(c);
        break;
      }
    }while(1);

    // #Vertice, #Faces, #Edges
    in >> nV >> nF >> nE;
    getline(in,line); // reads eol

    mesh.clear();
    mesh.reserve(nV, (std::max)(3*nV, nE), nF);

    // read vertices: pos [normal] [color] [texcoord]
    for (i=0; i<nV && in.good(); ++i)
    {
        // read line
        getline(in, line);
        if(line[0] == '#') // if the first column is a # we are looking at a comment line
        {
          --i;
          continue;
        }

        // position
        std::istringstream iss(line);
        iss >> iformat(buffer[0]) >> iformat(buffer[1]) >> iformat(buffer[2]);
        v = mesh.add_vertex(Point_3(buffer[0], buffer[1], buffer[2]));

        // normal
        if (has_normals)
        {
            iss >> iformat(buffer[0]) >> iformat(buffer[1]) >> iformat(buffer[2]);
        }

        // tex coord
        if (has_texcoords)
        {
            iss >> iformat(buffer[0]) >> iformat(buffer[1]);
            texcoords[v] = Vector_3(buffer[0], buffer[1], 0.0);
        }
    }

    // read faces: #N v[1] v[2] ... v[n-1]
    std::vector<typename Mesh::Vertex_index> vertices;
    for (i=0; i<nF; ++i)
    {
        // read line
        getline(in, line);
        if(line[0] == '#') // if the first column is a # we are looking at a comment line
        {
          --i;
          continue;
        }

        // #vertices
        std::istringstream iss(line);
        iss >> nV;
        vertices.resize(nV);

        // indices
        for (j=0; j<nV; ++j)
        {
           iss >> idx;
            vertices[j] = typename Mesh::Vertex_index(idx);
        }

        if(!mesh.add_face(vertices).is_valid()) {
            // adding a face did not succeed, stop reading the rest
            return false;
        }
    }
    return true;
}
}

#if 0
/// \addtogroup PkgSurfaceMeshIO
///
/// I/O functionality for `Surface_mesh`. The top-level functions
/// `read_mesh()` and `write_mesh()` dispatch on the available readers
/// according to the file extension. Currently only `OFF` files are
/// supported.
///
/// @{

/// This function reads an `OFF` file into a `Surface_mesh`. It
/// supports the `OFF` vertex properties normal, color, and vertex
/// coordinates. If a property is detected in the `OFF` file, it will be
/// read into the `mesh` as a vertex property with the name
/// `"v:normal"`, `"v:color"`, and `"v:texcoord"`, respectivly.
///
/// @param mesh The mesh that should contain the file contents.
/// @param filename The name of the file to be read.
///
/// @returns `true`, if reading succeeded, `false` otherwise
///
#endif
template <typename K>
bool read_off(Surface_mesh<K>& mesh, const std::string& filename)
{
  std::string  line;
    bool  has_texcoords = false;
    bool  has_normals   = false;
    bool  has_hcoords   = false;
    bool  has_dim       = false;
    bool  is_binary     = false;

    // open file (in ASCII mode)
    std::ifstream in(filename.c_str());
    if (!in) return false;

    // read header: [ST][C][N][4][n]OFF BINARY
    std::getline(in,line);
    const char *c = line.c_str();
    if (c[0] == 'S' && c[1] == 'T') { has_texcoords = true; c += 2; }
    if (c[0] == 'N') { has_normals = true; ++c; }
    if (c[0] == '4') { has_hcoords = true; ++c; }
    if (c[0] == 'n') { has_dim     = true; ++c; }
    if (strncmp(c, "OFF", 3) != 0) { in.close(); return false; } // no OFF
    if (strncmp(c+4, "BINARY", 6) == 0) is_binary = true;


    // homogeneous coords, and vertex dimension != 3 are not supported
    if (has_hcoords || has_dim)
    {
        in.close();
        return false;
    }


    // if binary: reopen file in binary mode
    if (is_binary)
    {
      in.close();
      in.open(filename.c_str(), std::ios::binary);
      std::getline(in,line);
    }

    // read as ASCII or binary
    bool ok = (is_binary ?
               internal::read_off_binary(mesh, in, has_normals, has_texcoords) :
               internal::read_off_ascii(mesh, in, has_normals, has_texcoords));

    in.close();
    return ok;
}

#if 0
/// This function writes a `Surface_mesh` into an ASCII `OFF`
/// file. It does not support properties.
///
/// @param mesh The mesh that should be written.
/// @param filename The name of the file to be written.
///
/// @returns `true`, if reading succeeded, `false` otherwise
///
#endif
template <typename K>
bool write_off(const Surface_mesh<K>& mesh, const std::string& filename)
{
    std::ofstream out(filename.c_str());
    if (out.fail())
        return false;

    out << std::setprecision(17);
    write_off(out, mesh);

    return !out.fail();
}

#if 0
/// Read a file into a `Surface_mesh`. The extension of the
/// filename determines which reader is used.
///
/// Mapping from extension to reader:
/// - off/OFF -> `read_off()`
///
/// @param mesh The mesh that should contain the input.
/// @param filename The name of the file to be read.
///
/// @return `true`, if reading succeeded, `false` otherwise
///
#endif
template <typename K>
bool read_mesh(Surface_mesh<K>& mesh, const std::string& filename) {
      // clear mesh before reading from file
    mesh.clear();

    // extract file extension
    std::string::size_type dot(filename.rfind("."));
    if (dot == std::string::npos) return false;
    std::string ext = filename.substr(dot+1, filename.length()-dot-1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    // extension determines reader
    if (ext == "off")
    {
        return read_off(mesh, filename);
    }

    // we didn't find a reader module
    return false;
}

#if 0
/// Write a `Surface_mesh` to a file. The extension of the
/// filename determines which writer is used.
///
/// Mapping from extension to writer:
/// - off/OFF -> `read_off()`
///
/// @param mesh The mesh to be written.
/// @param filename The name of the file to be written.
///
/// @return `true`, if writing succeeded, `false` otherwise
///
#endif
template <typename K>
bool write_mesh(const Surface_mesh<K>& mesh, const std::string& filename)
{
      // extract file extension
    std::string::size_type dot(filename.rfind("."));
    if (dot == std::string::npos) return false;
    std::string ext = filename.substr(dot+1, filename.length()-dot-1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    // extension determines reader
    if (ext == "off")
    {
        return write_off(mesh, filename);
    }

    // we didn't find a writer module
    return false;
}

/// group io
template <class P, class Writer>
void
generic_print_surface_mesh( std::ostream&   out,
                          const Surface_mesh<P>&       M,
                          Writer&           writer) {
  // writes M to `out' in the format provided by `writer'.
  typedef typename boost::graph_traits<Surface_mesh<P> >::vertex_iterator VCI;
  typedef typename boost::graph_traits<Surface_mesh<P> >::face_iterator   FCI;
  typedef typename Surface_mesh<P>::Halfedge_around_face_circulator            HFCC;
  typedef typename boost::property_map<Surface_mesh<P>,CGAL::vertex_point_t>::type VPmap;
  VPmap map = get(CGAL::vertex_point, M);
  // Print header.
  writer.write_header( out,
                       num_vertices(M),
                       num_halfedges(M),
                       num_faces(M));

  std::map<typename Surface_mesh<P>::vertex_index, std::size_t> index_map;
  typename std::map<typename Surface_mesh<P>::vertex_index, std::size_t>::iterator hint = index_map.begin();
  std::size_t id = 0;

  for( VCI vi = vertices(M).begin(); vi != vertices(M).end(); ++vi) {
    writer.write_vertex( ::CGAL::to_double( get(map, *vi).x()),
                         ::CGAL::to_double( get(map, *vi).y()),
                         ::CGAL::to_double( get(map, *vi).z()));

    hint = index_map.insert(hint, std::make_pair(*vi, id++));
  }
  typedef typename boost::property_traits<VPmap>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel K;
  typename Surface_mesh<P>::template Property_map<typename Surface_mesh<P>::Vertex_index, typename K::Vector_3 > vnormals;
  bool has_normals = false;
  boost::tie(vnormals, has_normals) = M.template property_map<typename Surface_mesh<P>::Vertex_index, typename K::Vector_3>("v:normal");
  if(has_normals)
  {
    for( VCI vi = vertices(M).begin(); vi != vertices(M).end(); ++vi) {
      writer.write_vertex_normal( ::CGAL::to_double( get(vnormals, *vi).x()),
                                  ::CGAL::to_double( get(vnormals, *vi).y()),
                                  ::CGAL::to_double( get(vnormals, *vi).z()));
    }
  }

  writer.write_facet_header();
  for( FCI fi = faces(M).begin(); fi != faces(M).end(); ++fi) {
    HFCC hc(halfedge(*fi, M), M);
    HFCC hc_end = hc;
    std::size_t n = circulator_size( hc);
    CGAL_assertion( n >= 3);
    writer.write_facet_begin( n);
    do {
      writer.write_facet_vertex_index(index_map[target(*hc, M)]);
      ++hc;
    } while( hc != hc_end);
    writer.write_facet_end();
  }
  writer.write_footer();
}
} // CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESH_IO_H

