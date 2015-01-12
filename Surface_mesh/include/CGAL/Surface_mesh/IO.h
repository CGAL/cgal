//=============================================================================
// Copyright (C) 2001-2005 by Computer Graphics Group, RWTH Aachen
// Copyright (C) 2011 by Graphics & Geometry Group, Bielefeld University
// Copyright (C) 2014 GeometryFactory
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

#ifndef CGAL_SURFACE_MESH_IO_H
#define CGAL_SURFACE_MESH_IO_H


//== INCLUDES =================================================================

#include <string>
#include <cstdio>
#include <cstring>


#include <boost/array.hpp>

#include <CGAL/assertions.h>
#include <CGAL/use.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Surface_mesh/Properties.h>

//=============================================================================

namespace CGAL {

namespace internal {
// helper function
template <typename T> void read(FILE* in, T& t)
{
  std::size_t err = 0;
    err = fread(&t, 1, sizeof(t), in);
    if(err != 0)
      throw std::runtime_error("fread error");
}

template <typename Point_3>
bool read_off_binary(Surface_mesh<Point_3>& mesh,
                     FILE* in,
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
    typename Mesh:: template Property_map<typename Mesh::Vertex_index, Normal>              normals;
    typename Mesh:: template Property_map<typename Mesh::Vertex_index, Texture_coordinate>  texcoords;
    if (has_normals)   normals   = mesh.template add_property_map<typename Mesh::Vertex_index, Normal>("v:normal").first;
    if (has_texcoords) texcoords = mesh.template add_property_map<typename Mesh::Vertex_index, Texture_coordinate>("v:texcoord").first;


    // #Vertice, #Faces, #Edges
    internal::read(in, nV);
    internal::read(in, nF);
    internal::read(in, nE);
    mesh.clear();
    mesh.reserve(nV, (std::max)(3*nV, nE), nF);


    // read vertices: pos [normal] [color] [texcoord]
    for (i=0; i<nV && !feof(in); ++i)
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
                    FILE* in,
                    const bool has_normals,
                    const bool has_texcoords)
{
    typedef Surface_mesh<Point_3> Mesh;
    typedef typename Kernel_traits<Point_3>::Kernel K;
    typedef typename K::Vector_3 Vector_3;
    typedef typename K::Vector_3 Normal;
    typedef typename K::Vector_3 Texture_coordinate;

    boost::array<double, 3> buffer;
    char                    line[100], *lp;
    unsigned int            i, j, items, idx;
    int                     nc;
    unsigned int            nV, nF, nE;
    typename Mesh::Vertex_index   v;

    // properties
    typename Mesh::template Property_map<typename Mesh::Vertex_index, Normal>                 normals;
    typename Mesh::template Property_map<typename Mesh::Vertex_index, Texture_coordinate>     texcoords;
    
    if (has_normals)   normals   = mesh.template add_property_map<typename Mesh::Vertex_index, Normal>("v:normal").first;
    if (has_texcoords) texcoords = mesh.template add_property_map<typename Mesh::Vertex_index, Texture_coordinate>("v:texcoord").first;

    int c;
    do {
      c = getc(in);
      if(c == '#'){
        fgets(line, 100, in);
      } else {
        ungetc(c,in);
        break;
      }
    }while(1);

    // #Vertice, #Faces, #Edges
    items = fscanf(in, "%d %d %d\n", (int*)&nV, (int*)&nF, (int*)&nE);

    mesh.clear();
    mesh.reserve(nV, (std::max)(3*nV, nE), nF);

    // read vertices: pos [normal] [color] [texcoord]
    for (i=0; i<nV && !feof(in); ++i)
    {
        // read line
        lp = fgets(line, 100, in);
        if(line[0] == '#') // if the first column is a # we are looking at a comment line
        {
          --i;
          continue;
        }

        // position
        items = sscanf(lp, "%lf %lf %lf%n", &(buffer[0]), &buffer[1], &buffer[2], &nc);
        CGAL_assertion(items==3);
        v = mesh.add_vertex(Point_3(buffer[0], buffer[1], buffer[2]));
        lp += nc;

        // normal
        if (has_normals)
        {
            if (sscanf(lp, "%lf %lf %lf%n", &buffer[0], &buffer[1], &buffer[2], &nc) == 3)
            {
              normals[v] = Vector_3(buffer[0], buffer[1], buffer[2]);
            }
            lp += nc;
        }

        // tex coord
        if (has_texcoords)
        {
            items = sscanf(lp, "%lf %lf%n", &buffer[0], &buffer[1], &nc);
            CGAL_assertion(items == 2);
            texcoords[v] = Vector_3(buffer[0], buffer[1], 0.0);
            lp += nc;
        }
    }

    // read faces: #N v[1] v[2] ... v[n-1]
    std::vector<typename Mesh::Vertex_index> vertices;
    for (i=0; i<nF; ++i)
    {
        // read line
        lp = fgets(line, 100, in);
        if(line[0] == '#') // if the first column is a # we are looking at a comment line
        {
          --i;
          continue;
        }

        // #vertices
        items = sscanf(lp, "%d%n", (int*)&nV, &nc);
        CGAL_assertion(items == 1);
        vertices.resize(nV);
        lp += nc;

        // indices
        for (j=0; j<nV; ++j)
        {
            items = sscanf(lp, "%d%n", (int*)&idx, &nc);
            CGAL_assertion(items == 1);
            vertices[j] = typename Mesh::Vertex_index(idx);
            lp += nc;
        }

        if(!mesh.add_face(vertices).is_valid()) {
            // adding a face did not succeed, stop reading the rest
            return false;
        }
    }
    CGAL_USE(items);
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
    char  line[100];
    bool  has_texcoords = false;
    bool  has_normals   = false;
    bool  has_hcoords   = false;
    bool  has_dim       = false;
    bool  is_binary     = false;

    // open file (in ASCII mode)
    FILE* in = std::fopen(filename.c_str(), "r");
    if (!in) return false;

    // read header: [ST][C][N][4][n]OFF BINARY
    char *c = std::fgets(line, 100, in);
    CGAL_assertion(c != NULL);
    c = line;
    if (c[0] == 'S' && c[1] == 'T') { has_texcoords = true; c += 2; }
    if (c[0] == 'N') { has_normals = true; ++c; }
    if (c[0] == '4') { has_hcoords = true; ++c; }
    if (c[0] == 'n') { has_dim     = true; ++c; }
    if (strncmp(c, "OFF", 3) != 0) { std::fclose(in); return false; } // no OFF
    if (strncmp(c+4, "BINARY", 6) == 0) is_binary = true;


    // homogeneous coords, and vertex dimension != 3 are not supported
    if (has_hcoords || has_dim)
    {
        std::fclose(in);
        return false;
    }


    // if binary: reopen file in binary mode
    if (is_binary)
    {
        std::fclose(in);
        in = std::fopen(filename.c_str(), "rb");
        c = std::fgets(line, 100, in);
        CGAL_assertion(c != NULL);
    }

    // read as ASCII or binary
    bool ok = (is_binary ?
               internal::read_off_binary(mesh, in, has_normals, has_texcoords) :
               internal::read_off_ascii(mesh, in, has_normals, has_texcoords));


    std::fclose(in);
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
    typedef Surface_mesh<K> Mesh;
    typedef typename Mesh::Point Point_3;

    FILE* out = fopen(filename.c_str(), "w");
    if (!out)
        return false;


    // header
    fprintf(out, "OFF\n%d %d 0\n", mesh.num_vertices(), mesh.num_faces());


    // vertices
    typename Mesh::template Property_map<typename Mesh::Vertex_index, Point_3> points 
      = mesh.template property_map<typename Mesh::Vertex_index, Point_3>("v:point").first;
    for (typename Mesh::Vertex_iterator vit=mesh.vertices_begin(); vit!=mesh.vertices_end(); ++vit)
    {
        const Point_3& p = points[*vit];
        fprintf(out, "%.10f %.10f %.10f\n", p[0], p[1], p[2]);
    }


    // faces
    for (typename Surface_mesh<K>::Face_iterator fit=mesh.faces_begin(); fit!=mesh.faces_end(); ++fit)
    {
        int nV = mesh.degree(*fit);
        fprintf(out, "%d", nV);
        typename Surface_mesh<K>::Vertex_around_face_circulator fvit(mesh.halfedge(*fit),mesh), fvend=fvit;
        do
        {
          typename Surface_mesh<K>::size_type idx = *fvit;
          fprintf(out, " %d", idx);
        }
        while (++fvit != fvend);
        fprintf(out, "\n");
    }

    fclose(out);
    return true;
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
    std::transform(ext.begin(), ext.end(), ext.begin(), tolower);

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
    std::transform(ext.begin(), ext.end(), ext.begin(), tolower);

    // extension determines reader
    if (ext == "off")
    {
        return write_off(mesh, filename);
    }

    // we didn't find a writer module
    return false;
}

/// group io
/// @}

} // CGAL


#endif // CGAL_SURFACE_MESH_IO_H

