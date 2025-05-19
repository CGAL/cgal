//=============================================================================
// Copyright (C) 2001-2005 by Computer Graphics Group, RWTH Aachen
// Copyright (C) 2011 by Graphics & Geometry Group, Bielefeld University
//
// SPDX-License-Identifier: GPL-2.0-only
//
//=============================================================================


//== INCLUDES =================================================================

#include "IO.h"
#include <stdio.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

//== IMPLEMENTATION ===========================================================


// helper function
template <typename T> void read(FILE* in, T& t)
{
    int err = 0;
    err = fread(&t, 1, sizeof(t), in);
}


//-----------------------------------------------------------------------------

template <typename NamedParameters>
bool read_off_ascii(Surface_mesh& mesh,
                    FILE* in,
                    const bool has_normals,
                    const bool has_texcoords,
                    const bool has_colors,
                    NamedParameters& np)
{
    char                 line[100], *lp;
    unsigned int         i, j, items, idx, nc;
    unsigned int         nV, nF, nE;
    Vec3f                p, n, c;
    Vec2f                t;
    Surface_mesh::Vertex v;
    typename CGAL::GetVertexPointMap<Surface_mesh, NamedParameters>::const_type
        vpm = choose_parameter(get_parameter(np, CGAL::internal_np::vertex_point),
                           get_const_property_map(CGAL::vertex_point, mesh));


    // properties
    Surface_mesh::Vertex_property<Normal>              normals;
    Surface_mesh::Vertex_property<Texture_coordinate>  texcoords;
    Surface_mesh::Vertex_property<Color>               colors;
    if (has_normals)   normals   = mesh.vertex_property<Normal>("v:normal");
    if (has_texcoords) texcoords = mesh.vertex_property<Texture_coordinate>("v:texcoord");
    if (has_colors)    colors    = mesh.vertex_property<Color>("v:color");


    // #Vertice, #Faces, #Edges
    items = fscanf(in, "%d %d %d\n", (int*)&nV, (int*)&nF, (int*)&nE);
    mesh.clear();
    mesh.reserve(nV, (std::max)(3*nV, nE), nF);


    // read vertices: pos [normal] [color] [texcoord]
    for (i=0; i<nV && !feof(in); ++i)
    {
        // read line
        lp = fgets(line, 100, in);

        // position
        items = sscanf(lp, "%f %f %f%n", &p[0], &p[1], &p[2], &nc);
        assert(items==3);
        Surface_mesh::Vertex v = mesh.add_vertex();
        put(vpm, v, (Point)p);
        lp += nc;

        // normal
        if (has_normals)
        {
            if (sscanf(lp, "%f %f %f%n", &n[0], &n[1], &n[2], &nc) == 3)
            {
                normals[v] = n;
            }
            lp += nc;
        }

        // color
        if (has_colors)
        {
            if (sscanf(lp, "%f %f %f%n", &c[0], &c[1], &c[2], &nc) == 3)
            {
                if (c[0]>1.0f || c[1]>1.0f || c[2]>1.0f) c *= (1.0/255.0);
                colors[v] = c;
            }
            lp += nc;
        }

        // tex coord
        if (has_texcoords)
        {
            items = sscanf(lp, "%f %f%n", &t[0], &t[1], &nc);
            assert(items == 2);
            texcoords[v][0] = t[0];
            texcoords[v][1] = t[1];
            lp += nc;
        }
    }



    // read faces: #N v[1] v[2] ... v[n-1]
    std::vector<Surface_mesh::Vertex> vertices;
    for (i=0; i<nF; ++i)
    {
        // read line
        lp = fgets(line, 100, in);

        // #vertices
        items = sscanf(lp, "%d%n", (int*)&nV, &nc);
        assert(items == 1);
        vertices.resize(nV);
        lp += nc;

        // indices
        for (j=0; j<nV; ++j)
        {
            items = sscanf(lp, "%d%n", (int*)&idx, &nc);
            assert(items == 1);
            vertices[j] = Surface_mesh::Vertex(idx);
            lp += nc;
        }
        mesh.add_face(vertices);
    }


    return true;
}


//-----------------------------------------------------------------------------

template<typename NamedParameters>
bool read_off_binary(Surface_mesh& mesh,
                     FILE* in,
                     const bool has_normals,
                     const bool has_texcoords,
                     const bool has_colors,
                     NamedParameters& np)
{
    unsigned int       i, j, idx;
    unsigned int       nV, nF, nE;
    Vec3f              p, n, c;
    Vec2f              t;
    Surface_mesh::Vertex  v;


    // binary cannot (yet) read colors
    if (has_colors) return false;


    // properties
    Surface_mesh::Vertex_property<Normal>              normals;
    Surface_mesh::Vertex_property<Texture_coordinate>  texcoords;
    if (has_normals)   normals   = mesh.vertex_property<Normal>("v:normal");
    if (has_texcoords) texcoords = mesh.vertex_property<Texture_coordinate>("v:texcoord");
    typename CGAL::GetVertexPointMap<Surface_mesh, NamedParameters>::const_type
        vpm = choose_parameter(get_parameter(np, CGAL::internal_np::vertex_point),
                           get_const_property_map(CGAL::vertex_point, mesh));


    // #Vertice, #Faces, #Edges
    read(in, nV);
    read(in, nF);
    read(in, nE);
    mesh.clear();
    mesh.reserve(nV, (std::max)(3*nV, nE), nF);


    // read vertices: pos [normal] [color] [texcoord]
    for (i=0; i<nV && !feof(in); ++i)
    {
        // position
        read(in, p);
        Surface_mesh::Vertex v = mesh.add_vertex();
        put(vpm, v, (Point)p);

        // normal
        if (has_normals)
        {
            read(in, n);
            normals[v] = n;
        }

        // tex coord
        if (has_texcoords)
        {
            read(in, t);
            texcoords[v][0] = t[0];
            texcoords[v][1] = t[1];
        }
    }


    // read faces: #N v[1] v[2] ... v[n-1]
    std::vector<Surface_mesh::Vertex> vertices;
    for (i=0; i<nF; ++i)
    {
        read(in, nV);
        vertices.resize(nV);
        for (j=0; j<nV; ++j)
        {
            read(in, idx);
            vertices[j] = Surface_mesh::Vertex(idx);
        }
        mesh.add_face(vertices);
    }


    return true;
}


//-----------------------------------------------------------------------------

template<typename NamedParameters>
bool read_off(Surface_mesh& mesh, const std::string& filename, NamedParameters& np)
{
    char  line[100];
    bool  has_texcoords = false;
    bool  has_normals   = false;
    bool  has_colors    = false;
    bool  has_hcoords   = false;
    bool  has_dim       = false;
    bool  is_binary     = false;


    // open file (in ASCII mode)
    FILE* in = fopen(filename.c_str(), "r");
    if (!in) return false;


    // read header: [ST][C][N][4][n]OFF BINARY
    char *c = fgets(line, 100, in);
    assert(c != NULL);
    c = line;
    if (c[0] == 'S' && c[1] == 'T') { has_texcoords = true; c += 2; }
    if (c[0] == 'C') { has_colors  = true; ++c; }
    if (c[0] == 'N') { has_normals = true; ++c; }
    if (c[0] == '4') { has_hcoords = true; ++c; }
    if (c[0] == 'n') { has_dim     = true; ++c; }
    if (strncmp(c, "OFF", 3) != 0) { fclose(in); return false; } // no OFF
    if (strncmp(c+4, "BINARY", 6) == 0) is_binary = true;


    // homogeneous coords, and vertex dimension != 3 are not supported
    if (has_hcoords || has_dim)
    {
        fclose(in);
        return false;
    }


    // if binary: reopen file in binary mode
    if (is_binary)
    {
        fclose(in);
        in = fopen(filename.c_str(), "rb");
        c = fgets(line, 100, in);
        assert(c != NULL);
    }


    // read as ASCII or binary
    bool ok = (is_binary ?
               read_off_binary(mesh, in, has_normals, has_texcoords, has_colors, np) :
               read_off_ascii(mesh, in, has_normals, has_texcoords, has_colors, np));


    fclose(in);
    return ok;
}


//-----------------------------------------------------------------------------


bool write_off(const Surface_mesh& mesh, const std::string& filename)
{
    FILE* out = fopen(filename.c_str(), "w");
    if (!out)
        return false;


    // header
    fprintf(out, "OFF\n%d %d 0\n", mesh.n_vertices(), mesh.n_faces());


    // vertices
    Surface_mesh::Vertex_property<Point> points = mesh.get_vertex_property<Point>("v:point");
    for (Surface_mesh::Vertex_iterator vit=mesh.vertices_begin(); vit!=mesh.vertices_end(); ++vit)
    {
        const Point& p = points[*vit];
        fprintf(out, "%.10f %.10f %.10f\n", p[0], p[1], p[2]);
    }


    // faces
    for (Surface_mesh::Face_iterator fit=mesh.faces_begin(); fit!=mesh.faces_end(); ++fit)
    {
        int nV = mesh.valence(*fit);
        fprintf(out, "%d", nV);
        Surface_mesh::Vertex_around_face_circulator fvit=mesh.vertices(*fit), fvend=fvit;
        do
        {
            fprintf(out, " %d", (*fvit).idx());
        }
        while (++fvit != fvend);
        fprintf(out, "\n");
    }

    fclose(out);
    return true;
}


//=============================================================================
