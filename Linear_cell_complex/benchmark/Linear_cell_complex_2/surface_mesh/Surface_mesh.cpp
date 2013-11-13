//=============================================================================
// Copyright (C) 2001-2005 by Computer Graphics Group, RWTH Aachen
// Copyright (C) 2011 by Graphics & Geometry Group, Bielefeld University
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public License
// as published by the Free Software Foundation, version 2.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//=============================================================================


//== INCLUDES =================================================================


#include "Surface_mesh.h"
#include "IO.h"


//== IMPLEMENTATION ===========================================================


Surface_mesh::
Surface_mesh()
{
    // allocate standard properties
    // same list is used in operator=() and assign()
    vconn_    = add_vertex_property<Vertex_connectivity>("v:connectivity");
    hconn_    = add_halfedge_property<Halfedge_connectivity>("h:connectivity");
    fconn_    = add_face_property<Face_connectivity>("f:connectivity");
    vpoint_   = add_vertex_property<Point>("v:point");
    vdeleted_ = add_vertex_property<bool>("v:deleted", false);
    edeleted_ = add_edge_property<bool>("e:deleted", false);
    fdeleted_ = add_face_property<bool>("f:deleted", false);

    deleted_vertices_ = deleted_edges_ = deleted_faces_ = 0;
    garbage_ = false;
}


//-----------------------------------------------------------------------------


Surface_mesh::
~Surface_mesh()
{
}


//-----------------------------------------------------------------------------


Surface_mesh&
Surface_mesh::
operator=(const Surface_mesh& rhs)
{
    if (this != &rhs)
    {
        // deep copy of property containers
        vprops_ = rhs.vprops_;
        hprops_ = rhs.hprops_;
        eprops_ = rhs.eprops_;
        fprops_ = rhs.fprops_;

        // property handles contain pointers, have to be reassigned
        vconn_    = vertex_property<Vertex_connectivity>("v:connectivity");
        hconn_    = halfedge_property<Halfedge_connectivity>("h:connectivity");
        fconn_    = face_property<Face_connectivity>("f:connectivity");
        vdeleted_ = vertex_property<bool>("v:deleted");
        edeleted_ = edge_property<bool>("e:deleted");
        fdeleted_ = face_property<bool>("f:deleted");
        vpoint_   = vertex_property<Point>("v:point");

        // how many elements are deleted?
        deleted_vertices_ = rhs.deleted_vertices_;
        deleted_edges_    = rhs.deleted_edges_;
        deleted_faces_    = rhs.deleted_faces_;
        garbage_          = rhs.garbage_;
    }

    return *this;
}


//-----------------------------------------------------------------------------


Surface_mesh&
Surface_mesh::
assign(const Surface_mesh& rhs)
{
    if (this != &rhs)
    {
        // clear properties
        vprops_.clear();
        hprops_.clear();
        eprops_.clear();
        fprops_.clear();

        // allocate standard properties
        vconn_    = add_vertex_property<Vertex_connectivity>("v:connectivity");
        hconn_    = add_halfedge_property<Halfedge_connectivity>("h:connectivity");
        fconn_    = add_face_property<Face_connectivity>("f:connectivity");
        vpoint_   = add_vertex_property<Point>("v:point");
        vdeleted_ = add_vertex_property<bool>("v:deleted", false);
        edeleted_ = add_edge_property<bool>("e:deleted", false);
        fdeleted_ = add_face_property<bool>("f:deleted", false);

        // copy properties from other mesh
        vconn_.array()     = rhs.vconn_.array();
        hconn_.array()     = rhs.hconn_.array();
        fconn_.array()     = rhs.fconn_.array();
        vpoint_.array()    = rhs.vpoint_.array();
        vdeleted_.array()  = rhs.vdeleted_.array();
        edeleted_.array()  = rhs.edeleted_.array();
        fdeleted_.array()  = rhs.fdeleted_.array();

        // resize (needed by property containers)
        vprops_.resize(rhs.vertices_size());
        hprops_.resize(rhs.halfedges_size());
        eprops_.resize(rhs.edges_size());
        fprops_.resize(rhs.faces_size());

        // how many elements are deleted?
        deleted_vertices_ = rhs.deleted_vertices_;
        deleted_edges_    = rhs.deleted_edges_;
        deleted_faces_    = rhs.deleted_faces_;
        garbage_          = rhs.garbage_;
    }

    return *this;
}


//-----------------------------------------------------------------------------


bool
Surface_mesh::
read(const std::string& filename)
{
    return read_mesh(*this, filename);
}


//-----------------------------------------------------------------------------


bool
Surface_mesh::
write(const std::string& filename) const
{
    return write_mesh(*this, filename);
}


//-----------------------------------------------------------------------------


void
Surface_mesh::
clear()
{
    vprops_.resize(0);
    hprops_.resize(0);
    eprops_.resize(0);
    fprops_.resize(0);

    vprops_.free_memory();
    hprops_.free_memory();
    eprops_.free_memory();
    fprops_.free_memory();

    deleted_vertices_ = deleted_edges_ = deleted_faces_ = 0;
    garbage_ = false;
}


//-----------------------------------------------------------------------------


void
Surface_mesh::
free_memory()
{
    vprops_.free_memory();
    hprops_.free_memory();
    eprops_.free_memory();
    fprops_.free_memory();
}


//-----------------------------------------------------------------------------


void
Surface_mesh::
reserve(unsigned int nvertices,
        unsigned int nedges,
        unsigned int nfaces )
{
    vprops_.reserve(nvertices);
    hprops_.reserve(2*nedges);
    eprops_.reserve(nedges);
    fprops_.reserve(nfaces);
}


//-----------------------------------------------------------------------------


void
Surface_mesh::
property_stats() const
{
    std::vector<std::string> props;

    std::cout << "vertex properties:\n";
    props = vertex_properties();
    for (unsigned int i=0; i<props.size(); ++i)
        std::cout << "\t" << props[i] << std::endl;

    std::cout << "halfedge properties:\n";
    props = halfedge_properties();
    for (unsigned int i=0; i<props.size(); ++i)
        std::cout << "\t" << props[i] << std::endl;

    std::cout << "edge properties:\n";
    props = edge_properties();
    for (unsigned int i=0; i<props.size(); ++i)
        std::cout << "\t" << props[i] << std::endl;

    std::cout << "face properties:\n";
    props = face_properties();
    for (unsigned int i=0; i<props.size(); ++i)
        std::cout << "\t" << props[i] << std::endl;
}


//-----------------------------------------------------------------------------


Surface_mesh::Vertex
Surface_mesh::
add_vertex(const Point& p)
{
    Vertex v = new_vertex();
    vpoint_[v] = p;
    return v;
}


//-----------------------------------------------------------------------------


Surface_mesh::Halfedge
Surface_mesh::
find_halfedge(Vertex start, Vertex end) const
{
    assert(is_valid(start) && is_valid(end));

    Halfedge h  = halfedge(start);
    const Halfedge hh = h;

    if (h.is_valid())
    {
        do
        {
            if (to_vertex(h) == end)
                return h;
            h = cw_rotated_halfedge(h);
        }
        while (h != hh);
    }

    return Halfedge();
}


//-----------------------------------------------------------------------------


void
Surface_mesh::
adjust_outgoing_halfedge(Vertex v)
{
    Halfedge h  = halfedge(v);
    const Halfedge hh = h;

    if (h.is_valid())
    {
        do
        {
            if (is_boundary(h))
            {
                set_halfedge(v, h);
                return;
            }
            h = cw_rotated_halfedge(h);
        }
        while (h != hh);
    }
}


//-----------------------------------------------------------------------------


Surface_mesh::Face
Surface_mesh::
add_triangle(Vertex v0, Vertex v1, Vertex v2)
{
    std::vector<Vertex> v(3);
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
    return add_face(v);
}


//-----------------------------------------------------------------------------


Surface_mesh::Face
Surface_mesh::
add_quad(Vertex v0, Vertex v1, Vertex v2, Vertex v3)
{
    std::vector<Vertex> v(4);
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
    v[3] = v3;
    return add_face(v);
}


//-----------------------------------------------------------------------------


Surface_mesh::Face
Surface_mesh::
add_face(const std::vector<Vertex>& vertices)
{
    Vertex                   v;
    unsigned int             i, ii, n((int)vertices.size()), id;
    std::vector<Halfedge>    halfedges(n);
    std::vector<bool>        is_new(n), needs_adjust(n, false);
    Halfedge                 inner_next, inner_prev,
    outer_next, outer_prev,
    boundary_next, boundary_prev,
    patch_start, patch_end;

    // cache for set_next_halfedge and vertex' set_halfedge
    typedef std::pair<Halfedge, Halfedge>  NextCacheEntry;
    typedef std::vector<NextCacheEntry>    NextCache;

    NextCache    next_cache;
    next_cache.reserve(3*n);


    // don't allow degenerated faces
    assert (n > 2);


    // test for topological errors
    for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
    {
        if ( !is_boundary(vertices[i]) )
        {
            std::cerr << "Surface_meshT::add_face: complex vertex\n";
            return Face();
        }

        halfedges[i] = find_halfedge(vertices[i], vertices[ii]);
        is_new[i]    = !halfedges[i].is_valid();

        if (!is_new[i] && !is_boundary(halfedges[i]))
        {
            std::cerr << "Surface_meshT::add_face: complex edge\n";
            return Face();
        }
    }


    // re-link patches if necessary
    for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
    {
        if (!is_new[i] && !is_new[ii])
        {
            inner_prev = halfedges[i];
            inner_next = halfedges[ii];

            if (next_halfedge(inner_prev) != inner_next)
            {
                // here comes the ugly part... we have to relink a whole patch

                // search a free gap
                // free gap will be between boundary_prev and boundary_next
                outer_prev = opposite_halfedge(inner_next);
                outer_next = opposite_halfedge(inner_prev);
                boundary_prev = outer_prev;
                do
                    boundary_prev = opposite_halfedge(next_halfedge(boundary_prev));
                while (!is_boundary(boundary_prev) || boundary_prev==inner_prev);
                boundary_next = next_halfedge(boundary_prev);
                assert(is_boundary(boundary_prev));
                assert(is_boundary(boundary_next));


                // ok ?
                if (boundary_next == inner_next)
                {
                    std::cerr << "Surface_meshT::add_face: patch re-linking failed\n";
                    return Face();
                }

                // other halfedges' handles
                patch_start = next_halfedge(inner_prev);
                patch_end   = prev_halfedge(inner_next);

                // relink
                next_cache.push_back(NextCacheEntry(boundary_prev, patch_start));
                next_cache.push_back(NextCacheEntry(patch_end, boundary_next));
                next_cache.push_back(NextCacheEntry(inner_prev, inner_next));
            }
        }
    }



    // create missing edges
    for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
        if (is_new[i])
            halfedges[i] = new_edge(vertices[i], vertices[ii]);



    // create the face
    Face f(new_face());
    set_halfedge(f, halfedges[n-1]);



    // setup halfedges
    for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
    {
        v          = vertices[ii];
        inner_prev = halfedges[i];
        inner_next = halfedges[ii];

        id = 0;
        if (is_new[i])  id |= 1;
        if (is_new[ii]) id |= 2;

        if (id)
        {
            outer_prev = opposite_halfedge(inner_next);
            outer_next = opposite_halfedge(inner_prev);

            // set outer links
            switch (id)
            {
                case 1: // prev is new, next is old
                    boundary_prev = prev_halfedge(inner_next);
                    next_cache.push_back(NextCacheEntry(boundary_prev, outer_next));
                    set_halfedge(v, outer_next);
                    break;

                case 2: // next is new, prev is old
                    boundary_next = next_halfedge(inner_prev);
                    next_cache.push_back(NextCacheEntry(outer_prev, boundary_next));
                    set_halfedge(v, boundary_next);
                    break;

                case 3: // both are new
                    if (!halfedge(v).is_valid())
                    {
                        set_halfedge(v, outer_next);
                        next_cache.push_back(NextCacheEntry(outer_prev, outer_next));
                    }
                    else
                    {
                        boundary_next = halfedge(v);
                        boundary_prev = prev_halfedge(boundary_next);
                        next_cache.push_back(NextCacheEntry(boundary_prev, outer_next));
                        next_cache.push_back(NextCacheEntry(outer_prev, boundary_next));
                    }
                    break;
            }

            // set inner link
            next_cache.push_back(NextCacheEntry(inner_prev, inner_next));
        }
        else needs_adjust[ii] = (halfedge(v) == inner_next);


        // set face handle
        set_face(halfedges[i], f);
    }



    // process next halfedge cache
    NextCache::const_iterator ncIt(next_cache.begin()), ncEnd(next_cache.end());
    for (; ncIt != ncEnd; ++ncIt)
        set_next_halfedge(ncIt->first, ncIt->second);



    // adjust vertices' halfedge handle
    for (i=0; i<n; ++i)
        if (needs_adjust[i])
            adjust_outgoing_halfedge(vertices[i]);


    return f;
}


//-----------------------------------------------------------------------------


unsigned int
Surface_mesh::
valence(Vertex v) const
{
    unsigned int count(0);

    Vertex_around_vertex_circulator vvit = vertices(v);
    Vertex_around_vertex_circulator vvend = vvit;
    if (vvit) do
    {
        ++count;
    } while (++vvit != vvend);

    return count;
}


//-----------------------------------------------------------------------------


unsigned int
Surface_mesh::
valence(Face f) const
{
    unsigned int count(0);

    Vertex_around_face_circulator fvit = vertices(f);
    Vertex_around_face_circulator fvend = fvit;
    do {
        ++count;
    } while (++fvit != fvend);

    return count;
}


//-----------------------------------------------------------------------------


bool
Surface_mesh::
is_triangle_mesh() const
{
    Face_iterator fit=faces_begin(), fend=faces_end();
    for (; fit!=fend; ++fit)
        if (valence(*fit) != 3)
            return false;

    return true;
}


//-----------------------------------------------------------------------------


void
Surface_mesh::
triangulate()
{
    /* The iterators will stay valid, even though new faces are added,
     because they are now implemented index-based instead of
     pointer-based.
     */
    Face_iterator fit=faces_begin(), fend=faces_end();
    for (; fit!=fend; ++fit)
        triangulate(*fit);
}


//-----------------------------------------------------------------------------


void
Surface_mesh::
triangulate(Face f)
{
    /*
     Split an arbitrary face into triangles by connecting
     each vertex of fh after its second to vh.

     - fh will remain valid (it will become one of the
     triangles)
     - the halfedge handles of the new triangles will
     point to the old halfedges
     */

    Halfedge base_h  = halfedge(f);
    Vertex   start_v = from_vertex(base_h);
    Halfedge next_h  = next_halfedge(base_h);

    while (to_vertex(next_halfedge(next_h)) != start_v)
    {
        Halfedge next_next_h(next_halfedge(next_h));

        Face new_f = new_face();
        set_halfedge(new_f, base_h);

        Halfedge new_h = new_edge(to_vertex(next_h), start_v);

        set_next_halfedge(base_h, next_h);
        set_next_halfedge(next_h, new_h);
        set_next_halfedge(new_h,  base_h);

        set_face(base_h, new_f);
        set_face(next_h, new_f);
        set_face(new_h,  new_f);

        base_h = opposite_halfedge(new_h);
        next_h = next_next_h;
    }
    set_halfedge(f, base_h);  //the last face takes the handle _fh

    set_next_halfedge(base_h, next_h);
    set_next_halfedge(next_halfedge(next_h), base_h);

    set_face(base_h, f);
}


//-----------------------------------------------------------------------------


void
Surface_mesh::
update_face_normals()
{
    if (!fnormal_)
        fnormal_ = face_property<Point>("f:normal");

    Face_iterator fit, fend=faces_end();

    for (fit=faces_begin(); fit!=fend; ++fit)
        fnormal_[*fit] = compute_face_normal(*fit);
}


//-----------------------------------------------------------------------------


Normal
Surface_mesh::
compute_face_normal(Face f) const
{
    Halfedge h = halfedge(f);
    Halfedge hend = h;

    Point p0 = vpoint_[to_vertex(h)];
    h = next_halfedge(h);
    Point p1 = vpoint_[to_vertex(h)];
    h = next_halfedge(h);
    Point p2 = vpoint_[to_vertex(h)];

    if (next_halfedge(h) == hend) // face is a triangle
    {
        return cross(p2-=p1, p0-=p1).normalize();
    }

    else // face is a general polygon
    {
        Normal n(0,0,0);

        hend = h;
        do
        {
            n += cross(p2-p1, p0-p1);
            h  = next_halfedge(h);
            p0 = p1;
            p1 = p2;
            p2 = vpoint_[to_vertex(h)];
        }
        while (h != hend);

        return n.normalize();
    }
}


//-----------------------------------------------------------------------------


void
Surface_mesh::
update_vertex_normals()
{
    if (!vnormal_)
        vnormal_ = vertex_property<Point>("v:normal");

    Vertex_iterator vit, vend=vertices_end();

    for (vit=vertices_begin(); vit!=vend; ++vit)
        vnormal_[*vit] = compute_vertex_normal(*vit);
}


//-----------------------------------------------------------------------------


Normal
Surface_mesh::
compute_vertex_normal(Vertex v) const
{
    Point     nn(0,0,0);
    Halfedge  h = halfedge(v);

    if (h.is_valid())
    {
        const Halfedge hend = h;
        const Point p0 = vpoint_[v];

        Point   n, p1, p2;
        Scalar  cosine, angle;

        do
        {
            if (!is_boundary(h))
            {
                p1 = vpoint_[to_vertex(h)];
                p1 -= p0;
                p1.normalize();

                p2 = vpoint_[from_vertex(prev_halfedge(h))];
                p2 -= p0;
                p2.normalize();

                cosine = dot(p1,p2) / sqrt(dot(p1,p1)*dot(p2,p2));
                if      (cosine < -1.0) cosine = -1.0;
                else if (cosine >  1.0) cosine =  1.0;
                angle = acos(cosine);

                n   = cross(p1,p2).normalize();
                n  *= angle;
                nn += n;
            }

            h  = cw_rotated_halfedge(h);
        }
        while (h != hend);

        nn.normalize();
    }

    return nn;
}


//-----------------------------------------------------------------------------


Scalar
Surface_mesh::
edge_length(Edge e) const
{
    return norm(vpoint_[vertex(e,0)] - vpoint_[vertex(e,1)]);
}


//-----------------------------------------------------------------------------


void
Surface_mesh::
split(Face f, Vertex v)
{
    /*
     Split an arbitrary face into triangles by connecting each vertex of fh to vh.
     - fh will remain valid (it will become one of the triangles)
     - the halfedge handles of the new triangles will point to the old halfeges
     */

    Halfedge hend = halfedge(f);
    Halfedge h    = next_halfedge(hend);

    Halfedge hold = new_edge(to_vertex(hend), v);

    set_next_halfedge(hend, hold);
    set_face(hold, f);

    hold = opposite_halfedge(hold);

    while (h != hend)
    {
        Halfedge hnext = next_halfedge(h);

        Face fnew = new_face();
        set_halfedge(fnew, h);

        Halfedge hnew = new_edge(to_vertex(h), v);

        set_next_halfedge(hnew, hold);
        set_next_halfedge(hold, h);
        set_next_halfedge(h,    hnew);

        set_face(hnew, fnew);
        set_face(hold, fnew);
        set_face(h,    fnew);

        hold = opposite_halfedge(hnew);

        h = hnext;
    }

    set_next_halfedge(hold, hend);
    set_next_halfedge(next_halfedge(hend), hold);

    set_face(hold, f);

    set_halfedge(v, hold);
}


//-----------------------------------------------------------------------------


void
Surface_mesh::
split(Edge e, Vertex v)
{
    Halfedge h0 = halfedge(e, 0);
    Halfedge o0 = halfedge(e, 1);

    Vertex   v2 = to_vertex(o0);

    Halfedge e1 = new_edge(v, v2);
    Halfedge t1 = opposite_halfedge(e1);

    Face     f0 = face(h0);
    Face     f3 = face(o0);

    set_halfedge(v, h0);
    set_vertex(o0, v);

    if (!is_boundary(h0))
    {
        Halfedge h1 = next_halfedge(h0);
        Halfedge h2 = next_halfedge(h1);

        Vertex   v1 = to_vertex(h1);

        Halfedge e0 = new_edge(v, v1);
        Halfedge t0 = opposite_halfedge(e0);

        Face f1 = new_face();
        set_halfedge(f0, h0);
        set_halfedge(f1, h2);

        set_face(h1, f0);
        set_face(t0, f0);
        set_face(h0, f0);

        set_face(h2, f1);
        set_face(t1, f1);
        set_face(e0, f1);

        set_next_halfedge(h0, h1);
        set_next_halfedge(h1, t0);
        set_next_halfedge(t0, h0);

        set_next_halfedge(e0, h2);
        set_next_halfedge(h2, t1);
        set_next_halfedge(t1, e0);
    }
    else
    {
        set_next_halfedge(prev_halfedge(h0), t1);
        set_next_halfedge(t1, h0);
        // halfedge handle of _vh already is h0
    }


    if (!is_boundary(o0))
    {
        Halfedge o1 = next_halfedge(o0);
        Halfedge o2 = next_halfedge(o1);

        Vertex v3 = to_vertex(o1);

        Halfedge e2 = new_edge(v, v3);
        Halfedge t2 = opposite_halfedge(e2);

        Face f2 = new_face();
        set_halfedge(f2, o1);
        set_halfedge(f3, o0);

        set_face(o1, f2);
        set_face(t2, f2);
        set_face(e1, f2);

        set_face(o2, f3);
        set_face(o0, f3);
        set_face(e2, f3);

        set_next_halfedge(e1, o1);
        set_next_halfedge(o1, t2);
        set_next_halfedge(t2, e1);

        set_next_halfedge(o0, e2);
        set_next_halfedge(e2, o2);
        set_next_halfedge(o2, o0);
    }
    else
    {
        set_next_halfedge(e1, next_halfedge(o0));
        set_next_halfedge(o0, e1);
        set_halfedge(v, e1);
    }

    if (halfedge(v2) == h0)
        set_halfedge(v2, t1);
}


//-----------------------------------------------------------------------------


void
Surface_mesh::
insert_vertex(Edge e, Vertex v)
{
    // before:
    //
    // v0      h0       v2
    //  o--------------->o
    //   <---------------
    //         o0
    //
    // after:
    //
    // v0  h0   v   h1   v2
    //  o------>o------->o
    //   <------ <-------
    //     o0       o1

    Halfedge h0 = halfedge(e, 0);
    Halfedge h2 = next_halfedge(h0);
    Halfedge o0 = opposite_halfedge(h0);
    Halfedge o2 = prev_halfedge(o0);
    Vertex   v0 = to_vertex(o0);
    Vertex   v2 = to_vertex(h0);
    Face     fh = face(h0);
    Face     fo = face(o0);

    Halfedge h1 = new_edge(v, v2);
    Halfedge o1 = opposite_halfedge(h1);

    // adjust halfedge connectivity
    set_next_halfedge(h1, h2);
    set_next_halfedge(h0, h1);
    set_vertex(h0, v);
    set_vertex(h1, v2);
    set_face(h1, fh);

    set_next_halfedge(o1, o0);
    set_next_halfedge(o2, o1);
    set_vertex(o1, v);
    set_face(o1, fo);

    // adjust vertex connectivity
    set_halfedge(v2, o1);
    adjust_outgoing_halfedge(v2);
    set_halfedge(v, h1);
    adjust_outgoing_halfedge(v);

    // adjust face connectivity
    if (fh.is_valid()) set_halfedge(fh, h0);
    if (fo.is_valid()) set_halfedge(fo, o1);
}


//-----------------------------------------------------------------------------


void
Surface_mesh::
insert_edge(Halfedge h0, Halfedge h1)
{
    assert(face(h0) == face(h1));
    assert(face(h0).is_valid());

    Vertex   v0 = to_vertex(h0);
    Vertex   v1 = to_vertex(h1);

    Halfedge h2 = next_halfedge(h0);
    Halfedge h3 = next_halfedge(h1);

    Halfedge h4 = new_edge(v0, v1);
    Halfedge h5 = opposite_halfedge(h4);

    Face     f0 = face(h0);
    Face     f1 = new_face();

    set_halfedge(f0, h0);
    set_halfedge(f1, h1);

    set_next_halfedge(h0, h4);
    set_next_halfedge(h4, h3);
    set_face(h4, f0);

    set_next_halfedge(h1, h5);
    set_next_halfedge(h5, h2);
    Halfedge h = h2;
    do
    {
        set_face(h, f1);
        h = next_halfedge(h);
    }
    while (h != h2);
}


//-----------------------------------------------------------------------------


bool
Surface_mesh::
is_flip_ok(Edge e) const
{
    // boundary edges cannot be flipped
    if (is_boundary(e)) return false;

    // check if the flipped edge is already present in the mesh

    Halfedge h0 = halfedge(e, 0);
    Halfedge h1 = halfedge(e, 1);

    Vertex v0 = to_vertex(next_halfedge(h0));
    Vertex v1 = to_vertex(next_halfedge(h1));

    if (v0 == v1)   // this is generally a bad sign !!!
        return false;

    if (find_halfedge(v0, v1).is_valid())
        return false;

    return true;
}


//-----------------------------------------------------------------------------


void
Surface_mesh::
flip(Edge e)
{
    // CAUTION : Flipping a halfedge may result in
    // a non-manifold mesh, hence check for yourself
    // whether this operation is allowed or not!

    //let's make it sure it is actually checked
    assert(is_flip_ok(e));

    Halfedge a0 = halfedge(e, 0);
    Halfedge b0 = halfedge(e, 1);

    Halfedge a1 = next_halfedge(a0);
    Halfedge a2 = next_halfedge(a1);

    Halfedge b1 = next_halfedge(b0);
    Halfedge b2 = next_halfedge(b1);

    Vertex   va0 = to_vertex(a0);
    Vertex   va1 = to_vertex(a1);

    Vertex   vb0 = to_vertex(b0);
    Vertex   vb1 = to_vertex(b1);

    Face     fa  = face(a0);
    Face     fb  = face(b0);

    set_vertex(a0, va1);
    set_vertex(b0, vb1);

    set_next_halfedge(a0, a2);
    set_next_halfedge(a2, b1);
    set_next_halfedge(b1, a0);

    set_next_halfedge(b0, b2);
    set_next_halfedge(b2, a1);
    set_next_halfedge(a1, b0);

    set_face(a1, fb);
    set_face(b1, fa);

    set_halfedge(fa, a0);
    set_halfedge(fb, b0);

    if (halfedge(va0) == b0)
        set_halfedge(va0, a1);
    if (halfedge(vb0) == a0)
        set_halfedge(vb0, b1);
}


//-----------------------------------------------------------------------------


bool
Surface_mesh::
is_collapse_ok(Halfedge v0v1)
{
    Halfedge  v1v0(opposite_halfedge(v0v1));
    Vertex    v0(to_vertex(v1v0));
    Vertex    v1(to_vertex(v0v1));
    Vertex    vv, vl, vr;
    Halfedge  h1, h2;


    // the edges v1-vl and vl-v0 must not be both boundary edges
    if (!is_boundary(v0v1))
    {
        vl = to_vertex(next_halfedge(v0v1));
        h1 = next_halfedge(v0v1);
        h2 = next_halfedge(h1);
        if (is_boundary(opposite_halfedge(h1)) && is_boundary(opposite_halfedge(h2)))
            return false;
    }


    // the edges v0-vr and vr-v1 must not be both boundary edges
    if (!is_boundary(v1v0))
    {
        vr = to_vertex(next_halfedge(v1v0));
        h1 = next_halfedge(v1v0);
        h2 = next_halfedge(h1);
        if (is_boundary(opposite_halfedge(h1)) && is_boundary(opposite_halfedge(h2)))
            return false;
    }


    // if vl and vr are equal or both invalid -> fail
    if (vl == vr) return false;


    // edge between two boundary vertices should be a boundary edge
    if ( is_boundary(v0) && is_boundary(v1) &&
        !is_boundary(v0v1) && !is_boundary(v1v0))
        return false;


    // test intersection of the one-rings of v0 and v1
    Vertex_around_vertex_circulator vv_it, vv_end;
    vv_it = vv_end = vertices(v0);
    do
    {
        vv = *vv_it;
        if (vv != v1 && vv != vl && vv != vr)
            if (find_halfedge(vv, v1).is_valid())
                return false;
    }
    while (++vv_it != vv_end);


    // passed all tests
    return true;
}


//-----------------------------------------------------------------------------


void
Surface_mesh::
collapse(Halfedge h)
{
    Halfedge h0 = h;
    Halfedge h1 = prev_halfedge(h0);
    Halfedge o0 = opposite_halfedge(h0);
    Halfedge o1 = next_halfedge(o0);

    // remove edge
    remove_edge(h0);

    // remove loops
    if (next_halfedge(next_halfedge(h1)) == h1)
        remove_loop(h1);
    if (next_halfedge(next_halfedge(o1)) == o1)
        remove_loop(o1);
}


//-----------------------------------------------------------------------------


void
Surface_mesh::
remove_edge(Halfedge h)
{
    Halfedge  hn = next_halfedge(h);
    Halfedge  hp = prev_halfedge(h);

    Halfedge  o  = opposite_halfedge(h);
    Halfedge  on = next_halfedge(o);
    Halfedge  op = prev_halfedge(o);

    Face      fh = face(h);
    Face      fo = face(o);

    Vertex    vh = to_vertex(h);
    Vertex    vo = to_vertex(o);



    // halfedge -> vertex
    Halfedge_around_vertex_circulator vh_it, vh_end;
    vh_it = vh_end = halfedges(vo);
    do
    {
        set_vertex(opposite_halfedge(*vh_it), vh);
    }
    while (++vh_it != vh_end);


    // halfedge -> halfedge
    set_next_halfedge(hp, hn);
    set_next_halfedge(op, on);


    // face -> halfedge
    if (fh.is_valid())  set_halfedge(fh, hn);
    if (fo.is_valid())  set_halfedge(fo, on);


    // vertex -> halfedge
    if (halfedge(vh) == o)  set_halfedge(vh, hn);
    adjust_outgoing_halfedge(vh);
    set_halfedge(vo, Halfedge());


    // delete stuff
    if (!vdeleted_) vdeleted_ = vertex_property<bool>("v:deleted", false);
    if (!edeleted_) edeleted_ = edge_property<bool>("e:deleted", false);
    vdeleted_[vo]      = true; ++deleted_vertices_;
    edeleted_[edge(h)] = true; ++deleted_edges_;
    garbage_ = true;
}


//-----------------------------------------------------------------------------


void
Surface_mesh::
remove_loop(Halfedge h)
{
    Halfedge  h0 = h;
    Halfedge  h1 = next_halfedge(h0);

    Halfedge  o0 = opposite_halfedge(h0);
    Halfedge  o1 = opposite_halfedge(h1);

    Vertex    v0 = to_vertex(h0);
    Vertex    v1 = to_vertex(h1);

    Face      fh = face(h0);
    Face      fo = face(o0);



    // is it a loop ?
    assert ((next_halfedge(h1) == h0) && (h1 != o0));


    // halfedge -> halfedge
    set_next_halfedge(h1, next_halfedge(o0));
    set_next_halfedge(prev_halfedge(o0), h1);


    // halfedge -> face
    set_face(h1, fo);


    // vertex -> halfedge
    set_halfedge(v0, h1);  adjust_outgoing_halfedge(v0);
    set_halfedge(v1, o1);  adjust_outgoing_halfedge(v1);


    // face -> halfedge
    if (fo.is_valid() && halfedge(fo) == o0)
        set_halfedge(fo, h1);


    // delete stuff
    if (!edeleted_) edeleted_ = edge_property<bool>("e:deleted", false);
    if (!fdeleted_) fdeleted_ = face_property<bool>("f:deleted", false);
    if (fh.is_valid()) { fdeleted_[fh] = true; ++deleted_faces_; }
    edeleted_[edge(h0)] = true; ++deleted_edges_;
    garbage_ = true;
}


//-----------------------------------------------------------------------------


void
Surface_mesh::
garbage_collection()
{
    int  i, i0, i1,
    nV(vertices_size()),
    nE(edges_size()),
    nH(halfedges_size()),
    nF(faces_size());

    Vertex    v;
    Halfedge  h;
    Face      f;


    // setup handle mapping
    Vertex_property<Vertex>      vmap = add_vertex_property<Vertex>("v:garbage-collection");
    Halfedge_property<Halfedge>  hmap = add_halfedge_property<Halfedge>("h:garbage-collection");
    Face_property<Face>          fmap = add_face_property<Face>("f:garbage-collection");
    for (i=0; i<nV; ++i)
        vmap[Vertex(i)] = Vertex(i);
    for (i=0; i<nH; ++i)
        hmap[Halfedge(i)] = Halfedge(i);
    for (i=0; i<nF; ++i)
        fmap[Face(i)] = Face(i);



    // remove deleted vertices
    if (nV > 0)
    {
        i0=0;  i1=nV-1;

        while (1)
        {
            // find first deleted and last un-deleted
            while (!vdeleted_[Vertex(i0)] && i0 < i1)  ++i0;
            while ( vdeleted_[Vertex(i1)] && i0 < i1)  --i1;
            if (i0 >= i1) break;

            // swap
            vprops_.swap(i0, i1);
        };

        // remember new size
        nV = vdeleted_[Vertex(i0)] ? i0 : i0+1;
    }


    // remove deleted edges
    if (nE > 0)
    {
        i0=0;  i1=nE-1;

        while (1)
        {
            // find first deleted and last un-deleted
            while (!edeleted_[Edge(i0)] && i0 < i1)  ++i0;
            while ( edeleted_[Edge(i1)] && i0 < i1)  --i1;
            if (i0 >= i1) break;

            // swap
            eprops_.swap(i0, i1);
            hprops_.swap(2*i0,   2*i1);
            hprops_.swap(2*i0+1, 2*i1+1);
        };

        // remember new size
        nE = edeleted_[Edge(i0)] ? i0 : i0+1;
        nH = 2*nE;
    }


    // remove deleted faces
    if (nF > 0)
    {
        i0=0;  i1=nF-1;

        while (1)
        {
            // find 1st deleted and last un-deleted
            while (!fdeleted_[Face(i0)] && i0 < i1)  ++i0;
            while ( fdeleted_[Face(i1)] && i0 < i1)  --i1;
            if (i0 >= i1) break;

            // swap
            fprops_.swap(i0, i1);
        };

        // remember new size
        nF = fdeleted_[Face(i0)] ? i0 : i0+1;
    }


    // update vertex connectivity
    for (i=0; i<nV; ++i)
    {
        v = Vertex(i);
        if (!is_isolated(v))
            set_halfedge(v, hmap[halfedge(v)]);
    }


    // update halfedge connectivity
    for (i=0; i<nH; ++i)
    {
        h = Halfedge(i);
        set_vertex(h, vmap[to_vertex(h)]);
        set_next_halfedge(h, hmap[next_halfedge(h)]);
        if (!is_boundary(h))
            set_face(h, fmap[face(h)]);
    }


    // update handles of faces
    for (i=0; i<nF; ++i)
    {
        f = Face(i);
        set_halfedge(f, hmap[halfedge(f)]);
    }


    // remove handle maps
    remove_vertex_property(vmap);
    remove_halfedge_property(hmap);
    remove_face_property(fmap);


    // finally resize arrays
    vprops_.resize(nV); vprops_.free_memory();
    hprops_.resize(nH); hprops_.free_memory();
    eprops_.resize(nE); eprops_.free_memory();
    fprops_.resize(nF); fprops_.free_memory();

    deleted_vertices_ = deleted_edges_ = deleted_faces_ = 0;
    garbage_ = false;
}

//=============================================================================
