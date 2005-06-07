// Copyright (c) 2005  INRIA (France).
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
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Laurent Saboret, Pierre Alliez


#ifndef CGAL_MESH_ADAPTOR_PATCH_VERTEX_H
#define CGAL_MESH_ADAPTOR_PATCH_VERTEX_H

#include <CGAL/parameterization_assertions.h>

#include <list>

CGAL_BEGIN_NAMESPACE


// Forward reference
template<class MeshAdaptor_3> class Mesh_adaptor_patch_vertex_const_handle;


// Class Mesh_adaptor_patch_vertex
// Represents a vertex of a Mesh_adaptor_patch_3<MeshAdaptor_3> mesh
//
// Implementation note:
// A Mesh_adaptor_patch_vertex object is basically a handle to a
// MeshAdaptor_3::Vertex + its position / seam.
// Mesh_adaptor_patch_vertex comparison methods basically compare
// the address of the MeshAdaptor_3::Vertex that it points to.
//
template<class MeshAdaptor_3>
class Mesh_adaptor_patch_vertex
{
// Public types
public:

    // Export template parameter type
    typedef MeshAdaptor_3                   Adaptor;

// Public operations
public:

    // Default constructor
    Mesh_adaptor_patch_vertex()
    {
        m_adaptor_vertex = NULL;
        m_prev_seam_vertex = NULL;
        m_next_seam_vertex = NULL;
    }

    // Constructor from an INNER adaptor vertex
    Mesh_adaptor_patch_vertex(typename Adaptor::Vertex_handle adaptor_vertex)
    {
        CGAL_parameterization_assertion(adaptor_vertex != NULL);

        m_adaptor_vertex   = adaptor_vertex;
        m_prev_seam_vertex = NULL;
        m_next_seam_vertex = NULL;
    }

    // Constructor from a SEAM/MAIN BORDER adaptor vertex
    // prev_seam_vertex/next_seam_vertex are the previous/next
    // vertices on the seam
    Mesh_adaptor_patch_vertex(typename Adaptor::Vertex_handle adaptor_vertex,
                              typename Adaptor::Vertex_handle prev_seam_vertex,
                              typename Adaptor::Vertex_handle next_seam_vertex)
    {
        CGAL_parameterization_assertion(adaptor_vertex != NULL);
        CGAL_parameterization_assertion(prev_seam_vertex != NULL);
        CGAL_parameterization_assertion(next_seam_vertex != NULL);

        m_adaptor_vertex   = adaptor_vertex;
        m_prev_seam_vertex = prev_seam_vertex;
        m_next_seam_vertex = next_seam_vertex;
    }

    // Default copy constructor and operator =() are fine

    // Comparison
    bool operator==(const Mesh_adaptor_patch_vertex& vertex) const {
        return m_adaptor_vertex   == vertex.m_adaptor_vertex
            && m_prev_seam_vertex == vertex.m_prev_seam_vertex
            && m_next_seam_vertex == vertex.m_next_seam_vertex;
    }
    bool operator!=(const Mesh_adaptor_patch_vertex& vertex) const {
        return ! (*this == vertex);
    }

    // Get content
    typename Adaptor::Vertex_handle get_adaptor_vertex() {
        return m_adaptor_vertex;
    }
    typename Adaptor::Vertex_const_handle get_adaptor_vertex() const {
        return m_adaptor_vertex;
    }
    typename Adaptor::Vertex_handle get_prev_seam_vertex() {
        return m_prev_seam_vertex;
    }
    typename Adaptor::Vertex_const_handle get_prev_seam_vertex() const {
        return m_prev_seam_vertex;
    }
    typename Adaptor::Vertex_handle get_next_seam_vertex() {
        return m_next_seam_vertex;
    }
    typename Adaptor::Vertex_const_handle get_next_seam_vertex() const {
        return m_next_seam_vertex;
    }

// Fields
private:
    // The decorated vertex
    typename Adaptor::Vertex_handle m_adaptor_vertex;

    // Previous and next vertices on the main boundary/seam (NULL for inner vertex)
    typename Adaptor::Vertex_handle m_prev_seam_vertex;
    typename Adaptor::Vertex_handle m_next_seam_vertex;

}; // Mesh_adaptor_patch_vertex


// Class Mesh_adaptor_patch_vertex_handle
//
// This class represents a handle to a Mesh_adaptor_patch_vertex object,
// thus has the same behavior as Mesh_adaptor_patch_vertex* pointer type.
//
// Design pattern:
// Mesh_adaptor_patch_vertex_handle is an Bridge (see [GOF95]).
//
// Implementation note:
// A Mesh_adaptor_patch_vertex_handle contains in fact a Mesh_adaptor_patch_vertex
// object, which is basically a handle to a MeshAdaptor_3::Vertex.
// Mesh_adaptor_patch_vertex_handle comparison methods basically compare
// the address of the MeshAdaptor_3::Vertex pointed by the handles.
//
template<class MeshAdaptor_3>
class Mesh_adaptor_patch_vertex_handle
{
// Private types
private:

    typedef Mesh_adaptor_patch_vertex_handle Self;
    typedef Mesh_adaptor_patch_vertex<MeshAdaptor_3>
                                            Vertex;

// Public types
public:

    // Export template parameter types
    typedef MeshAdaptor_3                   Adaptor;

    // Iterator types
    typedef Vertex                          value_type;
    typedef std::ptrdiff_t                  difference_type;
    typedef std::size_t                     size_type;
    typedef Vertex&                         reference;
    typedef Vertex*                         pointer;

// Public operations
public:

    // Constructor
    Mesh_adaptor_patch_vertex_handle(pointer ptr = NULL)
    {
        if (ptr == NULL)
        {
            m_vertex = Vertex();
            m_ptr = NULL;
        }
        else
        {
            m_vertex = *ptr;
            m_ptr = &m_vertex;
        }
    }

    // Copy constructor
    Mesh_adaptor_patch_vertex_handle(const Self& hdl)
    {
        if (hdl.m_ptr == NULL)
        {
            m_vertex = Vertex();
            m_ptr = NULL;
        }
        else
        {
            m_vertex = *hdl.m_ptr;
            m_ptr = &m_vertex;
        }
    }

    // operator =()
    Self& operator =(const Self& hdl)
    {
        if (hdl.m_ptr == NULL)
        {
            m_vertex = Vertex();
            m_ptr = NULL;
        }
        else
        {
            m_vertex = *hdl.m_ptr;
            m_ptr = &m_vertex;
        }

        return *this;
    }

    // Compare patch vertices instead of patch vertex pointers (2 patch
    // vertex handles are equal iff they point to the same adaptor vertex)
    bool operator==(const Self& hdl) const
    {
        if (m_ptr == NULL || hdl.m_ptr == NULL)
            return m_ptr == hdl.m_ptr;
        else
            return *m_ptr == *(hdl.m_ptr);
    }
    bool operator!=(const Self& hdl) const { return ! (*this == hdl); }

    // Comparison to NULL pointer
    bool operator==(CGAL_NULL_TYPE ptr) const {
        CGAL_parameterization_assertion(ptr == NULL);
        return m_ptr == NULL;
    }
    bool operator!=(CGAL_NULL_TYPE ptr) const { return ! (*this == ptr); }

    pointer operator->()  const { return  m_ptr; }
    reference operator*() const { return *m_ptr; }
    operator Mesh_adaptor_patch_vertex_const_handle<Adaptor>() const {
        return m_ptr;
    }

// Fields
private:
    // The actual pointer
    pointer m_ptr;

    // Internal MeshAdaptor_3 vertex pointed to by m_ptr (except if NULL)
    Vertex m_vertex;

}; // Mesh_adaptor_patch_vertex_handle


// Class Mesh_adaptor_patch_vertex_const_handle
//
// This class represents a handle to a Mesh_adaptor_patch_vertex object,
// thus has the same behavior as const Mesh_adaptor_patch_vertex* pointer type.
//
// Design pattern:
// Mesh_adaptor_patch_vertex_const_handle is an Bridge (see [GOF95]).
//
// Implementation note:
// A Mesh_adaptor_patch_vertex_const_handle contains in fact a Mesh_adaptor_patch_vertex
// object, which is basically a handle to a MeshAdaptor_3::Vertex.
// Mesh_adaptor_patch_vertex_const_handle comparison methods basically compare
// the address of the MeshAdaptor_3::Vertex pointed by the handles.
//
template<class MeshAdaptor_3>
class Mesh_adaptor_patch_vertex_const_handle
{
// Private types
private:

    typedef Mesh_adaptor_patch_vertex_const_handle
                                            Self;
    typedef Mesh_adaptor_patch_vertex<MeshAdaptor_3>
                                            Vertex;

// Public types
public:

    // Export template parameter types
    typedef MeshAdaptor_3                   Adaptor;

    // Iterator types
    typedef Vertex                          value_type;
    typedef std::ptrdiff_t                  difference_type;
    typedef std::size_t                     size_type;
    typedef const Vertex&                   reference;
    typedef const Vertex*                   pointer;

// Public operations
public:

    // Constructor
    Mesh_adaptor_patch_vertex_const_handle(pointer ptr = NULL)
    {
        if (ptr == NULL)
        {
            m_vertex = Vertex();
            m_ptr = NULL;
        }
        else
        {
            m_vertex = *ptr;
            m_ptr = &m_vertex;
        }
    }

    // Copy constructor
    Mesh_adaptor_patch_vertex_const_handle(const Self& hdl)
    {
        if (hdl.m_ptr == NULL)
        {
            m_vertex = Vertex();
            m_ptr = NULL;
        }
        else
        {
            m_vertex = *hdl.m_ptr;
            m_ptr = &m_vertex;
        }
    }

    // operator =()
    Self& operator =(const Self& hdl)
    {
        if (hdl.m_ptr == NULL)
        {
            m_vertex = Vertex();
            m_ptr = NULL;
        }
        else
        {
            m_vertex = *hdl.m_ptr;
            m_ptr = &m_vertex;
        }

        return *this;
    }

    // Compare patch vertices instead of patch vertex pointers (2 patch
    // vertex handles are equal iff they point to the same adaptor vertex)
    bool operator==(const Self& hdl) const
    {
        if (m_ptr == NULL || hdl.m_ptr == NULL)
            return m_ptr == hdl.m_ptr;
        else
            return *m_ptr == *(hdl.m_ptr);
    }
    bool operator!=(const Self& hdl) const { return ! (*this == hdl); }

    // Comparison to NULL pointer
    bool operator==(CGAL_NULL_TYPE ptr) const {
        CGAL_parameterization_assertion(ptr == NULL);
        return m_ptr == NULL;
    }
    bool operator!=(CGAL_NULL_TYPE ptr) const { return ! (*this == ptr); }

    pointer operator->()  const { return  m_ptr; }
    reference operator*() const { return *m_ptr; }

// Fields
private:
    // The actual pointer
    pointer m_ptr;

    // Internal MeshAdaptor_3 vertex pointed to by m_ptr (except if NULL)
    Vertex m_vertex;

}; // Mesh_adaptor_patch_vertex_const_handle


CGAL_END_NAMESPACE

#endif //CGAL_MESH_ADAPTOR_PATCH_VERTEX_H

