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
template<class PatchableMeshAdaptor_3> class Mesh_adaptor_patch_vertex_const_handle;


/// Class Mesh_adaptor_patch_vertex
/// Represents a vertex of a Mesh_adaptor_patch_3<PatchableMeshAdaptor_3> mesh
///
/// Implementation note:
/// A Mesh_adaptor_patch_vertex object is basically a handle to a
/// PatchableMeshAdaptor_3::Vertex + its position / seam.
/// Mesh_adaptor_patch_vertex comparison methods compare the pointers.
///
template<class PatchableMeshAdaptor_3>
class Mesh_adaptor_patch_vertex
{
// Private types
private:

    typedef Mesh_adaptor_patch_vertex       Self;

// Public types
public:

    /// Export template parameter type
    typedef PatchableMeshAdaptor_3           Adaptor;

// Public operations
public:

    /// Default constructor
    Mesh_adaptor_patch_vertex()
    {
        m_vertex            = NULL;
        m_last_cw_neighbor  = NULL;
        m_first_cw_neighbor = NULL;
    }

    /// Constructor:
    /// - for an INNER adaptor vertex, last_cw_neighbor and first_cw_neighbor
    ///   must be NULL
    /// - for a SEAM/MAIN BORDER vertex, [first_cw_neighbor, last_cw_neighbor]
    /// defines the range of the valid neighbors of adaptor_vertex (included).
    explicit Mesh_adaptor_patch_vertex(
        typename Adaptor::Vertex_handle adaptor_vertex,
        typename Adaptor::Vertex_handle last_cw_neighbor  = NULL,
        typename Adaptor::Vertex_handle first_cw_neighbor = NULL)
    {
        CGAL_parameterization_assertion(adaptor_vertex != NULL);
        CGAL_parameterization_assertion( (last_cw_neighbor == NULL) ==
                                         (first_cw_neighbor == NULL) );

        m_vertex            = adaptor_vertex;
        m_last_cw_neighbor  = last_cw_neighbor;
        m_first_cw_neighbor = first_cw_neighbor;
    }

    /// Copy constructor
    Mesh_adaptor_patch_vertex(const Self& hdl)
    {
        m_vertex            = hdl.m_vertex;
        m_last_cw_neighbor  = hdl.m_last_cw_neighbor;
        m_first_cw_neighbor = hdl.m_first_cw_neighbor;
    }

    /// operator =()
    Self& operator =(const Self& hdl)
    {
        m_vertex            = hdl.m_vertex;
        m_last_cw_neighbor  = hdl.m_last_cw_neighbor;
        m_first_cw_neighbor = hdl.m_first_cw_neighbor;

        return *this;
    }

    /// Comparison
    bool operator==(const Mesh_adaptor_patch_vertex& vertex) const {
        return m_vertex            == vertex.m_vertex
            && m_last_cw_neighbor  == vertex.m_last_cw_neighbor
            && m_first_cw_neighbor == vertex.m_first_cw_neighbor;
    }
    bool operator!=(const Mesh_adaptor_patch_vertex& vertex) const {
        return ! (*this == vertex);
    }

    /// Get content
    typename Adaptor::Vertex_handle vertex() {
        return m_vertex;
    }
    typename Adaptor::Vertex_const_handle vertex() const {
        return m_vertex;
    }
    typename Adaptor::Vertex_handle last_cw_neighbor() {
        return m_last_cw_neighbor;
    }
    typename Adaptor::Vertex_const_handle last_cw_neighbor() const {
        return m_last_cw_neighbor;
    }
    typename Adaptor::Vertex_handle first_cw_neighbor() {
        return m_first_cw_neighbor;
    }
    typename Adaptor::Vertex_const_handle first_cw_neighbor() const {
        return m_first_cw_neighbor;
    }

// Fields
private:
    /// The decorated vertex
    typename Adaptor::Vertex_handle m_vertex;

    /// [m_first_cw_neighbor, m_last_cw_neighbor] defines the range of the valid 
    /// neighbors of m_vertex (included) if m_vertex is on the main boundary/seam 
    /// (NULL if inner vertex)
    typename Adaptor::Vertex_handle m_last_cw_neighbor;
    typename Adaptor::Vertex_handle m_first_cw_neighbor;

}; // Mesh_adaptor_patch_vertex


/// Class Mesh_adaptor_patch_vertex_handle
///
/// This class represents a handle to a Mesh_adaptor_patch_vertex object,
/// thus has the same behavior as Mesh_adaptor_patch_vertex* pointer type.
///
/// Design pattern:
/// Mesh_adaptor_patch_vertex_handle is an Bridge (see [GOF95]).
///
/// Implementation note:
/// A Mesh_adaptor_patch_vertex_handle contains in fact a Mesh_adaptor_patch_vertex
/// object, which is basically a handle to a PatchableMeshAdaptor_3::Vertex.
/// Mesh_adaptor_patch_vertex_handle comparison methods basically compare
/// the address of the PatchableMeshAdaptor_3::Vertex pointed by the handles.
///
template<class PatchableMeshAdaptor_3>
class Mesh_adaptor_patch_vertex_handle
{
// Private types
private:

    typedef Mesh_adaptor_patch_vertex_handle Self;

// Public types
public:

    /// Export template parameter type
    typedef PatchableMeshAdaptor_3           Adaptor;

    // Iterator types
    typedef Mesh_adaptor_patch_vertex<PatchableMeshAdaptor_3>
                                            Vertex;
    typedef Vertex                          value_type;
    typedef std::ptrdiff_t                  difference_type;
    typedef std::size_t                     size_type;
    typedef Vertex&                         reference;
    typedef Vertex*                         pointer;

// Public operations
public:

    /// Constructor from PatchableMeshAdaptor_3::Vertex pointer
    Mesh_adaptor_patch_vertex_handle(Vertex* ptr = NULL)
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

        assert(m_ptr == NULL || m_ptr == &m_vertex);
    }

    /// Extra constructor that will create the Mesh_adaptor_patch_3<PatchableMeshAdaptor_3>::Vertex on the fly
    /// - for an INNER adaptor vertex, last_cw_neighbor and first_cw_neighbor
    ///   must be NULL
    /// - for a SEAM/MAIN BORDER vertex, [first_cw_neighbor, last_cw_neighbor]
    /// defines the range of the valid neighbors of adaptor_vertex (included).
    explicit Mesh_adaptor_patch_vertex_handle(
        typename Adaptor::Vertex_handle adaptor_vertex,
        typename Adaptor::Vertex_handle last_cw_neighbor  = NULL,
        typename Adaptor::Vertex_handle first_cw_neighbor = NULL)
    {
        CGAL_parameterization_assertion(adaptor_vertex != NULL);
        m_vertex = Vertex(adaptor_vertex, last_cw_neighbor, first_cw_neighbor);
        m_ptr = &m_vertex;
    }

    /// Copy constructor
    Mesh_adaptor_patch_vertex_handle(const Self& hdl)
    {
        if (hdl.m_ptr == NULL)
        {
            m_vertex = Vertex();
            m_ptr = NULL;
        }
        else
        {
            m_vertex = hdl.m_vertex;
            m_ptr = &m_vertex;
        }

        assert(m_ptr == NULL || m_ptr == &m_vertex);
    }

    /// operator =()
    Self& operator =(const Self& hdl)
    {
        if (hdl.m_ptr == NULL)
        {
            m_vertex = Vertex();
            m_ptr = NULL;
        }
        else
        {
            m_vertex = hdl.m_vertex;
            m_ptr = &m_vertex;
        }

        assert(m_ptr == NULL || m_ptr == &m_vertex);

        return *this;
    }

    /// Compare patch vertices instead of patch vertex pointers (2 patch
    /// vertex handles are equal iff they point to the same adaptor vertex)
    bool operator==(const Self& hdl) const
    {
        if (m_ptr == NULL || hdl.m_ptr == NULL)
            return m_ptr == hdl.m_ptr;
        else
            return *m_ptr == *(hdl.m_ptr);
    }
    bool operator!=(const Self& hdl) const { return ! (*this == hdl); }

    /// Comparison to NULL pointer
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
    /// The actual pointer
    pointer m_ptr;

    /// Internal Mesh_adaptor_patch_3<PatchableMeshAdaptor_3> vertex 
    /// pointed to by m_ptr (except if NULL)
    Vertex m_vertex;

}; // Mesh_adaptor_patch_vertex_handle


/// Class Mesh_adaptor_patch_vertex_const_handle
///
/// This class represents a handle to a Mesh_adaptor_patch_vertex object,
/// thus has the same behavior as const Mesh_adaptor_patch_vertex* pointer type.
///
/// Design pattern:
/// Mesh_adaptor_patch_vertex_const_handle is an Bridge (see [GOF95]).
///
/// Implementation note:
/// A Mesh_adaptor_patch_vertex_const_handle contains in fact a Mesh_adaptor_patch_vertex
/// object, which is basically a handle to a PatchableMeshAdaptor_3::Vertex.
/// Mesh_adaptor_patch_vertex_const_handle comparison methods basically compare
/// the address of the PatchableMeshAdaptor_3::Vertex pointed by the handles.
///
template<class PatchableMeshAdaptor_3>
class Mesh_adaptor_patch_vertex_const_handle
{
// Private types
private:

    typedef Mesh_adaptor_patch_vertex_const_handle
                                            Self;

// Public types
public:

    /// Export template parameter type
    typedef PatchableMeshAdaptor_3           Adaptor;

    // Iterator types
    typedef Mesh_adaptor_patch_vertex<PatchableMeshAdaptor_3>
                                            Vertex;
    typedef Vertex                          value_type;
    typedef std::ptrdiff_t                  difference_type;
    typedef std::size_t                     size_type;
    typedef const Vertex&                   reference;
    typedef const Vertex*                   pointer;

// Public operations
public:

    /// Constructor from PatchableMeshAdaptor_3::Vertex pointer
    Mesh_adaptor_patch_vertex_const_handle(const Vertex* ptr = NULL)
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

        assert(m_ptr == NULL || m_ptr == &m_vertex);
    }

    /// Extra constructor that will create the Mesh_adaptor_patch_3<PatchableMeshAdaptor_3>::Vertex on the fly
    /// - for an INNER adaptor vertex, last_cw_neighbor and first_cw_neighbor
    ///   must be NULL
    /// - for a SEAM/MAIN BORDER vertex, [first_cw_neighbor, last_cw_neighbor]
    /// defines the range of the valid neighbors of adaptor_vertex (included).
    explicit Mesh_adaptor_patch_vertex_const_handle(
        typename Adaptor::Vertex_const_handle adaptor_vertex,
        typename Adaptor::Vertex_const_handle last_cw_neighbor  = NULL,
        typename Adaptor::Vertex_const_handle first_cw_neighbor = NULL)
    {
        CGAL_parameterization_assertion(adaptor_vertex != NULL);
        m_vertex = Vertex((typename Adaptor::Vertex*)&*adaptor_vertex,
                          (typename Adaptor::Vertex*)&*last_cw_neighbor,
                          (typename Adaptor::Vertex*)&*first_cw_neighbor);
        m_ptr = &m_vertex;
    }

    /// Copy constructor
    Mesh_adaptor_patch_vertex_const_handle(const Self& hdl)
    {
        if (hdl.m_ptr == NULL)
        {
            m_vertex = Vertex();
            m_ptr = NULL;
        }
        else
        {
            m_vertex = hdl.m_vertex;
            m_ptr = &m_vertex;
        }

        assert(m_ptr == NULL || m_ptr == &m_vertex);
    }

    /// operator =()
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

        assert(m_ptr == NULL || m_ptr == &m_vertex);

        return *this;
    }

    /// Compare patch vertices instead of patch vertex pointers (2 patch
    /// vertex handles are equal iff they point to the same adaptor vertex)
    bool operator==(const Self& hdl) const
    {
        if (m_ptr == NULL || hdl.m_ptr == NULL)
            return m_ptr == hdl.m_ptr;
        else
            return *m_ptr == *(hdl.m_ptr);
    }
    bool operator!=(const Self& hdl) const { return ! (*this == hdl); }

    /// Comparison to NULL pointer
    bool operator==(CGAL_NULL_TYPE ptr) const {
        CGAL_parameterization_assertion(ptr == NULL);
        return m_ptr == NULL;
    }
    bool operator!=(CGAL_NULL_TYPE ptr) const { return ! (*this == ptr); }

    pointer operator->()  const { return  m_ptr; }
    reference operator*() const { return *m_ptr; }

// Fields
private:
    /// The actual pointer
    pointer m_ptr;

    /// Internal Mesh_adaptor_patch_3<PatchableMeshAdaptor_3> vertex 
    /// pointed to by m_ptr (except if NULL)
    Vertex m_vertex;

}; // Mesh_adaptor_patch_vertex_const_handle


CGAL_END_NAMESPACE

#endif //CGAL_MESH_ADAPTOR_PATCH_VERTEX_H

