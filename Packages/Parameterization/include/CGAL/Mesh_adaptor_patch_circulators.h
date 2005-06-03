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


#ifndef CGAL_MESH_ADAPTOR_PATCH_CIRCULATORS_H
#define CGAL_MESH_ADAPTOR_PATCH_CIRCULATORS_H

#include <CGAL/iterator.h>
#include <CGAL/circulator.h>
#include <CGAL/Convertible_iterator_project.h>
#include <CGAL/Convertible_circulator_project.h>
#include <CGAL/function_objects.h>
#include <CGAL/HalfedgeDS_iterator.h>
#include <CGAL/Mesh_adaptor_patch_vertex.h>

#include <CGAL/parameterization_assertions.h>

CGAL_BEGIN_NAMESPACE


// Class Mesh_patch_vertex_around_vertex_cir
// Represents a (clockwise) circulator around a vertex 
// of a Mesh_adaptor_patch_3<MeshAdaptor_3> mesh
template<class MeshPatchPtrType,    // = [const] Mesh_adaptor_patch_3*
         class VertexHandleType,    // = Mesh_adaptor_patch_3::Vertex_[const_]handle
         class AdaptorVertexAroundVertexCirculatorType>   
                                    // = Adaptor::Vertex_around_vertex_[const_]circulator
class Mesh_patch_vertex_around_vertex_cir 
    : public VertexHandleType
{
    typedef VertexHandleType                    Base;
    typedef Mesh_patch_vertex_around_vertex_cir Self;
    typedef typename std::iterator_traits<MeshPatchPtrType>::value_type 
                                                Mesh_patch;
    typedef typename Mesh_patch::Adaptor        Adaptor;

public:

// TYPES
// -----

    // Export template parameter types
    typedef MeshPatchPtrType                    Mesh_patch_ptr_t;
    typedef VertexHandleType                    Vertex_handle_t;
    typedef AdaptorVertexAroundVertexCirculatorType 
                                                Adaptor_vertex_around_vertex_cir_t;

    // Iterator types
    typedef typename Adaptor_vertex_around_vertex_cir_t::iterator_category 
                                                iterator_category;
    typedef typename Vertex_handle_t::value_type value_type;
    typedef std::ptrdiff_t                      difference_type;
    typedef std::size_t                         size_type;
    typedef typename Vertex_handle_t::reference reference;
    typedef typename Vertex_handle_t::pointer   pointer;

// CREATION
// --------

    // Circulator pointing to NULL
    Mesh_patch_vertex_around_vertex_cir() {}

    // Get circulator over the vertices incident to 'vertex'
    //
    // Implementation note:
    // If m_center_vertex is a seam/main boundary vertex, then [get_prev_seam_vertex(),
    // get_next_seam_vertex()] is the range of its inner neighbor vertices.
    Mesh_patch_vertex_around_vertex_cir(Mesh_patch_ptr_t mesh, 
                                        Vertex_handle_t vertex, 
                                        Adaptor_vertex_around_vertex_cir_t adaptor_circulator)
        : m_mesh_patch(mesh),
          m_center_vertex(vertex),
          m_adaptor_circulator(adaptor_circulator)
    {
        CGAL_parameterization_assertion(m_mesh_patch != NULL);
        CGAL_parameterization_assertion(m_mesh_patch->is_valid(vertex));
        CGAL_parameterization_assertion(adaptor_circulator != NULL);

        // Update the inherited vertex handle
        update_inherited_handle();
    }

    // Copy constructor
    Self(const Self& cir)
      : Base(cir),
        m_mesh_patch(cir.m_mesh_patch),
        m_center_vertex(cir.m_center_vertex),
        m_adaptor_circulator(cir.m_adaptor_circulator)
    {
        // Update the inherited vertex handle
        update_inherited_handle();
    }

    // operator =()
    Self& operator =(const Self& cir) 
    {
        Base::operator=()(cir);
        m_mesh_patch = cir.m_mesh_patch;
        m_center_vertex = cir.m_center_vertex;
        m_adaptor_circulator = cir.m_adaptor_circulator;

        // Update the inherited vertex handle
        update_inherited_handle();

        return *this;
    }

// OPERATIONS Forward Category
// ---------------------------

    bool operator==(const Self& cir)    const { return Base::operator==(cir); }
    bool operator!=(const Self& cir)    const { return !(*this == cir); }
    bool operator==(CGAL_NULL_TYPE ptr) const { return Base::operator==(ptr); }
    bool operator!=(CGAL_NULL_TYPE ptr) const { return !(*this == ptr); }

    //  operator*() and operator->() are inherited

    // Clockwise rotation
    Self& operator++()
    {
        // Increment m_adaptor_circulator
        //
        // If non seam vertex, m_adaptor_circulator simply
        // circulates over m_center_vertex's vertices
        if (m_center_vertex->get_next_seam_vertex() == NULL)
        {
            m_adaptor_circulator++;
        }
        else // if seam vertex, circulates only from get_next_seam_vertex()
        {    //                 to get_prev_seam_vertex() (included)

            // TODO: check circulator way
            assert(false);

            if (m_center_vertex->get_prev_seam_vertex() == m_adaptor_circulator)
            {
                do
                    m_adaptor_circulator++;
                while (m_center_vertex->get_next_seam_vertex() != m_adaptor_circulator);
            }
            else
            {
                m_adaptor_circulator++;
            }
        }

        // Update the inherited vertex handle
        update_inherited_handle();

        return *this;
    }
    //
    Self  operator++(int) {
        Self tmp = *this;
        ++*this;
        return tmp;
    }

// OPERATIONS Bidirectional Category
// ---------------------------------

    // Counterclockwise rotation
    Self& operator--()
    {
        // Decrement m_adaptor_circulator
        //
        // If non seam vertex, m_adaptor_circulator simply
        // circulates over m_center_vertex's vertices
        if (m_center_vertex->get_next_seam_vertex() == NULL)
        {
            m_adaptor_circulator--;
        }
        else // if seam vertex, circulates only from get_prev_seam_vertex()
        {    //                 to get_next_seam_vertex() (included)

            // TODO: check circulator way
            assert(false);

            if (m_center_vertex->get_next_seam_vertex() == m_adaptor_circulator)
            {
                do
                    m_adaptor_circulator--;
                while (m_center_vertex->get_prev_seam_vertex() != m_adaptor_circulator);
            }
            else
            {
                m_adaptor_circulator--;
            }
        }

        // Update the inherited vertex handle
        update_inherited_handle();

        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }

private:
    // Update the inherited vertex handle
    //
    // Precondition: m_adaptor_circulator is valid
    void update_inherited_handle()
    {
        // TODO: Support the case where the current vertex m_adaptor_circulator 
        //       and the center vertex m_center_vertex are on the same (non-oriented) seam.
        //       See Mesh_patch_vertex_around_facet_cir::update_inherited_handle() 
        //       which does it with 3 points. Use a rotation to get the 3rd point?
        assert(m_mesh_patch->m_mesh_adaptor->get_edge_seaming(
                m_adaptor_circulator, m_center_vertex->get_adaptor_vertex()) != Adaptor::BORDER 
            || m_mesh_patch->m_mesh_adaptor->get_edge_seaming(
                m_center_vertex->get_adaptor_vertex(), m_adaptor_circulator) != Adaptor::BORDER);

        // Get current decorated vertex
        Vertex_handle_t current_decorated_vertex 
            = m_mesh_patch->get_decorated_vertex(m_adaptor_circulator, 
                                                 m_center_vertex->get_adaptor_vertex());

        // Update the inherited vertex handle
        Base::operator=(current_decorated_vertex);
    }

private:
    // The mesh that we are circulating on:
    Mesh_patch_ptr_t m_mesh_patch;

    // Vertex center of the circulation (+ circulator range for a border vertex)
    Vertex_handle_t m_center_vertex;

    // Internal circulator
    Adaptor_vertex_around_vertex_cir_t m_adaptor_circulator;

}; // Mesh_patch_vertex_around_vertex_cir


// Class Mesh_patch_vertex_around_facet_cir
// Represents a (clockwise) circulator around a facet 
// of a Mesh_adaptor_patch_3<MeshAdaptor_3> mesh
template<class MeshPatchPtrType,    // = [const] Mesh_adaptor_patch_3*
         class VertexHandleType,    // = Mesh_adaptor_patch_3::Vertex_[const_]handle
         class AdaptorVertexAroundFacetCirculatorType>   
                                    // = Adaptor::Vertex_around_facet_[const_]circulator
class Mesh_patch_vertex_around_facet_cir 
    : public VertexHandleType
{
    typedef VertexHandleType                    Base;
    typedef Mesh_patch_vertex_around_facet_cir  Self;
    typedef typename std::iterator_traits<MeshPatchPtrType>::value_type 
                                                Mesh_patch;
    typedef typename Mesh_patch::Adaptor        Adaptor;
    typedef typename Mesh_patch::Vertex         Vertex;

public:

// TYPES
// -----

    // Export template parameter types
    typedef MeshPatchPtrType                    Mesh_patch_ptr_t;
    typedef VertexHandleType                    Vertex_handle_t;
    typedef AdaptorVertexAroundFacetCirculatorType 
                                                Adaptor_vertex_around_facet_cir_t;

    // Iterator types
    typedef typename Adaptor_vertex_around_facet_cir_t::iterator_category 
                                                iterator_category;
    typedef typename Vertex_handle_t::value_type value_type;
    typedef std::ptrdiff_t                      difference_type;
    typedef std::size_t                         size_type;
    typedef typename Vertex_handle_t::reference reference;
    typedef typename Vertex_handle_t::pointer   pointer;

// CREATION
// --------

    // Circulator pointing to NULL
    Mesh_patch_vertex_around_facet_cir() {}

    Mesh_patch_vertex_around_facet_cir(Mesh_patch_ptr_t mesh, 
                                       Adaptor_vertex_around_facet_cir_t adaptor_circulator)
        : m_mesh_patch(mesh),
          m_adaptor_circulator(adaptor_circulator)
    {
        CGAL_parameterization_assertion(m_mesh_patch != NULL);
        CGAL_parameterization_assertion(adaptor_circulator != NULL);

        // Update the inherited vertex handle
        update_inherited_handle();
    }

    // Copy constructor
    Self(const Self& cir)
      : Base(cir),
        m_mesh_patch(cir.m_mesh_patch),
        m_adaptor_circulator(cir.m_adaptor_circulator)
    {
        // Update the inherited vertex handle
        update_inherited_handle();
    }

    // operator =()
    Self& operator =(const Self& cir) 
    {
        Base::operator=(cir);
        m_mesh_patch = cir.m_mesh_patch;
        m_adaptor_circulator = cir.m_adaptor_circulator;

        // Update the inherited vertex handle
        update_inherited_handle();

        return *this;
    }

// OPERATIONS Forward Category
// ---------------------------

    bool operator==(CGAL_NULL_TYPE ptr) const { return Base::operator==(ptr); }
    bool operator!=(CGAL_NULL_TYPE ptr) const { return !(*this == ptr); }
    bool operator==(const Self& cir)    const { return Base::operator==(cir); }
    bool operator!=(const Self& cir)    const { return !(*this == cir); }

    //  operator*() and operator->() are inherited

    // Clockwise rotation
    Self& operator++()
    {
        // Increment m_adaptor_circulator
        m_adaptor_circulator++;

        // Update the inherited vertex handle
        update_inherited_handle();

        return *this;
    }
    //
    Self  operator++(int) {
        Self tmp = *this;
        ++*this;
        return tmp;
    }

// OPERATIONS Bidirectional Category
// ---------------------------------

    // Counterclockwise rotation
    Self& operator--()
    {
        // Decrement m_adaptor_circulator
        m_adaptor_circulator--;

        // Update the inherited vertex handle
        update_inherited_handle();

        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }

private:
    // Update the inherited vertex handle
    //
    // Precondition: m_adaptor_circulator is valid
    void update_inherited_handle()
    {
        Vertex_handle_t current_decorated_vertex = NULL;

        // Get neighbor adaptor vertices on facet
        Adaptor_vertex_around_facet_cir_t next_vertex = m_adaptor_circulator;
        next_vertex++;
        Adaptor_vertex_around_facet_cir_t prev_vertex = m_adaptor_circulator;
        prev_vertex--;

        // If (m_adaptor_circulator, next_vertex) isn't a seam (non-oriented) edge
        if (m_mesh_patch->m_mesh_adaptor->get_edge_seaming(
                        m_adaptor_circulator, next_vertex) != Adaptor::BORDER
         || m_mesh_patch->m_mesh_adaptor->get_edge_seaming(
                        next_vertex, m_adaptor_circulator) != Adaptor::BORDER)
        {
            current_decorated_vertex 
                    = m_mesh_patch->get_decorated_vertex(m_adaptor_circulator,
                                                         next_vertex);
        }
        // If (m_adaptor_circulator, prev_vertex) isn't a seam (non-oriented) edge
        else if (m_mesh_patch->m_mesh_adaptor->get_edge_seaming(
                        m_adaptor_circulator, prev_vertex) != Adaptor::BORDER
         || m_mesh_patch->m_mesh_adaptor->get_edge_seaming(
                        prev_vertex, m_adaptor_circulator) != Adaptor::BORDER)
        {
            current_decorated_vertex 
                    = m_mesh_patch->get_decorated_vertex(m_adaptor_circulator,
                                                         prev_vertex);
        }
        // Special case: Both edges belong to the seam
        else
        {
            Vertex vertex((typename Adaptor::Vertex*)&*m_adaptor_circulator, 
                          (typename Adaptor::Vertex*)&*next_vertex,     // order...
                          (typename Adaptor::Vertex*)&*prev_vertex);    // ...inverted!

            // Implementation note:
            // The next line seems to return a reference to a local Vertex variable. 
            // In fact, Vertex_[const_]handle constructor copies the Vertex object.
            // The purpose is to save the time of searching the Vertex in 
            // m_inner_and_border_vertices list.
            current_decorated_vertex = &vertex; 
        }

        // Update the inherited vertex handle
        assert(m_mesh_patch->is_valid(current_decorated_vertex));
        Base::operator=(current_decorated_vertex);
    }

private:
    // The mesh that we are circulating on:
    Mesh_patch_ptr_t m_mesh_patch;

    // Internal circulator
    Adaptor_vertex_around_facet_cir_t m_adaptor_circulator;

}; // Mesh_patch_vertex_around_facet_cir


CGAL_END_NAMESPACE

#endif //CGAL_MESH_ADAPTOR_PATCH_CIRCULATORS_H

