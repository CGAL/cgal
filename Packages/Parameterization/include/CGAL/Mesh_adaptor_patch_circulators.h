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
// represents a (clockwise) circulator around a vertex
// of a Mesh_adaptor_patch_3<MeshAdaptor_3> mesh
template<class MeshPatchPtrType,    // = [const] Mesh_adaptor_patch_3*
         class VertexHandleType,    // = Mesh_adaptor_patch_3::Vertex_[const_]handle
         class AdaptorVertexAroundVertexCirculatorType,
                                    // = Adaptor::Vertex_around_vertex_[const_]circulator
         class AdaptorVertexHandleType>
                                    // = Adaptor::Vertex_[const_]handle
class Mesh_patch_vertex_around_vertex_cir
    : public VertexHandleType
{
// PRIVATE TYPES
// -------------

    // Base and self classes
    typedef VertexHandleType                    Base;
    typedef Mesh_patch_vertex_around_vertex_cir Self;

    // Mesh_patch types
    typedef typename std::iterator_traits<MeshPatchPtrType>::value_type
                                                Mesh_patch;
    typedef typename Mesh_patch::Adaptor        Adaptor;

public:

// PUBLIC TYPES
// ------------

    // Export template parameter types
    typedef MeshPatchPtrType                    Mesh_patch_c_ptr;
    typedef VertexHandleType                    Vertex_c_handle;
    typedef AdaptorVertexAroundVertexCirculatorType
                                                Adaptor_vertex_around_vertex_c_cir;
    typedef AdaptorVertexHandleType             Adaptor_vertex_c_handle;

    // Iterator types
    typedef typename Adaptor_vertex_around_vertex_c_cir::iterator_category
                                                iterator_category;
    typedef typename Vertex_c_handle::value_type value_type;
    typedef std::ptrdiff_t                      difference_type;
    typedef std::size_t                         size_type;
    typedef typename Vertex_c_handle::reference reference;
    typedef typename Vertex_c_handle::pointer   pointer;

// CREATION
// --------

    // Circulator pointing to NULL
    Mesh_patch_vertex_around_vertex_cir() {}

    // Get circulator over the vertices incident to 'vertex'
    // 'start_position' defines the initial position of the circulator [required}
    Mesh_patch_vertex_around_vertex_cir(Mesh_patch_c_ptr mesh,
                                        Vertex_c_handle vertex,
                                        Vertex_c_handle start_position)
      : Base(start_position),
        m_mesh_patch(mesh),
        m_center(vertex)
    {
        CGAL_parameterization_assertion(m_mesh_patch != NULL);
        CGAL_parameterization_assertion(m_mesh_patch->m_mesh_adaptor != NULL);
        CGAL_parameterization_assertion(m_mesh_patch->is_valid(vertex));
        CGAL_parameterization_assertion(m_mesh_patch->is_valid(start_position));

//#ifdef DEBUG_TRACE
//        std::cerr << "    Mesh_patch_vertex_around_vertex_cir(";
//        std::cerr << "#" << m_mesh_patch->get_vertex_index(vertex) << ",";
//        std::cerr << "#" << m_mesh_patch->get_vertex_index(start_position) << ")\n";
//#endif

        // Construct an adaptor circulator over the vertices
        // incident to vertex->vertex()
        m_adaptor_circulator = m_mesh_patch->m_mesh_adaptor->vertices_around_vertex_begin(
            vertex->vertex(), start_position->vertex());
    }

    // Copy constructor
    Mesh_patch_vertex_around_vertex_cir(const Self& cir)
      : Base(cir),
        m_mesh_patch(cir.m_mesh_patch),
        m_center(cir.m_center),
        m_adaptor_circulator(cir.m_adaptor_circulator)
    {
    }

    // operator =()
    Self& operator =(const Self& cir)
    {
        Base::operator=(cir);
        m_mesh_patch = cir.m_mesh_patch;
        m_center = cir.m_center;
        m_adaptor_circulator = cir.m_adaptor_circulator;

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
//#ifdef DEBUG_TRACE
//        std::cerr << "    Mesh_patch_vertex_around_vertex_cir ++\n";
//#endif

        // Check that the inherited vertex handle is valid
        assert((*this)->vertex() == m_adaptor_circulator);

        // If m_center is an inner vertex, m_adaptor_circulator 
        // already circulates over m_center's vertices
        if (m_center->first_cw_neighbor() == NULL)
        {
            m_adaptor_circulator++;

            // Update the inherited vertex handle
            update_inherited_handle();

            return *this;
        }

        // If m_center is a border/seam vertex (except a seam extremity), 
        // *this must circulate only from
        // first_cw_neighbor() to last_cw_neighbor() (included)
        if (m_center->last_cw_neighbor() != m_center->first_cw_neighbor())
        {
            // if the previous position m_adaptor_circulator
            // is m_center's last clockwise neighbor
            if (m_center->last_cw_neighbor() == m_adaptor_circulator)
            {
                // Skip the outer part of the seam
                while (m_center->first_cw_neighbor() != m_adaptor_circulator)
                    m_adaptor_circulator++;
            }
            else // if the previous position is inner,
            {    // simply rotate the adaptor circulator
                m_adaptor_circulator++;
            }

            // Update the inherited vertex handle
            update_inherited_handle();

            return *this;
        }

        // Special case: m_center is the extremity of the seam.
        // Its seam adaptor neighbor vertex is "virtually" duplicated.
        assert(m_center->last_cw_neighbor() == m_center->first_cw_neighbor());

        // If the previous position is the last "virtual" clockwise neighbor, 
        // move to first "virtual" clockwise neighbor 
        // without rotating m_adaptor_circulator.
        if (m_center->last_cw_neighbor() == m_adaptor_circulator
         && (*this)->first_cw_neighbor() == m_center->vertex())
        {   
            // Update directly the inherited vertex handle
            // because this case is ambiguous for update_inherited_handle()
            Vertex_c_handle current_decorated_vertex((*this)->vertex(), 
                                                     (*this)->first_cw_neighbor(),  // order...
                                                     (*this)->last_cw_neighbor());  // ...inverted!
//#ifdef DEBUG_TRACE
//            std::cerr << "      Mesh_patch_vertex_around_vertex_cir = (";
//            std::cerr << "#" << m_mesh_patch->get_vertex_index(m_center) << ",";
//            std::cerr << "#" << m_mesh_patch->get_vertex_index(current_decorated_vertex);
//            std::cerr << ")\n";
//#endif
            assert(m_mesh_patch->is_valid(current_decorated_vertex));
            Base::operator=(current_decorated_vertex);

            return *this;
        }

        // If the previous position is NOT the last "virtual" clockwise neighbor, 
        // simply rotate the adaptor circulator
        m_adaptor_circulator++;

        // Update the inherited vertex handle
        //
        // If the new position is the last "virtual" clockwise neighbor, 
        // update directly the inherited vertex handle
        // because this case is ambiguous for update_inherited_handle()
        if (m_center->last_cw_neighbor() == m_adaptor_circulator)
        {   
            Vertex_c_handle current_decorated_vertex
                    = m_mesh_patch->get_decorated_border_vertex(m_adaptor_circulator,
                                                                NULL,
                                                                m_center->vertex());
//#ifdef DEBUG_TRACE
//            std::cerr << "      Mesh_patch_vertex_around_vertex_cir = (";
//            std::cerr << "#" << m_mesh_patch->get_vertex_index(m_center) << ",";
//            std::cerr << "#" << m_mesh_patch->get_vertex_index(current_decorated_vertex);
//            std::cerr << ")\n";
//#endif
            assert(m_mesh_patch->is_valid(current_decorated_vertex));
            Base::operator=(current_decorated_vertex);

            return *this;
        }
        else
        {
            update_inherited_handle();

            return *this;
        }

    }
    //
    Self  operator++(int) {
        Self tmp = *this;
        ++*this;
        return tmp;
    }

// OPERATIONS Bidirectional Category
// ---------------------------------

    // Counter-clockwise rotation
    Self& operator--()
    {
//#ifdef DEBUG_TRACE
//        std::cerr << "    Mesh_patch_vertex_around_vertex_cir --\n";
//#endif

        // Check that the inherited vertex handle is valid
        assert((*this)->vertex() == m_adaptor_circulator);

        // If m_center is an inner vertex, m_adaptor_circulator 
        // already circulates over m_center's vertices
        if (m_center->first_cw_neighbor() == NULL)
        {
            m_adaptor_circulator--;

            // Update the inherited vertex handle
            update_inherited_handle();

            return *this;
        }

        // If m_center is a border/seam vertex (except a seam extremity), 
        // *this must circulate only from
        // last_cw_neighbor() to first_cw_neighbor() (included)
        if (m_center->last_cw_neighbor() != m_center->first_cw_neighbor())
        {
            // if the previous position m_adaptor_circulator
            // is m_center's last counter-clockwise neighbor
            if (m_center->first_cw_neighbor() == m_adaptor_circulator)
            {
                // Skip the outer part of the seam
                while (m_center->last_cw_neighbor() != m_adaptor_circulator)
                    m_adaptor_circulator--;
            }
            else // if the previous position is inner,
            {    // simply rotate the adaptor circulator
                m_adaptor_circulator--;
            }

            // Update the inherited vertex handle
            update_inherited_handle();

            return *this;
        }

        // Special case: m_center is the extremity of the seam.
        // Its seam adaptor neighbor vertex is "virtually" duplicated.
        assert(m_center->last_cw_neighbor() == m_center->first_cw_neighbor());

        // If the previous position is on the last "virtual" counter-clockwise  
        // neighbor, move to first "virtual" counter-clockwise neighbor 
        // without rotating m_adaptor_circulator.
        if (m_center->last_cw_neighbor() == m_adaptor_circulator
         && (*this)->last_cw_neighbor() == m_center->vertex())
        {   
            // Update directly the inherited vertex handle
            // because this case is ambiguous for update_inherited_handle()
            Vertex_c_handle current_decorated_vertex((*this)->vertex(), 
                                                     (*this)->first_cw_neighbor(),  // order...
                                                     (*this)->last_cw_neighbor());  // ...inverted!
//#ifdef DEBUG_TRACE
//            std::cerr << "      Mesh_patch_vertex_around_vertex_cir = (";
//            std::cerr << "#" << m_mesh_patch->get_vertex_index(m_center) << ",";
//            std::cerr << "#" << m_mesh_patch->get_vertex_index(current_decorated_vertex);
//            std::cerr << ")\n";
//#endif
            assert(m_mesh_patch->is_valid(current_decorated_vertex));
            Base::operator=(current_decorated_vertex);

            return *this;
        }

        // If the previous position is NOT the last "virtual" counter-clockwise 
        //  neighbor, simply rotate the adaptor circulator
        m_adaptor_circulator--;

        // Update the inherited vertex handle
        //
        // If the new position is the last "virtual" counter-clockwise neighbor, 
        // update directly the inherited vertex handle
        // because this case is ambiguous for update_inherited_handle()
        if (m_center->last_cw_neighbor() == m_adaptor_circulator)
        {   
            Vertex_c_handle current_decorated_vertex
                    = m_mesh_patch->get_decorated_border_vertex(m_adaptor_circulator,
                                                                m_center->vertex(),
                                                                NULL);
//#ifdef DEBUG_TRACE
//            std::cerr << "      Mesh_patch_vertex_around_vertex_cir = (";
//            std::cerr << "#" << m_mesh_patch->get_vertex_index(m_center) << ",";
//            std::cerr << "#" << m_mesh_patch->get_vertex_index(current_decorated_vertex);
//            std::cerr << ")\n";
//#endif
            assert(m_mesh_patch->is_valid(current_decorated_vertex));
            Base::operator=(current_decorated_vertex);

            return *this;
        }
        else
        {
            update_inherited_handle();

            return *this;
        }

    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }

private:
    // Update the inherited vertex handle
    //
    // Precondition: m_adaptor_circulator and m_center are valid
    void update_inherited_handle()
    {
        Vertex_c_handle current_decorated_vertex = NULL;

        // Ambiguous case: m_center is the extremity of the seam
        //                 and m_adaptor_circulator is on the seam
        assert(m_center->last_cw_neighbor() != m_adaptor_circulator
            || m_center->first_cw_neighbor() != m_adaptor_circulator);

        // Get next vertex on facet
        Adaptor_vertex_c_handle next_vertex = NULL;
        Adaptor_vertex_around_vertex_c_cir ccw_neighbor = m_adaptor_circulator;
        ccw_neighbor--;
        if (m_center->first_cw_neighbor() == NULL) // if inner vertex
        {
            next_vertex = ccw_neighbor;
        }
        else // if border/seam vertex, circulates only from 
        {    // first_cw_neighbor() to last_cw_neighbor() (included)
            if (m_center->first_cw_neighbor() == m_adaptor_circulator)
                next_vertex = m_center->vertex();
            else
                next_vertex = ccw_neighbor;
        }

        // If (m_adaptor_circulator, next_vertex) isn't a seam (non-oriented) edge
        if (m_mesh_patch->m_mesh_adaptor->get_halfedge_seaming(
                        m_adaptor_circulator, next_vertex) != Mesh_patch::BORDER
         || m_mesh_patch->m_mesh_adaptor->get_halfedge_seaming(
                        next_vertex, m_adaptor_circulator) != Mesh_patch::BORDER)
        {
            current_decorated_vertex
                    = m_mesh_patch->get_decorated_inner_vertex(m_adaptor_circulator,
                                                               next_vertex);
        }
        // Special case: both vertices belong to the seam
        else 
        {
            current_decorated_vertex
                    = m_mesh_patch->get_decorated_border_vertex(m_adaptor_circulator,
                                                                next_vertex,
                                                                NULL);
        }

//#ifdef DEBUG_TRACE
//        std::cerr << "      Mesh_patch_vertex_around_vertex_cir = (";
//        std::cerr << "#" << m_mesh_patch->get_vertex_index(m_center) << ",";
//        std::cerr << "#" << m_mesh_patch->get_vertex_index(current_decorated_vertex);
//        std::cerr << ")\n";
//#endif

        // Update the inherited vertex handle
        assert(m_mesh_patch->is_valid(current_decorated_vertex));
        Base::operator=(current_decorated_vertex);
    }

private:
    // The mesh that we are circulating on:
    Mesh_patch_c_ptr m_mesh_patch;

    // Vertex center of the circulation (+ circulator range for a border vertex)
    Vertex_c_handle m_center;

    // Internal circulator
    Adaptor_vertex_around_vertex_c_cir m_adaptor_circulator;

}; // Mesh_patch_vertex_around_vertex_cir


// Class Mesh_patch_vertex_around_facet_cir
// represents a (clockwise) circulator around a facet
// of a Mesh_adaptor_patch_3<MeshAdaptor_3> mesh

template<class MeshPatchPtrType,    // = [const] Mesh_adaptor_patch_3*
         class VertexHandleType,    // = Mesh_adaptor_patch_3::Vertex_[const_]handle
         class AdaptorVertexAroundFacetCirculatorType>
                                    // = Adaptor::Vertex_around_facet_[const_]circulator
class Mesh_patch_vertex_around_facet_cir
    : public VertexHandleType
{
// PRIVATE TYPES
// -------------

    // Base and self classes
    typedef VertexHandleType                    Base;
    typedef Mesh_patch_vertex_around_facet_cir  Self;

    // Mesh_patch types
    typedef typename std::iterator_traits<MeshPatchPtrType>::value_type
                                                Mesh_patch;
    typedef typename Mesh_patch::Adaptor        Adaptor;

public:

// PUBLIC TYPES
// ------------

    // Export template parameter types
    typedef MeshPatchPtrType                    Mesh_patch_c_ptr;
    typedef VertexHandleType                    Vertex_c_handle;
    typedef AdaptorVertexAroundFacetCirculatorType
                                                Adaptor_vertex_around_facet_c_cir;

    // Iterator types
    typedef typename Adaptor_vertex_around_facet_c_cir::iterator_category
                                                iterator_category;
    typedef typename Vertex_c_handle::value_type value_type;
    typedef std::ptrdiff_t                      difference_type;
    typedef std::size_t                         size_type;
    typedef typename Vertex_c_handle::reference reference;
    typedef typename Vertex_c_handle::pointer   pointer;

// CREATION
// --------

    // Circulator pointing to NULL
    Mesh_patch_vertex_around_facet_cir() {}

    Mesh_patch_vertex_around_facet_cir(Mesh_patch_c_ptr mesh,
                                       Adaptor_vertex_around_facet_c_cir adaptor_circulator)
        : m_mesh_patch(mesh),
          m_adaptor_circulator(adaptor_circulator)
    {
        CGAL_parameterization_assertion(m_mesh_patch != NULL);
        CGAL_parameterization_assertion(adaptor_circulator != NULL);

        // Update the inherited vertex handle
        update_inherited_handle();
    }

    // Copy constructor
    Mesh_patch_vertex_around_facet_cir(const Self& cir)
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
        Vertex_c_handle current_decorated_vertex = NULL;

        // Get next vertex on facet
        Adaptor_vertex_around_facet_c_cir next_vertex = m_adaptor_circulator;
        next_vertex++;

        // If (m_adaptor_circulator, next_vertex) isn't a seam (non-oriented) edge
        if (m_mesh_patch->m_mesh_adaptor->get_halfedge_seaming(
                        m_adaptor_circulator, next_vertex) != Mesh_patch::BORDER
         || m_mesh_patch->m_mesh_adaptor->get_halfedge_seaming(
                        next_vertex, m_adaptor_circulator) != Mesh_patch::BORDER)
        {
            current_decorated_vertex
                    = m_mesh_patch->get_decorated_inner_vertex(m_adaptor_circulator,
                                                               next_vertex);
        }
        // Special case: both vertices belong to the seam
        else 
        {
            current_decorated_vertex
                    = m_mesh_patch->get_decorated_border_vertex(m_adaptor_circulator,
                                                                next_vertex,
                                                                NULL);
        }

        // Update the inherited vertex handle
        assert(m_mesh_patch->is_valid(current_decorated_vertex));
        Base::operator=(current_decorated_vertex);
    }

private:
    // The mesh that we are circulating on:
    Mesh_patch_c_ptr m_mesh_patch;

    // Internal circulator
    Adaptor_vertex_around_facet_c_cir m_adaptor_circulator;

}; // Mesh_patch_vertex_around_facet_cir


CGAL_END_NAMESPACE

#endif //CGAL_MESH_ADAPTOR_PATCH_CIRCULATORS_H

