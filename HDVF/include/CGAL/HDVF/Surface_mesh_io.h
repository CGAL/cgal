// Copyright (c) 2025 LIS Marseille (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alexandra Bac <alexandra.bac@univ-amu.fr>

#ifndef CGAL_HDVF_SURFACE_MESH_IO_H
#define CGAL_HDVF_SURFACE_MESH_IO_H

#include <CGAL/license/HDVF.h>

#include <iostream>
#include <cassert>
#include <vector>
#include <map>
#include <cmath>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/HDVF/Mesh_object_io.h>
#include <CGAL/HDVF/Abstract_simplicial_chain_complex.h>

namespace CGAL {
namespace Homological_discrete_vector_field {

/*!
 \ingroup PkgHDVFIOClasses

 The class `Surface_mesh_io` is an intermediate %IO class, used to load collections of 2-simplices from a triangle mesh.  These collections are then used as input to produce `Abstract_simplicial_chain_complex`.

 \tparam TriangleMesh a model of `FaceGraph` and `HalfedgeGraph` concepts, e.g., a `CGAL::Surface_mesh`.
 \tparam Traits a geometric traits class model of the `HDVFTraits` concept.
 */

template <typename TriangleMesh, typename Traits>
class Surface_mesh_io : public Mesh_object_io<Traits>
{
public:
    /** \brief Type of the surface mesh. */
    typedef TriangleMesh Surface_mesh;

    /** \brief Type of halfedges of the surface mesh. */
    typedef typename boost::graph_traits<Surface_mesh>::halfedge_descriptor halfedge_descriptor ;

    /** Type encoding chains of simplices (with coefficients in `CoefficientRing`, a model of `EuclideanRingWithoutDivision`). */
    template <typename CoefficientRing>
    using Simplex_chain = std::map<halfedge_descriptor, CoefficientRing>;

    /** \brief Iterator mode (determines if `iterator` enumerates homology/cohomology generators or both). */
    enum Iterator_mode {
        NONE = 0,
        HOMOLOGY = 0b0010,
        COHOMOLOGY = 0b0001,
        BOTH = 0b0011
    };

protected:
    typedef typename boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor ;
    typedef typename boost::graph_traits<Surface_mesh>::edge_descriptor edge_descriptor ;
    typedef typename boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor ;


    /* \brief Storage of vertices Io_cell_type <-> Vertex_index permutation.
     *
     * `_io_cell_to_he_index` maps a Io_cell_type to its associated Halfedge_index in the mesh
     *
     * `_he_index_to_io_cell.at(q)` stores, for each dimension `q`,  the map between Halfedge_index and Io_cell_type of cells of dimension `q`
     *
     * All cells, whatever their dimension, are identified by halfedges:
     * - a vertex is the `target` of the halfedge
     * - an edge is the `edge` of the halfedge
     * - a face is the `face` of the halfedge
     */
    std::map<Io_cell_type, halfedge_descriptor> _io_cell_to_he_index;
    std::vector<std::map<halfedge_descriptor, Io_cell_type> > _he_index_to_io_cell;

    // Underlying (triangular) surface mesh
    const TriangleMesh& _sm;

public:
    /** \brief Constructor from a triangular mesh.
     *
     * Loads vertices, edges and triangles of a surface mesh `mesh`. Vertices and simplices are identified by their index together with their dimension. In a given a dimension, indices coincide with the enumeration order of corresponding surface mesh iterators.
     *
     * As vertices are identified by their index, vertices with similar underlying points and different vertex indices will be considered as different.
     */
    Surface_mesh_io(const TriangleMesh& mesh) : Mesh_object_io<Traits>(2), _sm(mesh) {
        typedef typename Traits::Point Point;

        this->_nvertices = num_vertices(mesh);
        this->_ncells = num_vertices(mesh) + num_edges(mesh) + num_faces(mesh) ;
        _he_index_to_io_cell.resize(3);

        auto vpm = get(CGAL::vertex_point, mesh);
        // Load nodes
        for (vertex_descriptor v : vertices(mesh)) {
            Point p(get(vpm,v));
            this->_nodes.push_back(p);
        }
        // Load vertices
        size_t tmp_index(0);
        for (vertex_descriptor v : vertices(mesh)) {
            // Corresponding Io_cell_type
            Io_cell_type tmp_io_cell({tmp_index});
            // Add to cells
            this->_cells.push_back(tmp_io_cell);
            // Associated he
            halfedge_descriptor vertex_he(halfedge(v,mesh));
            // Store the permutation
            _io_cell_to_he_index[tmp_io_cell] = vertex_he;
            _he_index_to_io_cell.at(0)[vertex_he] = tmp_io_cell;
            ++tmp_index;
        }
        // Load edges
        for (edge_descriptor e : edges(mesh)) {
            // Associated he
            halfedge_descriptor edge_he(mesh.halfedge(e));
            // Compute corresponding Io_cell
            std::vector<size_t> tmp_cell;
            // Vertex1
            {
                // Get the halfedge "encoding" the target vertex
                halfedge_descriptor vert_he_ind(halfedge(target(edge_he, mesh), mesh));
                // Get corresponding io_cell
                Io_cell_type vert(_he_index_to_io_cell.at(0).at(vert_he_ind));
                assert(vert.size()==1); // vertex
                // Push_back the index
                tmp_cell.push_back(vert.at(0));
            }
            // Vertex2
            {
                // Get the halfedge "encoding" the source vertex
                halfedge_descriptor vert_he_ind(halfedge(target(opposite(edge_he, mesh), mesh), mesh));
                // Get the corresponding io_cell
                Io_cell_type vert(_he_index_to_io_cell.at(0).at(vert_he_ind));
                assert(vert.size()==1); // vertex
                // Push_back the index
                tmp_cell.push_back(vert.at(0));
            }
            std::sort(tmp_cell.begin(), tmp_cell.end());
            // Add to cells
            this->_cells.push_back(tmp_cell);
            // Set the permutation
            _io_cell_to_he_index[tmp_cell] = edge_he;
            _he_index_to_io_cell.at(1)[edge_he] = tmp_cell;
        }
        // Load faces
        for (face_descriptor f : faces(mesh)) {
            // Associated he
            halfedge_descriptor face_he(halfedge(f,mesh));
            // Compute corresponding Io_cell
            std::vector<size_t> tmp_cell;
            // Visit vertices around the face
            size_t cpt_verts(0);
            for(vertex_descriptor v : vertices_around_face(halfedge(f,mesh), mesh)){
                ++cpt_verts;
                // Get the halfedge stored in the vertex
                halfedge_descriptor vert_he_ind(halfedge(v,mesh));
                // Get the corresponding io_cell
                Io_cell_type vert(_he_index_to_io_cell.at(0).at(vert_he_ind));
                assert(vert.size()==1); // vertex
                // Push_back the index
                tmp_cell.push_back(vert.at(0));
            }
            // Faces must be triangles
            assert(cpt_verts == 3);
            std::sort(tmp_cell.begin(), tmp_cell.end());
            // Add to cells
            this->_cells.push_back(tmp_cell);
            // Set the permutation
            _io_cell_to_he_index[tmp_cell] = face_he;
            _he_index_to_io_cell.at(2)[face_he] = tmp_cell;
        }
    }

    Io_cell_type halfedge_handle_to_io_cell (halfedge_descriptor he, int q) const {
        return _he_index_to_io_cell.at(q).at(he);
    }

    halfedge_descriptor io_cell_to_halfedge_handle (Io_cell_type io_cell, int q) const {
        return _io_cell_to_he_index.at(io_cell);
    }

    /*! \brief Iterator over homology/cohomology generators on the underlying surface mesh encoded as `Simplex_chain`.
     *
     * Given a perfect HDVF, iterates over critical cells and associated homology/cohomlogy generators.
     */
    struct iterator
    {
        // Iterator tags
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = Hole_information;

        /*! \brief Iterator constructor
         *
         * \param hdvf Constant reference over the HDVF iterated.
         * \param dimension Dimension of iteration (if positive, iterate over a single dimension, if -1, iterate over all dimensions computed).
         * \param mode Iteration mode (iterate homology and/or cohomology generators). Default deduces the iteration mode from `hdvf` computation options.
         *
         * \exception Invalid_hdvf_options If computation options of the underlying `hdvf` are incompatible (no homology/cohomology generators computed), raises a `%std::runtime_error`
         */
        iterator(const Hdvf& hdvf, int dimension = -1, Iterator_mode mode = 0) : _i(0), _mode(mode), _hdvf(hdvf)
        {
            // Iteration mode
            if (_mode == 0) // Get mode from hdvf_opts
            {
                if (hdvf.hdvf_opts() == OPT_G)
                    _mode = HOMOLOGY;
                else if (hdvf.hdvf_opts() == OPT_F)
                    _mode = COHOMOLOGY;
                else if (hdvf.hdvf_opts() == OPT_FULL)
                    _mode = BOTH;
                else {
                    std::cerr << "Cannot iterate homology/cohomolgy generators: incompatible hdvf options" << std::endl;
                    throw std::runtime_error("Cannot iterate homology/cohomolgy generators: incompatible hdvf options");
                }
            }
            else { // Check the compatibility with HDVF options
                std::cout << "=== " << compatible_flags(_mode, hdvf.hdvf_opts()) << std::endl;
                if (!(compatible_flags(_mode, hdvf.hdvf_opts()))) { // _mode is incompatible with computed pieces of the reduction
                    std::cerr << "Cannot iterate homology/cohomolgy generators: incompatible hdvf options" << std::endl;
                    throw std::runtime_error("Cannot iterate homology/cohomolgy generators: incompatible hdvf options");
                }
            }

            // Dimension(s) and initialisation of criticals visited
            _criticals.resize(hdvf.complex().dimension()+1);
            if (dimension < 0)
            {
                for (int q = hdvf._min_dimension; q <= hdvf._max_dimension; ++q) {
                    _criticals[q] = hdvf.psc_flag(CRITICAL, q);
                }
                _dim = hdvf._min_dimension;
            }
            else {
                if ((dimension >= hdvf._min_dimension) && (dimension <= hdvf._max_dimension)) {
                    _criticals[dimension] = hdvf.psc_flag(CRITICAL, dimension);
                    _dim = dimension;
                }
                else
                    _dim = hdvf.complex().dimension()+1; // Past the end iterator
            }
        }

        // Operators
        /*! \brief Iterator dereference
         *
         * \returns A `Hole_information` structure containing all information related to the current persistence hole.
         */
        value_type operator*() const
        {
            Hole_information res ;
            const size_t sigma(_criticals.at(_dim).at(_i));
            res.cell_index = sigma;
            res.dimension = _dim;

            if (_mode & HOMOLOGY)
            {
                res.homology_chain = _hdvf.homology_chain(sigma, _dim);
            }
            if (_mode & COHOMOLOGY)
            {
                res.cohomology_chain = _hdvf.cohomology_chain(sigma, _dim);
            }
            return res ;
        }

        /**
         * \brief Prefix increment.
         *
         * \returns The reference to the current iterator.
         */
        iterator& operator++()
        {
            if (_dim <= _hdvf.complex().dimension()) {
                ++_i;
                if (_i >= _criticals.at(_dim).size()) {
                    _i = 0;
                    ++_dim;
                    // Skip dimensions without critical cells
                    while ((_dim <= _hdvf.complex().dimension()) && (_criticals.at(_dim).size() == 0))
                        ++_dim;
                }
            }
            return *this;
        }

        /**
         * \brief Postfix increment.
         * \returns The pre-increment iterator.
         */
        iterator operator++(int) { iterator tmp = *this; ++(*this); return tmp; }

        /**
         * \brief Equality check.
         * \returns `true` if the indices are equal.
         */
        friend bool operator== (const iterator& a, const iterator& b) { return ((a._i == b._i) && (a._dim == b._dim)); };

        /**
         * \brief Inequality check.
         * \returns `true` if the indices are different.
         */
        friend bool operator!= (const iterator& a, const iterator& b) { return ((a._i != b._i) || (a._dim != b._dim)); };

    private:
        std::vector<std::vector<size_t>> _criticals; // Criticals to visit
        size_t _i; // Current index along _criticals in dimension _dim
        int _dim; // Current dimension visited
        Iterator_mode _mode; // HDVF_opt considered
        const Hdvf& _hdvf; // HDVF iterated
    };

    /**
     * \brief Iterator to the beginning of homology/cohomology generators.
     *
     * \param q If positive, the iterator visits only holes of dimension q, if negative (default), visits holes in any dimension computed (depending on dimension restriction)
     * \param mode Defines which information is returned by the iterator (HOMOLOGY, COHOMOLOGY or BOTH). This argument must be coherent with the HDVF option chosen for computation. By default, copies the current HDVF options.
     *
     * \returns The iterator to the first hole visited.
     */
    iterator begin(int q = -1, Iterator_mode mode = Iterator_mode(0)) const { return iterator(*this, q, mode) ; }

    /**
     * \brief Iterator past-the-end of homology/cohomology generators.
     *
     * \returns The past-the-end iterator of  persistent intervals.
     */
    iterator end() const { return iterator(*this, this->complex().dimension()+1, Iterator_mode(0)) ; }
} ;

/** \relates Surface_mesh_io
 *
 * \brief Returns the index of the cell of dimension `q` (in an `AbstractSimplicialChainComplex` `K` built from the `Surface_mesh_io` `SM`) associated to a given halfedge handle (of the underlying surface mesh of `SM`).
 *
 * \tparam _TM type of the surface mesh `SM` (a model of `FaceGraph` and `HalfedgeGraph` concepts, e.g., a `CGAL::Surface_mesh`).
 * \tparam _Traits a geometric traits class model of the `HDVFTraits` concept.
 * \tparam _CT type of the complex `K` (a model of `AbstractSimplicialChainComplex`).
 *
 * \param he Halfedge handle considered.
 * \param q Dimension of the cell associated to the halfedge
 * \param SM `Surface_mesh_io`associated to the initial surface mesh
 * \param K Abstract simplicial chain complex associated to `SM`
 */
template <typename _TM, typename _Traits, typename _CT>
size_t halfedge_handle_to_cell_index (typename Surface_mesh_io<_TM, _Traits>::halfedge_descriptor he, int q, const Surface_mesh_io<_TM, _Traits>& SM, const _CT& K) {
    typedef Surface_mesh_io<_TM, _Traits> SM_io;
    // Get the io_cell associated by SM to he in dimension q
    Io_cell_type io_cell(SM.halfedge_handle_to_io_cell(he, q));
    // Build the corresponding simplex
    Simplex s(io_cell);
    // Get the index of s in K
    return K.cell_to_index(s);
}


/** \relates Surface_mesh_io
 *
 * \brief Returns the halfedge handle (in the underlying surface mesh of `SM`) associated to the cell of index `i` and dimension `q` (in an `AbstractSimplicialChainComplex` `K` built from the `Surface_mesh_io` `SM`).
 *
 * \tparam _TM type of the surface mesh `SM` (a model of `FaceGraph` and `HalfedgeGraph` concepts, e.g., a `CGAL::Surface_mesh`).
 * \tparam _Traits a geometric traits class model of the `HDVFTraits` concept.
 * \tparam _CT type of the complex `K` (a model of `AbstractSimplicialChainComplex`).
 *
 * \param i Index of the cell of `K`
 * \param q Dimension of the cell of `K`
 * \param SM `Surface_mesh_io`associated to the initial surface mesh
 * \param K Abstract simplicial chain complex associated to `SM`
 */
template <typename _TM, typename _Traits, typename _CT>
typename Surface_mesh_io<_TM, _Traits>::halfedge_descriptor cell_index_to_halfedge_handle (size_t i, int q, const Surface_mesh_io<_TM, _Traits>& SM, const _CT& K) {
    typedef Surface_mesh_io<_TM, _Traits> SM_io;
    // Get the cell associated to the index i in K
    Simplex simplex(K.index_to_cell(i,q));
    // Build the corresponding io_cell
    Io_cell_type io_cell(simplex.vertices());
    // Get the halfedge associated by SM to he in dimension q
    return SM.io_cell_to_halfedge_handle(io_cell, q);
}

/** \relates Surface_mesh_io
 *
 * \brief Returns the index in `K` of the cell of dimension `q` associated to a given halfedge handle in a complex `K` built from the `Surface_mesh_io` `SM`.
 *
 * \tparam _SCT Type of the the sparse chain `chain` (a model of `SparseChain`).
 * \tparam _SMIO Type of `SM` (an instance of `Surface_mesh_io`).
 * \tparam _CT Type of the complex `K` (a model of `AbstractSimplicialChainComplex`).
 *
 * \param chain Chain of `K`.
 * \param q Dimension of the chain
 * \param SM `Surface_mesh_io`associated to the initial surface mesh
 * \param K Abstract simplicial chain complex associated to `SM`
 */
template <typename _SCT, typename _SMIO, typename _CT>
_SMIO::template Simplex_chain<typename _CT::Coefficient_ring> chain_to_simplex_chain (const _SCT& chain, int q, const _SMIO& SM, const _CT& K) {
    typedef typename _CT::Coefficient_ring Coefficient_ring;
    typedef typename _SMIO::template Simplex_chain<Coefficient_ring> Simplex_chain;
    Simplex_chain res;
    // Visit chain and insert corresponding simplices/coefficients in res
    for (typename _SCT::const_iterator it = chain.cbegin(); it != chain.cend(); ++it) {
        const typename _SMIO::halfedge_descriptor he(cell_index_to_halfedge_handle(it->first, q, SM, K));
        res[he] = it->second;
    }
    return res;
}

} /* end namespace Homological_discrete_vector_field */
} /* end namespace CGAL */


#endif // CGAL_HDVF_SURFACE_MESH_IO_H
