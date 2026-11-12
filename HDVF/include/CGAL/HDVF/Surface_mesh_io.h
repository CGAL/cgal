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
#include <CGAL/HDVF/Hdvf.h>

namespace CGAL {
namespace Homological_discrete_vector_field {

/*!
 \ingroup PkgHDVFIOClasses

 The class `Surface_mesh_io` is an intermediate %IO class, used to load the collection of simplices of a triangle mesh.  These collections are then used as input to produce `Abstract_simplicial_chain_complex`.

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

    /*! \brief Hole information returned by the  iterator.
     *
     * Information comprises:
     * - `cell_index` index of the critical cell.
     * - `dimension` dimension of the critical cell.
     * - `homology_chain`, `cohomology_chain` associated homology and/or cohomolgy generators.
     */
    template <typename CoefficientRing>
    struct Hole_information {
        halfedge_descriptor critical_cell; // Index of the critical cell
        int dimension; // Dimension of the critical cell
        Simplex_chain<CoefficientRing> homology_chain, cohomology_chain; // Homology and cohomology generators
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
     * Loads vertices, edges and triangles of a surface mesh `mesh` (called the underlying triangle mesh). Vertices and simplices are identified by their index together with their dimension. In a given a dimension, indices coincide with the enumeration order of corresponding surface mesh iterators.
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

    /** \brief */
    const TriangleMesh& triangle_mesh () const { return _sm; }

    Io_cell_type halfedge_handle_to_io_cell (halfedge_descriptor he, int q) const {
        return _he_index_to_io_cell.at(q).at(he);
    }

    halfedge_descriptor io_cell_to_halfedge_handle (Io_cell_type io_cell, int q) const {
        return _io_cell_to_he_index.at(io_cell);
    }

    /** \brief Returns the index of the cell of dimension `q` (in an `AbstractSimplicialChainComplex` `K` built from the `Surface_mesh_io`) associated to a given halfedge handle (in the underlying triangle mesh).
     *
     * \tparam _CT type of the complex `K` (a model of `AbstractSimplicialChainComplex`).
     *
     * \param he Halfedge handle considered.
     * \param q Dimension of the cell associated to the halfedge
     * \param K Abstract simplicial chain complex built from the `Surface_mesh_io`
     */
    template <typename _CT>
    size_t halfedge_handle_to_cell_index (halfedge_descriptor he, int q, const _CT& K) const {
        // Get the io_cell associated by SM to he in dimension q
        Io_cell_type io_cell(halfedge_handle_to_io_cell(he, q));
        // Build the corresponding simplex
        Simplex s(io_cell);
        // Get the index of s in K
        return K.cell_to_index(s);
    }

    /** \brief Returns the halfedge handle (in the underlying triangle mesh) associated to the cell of index `i` and dimension `q` (in an `AbstractSimplicialChainComplex` `K` built from the `Surface_mesh_io`).
     *
     * \tparam _CT type of the complex `K` (a model of `AbstractSimplicialChainComplex`).
     *
     * \param i Index of the cell of `K`
     * \param q Dimension of the cell of `K`
     * \param K Abstract simplicial chain complex built from the `Surface_mesh_io`
     */
    template <typename _CT>
    halfedge_descriptor cell_index_to_halfedge_handle (size_t i, int q, const _CT& K) const {
        // Get the cell associated to the index i in K
        Simplex simplex(K.index_to_cell(i,q));
        // Build the corresponding io_cell
        Io_cell_type io_cell(simplex.vertices());
        // Get the halfedge associated by SM to he in dimension q
        return io_cell_to_halfedge_handle(io_cell, q);
    }

    /** \brief Returns the `Simplex_chain` associated to a given `Sparse_chain`.
     *
     * Homology/cohomolgy generators computed by HDVFs are represented by sparse chains (where cells are identified by their dimension and index in the abstract simplicial chain complex). A `Simplex_chain` maps such a chain to simplices of a surface mesh (identified by their corresponding halfedge).
     *
     * \tparam _SCT Type of the the sparse chain `chain` (a model of `SparseChain`).
     * \tparam _CT Type of the complex `K` (a model of `AbstractSimplicialChainComplex`).
     *
     * \param chain Chain of `K`.
     * \param q Dimension of the chain
     * \param K Abstract simplicial chain complex built from the `Surface_mesh_io`
     */
    template <typename _SCT, typename _CT>
    Simplex_chain<typename _CT::Coefficient_ring> chain_to_simplex_chain (const _SCT& chain, int q, const _CT& K) {
        typedef typename _CT::Coefficient_ring Coefficient_ring;
        typedef Simplex_chain<Coefficient_ring> SC;
        SC res;
        // Visit chain and insert corresponding simplices/coefficients in res
        for (typename _SCT::const_iterator it = chain.cbegin(); it != chain.cend(); ++it) {
            const halfedge_descriptor he(cell_index_to_halfedge_handle(it->first, q, K));
            res[he] = it->second;
        }
        return res;
    }

    /*! \brief Iterator over homology/cohomology generators on the underlying surface mesh encoded as `Simplex_chain`.
     *
     * Given a perfect HDVF, iterates over critical cells and associated homology/cohomlogy generators.
     */
    template <typename HDVFType>
    struct iterator
    {
        // Iterator tags
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = Hole_information<typename HDVFType::Coefficient_ring>;

        using Coefficient_ring = typename HDVFType::Coefficient_ring;

        typedef Hole_information<typename HDVFType::Coefficient_ring> Hole_information;

        /*! \brief Iterator constructor
         *
         * \param hdvf Constant reference over the HDVF iterated.
         * \param dimension Dimension of iteration (if positive, iterate over a single dimension, if -1, iterate over all dimensions computed).
         * \param mode Iteration mode (iterate homology and/or cohomology generators). Default deduces the iteration mode from `hdvf` computation options.
         *
         * \exception Invalid_hdvf_options If computation options of the underlying `hdvf` are incompatible (no homology/cohomology generators computed), raises a `%std::runtime_error`
         */
        iterator(const Surface_mesh_io& sm_io, const HDVFType& hdvf, int dimension = -1, Iterator_mode mode = IT_NONE) : _sm_io(sm_io), _hdvf(hdvf), _it(hdvf.begin(dimension, mode)) {}

        // Operators
        /*! \brief Iterator dereference
         *
         * \returns A `Hole_information` structure containing all information related to the current persistence hole.
         */
        value_type operator*() const
        {
            value_type hdvf_it_value(*_it);
            Hole_information res ;
            const int q(hdvf_it_value.dimension);
            res.critical_cell = _sm_io.cell_index_to_halfedge_handle(hdvf_it_value.cell_index, q, _hdvf.complex());
            res.dimension = q;

            if (_it.mode() & HOMOLOGY)
            {
                res.homology_chain = _sm_io.chain_to_simplex_chain(hdvf_it_value.homology_chain, q, _hdvf.complex());
            }
            if (_it.mode() & COHOMOLOGY)
            {
                res.cohomology_chain = _sm_io.chain_to_simplex_chain(hdvf_it_value.cohomology_chain, q, _hdvf.complex());
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
            ++_it;
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
        friend bool operator== (const iterator& a, const iterator& b) { return (a._it == b._it); };

        /**
         * \brief Inequality check.
         * \returns `true` if the indices are different.
         */
        friend bool operator!= (const iterator& a, const iterator& b) { return (a._it != b._it); };

    private:
        const HDVFType& _hdvf; // HDVF iterated
        const Surface_mesh_io& _sm_io; // Corresponding Surface_mesh_io
        HDVFType::iterator& _it;
    };

    /**
     * \brief Iterator to the beginning of homology/cohomology generators.
     *
     * \param q If positive, the iterator visits only holes of dimension q, if negative (default), visits holes in any dimension computed (depending on dimension restriction)
     * \param mode Defines which information is returned by the iterator (HOMOLOGY, COHOMOLOGY or BOTH). This argument must be coherent with the HDVF option chosen for computation. By default, copies the current HDVF options.
     *
     * \returns The iterator to the first hole visited.
     */
    template <typename HDVFType>
    iterator<HDVFType> begin(const HDVFType& hdvf, int q = -1, Iterator_mode mode = Iterator_mode(0)) const {
        return iterator<HDVFType>(*this, hdvf, q, mode) ;
    }

    /**
     * \brief Iterator past-the-end of homology/cohomology generators.
     *
     * \returns The past-the-end iterator of  persistent intervals.
     */
    template <typename HDVFType>
    iterator<HDVFType> end(const HDVFType& hdvf) const { return iterator<HDVFType>(*this, hdvf, this->complex().dimension()+1, Iterator_mode(0)) ; }
} ;



} /* end namespace Homological_discrete_vector_field */
} /* end namespace CGAL */


#endif // CGAL_HDVF_SURFACE_MESH_IO_H
