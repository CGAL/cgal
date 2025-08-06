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

#ifndef CGAL_HDVF_ABSTRACT_SIMPLICIAL_CHAIN_COMPLEX_H
#define CGAL_HDVF_ABSTRACT_SIMPLICIAL_CHAIN_COMPLEX_H

#include <CGAL/license/HDVF.h>

#include <vector>
#include <map>
#include <CGAL/HDVF/Simplex.h>
#include <CGAL/HDVF/Mesh_object_io.h>
#include <CGAL/OSM/OSM.h>

namespace CGAL {
namespace HDVF {

// Forward declaration of SimpComplexTools
template<typename CoefficientType> class Duality_simplicial_complex_tools ;

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Abstract_simplicial_chain_complex` represents (topological) chain complexes associated to abstract simplicial complexes.

 An abstract simplicial complex is a set of simplices, also called cells (class `Simplex`) such that: all the faces of a given simplex also belong to the complex and any two simplices intersect exactly along a common face.

 <img src="simplicial_complex.png" align="center" width=45%/>

 A simplex of dimension q contains exactly q+1 vertices (we will thus denote it by \f$\langle v_0, \ldots, v_q \rangle\f$). A 0-cell is thus a vertex, a 1-cell contains two vertices (edge), a 2-cell contains three vertices (triangle) while a 3-cell is a tetrahedron.

 The boundary map of the complex is computed by the constructor of the class using the standard formula:
 \f[ \partial_q\left( \langle v_0, \ldots, v_q \rangle\right) = \sum_{i=0}^q (-1)^i\cdot \langle v_0, \ldots, \widehat{v_i}, \cdots, v_q \rangle\f]
 where \f$\langle v_0, \ldots, \widehat{v_i}, \cdots, v_q \rangle\f$ denotes the \f$q-1\f$-simplex with \f$v_i\f$ omitted. Hence, matrices of boundary maps are stored in each dimension using sparse matrices (in column-major mode).

 Let us also point out that cells are indexed along each dimension, thus each simplex is uniquely determined by its dimension and its index in this dimension.

 \cgalModels{AbstractChainComplex}

 \tparam CoefficientType a model of the `Ring` concept.
 */


template<typename CoefficientType>
class Abstract_simplicial_chain_complex {
public:
    /**
     * \brief %Default constructor (empty simplicial complex of dimension `q`).
     *
     * Builds an empty abstract simplicial complex of dimension `q`.
     */
    Abstract_simplicial_chain_complex(int q = 0) ;

    /**
     * \brief Constructor from a `Mesh_object_io`.
     *
     * Builds the abstract simplicial complex associated to a triangular mesh (i.e. performs the down closure of cells and set the boundary matrices in any dimension).
     *
     * \param[in] mesh A `Mesh_object_io` containing a triangular mesh.
     */
    Abstract_simplicial_chain_complex(const Mesh_object_io& mesh);

    /** \brief Type of column-major chains */
    typedef CGAL::OSM::Sparse_chain<CoefficientType, CGAL::OSM::COLUMN> Col_chain;
    /** \brief Type of row-major chains */
    typedef CGAL::OSM::Sparse_chain<CoefficientType, CGAL::OSM::ROW> Row_chain ;
    /** \brief Type of column-major sparse matrices */
    typedef CGAL::OSM::Sparse_matrix<CoefficientType, CGAL::OSM::COLUMN> Col_matrix;


    /**
     * \brief Assignment operator for abstract simplicial chain complexes.
     *
     * Stores a copy of an abstract simplicial chain complex in `*this`.
     *
     * \param[in] complex The abstract simplicial chain complex which will be copied.
     */
    Abstract_simplicial_chain_complex& operator= (const Abstract_simplicial_chain_complex& complex)
    {
        _dim = complex._dim;
        _ind2simp = complex._ind2simp;
        _simp2ind = complex._simp2ind;
        _nb_cells = complex._nb_cells;
        _d = complex._d;
        return *this ;
    }

    // Methods of the AbstractChainComplex concept

    /**
     * \brief Returns the boundary of the cell id_cell in dimension q.
     *
     * Returns a copy of the column-major chain stored in the boundary matrix of dimension dim: boundary of the cell id_cell in dimension q.
     *
     * \param[in] id_cell %Index of the cell.
     * \param[in] q Dimension of the cell.
     *
     * \return The column-major chain containing the boundary of the cell id_cell in dimension q.
     */
    Col_chain d(size_t id_cell, int q) const
    {
        if (q > 0)
            return OSM::get_column(_d[q], id_cell);
        else
            return Col_chain(0) ;
    }

    /**
     * \brief Returns the co-boundary of the cell id_cell in dimension q.
     *
     * Returns a row-major chain containing the co-boundary of the cell id_cell in dimension q (so actually a row of the boundary matrix).
     *
     * \warning As the boundary matrix is stored column-major, this entails crossing the full matrix to extract the row coefficients (O(number of non empty columns))
     *
     * \param[in] id_cell %Index of the cell.
     * \param[in] q Dimension of the cell.
     *
     * \return The row-major chain containing the co-boundary of the cell id_cell in dimension q.
     */
    Row_chain cod(size_t id_cell, int q) const
    {
        if (q < _dim)
            return OSM::get_row(_d[q+1], id_cell);
        else
            return Row_chain(0) ;
    }

    /**
     * \brief Returns the dimension of the complex.
     *
     * Returns the dimension of the simplicial complex (i.e. largest dimension of cells).
     *
     * \return The dimension of the complex..
     */
    int dim() const { return _dim ;}

    /**
     * \brief Returns the number of cells in a given dimension.
     *
     * \param[in] q Dimension along which the number of cells is returned.
     *
     * \return Number of cells in dimension q.
     */
    size_t nb_cells(int q) const
    {
        if ((q >=0) && (q<=_dim))
            return _nb_cells[q] ;
        else
            return 0 ;
    }

    /**
     * \brief Returns a constant reference to the vector of boundary matrices (along each dimension).
     *
     * Returns a constant reference to the vector of boundary matrices along each dimension. The q-th element of this vector is a column-major sparse matrix containing the boundaries of q-cells (i.e. rows encode q-1 cells and columns q cells).
     *
     * \return Returns a constant reference to the vector of column-major boundary matrices along each dimension.
     */
    const vector<Col_matrix> & get_boundary_matrices() const
    {
        return _d ;
    }

    /**
     * \brief Returns a copy of the dim-th boundary matrix (i.e.\ column-major matrix of \f$\partial_q\f$).
     *
     * It is a column-major sparse matrix containing the boundaries of q-cells (i.e. rows encode q-1 cells and columns q cells).
     *
     * \param[in] q Dimension of the boundary matrix (i.e. columns will contain the boundary of dimension q cells).
     *
     * \return A column-major sparse matrix containing the matrix of the boundary operator of dimension q.
     */
    const Col_matrix & get_boundary_matrix(int q) const
    {
        return _d.at(q) ;
    }

    /**
     * \brief Returns dimension 0 simplices indexes included in the cell with index id_cell of dimension q.
     *
     * Returns the dimension 0 simplices indexes included in the cell with index id_cell of dimension q.
     *
     * \warning This does not come to return vertices indices, as dimension 0 simplices enumerate vertices in any order. For instance, if an abstract simplicial complex is built from 3 vertices {1,2,3} such that the enumeration of dimension 0 simplices is:
     *  id0: 3, id1 : 2, id2: 1
     * then the bottom_faces of the 1-simplex {1,2} are two 0-simplices with id 2 and 1.
     *
     * \param[in] id_cell %Index of the cell.
     * \param[in] q Dimension of the cell.
     *
     * \return A vector of 0-simplices indexes.
     */
    std::vector<size_t> bottom_faces(size_t id_cell, int q) const
    {
        std::set<size_t> verts(_ind2simp.at(q).at(id_cell).get_vertices()) ;
        std::vector<size_t> res ;
        // For each vertex in verts, compute the corresponding dimension 0 cell
        for (size_t vert_id : verts)
        {
            const size_t i(_simp2ind.at(0).at(Simplex(std::set<size_t>({vert_id})))) ;
            res.push_back(i) ;
        }
        return res ;
    }

    /*!
     * \brief Returns the cofaces of a given chain in dimension `q`.
     *
     * The resulting chain lies in dimension `q`+1 and is null if this dimension exceeds the dimension of the complex.
    */
    template <typename CoefficientT, int ChainTypeF>
    Col_chain cofaces_chain (OSM::Sparse_chain<CoefficientT, ChainTypeF> chain, int q) const
    {
        typedef OSM::Sparse_chain<CoefficientT, ChainTypeF> ChainType;
        // Compute the cofaces
        if (q < dim())
        {
            Col_chain fstar_cofaces(nb_cells(q+1)) ;
            for (typename ChainType::const_iterator it = chain.cbegin(); it != chain.cend(); ++it)
            {
                // Set the cofaces of it->first in dimension dim+1
                Row_chain cofaces(cod(it->first,q)) ;
                for (typename Row_chain::const_iterator it2 =  cofaces.cbegin(); it2 != cofaces.cend(); ++it2)
                    fstar_cofaces.set_coef(it2->first, 1) ;
            }
            return fstar_cofaces ;
        }
        else
            return Col_chain(0) ;
    }

protected:
    /*
     * \brief Prints informations on the complex.
     *
     * Displays the number of cells in each dimension and the boundary matrix in each dimension.
     */
    std::ostream& print_complex(std::ostream& out = std::cout) const;
public:
    /**
     * \brief Prints informations on the complex.
     *
     * Displays the number of cells in each dimension and the boundary matrix in each dimension.
     */
    template <typename _CT>
    friend ostream& operator<< (ostream& out, const Abstract_simplicial_chain_complex<_CT>& complex);

    /** \brief Get (unique) object Id.
     * For comparison of constant references to the complex.
     */
    size_t get_id () const { return _complex_id; }

protected:
    /* \brief Dimension of the complex */
    int _dim;
    /* \brief Vector of simplices along each dimension: _ind2simp.at(q) contains the vector of all simplices of dimension q */
    std::vector<std::vector<Simplex> > _ind2simp ;
    /* \brief Vector of maps associating indices to simplices in each dimension: _simp2ind.at(q) maps q-simplices to their index in _ind2simp  */
    std::vector<std::map<Simplex, size_t> > _simp2ind ;
    /* \brief Vector of number of cells in each dimension */
    std::vector<size_t> _nb_cells ;
    /* \brief Vector of column-major boundary matrices: _d.at(q) is the matrix of the boundary operator in dimension q. */
    mutable std::vector<Col_matrix>  _d;  // Boundary matrices

    // Protected methods

    /*
     * \brief Method filling the matrix _d.at(q), i.e. computing the boundaries of cells of dimension q.
     *
     * Compute the boundary with the standard simplicial complex definition:
     * \f[\partial (\{v_0,\ldots,v_q\}) = \sum_{i=0}^q (-1)^q \{v_0,\ldots, \hat{v_i},\ldots,v_q\}\f]
     *
     * \param[in] q Dimension considered for computation.
     */
    void  calculate_d(int q) const;

    /*
     * \brief Method inserting a simplex (and its faces if necessary) into the abstract simplicial complex.
     *
     * The method (recursively) inserts the simplex into _ind2simp and _simp2ind and its faces. If the simplex is already presend, the method does nothing.
     *
     * \warning insert_simplex does not update the boundary matrix.
     *
     * \param[in] tau The simplex inserted.
     */
    void insert_simplex(const Simplex& tau);

private:
    /* \brief Static counter for objects ids.
     * Initialized to 0.
     */
    static size_t _id_generator ; // Initialisation 0
    /* \brief Unique object id (for comparison of constant references to the complex). */
    const size_t _complex_id ;
};

// Initialization of _id_generator
template <typename CoefficientType>
size_t Abstract_simplicial_chain_complex<CoefficientType>::_id_generator(0);

//constructors
template<typename CoefficientType>
Abstract_simplicial_chain_complex<CoefficientType>::Abstract_simplicial_chain_complex(int q) : _complex_id(_id_generator++) {
    // Initialize attributes
    _dim = q;

    // Initialize vectors of Simplices and cell counts
    _ind2simp.resize(_dim + 1);
    _simp2ind.resize(_dim + 1);
    _nb_cells.resize(_dim + 1, 0);
    // Initialize boundary matrix
    _d.resize(_dim+1);
}

template<typename CoefficientType>
Abstract_simplicial_chain_complex<CoefficientType>::Abstract_simplicial_chain_complex(const Mesh_object_io& mesh) : _complex_id(_id_generator++) {
    // Initialize attributes

    _dim = abs(mesh.dim);

    // Initialize vectors of Simplices and cell counts
    _ind2simp.resize(_dim + 1);
    _simp2ind.resize(_dim + 1);
    _nb_cells.resize(_dim + 1, 0);

    // Iterate through the mesh cells and add them to the complex
    for (const auto& cell : mesh.cells) {
        //        Simplex simplex(std::set<size_t>(cell.begin(), cell.end()));
        //        insert_simplex(simplex);
        insert_simplex(cell) ;
    }

    _d.resize(_dim+1);
    for (int dim = 0; dim <= _dim; ++dim) {
        calculate_d(dim);
    }
}

// Recursive insertion method
template<typename CoefficientType>
void Abstract_simplicial_chain_complex<CoefficientType>::insert_simplex(const Simplex& tau) {
    int q = tau.dimension();
    if (q == -1) return;

    if (_simp2ind[q].find(tau) == _simp2ind[q].end()) {
        size_t i = _ind2simp[q].size();
        _ind2simp[q].push_back(tau);
        _simp2ind[q][tau] = i;
        _nb_cells[q]++;

        std::vector<Simplex> bord = tau.boundary();
        for (const auto& sigma : bord) {
            insert_simplex(sigma);
        }
    }
}

// calculate _d boundary matrix
template<typename CoefficientType>
void Abstract_simplicial_chain_complex<CoefficientType>::calculate_d(int dim) const {
    size_t nb_lignes = (dim == 0) ? 0 : _nb_cells[dim - 1];
    _d[dim] = Col_matrix(nb_lignes, _nb_cells[dim]);

    // Iterate through the cells of dimension dim
    if (dim>0)
    {
        for (size_t i = 0; i < _nb_cells[dim]; ++i) {
            // Boundary of the i-th simplex of dimension dim
            const Simplex& s = _ind2simp[dim][i];
            std::vector<Simplex> bord = s.boundary();

            // Create a chain with the correct size
            Col_chain chain(nb_lignes);

            // For each element of the boundary
            for (size_t j = 0; j < bord.size(); ++j) {
                auto it = _simp2ind[dim - 1].find(bord[j]); // Find the index of Simplex j
                if (it != _simp2ind[dim - 1].end()) { // If Simplex j is found
                    size_t ind_j = it->second; // Retrieve the index of Simplex j
                    chain.set_coef(ind_j, (j % 2 == 0) ? 1 : -1);
                }
                else
                    throw "calculate_d boundary simplex not found!";
            }

            // Insert the chain into the corresponding column of the delta matrix
            OSM::set_column(_d[dim], i, chain);
        }
    }
}

// Method to display the complex's information
template<typename CoefficientType>
std::ostream& Abstract_simplicial_chain_complex<CoefficientType>::print_complex(std::ostream& out) const {
    out << "Complex dimension: " << _dim << std::endl;

    // Total number of cells
    size_t nb_total_cells = 0;
    for (size_t i = 0; i <= _dim; ++i) {
        nb_total_cells += _nb_cells[i];
    }
    out << "Total number of cells: " << nb_total_cells << std::endl;

    // Cells per dimension
    for (int q = 0; q <= _dim; ++q) {
        out << "--- dimension " << q << std::endl;
        out << nb_cells(q) << " cells" << std::endl ;
        //        for (size_t j = 0; j < _nb_cells.at(q); ++j) {
        //            Simplex s(_ind2Simp.at(q).at(j));
        //            std::cout << j << " -> " << s << " -> " << _Simp2ind.at(q).at(s) << std::std::endl;
        //        }
    }

    // Boundary matrices
    out << "---------------------------" << std::endl << "Boundary matrices" << std::endl;
    for (int q = 1; q <= _dim; ++q)
        out << "_d[" << q << "] : " << _d[q].dimensions().first << "x" << _d[q].dimensions().second << std::endl <<  _d[q] << std::endl;
    return out ;
}

template <typename CoefficientType>
std::ostream& operator<< (std::ostream& out, const Abstract_simplicial_chain_complex<CoefficientType>& complex)
{
    return complex.print_complex(out);
}

} /* end namespace HDVF */
} /* end namespace CGAL */

#endif // CGAL_HDVF_ABSTRACT_SIMPLICIAL_CHAIN_COMPLEX_H
