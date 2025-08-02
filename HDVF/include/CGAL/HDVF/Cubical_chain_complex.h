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

#ifndef CGAL_HDVF_CUBICAL_CHAIN_COMPLEX_H
#define CGAL_HDVF_CUBICAL_CHAIN_COMPLEX_H

#include <vector>
#include <map>
#include <stdexcept>
#include <unordered_set>
#include <CGAL/HDVF/Cub_object_io.h>
#include <CGAL/OSM/OSM.h>

namespace CGAL {
namespace HDVF {

// Forward declaration of Duality_cubical_complex_tools
template<typename T> class Duality_cubical_complex_tools ;

// TOOLS

ostream & operator<<(ostream & out, std::vector<size_t> c)
{
    for (size_t i:c)
        out << i << " " ;
    out << endl ;
    return out ;
}

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Cubical_chain_complex` represents (topological) chain complexes associated to cubical complexes.

 \section Description Description

 A cubical complex is a set of "square" cells such that: all the faces of a given cell also belong to the complex and any two cells intersect exactly along a common face.

 <img src="cubical_complex.png" align="center" width=20%/>

 A cell in a cubical complex of dimension q is the %cartesian product of \f$q\f$ intervals \f$[n_1,n_1+\delta_1]\times [n_q,n_q+\delta_q] \f$ with \f$n_i\in \mathbb N \f$ and \f$\delta_i = 0,1\f$.
 The dimension of such a cell is \f$\sum_{i=1}^q \delta_i\f$.
 A 0-cell is thus a vertex, a 1-cell contains one non zero delta (edge along one of the axes), a 2-cell contains two non zero deltas (square along two of the axes), while a 3-cell is a cube...

 For instance, \f$v = [2,2]\times[1,1]\times[0,0]\f$, \f$e = [0,1]\times[1,1]\times[0,0]\f$ and \f$f = [2,3]\times[0,1]\times[1,1]\f$.

 Given a cubical complex of dimension \f$q\f$ included in \f$[0,B_1]\times\cdots \times[0,B_q]\f$, Khalimsky coordinates associate to each cell of the complex a unique coordinate in \f$[0,2B_1]\times\cdots \times[0,2B_q]\f$:
 \f[
 \begin{array}{rcl}
 \sigma & \mapsto & \mathrm{khal}(\sigma) \\
 [n_1,n_1+\delta_1]\times\cdots \times [n_q,n_q+\delta_q] & \,\to\, & (2n_1+\delta_1,\ldots,2n_1+\delta_1) \\
 \end{array}\f]
 we set \f$N_i = 2B_i+1\f$ the size of the Khalimsky bounding box along the \f$i\f$th axis.

 The dimension of a cell is the given by the number of odd coordinates.
 In our previous example, \f$B_1=3, B_2=2, B_3=1\f$, hence \f$N_1=7, N_2=5, N_3=3\f$ and we get the following Khalimsky coordinates: \f$\mathrm{khal}(v) = (4,2,0)\f$, \f$\mathrm{khal}(e) = (1,2,0)\f$ and \f$\mathrm{khal}(f) = (5,1,2)\f$.

 The boundary map of the complex is computed by the constructor of the class using the standard formula with Khalimsky coordinates. Given a cell \f$\sigma\f$ with \f$\mathrm{khal}(\sigma) = (x_1,\ldots,x_q)\f$
 \f[ \partial(x_1,\ldots,x_q) = \sum_{\substack{i=0\\x_i\text{ odd}}}^q (-1)^{j(i)} \left( (x_1,\ldots,x_i+1,\ldots, x_q) - (x_1,\ldots,x_i-1,\ldots, x_q)\right)\f]
 where \f$j(i)\f$ is the number of odd coordinates between \f$x_1\f$ and \f$x_i\f$.



 Let us also point out that besides Khalimsky coordinates, cells are indexed along each dimension, thus each cell is uniquely determined by its dimension and its index in this dimension (called "base index").

 \section SecHdvfImplementation Implementation Details

 As described above, Khalimsky coordinates provide a convenient tool to identify cells of any dimension. Hence a complex of dimension \f$q\f$ can be encoded by a boolean matrix of size \f$N_1\times\cdots\times N_q\f$ (with previous notations). For convenience, this matrix is vectorized and thus the complex is stored in a boolean vector (denoted by  `_cells`) of size \f$N_1\times \cdots \times N_q\f$.

 Cells of any dimension are thus repesented by a given element of this boolean vector and the corresponding index is called their *boolean index*.

 As stated above, besides this boolean representation, topological computations require to identify the bases of cells in any dimension (bases of the free chain groups). Hence, a cell of dimension \f$d\f$ is also identified by its index in the basis of \f$d\f$-dimensional cells. This index is called its *basis index*. The vector `_bool2base` stores, for each dimension, the map between boolean and base indices, while the `_base2bool` stores, for each dimension, the permutation between base and boolean indices.

 <img src="cubical_implementation.png" align="center" width=50%/>

 \cgalModels{GeometricChainComplex}

 \tparam CoefficientType a model of the `Ring` concept.
 */


template<typename CoefficientType>
class Cubical_chain_complex {
public:
    /** \brief Type of vertices coordinates */
    typedef std::vector<double> Point ;
private:
    /* \brief Vector of VTK types associated to cells in each dimension
     e.g. {1, 3, 8, 11} */
    static const std::vector<int> VTK_cubtypes ;

public:
    /** \brief Type used to encode primal or dual construction. */
    enum typeComplexCube {PRIMAL, DUAL};

    /**
     * \brief Default constructor (empty cubical complex).
     *
     * Builds an empty cubical complex.
     */
    Cubical_chain_complex() : _complex_id(_id_generator++) {} ;

    /**
     * \brief Constructor from a Cub_object_io (builds PRIMAL or DUAL associated complex depending on `type`).
     *
     * Builds the cubical complex associated to a a set of cells (vertices, edges, squares, cubes...), ie. performs the down closure of cells and set the boundary matrices in any dimension. Given a set of cells:
     *
     * - if the `type` is PRIMAL, the constructor builds the associated complex as such (see below middle), which comes to encode \f$3^q-1\f$ connectivity (with \f$q\f$ the dimension of the complex)
     * - if the `type` is DUAL and the Cub_object_io contains only cells of maximal dimension (ie. binary object), the constructor build the dual associated complex (see below right), which comes to encode \f$2q\f$ connectivity (with \f$q\f$ the dimension of the complex)
     *
     *<img src="primal_dual.png" align="center" width=20%/>
     *
     * \param[in] cub A Cub_object_io containing a set of "cubical" cells.
     * \param[in] type Type of construction used (PRIMAL or DUAL).
     */
    Cubical_chain_complex(const Cub_object_io& cub,typeComplexCube type);

    /** \brief Friend class `Duality_cubical_complex_tools` provides tools for Alexander duality. */
    friend Duality_cubical_complex_tools<CoefficientType> ;

    /** \brief Type of column-major chains */
    typedef OSM::Sparse_chain<CoefficientType, OSM::COLUMN> CChain;
    /** \brief Type of row-major chains */
    typedef OSM::Sparse_chain<CoefficientType, OSM::ROW> RChain ;
    /** \brief Type of column-major sparse matrices */
    typedef OSM::Sparse_matrix<CoefficientType, OSM::COLUMN> CMatrix;

    /**
     * \brief Affectation operator for cubical chain complexes.
     *
     * Stores a copy of an cubical chain complex in *this.
     *
     * \param[in] complex The cubical chain complex which will be copied.
     */
    Cubical_chain_complex& operator=(const Cubical_chain_complex& complex)
    {
        _dim = complex._dim;
        _size_bb = complex._size_bb;
        _P = complex._P;
        _cells = complex._cells;
        _base2bool = complex._base2bool;
        _bool2base = complex._bool2base;
        _visited_cells = complex._visited_cells;
        return *this;
    }

    /// Methods of the CubComplex concept

    /**
     * \brief Method returning the boundary of the cell id_cell in dimension q.
     *
     * Returns a copy of the column-major chain stored in the boundary matrix of dimension q: boundary of the cell id_cell in dimension q.
     *
     * \param[in] id_cell Index of the cell.
     * \param[in] q Dimension of the cell.
     *
     * \return The column-major chain containing the boundary of the cell id_cell in dimension q.
     */
    CChain d(size_t id_cell, int q) const
    {
        if (q > 0)
            return OSM::cget_column(_d[q], id_cell);
        else
            return CChain(0) ;
    }

    /**
     * \brief Method returning the co-boundary of the cell id_cell in dimension q.
     *
     * Returns a row-major chain containing the co-boundary of the cell id_cell in dimension q (so actually a row of the boundary matrix).
     *
     * \warning As the boundary matrix is stored column-major, this entails crossing the full matrix to extract the row coefficients (O(number of non empty columns))
     *
     * \param[in] id_cell Index of the cell.
     * \param[in] q Dimension of the cell.
     *
     * \return The row-major chain containing the co-boundary of the cell id_cell in dimension q.
     */
    RChain cod(size_t id_cell, int q) const
    {
        if (q < _dim)
            return OSM::get_row(_d[q+1], id_cell);
        else
            return RChain(0) ;
    }

    /**
     * \brief Method returning the dimension of the complex.
     *
     * Returns the dimension of the cubical complex.
     *
     * \return The dimension of the complex..
     */
    int dim() const
    {
        return _dim ;
    }

    /**
     * \brief Method returning the number of cells in a given dimension.
     *
     * \param[in] q Dimension along which the number of cells is returned.
     *
     * \return Number of cells in dimension q.
     */
    size_t nb_cells(int q) const
    {
        if ((q >=0) && (q <= _dim))
            return _base2bool.at(q).size() ;
        else
            return 0 ;
    }

    /**
     * \brief Method returning a constant reference to the vector of boundary matrices (along each dimension).
     *
     * Returns a constant reference to the vector of boundary matrices along each dimension. The q-th element of this vector is a column-major sparse matrix containing the boundaries of q-cells (ie. rows encode q-1 cells and columns q cells).
     *
     * \return Returns a constant reference to the vector of column-major boundary matrices along each dimension.
     */
    const vector<CMatrix> & get_bnd_matrices() const
    {
        return _d ;
    }

    /**
     * \brief Method returning a copy of the dim-th boundary matrix (ie. column-major matrix of \f$\partial_q\f$).
     *
     * It is a column-major sparse matrix containing the boundaries of dim-cells (ie. rows encode q-1 cells and columns q cells).
     *
     * \param[in] q Dimension of the boundary matrix (ie. columns will contain the boundary of dimension q cells).
     *
     * \return A column-major sparse matrix containing the matrix of the boundary operator of dimension q.
     */
    const CMatrix & get_bnd_matrix(int q) const
    {
        return _d.at(q) ;
    }


    /**
     * \brief Method returning dimension 0 cells indexes included in the cell with index id_cell of dimension q.
     *
     * Returns the dimension 0 vertices indexes included in the cell with index id_cell of dimension q.
     *
     * \param[in] id_cell Index of the cell.
     * \param[in] q Dimension of the cell.
     *
     * \return A vector of 0-cells indices.
     */
    std::vector<size_t> bottom_faces(size_t id_cell, int q) const
    {
        // Khalimsky coordinates of the cell
        const vector<size_t> coords(ind2khal(_base2bool.at(q).at(id_cell))) ;
        return khal_to_verts(coords) ;
    }

    /**
     * \brief Returns the cofaces of a given chain in dimension `q`.
     *
     * The resulting chain lies in dimension `q`+1 and is null if this dimension exceeds the dimension of the complex.
    */
    template <typename CoefficientT, int ChainTypeF>
    CChain cofaces_chain (OSM::Sparse_chain<CoefficientT, ChainTypeF> chain, int q) const
    {
        typedef OSM::Sparse_chain<CoefficientT, ChainTypeF> ChainType;
        // Compute the cofaces
        if (q < dim())
        {
            CChain fstar_cofaces(nb_cells(q+1)) ;
            for (typename ChainType::const_iterator it = chain.cbegin(); it != chain.cend(); ++it)
            {
                // Set the cofaces of it->first in dimension dim+1
                RChain cofaces(cod(it->first,q)) ;
                for (typename RChain::const_iterator it2 =  cofaces.cbegin(); it2 != cofaces.cend(); ++it2)
                    fstar_cofaces.set_coef(it2->first, 1) ;
            }
            return fstar_cofaces ;
        }
        else
            return CChain(0) ;
    }

    /**
     * \brief Method printing informations of the complex.
     *
     * Displays the number of cells in each dimension and the boundary matrix in each dimension.
     */
    void print_complex() const {
        for (int q = 0; q <= _dim; ++q) {
            std::cout << "-------- dimension " << q << std::endl;
            std::cout << "cellules de dimension " << q << " : " << _base2bool.at(q).size() << std::endl;
            for (size_t id_base = 0; id_base < _base2bool.at(q).size(); ++id_base) {
                size_t id_bool = _base2bool.at(q).at(id_base);
                std::vector<size_t> khal = ind2khal(id_bool);
                std::cout << id_base << " -> " << id_bool << " -> " << _bool2base.at(q).at(id_bool) << " -> ";
                for (size_t k : khal) std::cout << k << " ";
                std::cout << std::endl;
            }

            if (_base2bool[q].size() > 0 && q <= _dim) {
                std::cout << "matrice de bord de dimension " << q << std::endl;
                std::cout << _d[q] << std::endl;
            }
        }
    }

    /** \brief Get (unique) object Id.
     * For comparison of constant references to the complex.
     */
    size_t get_id () const { return _complex_id; }

    /** \brief Get the coordinates of the ith vertex */

    Point get_vertex_coords(size_t i) const
    {
        const vector<size_t> coords(ind2khal(_base2bool.at(0).at(i))) ;
        vector<double> res ;
        for (size_t c : coords)
            res.push_back(c/2. + .5) ;
        for (size_t i = coords.size(); i<3; ++i) // points must be 3D
            res.push_back(0) ;
        return res ;
    }

    /** \brief Get the vector of vertices coordinates  */
    const std::vector<Point>& get_vertices_coords() const
    {
        std::vector<Point> res ;
        for (size_t i=0; i<nb_cells(0); ++i)
            res.push_back(get_vertex_coords(i)) ;
        return res ;
    }

    /// END Methods of the Cubical_chain_complex concept

    // VTK export

    /**
     * \brief Method exporting a cubical complex (plus, optionally, labels) to a VTK file.
     *
     * The method generates legacy text VTK files. Labels are exported as such in a VTK property, together with CellID property, containing the index of each cell.
     *
     * \tparam LabelType Type of labels provided (default: int).
     *
     * \param[in] K Cubical complex exported.
     * \param[in] filename Output file root (output filenames will be built from this root).
     * \param[in] labels Pointer to a vector of labels in each dimension. (*labels).at(q) is the set of integer labels of cells of dimension q. If labels is NULL, only CellID property is exported.
     * \param[in] label_type_name Typename used in vtk export (e.g. "int" or "unsigned_long", see <a href = "https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html">VTK manual </a>).
     */
    template <typename LabelType = int>
    static void chain_complex_to_vtk(const Cubical_chain_complex<CoefficientType> &K, const std::string &filename, const std::vector<std::vector<LabelType> > *labels=NULL, string label_type_name = "int")
    {
        bool with_scalars = (labels != NULL) ;

        // Load out file...
        std::ofstream out ( filename, std::ios::out | std::ios::trunc);

        if ( not out . good () ) {
            std::cerr << "CubComplex_to_vtk. Fatal Error:\n  " << filename << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }

        // Header
        out << "# vtk DataFile Version 2.0" << endl ;
        out << "generators" << endl ;
        out << "ASCII" << endl ;
        out << "DATASET  UNSTRUCTURED_GRID" << endl ;

        // Points
        size_t nnodes = K.nb_cells(0) ;
        out << "POINTS " << nnodes << " double" << endl ;
        for (size_t n = 0; n < nnodes; ++n)
        {
            vector<double> p(K.get_vertex_coords(n)) ;
            for (double x : p)
                out << x << " " ;
            for (size_t i = p.size(); i<3; ++i) // points must be 3D -> add zeros
                out << "0 " ;
            out << endl ;
        }

        size_t ncells_tot = 0, size_cells_tot = 0 ;
        std::vector<int> types ;
        std::vector<LabelType> scalars ;
        std::vector<size_t> ids ;
        // all cells must be printed
        {
            // Cells up to dimension 3
            // Size : size of a cell of dimension q : 2^q
            for (int q=0; q<=K.dim(); ++q)
            {
                ncells_tot += K.nb_cells(q) ;
                const size_t size_cell = 1<<q ;
                size_cells_tot += (size_cell+1)*K.nb_cells(q) ;
            }
            out << "CELLS " << ncells_tot << " " << size_cells_tot << endl ;
            // Output cells by increasing dimension

            // Vertices
            for (size_t i = 0; i<K.nb_cells(0); ++i)
            {
                out << "1 " << i << endl ;
                types.push_back(Cubical_chain_complex<CoefficientType>::VTK_cubtypes.at(0)) ;
                if (with_scalars)
                {
                    scalars.push_back((*labels).at(0).at(i)) ;
                    ids.push_back(i) ;
                }
            }
            // Cells of higher dimension
            for (int q=1; q<=K.dim(); ++q)
            {
                const size_t size_cell = 1<<q ; //int_exp(2, q) ;
                for (size_t id =0; id < K.nb_cells(q); ++id)
                {
                    vector<size_t> verts(K.khal_to_verts(K.ind2khal(K._base2bool.at(q).at(id)))) ;
                    out << size_cell << " " ;
                    for (size_t i : verts)
                        out << i << " " ;
                    out << endl ;
                    types.push_back(Cubical_chain_complex<CoefficientType>::VTK_cubtypes.at(q)) ;
                    if (with_scalars)
                    {
                        scalars.push_back((*labels).at(q).at(id)) ;
                        ids.push_back(id) ;
                    }
                }
            }

            // CELL_TYPES
            out << "CELL_TYPES " << ncells_tot << endl ;
            for (int t : types)
                out << t << " " ;
            out << endl ;
        }

        if (with_scalars)
        {
            // CELL_LABEL
            out << "CELL_DATA " << ncells_tot << endl ;
            out << "SCALARS Label " << label_type_name << " 1" << endl ;
            out << "LOOKUP_TABLE default" << endl ;
            for (LabelType s : scalars)
                out << s << " " ;
            out << endl ;
            // CELL_IDs
            out << "SCALARS CellId " << "int" << " 1" << endl ;
            out << "LOOKUP_TABLE default" << endl ;
            for (size_t i : ids)
                out << i << " " ;
            out << endl ;
        }
        out.close() ;
    }

    /**
     * \brief Method exporting a chain over a cubical complex to a VTK file.
     *
     * The method generates legacy text VTK files. All the cells of the chain with non zero coefficient are exported. If a cellId is provided, labels are exported in a VTK property (2 for all cells, 0 for cell of index cellId).  The index of each cell is exported in a CellID property.
     *
     * \param[in] K Cubical complex exported.
     * \param[in] filename Output file root (output filenames will be built from this root).
     * \param[in] chain Sparse_chain exported (all the cells with non-zero coefficients in the chain are exported to vtk).
     * \param[in] q Dimension of the cells of the chain.
     * \param[in] cellId If different from MAX_SIZE_T, labels are exported to distinguish cells of the chain (label 2) from cellId cell (label 0).
     */
    static void chain_complex_chain_to_vtk(const Cubical_chain_complex<CoefficientType> &K, const std::string &filename, const OSM::Sparse_chain<CoefficientType, OSM::COLUMN>& chain, int q, size_t cellId = -1) ;

protected:
    // Methods to access data
    /* \brief Get the size of the Khalimsky bounding box \f$(N_1,\ldots,N_{\mathrm{\_dim}})\f$ */
    std::vector<size_t> get_size_bb() const {
        return _size_bb;
    }

    /* \brief Get _P (coefficients used for the vectorisation of the Khalimsky boolean matrix) */
    std::vector<size_t> get_P() const {
        return _P;
    }

    /* \brief Get the boolean vector of the complex in Khalimsky coordinates */
    std::vector<bool> get_cells() const {
        return _cells;
    }

    /* \brief Get the vector of vectors _base2bool */
    std::vector<std::vector<size_t> > get_base2bool() const {
        return _base2bool;
    }

    /* \brief Get the vector of maps _bool2base */
    std::vector<std::map<size_t, size_t> > get_bool2base() const {
        return _bool2base;
    }

    /* \brief Method computing the boundary of a given cell */
    CChain boundary_cell(size_t index_base, int dim) const {
        // Ensure dim is within valid range
        if (dim < 0 || dim >= _base2bool.size()) {
            throw std::out_of_range("Invalid dimension: " + std::to_string(dim));
        }

        // Check if index_base is within bounds of _base2bool[dim]
        if (index_base < 0 || index_base >= _base2bool[dim].size()) {
            throw std::out_of_range("index_base " + std::to_string(index_base) + " not found in _base2bool[" + std::to_string(dim) + "]");
        }

        size_t nb_lignes = (dim == 0) ? 0 : nb_cells(dim - 1);
        CChain boundary(nb_lignes);

        size_t index_bool = _base2bool[dim][index_base];
        std::vector<size_t> c = ind2khal(index_bool);

        CoefficientType sign = 1;
        for (size_t i = 0; i < _dim; ++i) {
            if (c[i] % 2 == 1) {
                // Calculate the coefficient based on the number of odd entries in c from 0 to i-1
                sign *= -1;

                size_t cell1 = index_bool + _P[i];
                size_t cell2 = index_bool - _P[i];

                // Insert the coefficient in the corresponding row of the boundary matrix
                if (cell1 >= 0 && cell1 < _cells.size()) {
                    size_t index = _bool2base[dim - 1].at(cell1);
                    boundary.set_coef(index, sign);
                }
                if (cell2 >= 0 && cell2 < _cells.size()) {
                    size_t index = _bool2base[dim - 1].at(cell2);
                    boundary.set_coef(index, -sign);
                }
            }
        }

        return boundary;
    }

    /* \brief Check if a cell (given in Khalimsky coordinates) is valid */
    bool is_valid_cell(const std::vector<size_t>& cell) const ;

    /* \brief Check if a cell (given by its boolean index) is valid */
    bool is_valid_cell(size_t id_cell) const ;

    /* \brief Compute vertices of a cell given by its Khalimsky coordinates */
    vector<size_t> khal_to_verts(vector<size_t> c) const
    {
        vector<vector<size_t> > vertices, vertices_tmp ;
        // Vertices are obtained by cartesian products
        for (size_t i=0; i<_dim; ++i)
        {
            if ((c[i]%2) == 1)
            {
                if (vertices.size()==0)
                {
                    vertices.push_back(vector<size_t>(1,c[i]-1)) ;
                    vertices.push_back(vector<size_t>(1,c[i]+1)) ;
                }
                else
                {
                    vertices_tmp.clear() ;
                    for (vector<size_t> vert : vertices)
                    {
                        vector<size_t> tmp(vert) ;
                        tmp.push_back(c[i]-1) ;
                        vertices_tmp.push_back(tmp) ;
                    }
                    for (vector<size_t> vert : vertices)
                    {
                        vector<size_t> tmp(vert) ;
                        tmp.push_back(c[i]+1) ;
                        vertices_tmp.push_back(tmp) ;
                    }
                    vertices = vertices_tmp ;
                }
            }
            else
            {
                if (vertices.size() == 0)
                    vertices.push_back(vector<size_t>(1,c[i])) ;
                else
                {
                    vertices_tmp.clear() ;
                    for (vector<size_t> vert : vertices)
                    {
                        vector<size_t> tmp(vert) ;
                        tmp.push_back(c[i]) ;
                        vertices_tmp.push_back(tmp) ;
                    }
                    vertices = vertices_tmp ;
                }
            }
        }
        vector<size_t> vertices_id ;
        for (vector<size_t> vert : vertices)
        {
            vertices_id.push_back(_bool2base[0].at(khal2ind(vert))) ;
        }
        return vertices_id ;
    }

    /// Member data
protected:
    /* \brief Dimension of the complex */
    int _dim;
    /* \brief Size of the Khalimsky bounding box \f$(N_1,\ldots,N_{\mathrm{\_dim}})\f$*/
    std::vector<size_t> _size_bb;
    /* \brief Vector of coefficients used for vectorization of the boolean representation */
    std::vector<size_t> _P;
    /* \brief Vectorized boolean representation of the complex (true if a cell is present, false otherwise) */
    std::vector<bool> _cells;
    /* \brief Maps from base indices to boolean indices (ie. indices in _cell) in each dimension */
    std::vector<std::vector<size_t>> _base2bool;
    /* \brief Maps from boolean indices (ie. indices in _cells) to base indices in each dimension */
    std::vector<std::map<size_t, size_t>> _bool2base;
    /* \brief Vector of boundary matrices in each dimension */
    std::vector<CMatrix>  _d;
private:
    std::vector<bool> _visited_cells; // Internal flag
    /* \brief Static counter for objects ids.
     * Initialized to 0.
     */
    static size_t _id_generator ;
    /* \brief Unique object id (for comparison of constant references to the complex). */
    const size_t _complex_id ;

    /// Protected methods
protected:
    /* Initialize _cells, _base2bool and _bool2base */
    void initialize_cells(const Cub_object_io& cub,typeComplexCube type);

    /** \brief Computes Khalimsky coordinates from boolean index */
    std::vector<size_t> ind2khal(size_t index) const {
        if (index > _P[_dim])
            throw std::invalid_argument("Index exceeds the size of boolean vector");
        std::vector<size_t> khal(_dim);
        for (size_t k = 0; k < _dim; ++k) {
            khal[k] = index % _size_bb[k];
            index /= _size_bb[k];
        }
        return khal;
    }

    /* \brief Computes boolean index from Khalimsky coordinates */
    size_t khal2ind(const std::vector<size_t>& base_indices) const {
        if (base_indices.size() != _dim) {
            throw std::invalid_argument("Dimension of base_indices does not match _dim");
        }

        size_t cell_index = 0;
        for (size_t i = 0; i < _dim; ++i) {
            cell_index += base_indices[i] * _P[i];
        }

        // verify if cell_index is within bounds of _cells
        if ((cell_index >= 0) && (cell_index < _P.at(_dim))) {
            return cell_index;
        } else {
            std::cerr << "Invalid cell index in _cells: " << cell_index << std::endl;
            throw std::out_of_range("Cell index out of bounds in _cells");
        }
    }

    /* \brief Computes voxel coordinates (in a binary object) from an index in a boolean vector
     *
     * This function is used ONLY for DUAL construction from a binary object (binary image in 2D, binary volume in 3D...)
     * Hence we consider a set of cells of maximal dimension vectorized to a boolean vector.
     */
    std::vector<size_t> ind2vox(size_t index, vector<size_t> B, size_t max_size) const
    {
        if (index > max_size)
            throw std::invalid_argument("ind2vox : index exceeding size of boolean vector");
        std::vector<size_t> coords(_dim);
        for (size_t k = 0; k < _dim; ++k) {
            coords[k] = index % B[k];
            index /= B[k];
        }
        return coords;
    }

    /* \brief Computes an index in a boolean vector from voxel coordinates (in a binary object)
     *
     * This function is used ONLY for DUAL construction from a binary object (binary image in 2D, binary volume in 3D...)
     * Hence we consider a set of cells of maximal dimension vectorized to a boolean vector.
     */
    size_t vox2ind(const std::vector<size_t>& base_indices, vector<size_t> B, size_t max_size) const {
        if (base_indices.size() != _dim) {
            throw std::invalid_argument("Dimension of base_indices does not match _dim");
        }

        size_t cell_index = 0;
        for (size_t i = 0; i < _dim; ++i) {
            cell_index += base_indices[i] * B[i];
        }

        // Verfiy if cell_index is within bounds of _cells
        if (cell_index >= 0 && cell_index < max_size)
            return cell_index;
        else {
            std::cerr << "Invalid cell index in _cells: " << cell_index << std::endl;
            throw std::out_of_range("Invalid cell index in _cells");
        }
    }


    /* \brief Computes the boundary matrix of dimension q */
    void  calculate_d(int q) ;

    /* \brief Insert a cell into the complex (and its faces if necessary) */
    void insert_cell(size_t cell);

    /* \brief Calculate the dimension of a cell (given in Khalimsky coordinates) */
    int calculate_dimension(const std::vector<size_t>& cell) const;
    /* \brief Calculate the dimension of a cell (given by its boolean index) */
    int calculate_dimension(size_t cell_index) const
    {
        return calculate_dimension(ind2khal(cell_index)) ;
    }

    /* \brief Compute (the boolean indices of) cells belonging to the boundary of `cell` (given by its boolean index) */
    std::vector<size_t> calculate_boundaries(size_t cell) const;

};

// Initialization of static VTK_cubtypes
template <typename CoefficientType> const
std::vector<int> Cubical_chain_complex<CoefficientType>::VTK_cubtypes({1, 3, 8, 11}) ;

// Initialization of _id_generator
template <typename CoefficientType>
size_t Cubical_chain_complex<CoefficientType>::_id_generator(0) ;

// Constructor implementation
template<typename CoefficientType>
Cubical_chain_complex<CoefficientType>::Cubical_chain_complex(const Cub_object_io& cub,typeComplexCube type) : _dim(cub.dim), _size_bb(_dim+1), _P(_dim+1,1), _base2bool(_dim+1), _bool2base(_dim+1), _complex_id(_id_generator++)

{
    // Initialize _size_bb and _P
    if (type==PRIMAL)
        _size_bb = cub.N;
    else
    {
        for (int q=0; q<_dim; ++q)
            _size_bb.at(q) = 2*cub.N.at(q)+1 ;
    }

    for (size_t i = 1; i <= _dim; ++i) {
        _P[i] = _P[i - 1] * _size_bb[i - 1];
    }

    _cells.resize(_P[_dim], false) ;

    // Initialize _visited_cells
    _visited_cells.resize(_P[_dim], false) ;

    // Initialize _cells, _base2bool, and _bool2base
    initialize_cells(cub,type);

    // Initialize _d
    _d.resize(_dim+1);
    for (int q = 0; q <= _dim; ++q) {
        calculate_d(q);
    }
}

// initialize_cells implementation
template<typename CoefficientType>
void Cubical_chain_complex<CoefficientType>::initialize_cells(const Cub_object_io& cub, typeComplexCube type)
{
    if (type == PRIMAL)
    {
        for (size_t i=0; i<cub.cubs.size(); ++i)
        {
            const size_t id(khal2ind(cub.cubs.at(i))) ;
            insert_cell(id);
        }
    }
    else if (type == DUAL)
    {
        size_t max_size(1) ;
        for (int q=0; q<_dim; ++q)
            max_size *= cub.N.at(q) ;

        //We iterate over all the voxels via indices
        for (size_t i=0; i<cub.cubs.size(); ++i)
        {
            vector<size_t> coords(cub.cubs.at(i)) ;
            // Calculate the coordinates of the voxel in the dual complex
            for (size_t i=0; i<_dim; ++i)
                coords.at(i)*=2 ;
            const size_t cell_index(khal2ind(coords)) ;

            _cells.at(cell_index) = true ;
            // Add the cell in _base2bool and _bool2base
            const int dim(calculate_dimension(coords)) ;
            const size_t n(_base2bool.at(dim).size()) ;
            _base2bool.at(dim).push_back(cell_index) ;
            _bool2base.at(dim)[cell_index] = n ;
        }

        // for the cells with dim>0
        for (int q = 1; q <= _dim; ++q) {

            for (size_t i = 0; i < _P[_dim]; ++i) {

                if (calculate_dimension(ind2khal(i)) == q) {
                    std::vector<size_t> boundaries = calculate_boundaries(i);
                    bool all_boundaries_present = true;
                    for (const auto& boundary_cell : boundaries) {
                        if (!_cells[boundary_cell]) {
                            all_boundaries_present = false;
                            break;
                        }
                    }
                    if (all_boundaries_present) {
                        _cells.at(i)=  true ;
                        // Add the cell in _base2bool and _bool2base
                        const int dim(calculate_dimension(i)) ;
                        const size_t n(_base2bool.at(dim).size()) ;
                        _base2bool.at(dim).push_back(i) ;
                        _bool2base.at(dim)[i] = n ;
                    }
                }
            }
        }
    }
    else
    {
        throw std::invalid_argument("Invalid typeComplexCube");
    }
}


// is_valid_cell implementation
template<typename CoefficientType>
bool Cubical_chain_complex<CoefficientType>::is_valid_cell(const std::vector<size_t>& cell) const {
    for (size_t i=0; i<_dim; ++i) {
        if (cell[i] < 0 || cell[i] >= _size_bb[i]) {
            return false;
        }
    }
    return true;
}

template<typename CoefficientType>
bool Cubical_chain_complex<CoefficientType>::is_valid_cell(size_t id_cell) const
{
    if ((id_cell < 0) || (id_cell > _P[_dim]))
        return false;
    else
        return true;
}


// insert_cell implementation
template<typename CoefficientType>
void Cubical_chain_complex<CoefficientType>::insert_cell(size_t cell) {
    if (!is_valid_cell(cell))
        throw std::out_of_range("insert_cell: trying to insert cell with invalid index");

    // verify if the cell has already been visited
    if (_visited_cells.at(cell)) {
        return; // cell has already been visited
    }

    std::vector<size_t> cell_coords(ind2khal(cell));
    int dim = calculate_dimension(cell_coords);

    _cells[cell] = true;
    size_t cell_base_index = _base2bool[dim].size();
    _base2bool[dim].push_back(cell);
    _bool2base[dim][cell] = cell_base_index;

    std::vector<size_t> boundaries(calculate_boundaries(cell));
    for (const auto& boundary : boundaries) {
        if (!_cells[boundary]) {
            insert_cell(boundary);
        }
    }

    // mark cell as visited
    _visited_cells.at(cell) = true ;
}


// calculate_d implementation
template<typename CoefficientType>
void Cubical_chain_complex<CoefficientType>::calculate_d(int dim)  {
    size_t nb_lignes = (dim == 0) ? 0 : nb_cells(dim - 1);

    _d[dim] = CMatrix(nb_lignes, nb_cells(dim));

    // Iterate through the cells of dimension dim
    for (size_t i = 0; i < nb_cells(dim); ++i) {
        // Boundary of the i-th cell of dimension dim
        CChain boundary = boundary_cell(i, dim);

        // Insert the chain into the corresponding column of the boundary matrix
        OSM::set_column(_d[dim], i, boundary);
    }
}

// calculate_dimension implementation
template<typename CoefficientType>
int Cubical_chain_complex<CoefficientType>::calculate_dimension(const std::vector<size_t>& cell) const {
    int dimension = 0;
    for (size_t index : cell) {
        if (index % 2 == 1) { // Un index impair indique une dimension plus élevée
            dimension++;
        }
    }
    return dimension;
}

// calculate_boundaries implementation
template<typename CoefficientType>
std::vector<size_t> Cubical_chain_complex<CoefficientType>::calculate_boundaries(size_t idcell) const {
    std::vector<size_t> boundaries;
    std::vector<size_t> c = ind2khal(idcell);

    for (size_t i = 0; i < _dim; ++i) {
        if (c[i] % 2 == 1)
        {
            // Calculate the coefficient based on the number of odd entries in c from 0 to i-1
            size_t cell1 = idcell + _P[i];
            if (is_valid_cell(cell1))
                boundaries.push_back(cell1) ;

            size_t cell2 = idcell - _P[i];
            if (is_valid_cell(cell2))
                boundaries.push_back(cell2) ;
        }
    }
    return boundaries;
}

/* \brief Export chain of dimension q to vtk.
 *           -> Only cells of the chain are exported
 *           -> If a cell Id is provided scalars are exported (0 for the given cellId / 2 for other cells)
 */

template <typename CoefficientType>
void Cubical_chain_complex<CoefficientType>::chain_complex_chain_to_vtk(const Cubical_chain_complex<CoefficientType> &K, const std::string &filename, const OSM::Sparse_chain<CoefficientType, OSM::COLUMN>& chain, int q, size_t cellId)
{
    bool with_scalars = (cellId != -1) ;

    // Load out file...
    std::ofstream out ( filename, std::ios::out | std::ios::trunc);

    if ( not out . good () ) {
        std::cerr << "SimpComplex_chain_to_vtk. Fatal Error:\n  " << filename << " not found.\n";
        throw std::runtime_error("File Parsing Error: File not found");
    }

    // Header
    out << "# vtk DataFile Version 2.0" << endl ;
    out << "generators" << endl ;
    out << "ASCII" << endl ;
    out << "DATASET  UNSTRUCTURED_GRID" << endl ;

    // Points
    size_t nnodes = K.nb_cells(0) ;
    out << "POINTS " << nnodes << " double" << endl ;
    for (size_t n = 0; n < nnodes; ++n)
    {
        vector<double> p(K.get_vertex_coords(n)) ;
        for (double x : p)
            out << x << " " ;
        for (size_t i = p.size(); i<3; ++i) // points must be 3D -> add zeros
            out << "0 " ;
        out << endl ;
    }

    size_t ncells_tot = 0, size_cells_tot = 0 ;
    std::vector<int> types ;
    std::vector<int> scalars ;
    std::vector<size_t> ids ;

    // output only cells of the chain (dimension q)
    {
        // 1 - Compute the number of cells / size of encoding
        {
            const size_t size_cell = 1<<q ;
            for (size_t id =0; id < K.nb_cells(q); ++id)
            {
                if (!chain.is_null(id))
                {
                    ++ncells_tot;
                    size_cells_tot += (size_cell+1) ;
                }
            }
        }
        // 2 - Output cells
        out << "CELLS " << ncells_tot << " " << size_cells_tot << endl ;

        {
            const size_t size_cell = 1<<q ;
            for (size_t id =0; id < K.nb_cells(q); ++id)
            {
                if (!chain.is_null(id))
                {
                    vector<size_t> khal(K.ind2khal(K._base2bool.at(q).at(id))) ;
                    vector<size_t> verts(K.khal_to_verts(K.ind2khal(K._base2bool.at(q).at(id)))) ;
                    out << size_cell << " " ;
                    for (size_t i : verts)
                        out << i << " " ;
                    out << endl ;
                    types.push_back(Cubical_chain_complex<CoefficientType>::VTK_cubtypes.at(q)) ;

                    if (with_scalars)
                    {
                        ids.push_back(id) ;
                        if (id != cellId)
                            scalars.push_back(2) ;
                        else
                            scalars.push_back(0) ;
                    }
                }
            }
        }
        // CELL_TYPES
        out << "CELL_TYPES " << ncells_tot << endl ;
        for (int t : types)
            out << t << " " ;
        out << endl ;
    }

    if (with_scalars)
    {
        // CELL_TYPES
        out << "CELL_DATA " << ncells_tot << endl ;
        out << "SCALARS Label " << "int" << " 1" << endl ;
        out << "LOOKUP_TABLE default" << endl ;
        for (int s : scalars)
            out << s << " " ;
        out << endl ;
        // CELL_IDs
        out << "SCALARS CellId " << "int" << " 1" << endl ;
        out << "LOOKUP_TABLE default" << endl ;
        for (size_t i : ids)
            out << i << " " ;
        out << endl ;
    }
    out.close() ;
}

} /* end namespace HDVF */
} /* end namespace CGAL */

#endif // CGAL_HDVF_CUBICAL_CHAIN_COMPLEX_H
