#ifndef SIMPCOMPLEX_HPP
#define SIMPCOMPLEX_HPP

#include <vector>
#include <map>
//#include "hdvf_common.hpp"
#include "tools_io.hpp"
#include "Simplex.hpp"
#include "CGAL/OSM/OSM.hpp"

// Forward declaration of SimpComplexTools
template<typename CoefficientType, typename VertexIdType> class Duality_SimpComplexTools ;

/**
 * \class Abstract_simplicial_chain_complex
 * \brief Implementation of abstract simplicial complexes.
 *
 * The Abstract_simplicial_chain_complex class contains constructors and functions associated to abstract simplicial complexes for HDVFs. It implements the Complex concept.
 *
 * \tparam CoefficientType The coefficient types used for boundary matrix computation (default is int)
 * \tparam VertexIdType Type used to store vertices indices (default is int)
 *
 * \author Alexandra Bac
 */
template<typename CoefficientType, typename VertexIdType = int>
class Abstract_simplicial_chain_complex {
public:
    /**
     * \brief Default constructor.
     *
     * Builds an empty abstract simplicial complex.
     *
     * \see \link Abstract_simplicial_chain_complex \endlink
     *
     * \author Alexandra Bac
     * */
    Abstract_simplicial_chain_complex() ;
    
    /**
     * \brief Constructor from a Mesh_object.
     *
     * Builds the abstract simplicial complex associated to a triangular mesh (ie. performs the down closure of cells and set the boundary matrices in any dimension).
     *
     * \param[in] mesh A Mesh_object containing a triangular mesh.
     *
     * \see \link Abstract_simplicial_chain_complex \endlink
     *
     * \author Alexandra Bac
     * */
    Abstract_simplicial_chain_complex(const Mesh_object& mesh);
    
    /** \brief Type of column-major chains */
    typedef OSM::Chain<CoefficientType, OSM::COLUMN> CChain;
    /** \brief Type of row-major chains */
    typedef OSM::Chain<CoefficientType, OSM::ROW> RChain ;
    /** \brief Type of column-major sparse matrices */
    typedef OSM::SparseMatrix<CoefficientType, OSM::COLUMN> CMatrix;
    
    
    /**
     * \brief Affectation operator for abstract simplicial complexes.
     *
     * Stores a copy of an abstract simplicial complex in *this.
     *
     * \param[in] complex The abstract simplicial complex which will be copied.
     *
     * \see \link Abstract_simplicial_chain_complex \endlink
     *
     * \author Alexandra Bac
     * */
    Abstract_simplicial_chain_complex& operator= (const Abstract_simplicial_chain_complex& complex)
    {
        _dim = complex._dim;
        _ind2simp = complex._ind2simp;
        _simp2ind = complex._simp2ind;
        _nb_cells = complex._nb_cells;
        _d = complex._d;
        return *this ;
    }
    
    /// Methods of the SimpComplex concept

    /**
     * \brief Method returning the boundary of the cell id_cell in dimension dim.
     *
     * Returns a copy of the column-major chain stored in the boundary matrix of dimension dim: boundary of the cell id_cell in dimension dim..
     *
     * \param[in] id_cell Index of the cell.
     * \param[in] dim Dimension of the cell.
     *
     * \return The column-major chain containing the boundary of the cell id_cell in dimension dim.
     *
     * \see \link Abstract_simplicial_chain_complex \endlink
     *
     * \author Alexandra Bac
     * */
     CChain d(int id_cell, int dim) const
    {
        if (dim > 0)
            return OSM::getColumn(_d[dim], id_cell);
        else
            return CChain(0) ;
    }
    
    /**
     * \brief Method returning the co-boundary of the cell id_cell in dimension dim.
     *
     * Returns a row-major chain containing the co-boundary of the cell id_cell in dimension dim (so actually a row of the boundary matrix).
     *
     * \warning As the boundary matrix is stored column-major, this entails crossing the full matrix to extract the row coefficients (O(number of non empty columns))
     *
     * \param[in] id_cell Index of the cell.
     * \param[in] dim Dimension of the cell.
     *
     * \return The row-major chain containing the co-boundary of the cell id_cell in dimension dim.
     *
     * \see \link Abstract_simplicial_chain_complex \endlink
     *
     * \author Alexandra Bac
     * */
    RChain cod(int id_cell, int dim) const
    {
        if (dim < _dim)
            return OSM::getRow(_d[dim+1], id_cell);
        else
            return RChain(0) ;
    }
    
    /**
     * \brief Method returning the dimension of the complex.
     *
     * Returns the dimension of the simplicial complex (ie. largest dimension of cells).
     *
     * \return The dimension of the complex..
     *
     * \see \link Abstract_simplicial_chain_complex \endlink
     *
     * \author Alexandra Bac
     * */
    int dim() const { return _dim ;}
    
    /**
     * \brief Method returning the number of cells in a given dimension.
     *
     * \param[in] dim Dimension along which the number of cells is returned.
     *
     * \return Number of cells in dimension dim.
     *
     * \see \link Abstract_simplicial_chain_complex \endlink
     *
     * \author Alexandra Bac
     * */
    int nb_cells(int dim) const
    {
        if ((dim >=0) && (dim<=_dim))
            return _nb_cells[dim] ;
        else
            return 0 ;
    }
    
    /**
     * \brief Method returning a constant reference to the vector of boundary matrices (along each dimension).
     *
     * Returns a constant reference to the vector of boundary matrices along each dimension. The q-th element of this vector is a column-major sparse matrix containing the boundaries of q-cells (ie. rows encode q-1 cells and columns q cells).
     *
     * \return Returns a constant reference to the vector of column-major boundary matrices along each dimension.
     *
     * \see \link Abstract_simplicial_chain_complex \endlink
     *
     * \author Alexandra Bac
     * */
    const vector<CMatrix> & get_bnd_matrices() const
    {
        return _d ;
    }
    
    /**
     * \brief Method returning a copy of the dim-th boundary matrix (ie. column-major matrix of d_dim).
     *
     * It is a column-major sparse matrix containing the boundaries of dim-cells (ie. rows encode dim-1 cells and columns dim cells).
     *
     * \param[in] dim Dimension of the boundary matrix (ie. columns will contain the boundary of dimension dim cells).
     *
     * \return A column-major sparse matrix containing the matrix of the boundary operator of dimension dim.
     *
     * \see \link Abstract_simplicial_chain_complex \endlink
     *
     * \author Alexandra Bac
     * */
    const CMatrix & get_bnd_matrix(int dim) const
    {
        return _d.at(dim) ;
    }

    /**
     * \brief Method returning dimension 0 simplicies indexes included in the cell with index id_cell of dimension dim.
     *
     * Returns the dimension 0 simplicies indexes included in the cell with index id_cell of dimension dim.
     *
     * \warning This does not come to return vertices indices, as dimension 0 simplices enumerate vertices in any order. For instance, if an abstract simplicial complex is build from 3 vertices {1,2,3} such that the enumeration of dimension 0 simplicies is:
     *  id0: 3, id1 : 2, id2: 1
     * then the bottom_faces of the 1-simplex {1,2} are two 0-simplices with id 2 and 1.
     *
     * \param[in] id_cell Index of the cell.
     * \param[in] dim Dimension of the cell.
     *
     * \return A vector of 0-simplices indices.
     *
     * \see \link Abstract_simplicial_chain_complex \endlink
     *
     * \author Alexandra Bac
     * */
    std::vector<int> bottom_faces(int id_cell, int dim) const
    {
        std::set<int> verts(_ind2simp.at(dim).at(id_cell).getVertices()) ;
        std::vector<int> res ;
        // For each vertex in verts, compute the corresponding dimension 0 cell
        for (int vert_id : verts)
        {
            const int i(_simp2ind.at(0).at(Simplex(std::set<int>({vert_id})))) ;
            res.push_back(i) ;
        }
        return res ;
    }
    
    /**
     * \brief Method printing informations of the complex.
     *
     * Displays the number of cells in each dimension and the boundary matrix in each dimension.
     *
     * \see \link Abstract_simplicial_chain_complex \endlink
     *
     * \author Alexandra Bac
     * */
    /** \brief   */
    void print_complex() const;

   
protected:
    /** \brief Dimension of the complex */
    int _dim;
    /** \brief Vector of simplices along each dimension: _ind2simp.at(q) contains the vector of all simplices of dimension q */
    std::vector<std::vector<Simplex<VertexIdType> > > _ind2simp ;
    /** \brief Vector of maps associating indices to simplices in each dimension: _simp2ind.at(q) maps q-simplices to their index in _ind2simp  */
    std::vector<std::map<Simplex<VertexIdType>, int> > _simp2ind ;
    /** \brief Vector of number of cells in each dimension */
    std::vector<int> _nb_cells ;
    /** \brief Vector of column-major boundary matrices: _d.at(q) is the matrix of the boundary operator in dimension q. */
    mutable std::vector<CMatrix>  _d;  // Boundary matrices

    // Protected methods
    
    /**
     * \brief Method filling the matrix _d.at(dim), ie. computing the boundaries of cells of dimension dim.
     *
     * Compute the boundary with the standard simplicial complex definition:
     * \f[\partial (\{v_0,\ldots,v_q\}) = \sum_{i=0}^q (-1)^q \{v_0,\ldots, \hat{v_i},\ldots,v_q\}\f]
     *
     * \param[in] dim Dimension considered for computation.
     *
     * \see \link Abstract_simplicial_chain_complex \endlink
     *
     * \author Alexandra Bac
     * */
    void  calculate_d(int dim) const;
    
    /**
     * \brief Method inserting a simplex (and its faces if necessary) into the abstract simplicial complex.
     *
     * The method (recursively) inserts the simplex into _ind2simp and _simp2ind and its faces. If the simplex is already presend, the method does nothing.
     *
     * \warning insert_simplex does not update the boundary matrix.
     *
     * \param[in] tau The simplex inserted.
     *
     * \see \link Abstract_simplicial_chain_complex \endlink
     *
     * \author Alexandra Bac
     * */
    void insert_simplex(const Simplex<VertexIdType>& tau);
};

//constructor
template<typename CoefficientType, typename VertexIdType>
Abstract_simplicial_chain_complex<CoefficientType, VertexIdType>::Abstract_simplicial_chain_complex(const Mesh_object& mesh) {
    // Initialize attributes
    
    _dim = abs(mesh.dim);
    
    // Initialize vectors of Simplices and cell counts
    _ind2simp.resize(_dim + 1);
    _simp2ind.resize(_dim + 1);
    _nb_cells.resize(_dim + 1, 0);

    // Iterate through the mesh cells and add them to the complex
    for (const auto& cell : mesh.cells) {
//        Simplex simplex(std::set<int>(cell.begin(), cell.end()));
//        insert_simplex(simplex);
        insert_simplex(cell) ;
    }

    _d.resize(_dim+1);
    for (int dim = 0; dim <= _dim; ++dim) {
        calculate_d(dim);
    }
}

// Recursive insertion method
template<typename CoefficientType, typename VertexIdType>
void Abstract_simplicial_chain_complex<CoefficientType, VertexIdType>::insert_simplex(const Simplex<VertexIdType>& tau) {
    int q = tau.dimension();
    if (q == -1) return;

    if (_simp2ind[q].find(tau) == _simp2ind[q].end()) {
        int i = _ind2simp[q].size();
        _ind2simp[q].push_back(tau);
        _simp2ind[q][tau] = i;
        _nb_cells[q]++;

        std::vector<Simplex<VertexIdType> > bord = tau.boundary();
        for (const auto& sigma : bord) {
            insert_simplex(sigma);
        }
    }
}

// calculate _d boundary matrix
template<typename CoefficientType, typename VertexIdType>
void Abstract_simplicial_chain_complex<CoefficientType,VertexIdType>::calculate_d(int dim) const {
    int nb_lignes = (dim == 0) ? 0 : _nb_cells[dim - 1];
    _d[dim] = CMatrix(nb_lignes, _nb_cells[dim]);

    // Iterate through the cells of dimension dim
    if (dim>0)
    {
        for (int i = 0; i < _nb_cells[dim]; ++i) {
            // Boundary of the i-th simplex of dimension dim
            const Simplex<VertexIdType>& s = _ind2simp[dim][i];
            std::vector<Simplex<VertexIdType> > bord = s.boundary();
            
            // Create a chain with the correct size
            CChain chain(nb_lignes);
            
            // For each element of the boundary
            for (int j = 0; j < bord.size(); ++j) {
                auto it = _simp2ind[dim - 1].find(bord[j]); // Find the index of Simplex j
                if (it != _simp2ind[dim - 1].end()) { // If Simplex j is found
                    int ind_j = it->second; // Retrieve the index of Simplex j
                    chain[ind_j] = (j % 2 == 0) ? 1 : -1;
                }
                else
                    throw "calculate_d boundary simplex not found!";
            }
            
            // Insert the chain into the corresponding column of the delta matrix
            OSM::setColumn(_d[dim], i, chain);
        }
    }
}

// Method to display the complex's information
template<typename CoefficientType, typename VertexIdType>
void Abstract_simplicial_chain_complex<CoefficientType, VertexIdType>::print_complex() const {
    std::cout << "Complex dimension: " << _dim << std::endl;

    // Total number of cells
    int nb_total_cells = 0;
    for (int i = 0; i <= _dim; ++i) {
        nb_total_cells += _nb_cells[i];
    }
    std::cout << "Total number of cells: " << nb_total_cells << std::endl;
    
    // Cells per dimension
    for (int q = 0; q <= _dim; ++q) {
        std::cout << "--- dimension " << q << std::endl;
        cout << nb_cells(q) << " cells" << endl ;
//        for (int j = 0; j < _nb_cells.at(q); ++j) {
//            Simplex s(_ind2Simp.at(q).at(j));
//            std::cout << j << " -> " << s << " -> " << _Simp2ind.at(q).at(s) << std::endl;
//        }
    }
    
    // Boundary matrices
    std::cout << "---------------------------" << std::endl << "Boundary matrices" << std::endl;
    for (int q = 1; q <= _dim; ++q)
        std::cout << "_d[" << q << "] : " << _d[q].dimensions().first << "x" << _d[q].dimensions().second << std::endl <<  _d[q] << std::endl;
}

/**
 * \class Simplicial_chain_complex
 * \brief Implementation of simplicial complexes.
 *
 * The Simplicial_chain_complex class contains constructors and functions associated to simplicial complexes for HDVFs. It implements the GeometricComplex concept. Hence, vtk output is available.
 *
 * \tparam _CoefficientType The coefficient types used for boundary matrix computation (default is int)
 *
 * \author Alexandra Bac
 */

template<typename CoefficientType, typename VertexIdType = int>
class Simplicial_chain_complex : public Abstract_simplicial_chain_complex<CoefficientType, VertexIdType> {
protected:
    std::vector<std::vector<double> > _coords ;
    
private:
    /// Associated to any dimension the type number of associated VTK cells
    /// {1, 3, 5, 10}
    static const std::vector<int> VTK_simptypes ;
    
public:
    
    /**
     * \brief Default constructor: builds an empty  simplicial complex
     *
     * \see \link Simplicial_chain_complex \endlink
     *
     * \author Alexandra Bac
     */
    Simplicial_chain_complex() {} ;
    // Constructor from a simplicial complex stored in a .simp file and read into a MeshObject
    /**
     * \brief Constructor from a Mesh_object.
     *
     * Builds a simplicial complex from a Mesh_object and a vector of vertices coordinates (
     *
     * \see \link Simplicial_chain_complex \endlink
     *
     * \author Alexandra Bac
     */
    Simplicial_chain_complex(const Mesh_object& mesh, std::vector<std::vector<double> > coords) : Abstract_simplicial_chain_complex<CoefficientType>(mesh), _coords(coords) {} ;
    
    /// Operator=
    Simplicial_chain_complex& operator= (const Simplicial_chain_complex& complex)
    {
        this->Abstract_simplicial_chain_complex<CoefficientType>::operator=(complex) ;
        _coords = complex._coords ;
        return *this ;
    }

    // Friend class : SimpObjectTools
    friend Duality_SimpComplexTools<CoefficientType, VertexIdType> ;
    
    /** \brief Get the vector of vertices coordinates  */
    const std::vector<std::vector<double> >& get_nodes_coords() const
    {
        return _coords ;
    }
    
    /** \brief Get the coordinates of the ith dimension 0 simplex  */
    std::vector<double> get_vertex_coords (int i) const
    {
        const Simplex simpl(this->_ind2simp.at(0).at(i)) ;
        const std::set<int> verts(simpl.getVertices()) ;
        const int id(*(verts.cbegin())) ;
        return _coords.at(id);
//        return _coords.at(i);
    }
    
    // VTK export

    /**
     * \brief Method exporting a simplicial complex (plus, optionally, labels) to a VTK file.
     *
     * The method generates legacy text VTK files. Labels are exported as such in a VTK property, together with CellID property, containing the index of each cell.
     *
     * \param[in] K Simplicial complex exported.
     * \param[in] filename Output file root (output filenames will be built from this root).
     * \param[in] labels Pointer to a vector of labels in each dimension. (*labels).at(q) is the set of integer labels of cells of dimension q. If labels is NULL, only CellID property is exported.
     *
     * \see \link Abstract_simplicial_chain_complex \endlink
     *
     * \author Alexandra Bac
     * */
    static void Simplicial_chain_complex_to_vtk(const Simplicial_chain_complex &K, const std::string &filename, const std::vector<std::vector<int> > *labels=NULL) ;

    /**
     * \brief Method exporting a chain over a simplicial complex to a VTK file.
     *
     * The method generates legacy text VTK files. All the cells of the chain with non zero coefficient are exported. If a cellId is provided, labels are exported in a VTK property (2 for all cells, 0 for cell of index cellId).  The index of each cell is exported in a CellID property.
     *
     * \param[in] K Simplicial complex exported.
     * \param[in] filename Output file root (output filenames will be built from this root).
     * \param[in] chain Chain exported (all the cells with non-zero coefficients in the chain are exported to vtk).
     * \param[in] q Dimension of the cells of the chain.
     * \param[in] cellId If a cellID is provided, labels are exported to distinguish cells of the chain (label 2) from cellId cell (label 0).
     *
     * \see \link Abstract_simplicial_chain_complex \endlink
     *
     * \author Alexandra Bac
     * */
    static void Simplicial_chain_complex_chain_to_vtk(const Simplicial_chain_complex &K, const std::string &filename, const OSM::Chain<CoefficientType, OSM::COLUMN>& chain, int q, int cellId = -1) ;
};

// Initialization of static VTK_simptypes
template <typename CoefficientType, typename VertexIdType> const
std::vector<int> Simplicial_chain_complex<CoefficientType, VertexIdType>::VTK_simptypes({1, 3, 5, 10});

template <typename CoefficientType, typename VertexIdType>
void Simplicial_chain_complex<CoefficientType, VertexIdType>::Simplicial_chain_complex_to_vtk(const Simplicial_chain_complex &K, const std::string &filename, const std::vector<std::vector<int> > *labels)
{
    typedef Simplicial_chain_complex<CoefficientType, VertexIdType> ComplexType;
    if (K._coords.size() != K.nb_cells(0))
    {
        std::cerr << "SimpComplex_to_vtk. Error, wrong number of points provided.\n";
        throw std::runtime_error("Geometry of points invalid.");
    }
    
    bool with_scalars = (labels != NULL) ;
    
    // Load out file...
    std::ofstream out ( filename, std::ios::out | std::ios::trunc);

    if ( not out . good () ) {
        std::cerr << "SimpComplex_to_vtk. Fatal Error:\n  " << filename << " not found.\n";
        throw std::runtime_error("File Parsing Error: File not found");
    }

    // Header
    out << "# vtk DataFile Version 2.0" << endl ;
    out << "generators" << endl ;
    out << "ASCII" << endl ;
    out << "DATASET  UNSTRUCTURED_GRID" << endl ;

    // Points
    size_t nnodes = K._coords.size() ;
    out << "POINTS " << nnodes << " double" << endl ;
    const std::vector<std::vector<double> >& coords(K.get_nodes_coords()) ;
    for (int n = 0; n < nnodes; ++n)
    {
        vector<double> p(coords.at(n)) ;
        for (double x : p)
            out << x << " " ;
        for (int i = p.size(); i<3; ++i) // points must be 3D -> add zeros
            out << "0 " ;
        out << endl ;
    }

    int ncells_tot = 0, size_cells_tot = 0 ;
    std::vector<int> types ;
    std::vector<int> scalars ;
    std::vector<int> ids ;
    // all cells must be printed
    {
        // Cells up to dimension 3
        // Size : size of a cell of dimension q : q+1
        for (int q=0; q<=K._dim; ++q)
        {
            ncells_tot += K.nb_cells(q) ;
            const int size_cell = q+1 ;
            size_cells_tot += (size_cell+1)*K.nb_cells(q) ;
        }
        out << "CELLS " << ncells_tot << " " << size_cells_tot << endl ;
        // Output cells by increasing dimension
        
//        // Vertices
//        for (int i = 0; i<K.nb_cells(0); ++i)
//        {
//            out << "1 " << i << endl ;
//            types.push_back(VTK_simptypes.at(0)) ;
//            if (with_scalars)
//            {
//                scalars.push_back((*labels).at(0).at(i)) ;
//                ids.push_back(i) ;
//            }
//        }
        // Cells of higher dimension
        for (int q=0; q<=K._dim; ++q)
        {
            const int size_cell = q+1 ;
            for (int id =0; id < K.nb_cells(q); ++id)
            {
                Simplex verts(K._ind2simp.at(q).at(id)) ;
                out << size_cell << " " ;
                for (typename Simplex<VertexIdType>::const_iterator it = verts.cbegin(); it != verts.cend(); ++it)
                    out << *it << " " ;
                out << endl ;
                types.push_back(ComplexType::VTK_simptypes.at(q)) ;
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
        out << "SCALARS Label " << "int" << " 1" << endl ;
        out << "LOOKUP_TABLE default" << endl ;
        for (int s : scalars)
            out << s << " " ;
        out << endl ;
        // CELL_IDs
        out << "SCALARS CellId " << "int" << " 1" << endl ;
        out << "LOOKUP_TABLE default" << endl ;
        for (int i : ids)
            out << i << " " ;
        out << endl ;
    }
    out.close() ;
}

/** \brief Export chain of dimension q to vtk.
*           -> Only cells of the chain are exported
*           -> If a cell Id is provided scalars are exported (0 for the given cellId / 2 for other cells)
*/

template <typename CoefficientType, typename VertexIdType>
void Simplicial_chain_complex<CoefficientType, VertexIdType>::Simplicial_chain_complex_chain_to_vtk(const Simplicial_chain_complex &K, const std::string &filename, const OSM::Chain<CoefficientType, OSM::COLUMN>& chain, int q, int cellId)
{
    typedef Simplicial_chain_complex<CoefficientType, VertexIdType> ComplexType ;
    if (K._coords.size() != K.nb_cells(0))
    {
        std::cerr << "SimpComplex_chain_to_vtk. Error, wrong number of points provided.\n";
        throw std::runtime_error("Geometry of points invalid.");
    }
    
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
    int nnodes = K._coords.size() ;
    out << "POINTS " << nnodes << " double" << endl ;
    const std::vector<std::vector<double> >& coords(K.get_nodes_coords()) ;
    for (int n = 0; n < nnodes; ++n)
    {
        vector<double> p(coords.at(n)) ;
        for (double x : p)
            out << x << " " ;
        for (int i = p.size(); i<3; ++i) // points must be 3D -> add zeros
            out << "0 " ;
        out << endl ;
    }

    int ncells_tot = 0, size_cells_tot = 0 ;
    std::vector<int> types ;
    std::vector<int> scalars ;
    std::vector<int> ids ;
    
    // output only cells of the chain (dimension q)
    {
        // 1 - Compute the number of cells / size of encoding
        {
            const int size_cell = q+1 ;
            for (int id =0; id < K.nb_cells(q); ++id)
            {
                if (!chain.isNull(id))
                {
                    ++ncells_tot;
                    size_cells_tot += (size_cell+1) ;
                }
            }
        }
        // 2 - Output cells
        out << "CELLS " << ncells_tot << " " << size_cells_tot << endl ;
        
        {
            const int size_cell = q+1 ;
            for (int id =0; id < K.nb_cells(q); ++id)
            {
                if (!chain.isNull(id))
                {
                    Simplex verts(K._ind2simp.at(q).at(id)) ;
                    out << size_cell << " " ;
                    for (typename Simplex<VertexIdType>::const_iterator it = verts.cbegin(); it != verts.cend(); ++it)
                        out << *it << " " ;
                    out << endl ;
                    types.push_back(ComplexType::VTK_simptypes.at(q)) ;
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
        for (int i : ids)
            out << i << " " ;
        out << endl ;
    }
    out.close() ;
}

#endif // SIMPCOMPLEX_HPP
