#ifndef CUBCOMPLEX_HPP
#define CUBCOMPLEX_HPP

#include <vector>
#include <map>
#include <stdexcept>
#include <unordered_set>
//#include "hdvf_common.hpp"
#include "tools_io.hpp"
#include "CGAL/OSM/OSM.hpp"

// Forward declaration of CubComplexTools
template<typename T> class CubComplexTools ;

// TOOLS

ostream & operator<<(ostream & out, std::vector<int> c)
{
    for (int i:c)
        out << i << " " ;
    out << endl ;
    return out ;
}

/**
 * \class Cubical_complex
 * \brief Implementation of cubical complexes.
 *
 * The Cubical_complex class contains constructors and functions associated to cubical complexes for HDVFs. It implements the GeometricComplex concept. Hence, vtk output is available.
 *
 * \tparam _CoefficientType The coefficient types used for boundary matrix computation (default is int)
 *
 * \author Alexandra Bac
 */

template<typename CoefficientType> // !!!! Substitute! -> coefficienttype everywhere...
class Cubical_complex {
private:
    /// Associate to any dimension the corresponding VTK type
    /// {1, 3, 8, 11}
    static const std::vector<int> VTK_cubtypes ;
    
public:
    /// Type of construction applied by the constructor
    enum typeComplexCube {PRIMAL, DUAL};
    
    /// Constructors of CubObject
    Cubical_complex() {} ;
    Cubical_complex(const CubObject& cub,typeComplexCube type);
    
    /// Friend class : CubObjectTools
    friend CubComplexTools<CoefficientType> ;
    
    ///typedef related to boundary matrices
    typedef OSM::Chain<CoefficientType, OSM::COLUMN> CChain;
    typedef OSM::Chain<CoefficientType, OSM::ROW> RChain ;
    typedef OSM::SparseMatrix<CoefficientType, OSM::COLUMN> CMatrix;
    
    /// operator=
    Cubical_complex& operator=(const Cubical_complex& complex)
    {
        _dim = complex._dim;
        _size_bb = complex._size_bb;
        _P = complex._P;
        _cells = complex._cells;
        _base2bool = complex._base2bool;
        _bool2base = complex._bool2base;
        visited_cells = complex.visited_cells;
        return *this;
    }
    
    /// Methods of the CubComplex concept
    
    CChain d(int id_cell, int dim) const //border of the cell of index id_cell of dimension dim
    {
        if (dim > 0)
            return OSM::cgetColumn(_d[dim], id_cell);
        else
            return CChain(0) ;
    }
    
    RChain cod(int id_cell, int dim) const // cobord of the cell of index id_cell of dimension dim
    {
        if (dim < _dim)
            return OSM::getRow(_d[dim+1], id_cell);
        else
            return RChain(0) ;
    }
    
    const vector<CMatrix> & get_bnd_matrices() const
    {
        return _d ;
    }
    
    const CMatrix & get_bnd_matrix(int q) const
    {
        return _d.at(q) ;
    }
    
    int dim() const
    {
        return _dim ;
    }
    
    int nb_cells(int dim) const
    {
        if ((dim >=0) && (dim <= _dim))
            return _base2bool.at(dim).size() ;
        else
            return 0 ;
    }
    
    void print_complex() const {
        for (int q = 0; q <= _dim; ++q) {
            std::cout << "-------- dimension " << q << std::endl;
            std::cout << "cellules de dimension " << q << " : " << _base2bool.at(q).size() << std::endl;
            for (size_t id_base = 0; id_base < _base2bool.at(q).size(); ++id_base) {
                int id_bool = _base2bool.at(q).at(id_base);
                std::vector<int> khal = ind2khal(id_bool);
                std::cout << id_base << " -> " << id_bool << " -> " << _bool2base.at(q).at(id_bool) << " -> ";
                for (int k : khal) std::cout << k << " ";
                std::cout << std::endl;
            }
            
            if (_base2bool[q].size() > 0 && q <= _dim) {
                std::cout << "matrice de bord de dimension " << q << std::endl;
                std::cout << _d[q] << std::endl;
            }
        }
    }
    
    std::vector<int> bottom_faces(int id_cell, int dim) const // Vertices included in the cell with index id_cell of dimension dim
    {
        // Khalimsky coordinates of the cell
        const vector<int> coords(ind2khal(_base2bool.at(dim).at(id_cell))) ;
        return khal_to_verts(coords) ;
    }
    
    /// END Methods of the SimpComplex concept
    
    // Methods to access data
    std::vector<int> get_size_bb() const {
        return _size_bb;
    }
    
    std::vector<int> get_P() const {
        return _P;
    }
    
    std::vector<bool> get_cells() const {
        return _cells;
    }
    
    std::vector<std::vector<int> > get_base2bool() const { 
        return _base2bool;
    }
    
    std::vector<std::map<int, int> > get_bool2base() const { 
        return _bool2base;
    }
    
    // Method to get number of cells for a given dimension
    int get_N(int dim) const {
        if (dim >= 0 && dim <= _dim) {
            return _base2bool[dim].size();
        } else {
            throw std::out_of_range("Invalid dimension");
        }
    }
    
    std::vector<double> get_vertex_coords(int i) const
    {
        const vector<int> coords(ind2khal(_base2bool.at(0).at(i))) ;
        vector<double> res ;
        for (int c : coords)
            res.push_back(c/2.) ;
        for (int i = coords.size(); i<3; ++i) // points must be 3D
            res.push_back(0) ;
        return res ;
    }
    
    CChain boundary_cell(int index_base, int dim) const {
        // Ensure dim is within valid range
        if (dim < 0 || dim >= _base2bool.size()) {
            throw std::out_of_range("Invalid dimension: " + std::to_string(dim));
        }
        
        // Check if index_base is within bounds of _base2bool[dim]
        if (index_base < 0 || index_base >= _base2bool[dim].size()) {
            throw std::out_of_range("index_base " + std::to_string(index_base) + " not found in _base2bool[" + std::to_string(dim) + "]");
        }
        
        int nb_lignes = (dim == 0) ? 0 : get_N(dim - 1);
        CChain boundary(nb_lignes);
        
        int index_bool = _base2bool[dim][index_base];
        std::vector<int> c = ind2khal(index_bool);
        
        CoefficientType sign = 1;
        for (int i = 0; i < _dim; ++i) {
            if (c[i] % 2 == 1) {
                // Calculate the coefficient based on the number of odd entries in c from 0 to i-1
                sign *= -1;
                
                int cell1 = index_bool + _P[i];
                int cell2 = index_bool - _P[i];
                
                // Insert the coefficient in the corresponding row of the boundary matrix
                if (cell1 >= 0 && cell1 < _cells.size()) {
                    int index = _bool2base[dim - 1].at(cell1);
                    boundary[index] = sign;
                }
                if (cell2 >= 0 && cell2 < _cells.size()) {
                    int index = _bool2base[dim - 1].at(cell2);
                    boundary[index] = -sign;
                }
            }
        }
        
        return boundary;
    }

    /// Check if a cell (given in Khalimsky coordinates) is valid
    bool is_valid_cell(const std::vector<int>& cell) const ;
    
    /// Check if a cell (given by its boolean index) is valid
    bool is_valid_cell(int id_cell) const ;
    
    vector<int> khal_to_verts(vector<int> c) const
    {
        vector<vector<int> > vertices, vertices_tmp ;
        // Vertices are obtained by cartesian products
        for (int i=0; i<_dim; ++i)
        {
            if ((c[i]%2) == 1)
            {
                if (vertices.size()==0)
                {
                    vertices.push_back(vector<int>(1,c[i]-1)) ;
                    vertices.push_back(vector<int>(1,c[i]+1)) ;
                }
                else
                {
                    vertices_tmp.clear() ;
                    for (vector<int> vert : vertices)
                    {
                        vector<int> tmp(vert) ;
                        tmp.push_back(c[i]-1) ;
                        vertices_tmp.push_back(tmp) ;
                    }
                    for (vector<int> vert : vertices)
                    {
                        vector<int> tmp(vert) ;
                        tmp.push_back(c[i]+1) ;
                        vertices_tmp.push_back(tmp) ;
                    }
                    vertices = vertices_tmp ;
                }
            }
            else
            {
                if (vertices.size() == 0)
                    vertices.push_back(vector<int>(1,c[i])) ;
                else
                {
                    vertices_tmp.clear() ;
                    for (vector<int> vert : vertices)
                    {
                        vector<int> tmp(vert) ;
                        tmp.push_back(c[i]) ;
                        vertices_tmp.push_back(tmp) ;
                    }
                    vertices = vertices_tmp ;
                }
            }
        }
        vector<int> vertices_id ;
        for (vector<int> vert : vertices)
        {
            vertices_id.push_back(_bool2base[0].at(khal2ind(vert))) ;
        }
        return vertices_id ;
    }
    
//    template <typename TCoefs>
//    friend void cubical_complex_write_to_vtk(const Cubical_complex<TCoefs> &K, const std::string &filename, const std::vector<std::vector<int> > *labels, TDataLabels tdata, const std::string scalar_type) ;
    
//    template <typename TCoefs>
    static void cubical_complex_to_vtk(const Cubical_complex<CoefficientType> &K, const std::string &filename, const std::vector<std::vector<int> > *labels=NULL) ;
    
//    template <typename TCoefs>
    static void cubical_complex_chain_to_vtk(const Cubical_complex<CoefficientType> &K, const std::string &filename, const OSM::Chain<CoefficientType, OSM::COLUMN>& chain, int q, int cellId = -1) ;
    
protected:
    // Protected data
    int _dim; // Complex dimension
    std::vector<int> _size_bb; // Size of the bounding box Ni
    std::vector<int> _P; // P vector
    std::vector<bool> _cells; // Cells of the complex (true if the cell is present, false otherwise)
    std::vector<std::vector<int>> _base2bool; // Maps from base indices to indices in _cells
    std::vector<std::map<int, int>> _bool2base; // Maps from indices in _cells to base indices

    std::set<int> visited_cells;

    // Protected methods
    void initialize_cells(const CubObject& cub,typeComplexCube type); // Initialize _cells and _base2bool and _bool2base
    
    // khal coordinates
    std::vector<int> ind2khal(int index) const {
        if (index > _P[_dim])
            throw std::invalid_argument("Index exceeds the size of boolean vector");
        std::vector<int> khal(_dim);
        for (int k = 0; k < _dim; ++k) {
            khal[k] = index % _size_bb[k];
            index /= _size_bb[k];
        }
        return khal;
    }

    
    int khal2ind(const std::vector<int>& base_indices) const {
        if (base_indices.size() != _dim) {
            throw std::invalid_argument("Dimension of base_indices does not match _dim");
        }
        
        int cell_index = 0;
        for (int i = 0; i < _dim; ++i) {
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

    // vox coordinates (for DUAL construction)
    std::vector<int> ind2vox(int index, vector<int> N, int max_size) const
    {
        if (index > max_size)
            throw std::invalid_argument("ind2vox : index exceeding size of boolean vector");
        std::vector<int> coords(_dim);
        for (int k = 0; k < _dim; ++k) {
            coords[k] = index % N[k];
            index /= N[k];
        }
        return coords;
    }
    
    int vox2ind(const std::vector<int>& base_indices, vector<int> N, int max_size) const {
        if (base_indices.size() != _dim) {
            throw std::invalid_argument("Dimension of base_indices does not match _dim");
        }
        
        int cell_index = 0;
        for (int i = 0; i < _dim; ++i) {
            cell_index += base_indices[i] * N[i];
        }
        
        // Verfiy if cell_index is within bounds of _cells
        if (cell_index >= 0 && cell_index < max_size)
            return cell_index;
        else {
            std::cerr << "Invalid cell index in _cells: " << cell_index << std::endl;
            throw std::out_of_range("Invalid cell index in _cells");
        }
    }
    
    mutable std::vector<CMatrix>  _d; // Boundary matrix
    
    void  Calculate_d(int dim) const; // Calculate the boundary matrix of dimension dim
    
    void insert_cell(int cell); // Insert a cell into the complex
    
    int calculate_dimension(const std::vector<int>& cell) const; // Calculate the dimension of a cell
    int calculate_dimension(int cell_index) const { return calculate_dimension(ind2khal(cell_index)) ; }
    
    std::vector<int> calculate_boundaries(int cell) const; // Calculate the boundaries of a cell
    
};

// Initialization of static VTK_cubtypes
template <typename CoefficientType> const
std::vector<int> Cubical_complex<CoefficientType>::VTK_cubtypes({1, 3, 8, 11}) ;

// Constructor implementation
template<typename CoefficientType>
Cubical_complex<CoefficientType>::Cubical_complex(const CubObject& cub,typeComplexCube type) : _dim(cub.dim), _size_bb(_dim+1), _P(_dim+1,1), _base2bool(_dim+1), _bool2base(_dim+1)

{
    // Initialize _size_bb and _P
    if (type==PRIMAL)
        _size_bb = cub.N;
    else
    {
        for (int q=0; q<_dim; ++q)
            _size_bb.at(q) = 2*cub.N.at(q)+1 ;
    }
    
    for (int i = 1; i <= _dim; ++i) {
        _P[i] = _P[i - 1] * _size_bb[i - 1];
    }
    
    _cells.resize(_P[_dim], false) ;
    
    // Initialize _cells, _base2bool, and _bool2base
    initialize_cells(cub,type);
    
    // Initialize _d
    _d.resize(_dim+1);
    for (int q = 0; q <= _dim; ++q) {
        Calculate_d(q);
    }
    
}

template<typename CoefficientType>
void Cubical_complex<CoefficientType>::initialize_cells(const CubObject& cub, typeComplexCube type)
{

    if (type == PRIMAL)
    {
        
        for (int i=0; i<cub.cubs.size(); ++i)
        {
            const int id(khal2ind(cub.cubs.at(i))) ;
            insert_cell(id);
        }
        
        
    }
    else if (type == DUAL)
    {
        
        
        
        int max_size(1) ;
        for (int q=0; q<_dim; ++q)
            max_size *= cub.N.at(q) ;
        
        
        //We iterate over all the voxels via indices 
        for (int i=0; i<cub.cubs.size(); ++i)
        {
            vector<int> coords(cub.cubs.at(i)) ;
            // Calculate the coordinates of the voxel in the dual complex
            for (int i=0; i<_dim; ++i)
                coords.at(i)*=2 ;
            const int cell_index(khal2ind(coords)) ;
            
            _cells.at(cell_index) = true ;
            // Add the cell in _base2bool and _bool2base
            const int dim(calculate_dimension(coords)) ;
            const int n(_base2bool.at(dim).size()) ;
            _base2bool.at(dim).push_back(cell_index) ;
            _bool2base.at(dim)[cell_index] = n ;
        }
        
        // for the cells with dim>0
        for (int q = 1; q <= _dim; ++q) {
            
            for (int i = 0; i < _P[_dim]; ++i) {
               
                if (calculate_dimension(ind2khal(i)) == q) {
                    std::vector<int> boundaries = calculate_boundaries(i);
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
                        const int n(_base2bool.at(dim).size()) ;
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


template<typename CoefficientType>
bool Cubical_complex<CoefficientType>::is_valid_cell(const std::vector<int>& cell) const {
    for (int i=0; i<_dim; ++i) {
        if (cell[i] < 0 || cell[i] >= _size_bb[i]) {
            return false;
        }
    }
    return true;
}

template<typename CoefficientType>
bool Cubical_complex<CoefficientType>::is_valid_cell(int id_cell) const
{
    if ((id_cell < 0) || (id_cell > _P[_dim]))
        return false;
    else
        return true;
}



template<typename CoefficientType>
void Cubical_complex<CoefficientType>::insert_cell(int cell) {
    if (!is_valid_cell(cell))
        throw std::out_of_range("insert_cell: trying to insert cell with invalid index");
    
    // verify if the cell has already been visited
    if (visited_cells.find(cell) != visited_cells.end()) {
        return; // cell has already been visited
    }
    
    std::vector<int> cell_coords(ind2khal(cell));
    int dim = calculate_dimension(cell_coords);

    _cells[cell] = true;
    int cell_base_index = _base2bool[dim].size();
    _base2bool[dim].push_back(cell);
    _bool2base[dim][cell] = cell_base_index;

    std::vector<int> boundaries(calculate_boundaries(cell));
    for (const auto& boundary : boundaries) {
        if (!_cells[boundary]) {
            insert_cell(boundary);
        }
    }

    // mark cell as visited
    visited_cells.insert(cell);
}



template<typename CoefficientType>
void Cubical_complex<CoefficientType>::Calculate_d(int dim) const {
    int nb_lignes = (dim == 0) ? 0 : get_N(dim - 1);
    
    _d[dim] = CMatrix(nb_lignes, get_N(dim));
    
    // Iterate through the cells of dimension dim
    for (int i = 0; i < get_N(dim); ++i) {
        // Boundary of the i-th cell of dimension dim
        CChain boundary = boundary_cell(i, dim);
        
        // Insert the chain into the corresponding column of the boundary matrix
        OSM::setColumn(_d[dim], i, boundary);
    }
    
    
    
}

template<typename CoefficientType>
int Cubical_complex<CoefficientType>::calculate_dimension(const std::vector<int>& cell) const {
    int dimension = 0;
    for (int index : cell) {
        if (index % 2 == 1) { // Un index impair indique une dimension plus élevée
            dimension++;
        }
    }
    return dimension;
}

template<typename CoefficientType>
std::vector<int> Cubical_complex<CoefficientType>::calculate_boundaries(int idcell) const {
    std::vector<int> boundaries;
    std::vector<int> c = ind2khal(idcell);
    

    
    for (int i = 0; i < _dim; ++i) {
        if (c[i] % 2 == 1)
        {
            // Calculate the coefficient based on the number of odd entries in c from 0 to i-1
            int cell1 = idcell + _P[i];
            if (is_valid_cell(cell1))
                boundaries.push_back(cell1) ;
            
            
            int cell2 = idcell - _P[i];
            if (is_valid_cell(cell2))
                boundaries.push_back(cell2) ;
            
        }
    }

    return boundaries;
}

// tdata: labels encode boolean (for FSTAR and G) or int data (for PSC)
// labels: data to mark cells (e.g. (co)homology generators labels
//      -> for P,S,C labelling, the type is int => export all cells with their label
//      -> for generator, the type is bool => export only cells with label true
// scalar_type: type of labels (e.g. T)


//template <typename CoefficientType>
//void cubical_complex_write_to_vtk(const Cubical_complex<CoefficientType> &K, const std::string &filename, const std::vector<std::vector<int> > *labels=NULL, TDataLabels tdata=INT_LABELS, const std::string scalar_type="none")
//{
//    bool with_scalars = (labels != NULL) ;
//    int type_data(0) ; // 1 for int (ie. P,S,C labelling), 2 for bool (generators), 0 when no labels
//    if (with_scalars)
//    {
//        if (tdata == INT_LABELS)
//            type_data = 1 ;
//        else
//            type_data = 2 ;
//    }
//    // Load out file...
//    std::ofstream out ( filename, std::ios::out | std::ios::trunc);
//
//    if ( not out . good () ) {
//        std::cerr << "SimplicialComplex::loadFromFile. Fatal Error:\n  " << filename << " not found.\n";
//        throw std::runtime_error("File Parsing Error: File not found");
//    }
//
//    // Header
//    out << "# vtk DataFile Version 2.0" << endl ;
//    out << "generators" << endl ;
//    out << "ASCII" << endl ;
//    out << "DATASET  UNSTRUCTURED_GRID" << endl ;
//
//    // Points
//    int nnodes = K.nb_cells(0) ;
//    out << "POINTS " << nnodes << " double" << endl ;
//    for (int i=0; i<K.nb_cells(0); ++i)
//    {
//        const vector<int> coords(K.ind2khal(K._base2bool[0].at(i))) ;
//        for (int c : coords)
//            out << c/2 << " " ;
//        for (int i = coords.size(); i<3; ++i) // points must be 3D
//            out << "0 " ;
//        out << endl ;
//    }
//
//    int ncells_tot = 0, size_cells_tot = 0 ;
//    std::vector<int> types ;
//    std::vector<int> scalars ;
//    std::vector<int> ids ;
//    if ((type_data==0) || (type_data==1)) // all cells must be printed
//    {
//        // Cells up to dimension 3
//        // Number of cells : for each chain, each number of vertices for each cell
//        // Size : size of a cell of dimension q : 2^q
//        for (int q=0; q<=K._dim; ++q)
//        {
//            ncells_tot += K.nb_cells(q) ;
//            const int size_cell = 1<<q ;
//            size_cells_tot += (size_cell+1)*K.nb_cells(q) ;
//        }
//        out << "CELLS " << ncells_tot << " " << size_cells_tot << endl ;
//        // Output cells by increasing dimension
//        
//        // Vertices
//        for (int i = 0; i<K.nb_cells(0); ++i)
//        {
//            out << "1 " << i << endl ;
//            types.push_back(Cubical_complex<CoefficientType>::VTK_cubtypes.at(0)) ;
//            if (with_scalars)
//            {
//                scalars.push_back((*labels).at(0).at(i)) ;
//                ids.push_back(i) ;
//            }
//        }
//        // Cells of higher dimension
//        for (int q=1; q<=K._dim; ++q)
//        {
//            const int size_cell = 1<<q ;
//            for (int id =0; id < K.nb_cells(q); ++id)
//            {
//                vector<int> verts(K.khal_to_verts(K.ind2khal(K._base2bool.at(q).at(id)))) ;
//                out << size_cell << " " ;
//                for (int i : verts)
//                    out << i << " " ;
//                out << endl ;
//                types.push_back(Cubical_complex<CoefficientType>::VTK_cubtypes.at(q)) ;
//                if (with_scalars)
//                {
//                    scalars.push_back((*labels).at(q).at(id)) ;
//                    ids.push_back(id) ;
//                }
//            }
//        }
//        
//        
//        // CELL_TYPES
//        out << "CELL_TYPES " << ncells_tot << endl ;
//        for (int t : types)
//            out << t << " " ;
//        out << endl ;
//    }
//    else // output only cells labelled true
//    {
//        // 1 - Compute the number of cells / size of encoding
//        // Vertices
//        for (int i = 0; i<K.nb_cells(0); ++i)
//        {
//            if ((*labels).at(0).at(i) == 1)
//            {
//                ++ncells_tot;
//                size_cells_tot += 2;
//            }
//        }
//        // Cells of higher dimension
//        for (int q=1; q<=K._dim; ++q)
//        {
//            const int size_cell = 1<<q ;
//            for (int id =0; id < K.nb_cells(q); ++id)
//            {
//                if ((*labels).at(q).at(id) == 1)
//                {
//                    ++ncells_tot;
//                    size_cells_tot += (size_cell+1) ;
//                }
//            }
//        }
//        // 2 - Output cells
//        out << "CELLS " << ncells_tot << " " << size_cells_tot << endl ;
//        // Vertices
//        for (int i = 0; i<K.nb_cells(0); ++i)
//        {
//            if ((*labels).at(0).at(i) == 1)
//            {
//                out << "1 " << i << endl ;
//                types.push_back(Cubical_complex<CoefficientType>::VTK_cubtypes.at(0)) ;
//            }
//        }
//        // Cells of higher dimension
//        for (int q=1; q<=K._dim; ++q)
//        {
//            const int size_cell = 1<<q ;
//            for (int id =0; id < K.nb_cells(q); ++id)
//            {
//                if ((*labels).at(q).at(id) == 1)
//                {
//                    vector<int> khal(K.ind2khal(K._base2bool.at(q).at(id))) ;
//                    vector<int> verts(K.khal_to_verts(K.ind2khal(K._base2bool.at(q).at(id)))) ;
//                    out << size_cell << " " ;
//                    for (int i : verts)
//                        out << i << " " ;
//                    out << endl ;
//                    types.push_back(Cubical_complex<CoefficientType>::VTK_cubtypes.at(q)) ;
//                }
//            }
//        }
//        // CELL_TYPES
//        out << "CELL_TYPES " << ncells_tot << endl ;
//        for (int t : types)
//            out << t << " " ;
//        out << endl ;
//    }
//
//    if (with_scalars && (type_data==1))
//    {
//        // CELL_TYPES
//        out << "CELL_DATA " << ncells_tot << endl ;
//        out << "SCALARS Label " << scalar_type << " 1" << endl ;
//        out << "LOOKUP_TABLE default" << endl ;
//        for (int s : scalars)
//            out << s << " " ;
//        out << endl ;
//        // CELL_IDs
//        out << "SCALARS CellId " << scalar_type << " 1" << endl ;
//        out << "LOOKUP_TABLE default" << endl ;
//        for (int i : ids)
//            out << i << " " ;
//        out << endl ;
//    }
//    out.close() ;
//}

/** \brief Export complex to vtk (with int labels if provided).
*           -> All cells are exported */

template <typename CoefficientType>
void Cubical_complex<CoefficientType>::cubical_complex_to_vtk(const Cubical_complex<CoefficientType> &K, const std::string &filename, const std::vector<std::vector<int> > *labels)
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
    for (int n = 0; n < nnodes; ++n)
    {
        vector<double> p(K.get_vertex_coords(n)) ;
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
        // Size : size of a cell of dimension q : 2^q
        for (int q=0; q<=K.dim(); ++q)
        {
            ncells_tot += K.nb_cells(q) ;
            const int size_cell = 1<<q ;
            size_cells_tot += (size_cell+1)*K.nb_cells(q) ;
        }
        out << "CELLS " << ncells_tot << " " << size_cells_tot << endl ;
        // Output cells by increasing dimension
        
        // Vertices
        for (int i = 0; i<K.nb_cells(0); ++i)
        {
            out << "1 " << i << endl ;
            types.push_back(Cubical_complex<CoefficientType>::VTK_cubtypes.at(0)) ;
            if (with_scalars)
            {
                scalars.push_back((*labels).at(0).at(i)) ;
                ids.push_back(i) ;
            }
        }
        // Cells of higher dimension
        for (int q=1; q<=K.dim(); ++q)
        {
            const int size_cell = 1<<q ; //int_exp(2, q) ;
            for (int id =0; id < K.nb_cells(q); ++id)
            {
                vector<int> verts(K.khal_to_verts(K.ind2khal(K._base2bool.at(q).at(id)))) ;
                out << size_cell << " " ;
                for (int i : verts)
                    out << i << " " ;
                out << endl ;
                types.push_back(Cubical_complex<CoefficientType>::VTK_cubtypes.at(q)) ;
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

template <typename CoefficientType>
void Cubical_complex<CoefficientType>::cubical_complex_chain_to_vtk(const Cubical_complex<CoefficientType> &K, const std::string &filename, const OSM::Chain<CoefficientType, OSM::COLUMN>& chain, int q, int cellId)
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
    int nnodes = K.nb_cells(0) ;
    out << "POINTS " << nnodes << " double" << endl ;
    for (int n = 0; n < nnodes; ++n)
    {
        vector<double> p(K.get_vertex_coords(n)) ;
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
            const int size_cell = 1<<q ;
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
                    vector<int> khal(K.ind2khal(K._base2bool.at(q).at(id))) ;
                    vector<int> verts(K.khal_to_verts(K.ind2khal(K._base2bool.at(q).at(id)))) ;
                    out << size_cell << " " ;
                    for (int i : verts)
                        out << i << " " ;
                    out << endl ;
                    types.push_back(Cubical_complex<CoefficientType>::VTK_cubtypes.at(q)) ;
                    
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

    
#endif // CUBCOMPLEX_HPP
