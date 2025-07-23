//
//  tools_io.hpp
//  duality
//
//  Created by umenohana on 09/03/2023.
//

#ifndef tools_io_hpp
#define tools_io_hpp

#include <iostream>
#include <cassert>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>

namespace CGAL {
namespace HDVF {

using namespace std ;

// ------ For simplicial complexes
// List of vertices of a cell (set of ids)
typedef std::set<size_t> IOCellType ;
// List of cells
typedef std::vector<IOCellType> IOChainType ;

// ------ For cubical complexes
// Khalimsky coordinates of a cube
typedef std::vector<size_t> IOCubCellType ;
// List of cells
typedef std::vector<IOCellType> IOCubChainType ;


// Type of points (for vertices coordinates in R^d)
class IONodeType {
    vector<double> _coords;
    
public:
    IONodeType(size_t d = 3, double x = 0.) : _coords(d,x) {}
    IONodeType(std::vector<double> v) : _coords(v) {}
    IONodeType(const IONodeType& v) : _coords(v._coords) {}
    
    size_t size() const { return _coords.size(); }
    double at(size_t i) const { return _coords.at(i) ;}
    double& operator[](size_t i) { return _coords.at(i) ;}
    IONodeType& operator=(const IONodeType& v) { _coords = v._coords ; return *this ; }
    void push_back(double x) { _coords.push_back(x) ; }
    std::vector<double> get_coords() const { return _coords; }
    
    friend IONodeType & operator+(const IONodeType &v1, const IONodeType &v2)
    {
        assert(v1.size() == v2.size()) ;
        IONodeType &tmp = *new(IONodeType)(v1.size(),0) ;
        for (size_t i =0; i<v1.size(); ++i)
            tmp[i] = (v1.at(i) + v2.at(i)) ;
        return tmp ;
    }

    friend IONodeType & operator-(const IONodeType &v1, const IONodeType &v2)
    {
        assert(v1.size() == v2.size()) ;
        IONodeType &tmp = *new(IONodeType)(v1.size(),0) ;
        for (size_t i =0; i<v1.size(); ++i)
            tmp[i] = (v1.at(i) - v2.at(i)) ;
        return tmp ;
    }

    friend IONodeType & operator+= (IONodeType &v1, const IONodeType &v2)
    {
        assert(v1.size() == v2.size()) ;
        for (size_t i =0; i<v1.size(); ++i)
            v1[i] += v2.at(i) ;
        return v1 ;
    }

    friend IONodeType operator/(const IONodeType &v1, double d)
    {
        IONodeType &tmp = *new(IONodeType)(v1.size(),0) ;
        for (size_t i =0; i<v1.size(); ++i)
            tmp[i] = v1.at(i)/d ;
        return tmp ;
    }

    friend IONodeType & operator/=(IONodeType &v1, double d)
    {
        for (size_t i =0; i<v1.size(); ++i)
            v1[i] /= d ;
        return v1 ;
    }

    friend IONodeType & operator*=(IONodeType &v, double d)
    {
        for (size_t i =0; i<v.size(); ++i)
            v[i] *= d ;
        return v ;
    }

    friend IONodeType & operator*=(IONodeType &v, IONodeType &lambda)
    {
        for (size_t i =0; i<v.size(); ++i)
            v[i] *= lambda.at(i) ;
        return v ;
    }

    friend double dist(const IONodeType &v1, const IONodeType &v2)
    {
        assert(v1.size() == v2.size()) ;
        double x=0., tmp ;
        for (size_t i =0; i<v1.size(); ++i)
        {
            tmp = v2.at(i) - v1.at(i) ;
            x += tmp*tmp ;
        }
        
        return sqrt(x) ;
    }

    friend IONodeType max(const IONodeType &v1, const IONodeType &v2)
    {
        assert(v1.size() == v2.size()) ;
        IONodeType tmp(v1) ;
        for (size_t i=0; i<v1.size(); ++i)
            tmp[i] = std::max(tmp.at(i), v2.at(i)) ;
        return tmp ;
    }

    friend IONodeType min(const IONodeType &v1, const IONodeType &v2)
    {
        assert(v1.size() == v2.size()) ;
        IONodeType tmp(v1) ;
        for (size_t i=0; i<v1.size(); ++i)
            tmp[i] = std::min(tmp.at(i), v2.at(i)) ;
        return tmp ;
    }

    friend void normalize(IONodeType &v)
    {
        const IONodeType zero = IONodeType(v.size(),0) ;
        double n = dist(zero, v) ;
        for (size_t i =0; i<v.size(); ++i)
            v[i] /= n ;
    }
};




// ----- VTK format -----

//// Associated to any dimension the type number of associated VTK cells
static std::vector<size_t> VTK_types_IO = {1, 3, 5, 10} ;

// Write vtk file
template <typename CoefficientType>
void write_vtk(const std::string &filename, const std::vector<IONodeType> &nodes, const std::vector<IOChainType> &chains, const std::vector<CoefficientType> *labels=NULL, const std::string scalar_type="none")
{
    bool with_scalars = (labels != NULL) ;
    // Load ...
    std::ofstream out ( filename, std::ios::out | std::ios::trunc);
    
    if ( not out . good () ) {
        std::cerr << "SimplicialComplex::loadFromFile. Fatal Error:\n  " << filename << " not found.\n";
        throw std::runtime_error("File Parsing Error: File not found");
    }
    
    // Header
    out << "# vtk DataFile Version 2.0" << endl ;
    out << "generators" << endl ;
    out << "ASCII" << endl ;
    out << "DATASET  UNSTRUCTURED_GRID" << endl ;
    
    // Points
    size_t nnodes = nodes.size() ;
    out << "POINTS " << nnodes << " double" << endl ;
    for (IONodeType n : nodes)
        out << n.at(0) << " " << n.at(1) << " " << n.at(2) << endl ;
    
    // Cells
    // Number of cells : for each chain, each number of vertices for each cell
    // Size : size of a cell defined by d vertices : d+1
    size_t ncells_tot = 0, size_cells_tot = 0 ;
    for (size_t i = 0; i<chains.size(); ++i)
    {
        // for each cells in the ith chain
        for (IOCellType c : chains.at(i))
        {
            ncells_tot += 1 ;
            size_cells_tot += (c.size()+1) ;
        }
        
    }
    out << "CELLS " << ncells_tot << " " << size_cells_tot << endl ;
    // Output cells
    std::vector<size_t> types ;
    std::vector<CoefficientType> scalars ;
    for (size_t i = 0; i<chains.size(); ++i)
    {
        const IOChainType cc = chains.at(i) ;
        
        for (IOCellType c : cc)
        {
            const int d = c.size() ;
            const size_t cell_type = VTK_types_IO.at(d-1) ;
            
            out << d << " " ;
            for (size_t j : c)
                out << j << " " ;
            out << endl ;
            types.push_back(cell_type) ;
            if (with_scalars)
                scalars.push_back(labels->at(i)) ;
        }
    }
    
    // CELL_TYPES
    out << "CELL_TYPES " << ncells_tot << endl ;
    for (size_t t : types)
        out << t << " " ;
    out << endl ;
    
    if (with_scalars)
    {
        // CELL_TYPES
        out << "CELL_DATA " << ncells_tot << endl ;
        out << "SCALARS CriticalCellsId " << scalar_type << " 1" << endl ;
        out << "LOOKUP_TABLE default" << endl ;
        for (CoefficientType s : scalars)
            out << s << " " ;
        out << endl ;
    }
    out.close() ;
}
//}


// ///////////////////

// ----- Generic -----

inline bool check_sanity_line(const std::string &line, const std::string &file)
{
    // Check that line is sanitized. If not, throw.
    for ( size_t i = 0; i < line.size(); ++ i ) {
        if ( ! ( std::isspace(line[i]) || std::isdigit(line[i]) || (line[i] == '-') || (line[i] == '.') || (line[i] == 'e')) ) {
            std::cerr << "Error:\n  Cannot parse file " << file << std::endl ;
            std::cerr << "--- " << line[i] << std::endl ;
            return false;
        }
    }
    return true ;
}

inline bool get_next_uncommented_line(std::ifstream &infile, std::string &result) {
    while(getline(infile,result)) {
        if(result.length() > 1 && result[0] != '#') {
            return true;
        }
    }
    return false;
}

// Generic Mesh_object class - for 3D triangular meshes
class Mesh_object
{
public:
    // Dimension :
    // - if dim > 0 : Mesh_object encodes a mesh and all cells have dimension d
    // - if dim < 0 : Mesh_object encodes a complex (possibly incomplete) of dimension d
    int dim = 0 ;
    size_t nvertices, ncells, nedges ;
    std::vector<IONodeType> nodes ; // Coordinates of vertices (optional)
    std::vector<IOCellType> cells ;
    
    Mesh_object() : dim(0), nvertices(0), ncells(0), nedges(0) {}
    
    Mesh_object(int d, const std::vector<IONodeType> &vnodes, const std::vector<IOCellType> &vcells) : dim(d), nvertices(vnodes.size()), ncells(vcells.size()), nedges(0), nodes(vnodes), cells(vcells) { check_dimension() ;}

    Mesh_object(int d, const std::vector<std::vector<double> > &vnodes, const std::vector<IOCellType> &vcells) : dim(d), nvertices(vnodes.size()), ncells(vcells.size()), nedges(0), cells(vcells)
    {
        for (std::vector<double> v : vnodes)
        {
            nodes.push_back(IONodeType(v)) ;
        }
        check_dimension() ;
    }
    
    Mesh_object(const Mesh_object &m) : dim(m.dim), nvertices(m.nvertices), ncells(m.ncells), nedges(m.nedges), nodes(m.nodes), cells(m.cells) {}
    
    std::vector<std::vector<double> > get_nodes ()
    {
        std::vector<std::vector<double> > res ;
        for (IONodeType v : nodes)
            res.push_back(v.get_coords()) ;
        return res ;
    }
    
    // Mesh operations
    void push_back(const Mesh_object &mesh)
    {
        size_t off = nvertices ; // The index of all the cells of mesh has to be incremented by off
        nvertices += mesh.nvertices ;
        ncells += mesh.ncells ;
        // Append all the vertices
        for (size_t i=0; i<mesh.nodes.size(); ++i)
            nodes.push_back(mesh.nodes.at(i)) ;
        // Append all the cells (and increment their indices by off)
        for (size_t i=0; i<mesh.cells.size(); ++i)
        {
            IOCellType tmp ;
            for (size_t c : mesh.cells.at(i))
                tmp.insert(c+off) ;
            cells.push_back(tmp) ;
        }
    }
    
    void add_node(const IONodeType &v) {nodes.push_back(v); ++nvertices ;}
    
    void clear_cells() { cells.clear() ; ncells = 0 ; }
    
    void clear_nodes() { nodes.clear() ; nvertices = 0 ; }
    
    void clear() { clear_nodes() ; clear_cells() ;}
    
    void add_cell(const IOCellType &c) {cells.push_back(c); ++ncells ;}
    
    size_t cells_of_dim (int q) const
    {
        size_t n = 0 ;
        for (IOCellType c : cells)
        {
            if (c.size() == (q+1))
                ++n ;
        }
        return n ;
    }
    // OFF
    bool read_off(const std::string &filename)
    {
        dim = 0 ;
        // 0 - open input file
        std::ifstream infile(filename);
        if(!infile.is_open()) {
            // failed to open the file
            cerr << "Warning: file " << filename << " does not exist" << endl ;
            return false;
        }
        // 1 - header
        std::string header;
        if (!get_next_uncommented_line(infile, header)) {
            return false;
        }
        // todo : check for header == "off"
        
        // 2 - number of vertices, number of faces, number of edges (can be ignored)
        std::string info;
        if (!get_next_uncommented_line(infile, info)) {
            return false;
        }
        std::istringstream info_stream;
        info_stream.str(info);
        
        info_stream >> nvertices >> ncells >> nedges;
        
        nodes.resize(nvertices) ;
        for(auto i=0; i < nvertices; ++i) {
            if (!get_next_uncommented_line(infile,info)) {
                return false;
            }
            std::istringstream info_stream(info);
            IONodeType p(3) ;
            info_stream >> p[0] >> p[1] >> p[2] ;
            nodes[i] = p ;
        }
        
        // 4 - the actual faces
        cells.resize(ncells);
        for(auto i=0; i < ncells; ++i) {
            if (!get_next_uncommented_line(infile,info)) {
                return false;
            }
            std::istringstream info_stream(info);
            unsigned long n;
            unsigned long index;
            info_stream >> n;
            IOCellType c ;
            for (auto j = 0; j < n; ++j) {
                info_stream >> index;
                c.insert(index) ;
            }
            cells[i] = c ;
            dim = (c.size()-1>dim)?c.size()-1:dim ;
        }
        
        infile.close();
        return true;
    }
    
    bool write_off(const std::string &filename)
    {
        // 0 - open input file
        std::ofstream outfile(filename);
        if(!outfile.is_open()) {
            cerr << "Warning: cannot open file " << filename << endl ;
            // failed to open the file
            return false;
        }
        // 1 - header
        outfile << "OFF" << endl ;
        outfile << nvertices << " " << cells_of_dim(2) << " " << nedges << endl ;
        // 2 - nodes
        for (size_t i=0; i<nvertices; ++i)
        {
            outfile << nodes.at(i).at(0) << " " << nodes.at(i).at(1) << " " << nodes.at(i).at(2) << endl ;
        }
        // 3 - cells (export only triangles)
        for (size_t i=0; i<ncells; ++i)
        {
            const size_t ni = cells.at(i).size() ;
            if (ni == 3)
            {
                outfile << ni << " " ;
                for (size_t k : cells.at(i))
                    outfile << k << " " ;
                outfile << endl ;
            }
        }
        outfile.close() ;
        return true;
    }
    
    // VTK
    void write_to_vtk(const std::string &filename)
    {
        std::vector<IOChainType> chains{cells} ;
        write_vtk<int>(filename, nodes, chains) ;
    }
    
    // SIMP
    bool write_simp(const std::string &filename)
    {
        // 0 - open input file
        std::ofstream outfile(filename);
        if(!outfile.is_open()) {
            cerr << "Warning: cannot open file " << filename << endl ;
            // failed to open the file
            return false;
        }
        // 1 - write cells
        for (size_t i=0; i<ncells; ++i)
        {
            const IOCellType cell = cells.at(i) ;
            for (size_t c : cell)
                outfile << c << " " ;
            outfile << endl ;
        }
        outfile.close() ;
        return true;
    }
    
    bool read_simp(const std::string &filename)
    {
        int d = 0 ;
        // 0 - open input file
        std::ifstream infile(filename);
        if(!infile.is_open()) {
            // failed to open the file
            cerr << "Warning: file " << filename << " does not exist" << endl ;
            return false;
        }
        std::size_t line_number = 0;
        while ( !(infile.eof()) )
        {
            ++line_number;
            std::string line;
            getline( infile, line );
            // Check that line is sanitized. If not, throw.
            for ( size_t i = 0; i < line.size(); ++ i ) {
                if ( ! ( std::isspace(line[i]) || std::isdigit(line[i]) ) ) {
                    std::cerr << "Fatal Error:\n  Cannot parse line #" << line_number << " of " << filename << "\n";
                    std::cerr << " --> " << line << "\n";
                    throw std::runtime_error("File Parsing Error: Invalid file");
                }
            }
            IOCellType cell ;
            std::istringstream is( line );
            size_t v;
            while ( is >> v )
                cell.insert(v);
            // Add this simplex to cells
            if (!(cell.empty()))
            {
                ++ncells ;
                cells.push_back(cell) ;
                const int dcell = cell.size()-1 ;
                if (dcell > d)
                    d = dcell ;
            }
        }
        dim = d ;
        infile.close() ;
        return true ;
    }
    
    bool read_nodes_file(const std::string &filename) ;
    
    void print_infos () const
    {
        cout << "Mesh_object infos - dim : "<< dim << ", nodes : " << nodes.size() << ", cells : " << cells.size() << endl ;
        for (int q = 0; q <= dim; ++q)
            cout << "cells of dim " << q << " : " << cells_of_dim(q) << endl ;
    }
    
    // Mesh computations
    IONodeType barycenter()
    {
        // Init the barycenter
        IONodeType bary(3,0) ;
        // Compute the barycenter
        for (IONodeType v : nodes)
        {
            bary += v ;
        }
        bary /= nodes.size() ;
        return bary;
    }
    
    double radius(const IONodeType &bary)
    {
        double r = 0 ;
        for (IONodeType v : nodes)
        {
            r = std::max (r, dist(v,bary)) ;
        }
        return r ;
    }
    pair<IONodeType, IONodeType> BB(double ratio=1.)
    {
        IONodeType minBB(nodes.at(0)), maxBB(nodes.at(0)) ;
        for (size_t i=1; i<nodes.size(); ++i)
        {
            minBB = min(minBB, nodes.at(i)) ;
            maxBB = max(maxBB, nodes.at(i)) ;
        }
        IONodeType c = (minBB+maxBB)/2. ;
        IONodeType rad = (maxBB-minBB)/2. ;
        rad *= ratio ;
        return std::pair<IONodeType, IONodeType>(c-rad, c+rad) ;
    }
private:
    void check_dimension()
    {
        if (dim > 0) // exact dimension
        {
            for (IOCellType cell : cells)
            {
                if (cell.size() != (dim+1))
                    throw "Mesh has a cell of inconsistent dimension" ;
            }
        }
        else // max dim
        {
            for (IOCellType cell : cells)
            {
                if ((cell.size() > (-dim+1)))
                    throw "Mesh has a cell of inconsistent dimension" ;
            }
        }
    }
} ;

inline Mesh_object mesh_BB(const IONodeType &BBmin, const IONodeType &BBmax)
{
    Mesh_object m ;
    IONodeType delta = BBmax-BBmin ;
    m.nvertices = 8 ;
    m.ncells = 12 ;
    m.nodes.resize(8) ;
    m.nodes[0] = IONodeType({0, 0, 0}) ;
    m.nodes[1] = IONodeType({1, 0, 0}) ;
    m.nodes[2] = IONodeType({1, 1, 0}) ;
    m.nodes[3] = IONodeType({0, 1, 0}) ;
    m.nodes[4] = IONodeType({0, 0, 1}) ;
    m.nodes[5] = IONodeType({1, 0, 1}) ;
    m.nodes[6] = IONodeType({1, 1, 1}) ;
    m.nodes[7] = IONodeType({0, 1, 1}) ;
    for (size_t i=0; i<8; ++i)
    {
        m.nodes[i] *= delta ;
        m.nodes[i] += BBmin ;
    }
    
    m.cells.resize(12) ;
    m.cells[0] = IOCellType({0, 1, 4}) ;
    m.cells[1] = IOCellType({1, 5, 4}) ;
    m.cells[2] = IOCellType({1, 2, 6}) ;
    m.cells[3] = IOCellType({1, 6, 5}) ;
    m.cells[4] = IOCellType({3, 1, 0}) ;
    m.cells[5] = IOCellType({2, 1, 3}) ;
    m.cells[6] = IOCellType({2, 3, 6}) ;
    m.cells[7] = IOCellType({6, 3, 7}) ;
    m.cells[8] = IOCellType({3, 0, 4}) ;
    m.cells[9] = IOCellType({7, 3, 4}) ;
    m.cells[10] = IOCellType({6, 4, 5}) ;
    m.cells[11] = IOCellType({6, 7, 4}) ;
    return m;
}

// TODO : check that de dimensions of nodes are consistent...



// Tetgen related

inline size_t read_nodes(const std::string &node_file, bool load_nodes, std::vector<IONodeType> *nodes)
{
    std::ifstream in_file (node_file) ;
    if ( not in_file . good () ) {
        std::cerr << "SimplicialComplex::loadFromFile. Fatal Error:\n  " << node_file << " not found.\n";
        throw std::runtime_error("File Parsing Error: File not found");
    }
    
    // First line is the number of nodes
    size_t nnodes, nnodes_tmp ;
    if (not in_file.eof())
    {
        std::string line;
        getline( in_file, line );
        check_sanity_line(line, node_file) ;
        // First number is the number of nodes
        std::istringstream is (line);
        is >> nnodes ;
    }
    nnodes_tmp =  nnodes ;
    while ( !(in_file.eof()) && (nnodes_tmp>0))
    {
        size_t trash ;
        double x ;
        IONodeType node ;
        --nnodes_tmp ;
        std::string line;
        getline( in_file, line );
        check_sanity_line(line, node_file) ;
        std::istringstream is (line);
        is >> trash ;
        for (int i = 0; i<3; ++i)
        {
            is >> x ;
            node.push_back(x) ;
        }
        if (load_nodes)
            nodes->push_back(node) ;
    }
    in_file.close() ;
    return nnodes;
}

inline bool Mesh_object::read_nodes_file(const std::string &filename)
{
    nvertices = read_nodes(filename, true, &nodes) ;
}

// Tetgen

class TetObject : public Mesh_object
{
public:
    TetObject(const std::string & prefix) : Mesh_object(), _prefix(prefix)
    {
        dim = -3 ;
        add_nodes() ;
        create_nodes() ;
        add_edges() ;
        add_faces() ;
        add_tets() ;
        ncells = cells.size() ;
    }
    
    void add_nodes()
    {
        const std::string file_node = fnodes_from_prefix(_prefix) ;
        cout << "file_node : " << file_node << endl ;
        read_nodes_file(file_node) ;
    }
    
    void create_nodes()
    {
        for (size_t i=0; i<nvertices; ++i)
        {
            IOCellType cell({i}) ;
            cells.push_back(cell) ;
        }
        cout << "--- " << nvertices << "vert" << endl ;
    }
    
    void add_edges()
    {
        const std::string file_edge = fedges_from_prefix(_prefix) ;
        std::ifstream input_file ( file_edge );
        size_t f_nedges ;
        // Open the input file
        if ( not input_file . good () )
        {
            std::cerr << "File " << file_edge << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }
        // First line is the number of edges
        std::size_t line_number = 0;
        if (not input_file.eof())
        {
            ++line_number ;
            std::string line;
            getline( input_file, line );
            check_sanity_line(line, file_edge) ;
            // First number is the number of nodes
            std::istringstream is (line);
            is >> f_nedges ;
        }
        while ( !(input_file.eof()) && (f_nedges>0))
        {
            size_t trash, i, j ;
            ++ line_number;
            --f_nedges ;
            std::string line;
            getline( input_file, line );
            check_sanity_line(line, file_edge) ;
            std::istringstream is (line);
            is >> trash ;
            is >> i ;
            is >> j ;
            IOCellType cell({i,j}) ;
            add_cell(cell) ;
        }
        input_file.close() ;
        cout << "--- " << f_nedges << " edges" << endl ;
    }
    
    void add_faces()
    {
        const std::string file_face = ffaces_from_prefix(_prefix) ;
        std::ifstream input_file ( file_face );
        size_t f_nfaces ;
        // Open the input file
        if ( not input_file . good () )
        {
            std::cerr << "File " << file_face << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }
        // First line is the number of edges
        std::size_t line_number = 0;
        if (not input_file.eof())
        {
            ++line_number ;
            std::string line;
            getline( input_file, line );
            check_sanity_line(line, file_face) ;
            // First number is the number of nodes
            std::istringstream is (line);
            is >> f_nfaces ;
        }
        while ( !(input_file.eof()) && (f_nfaces>0))
        {
            size_t trash, i, j, k ;
            ++ line_number;
            --f_nfaces ;
            std::string line;
            getline( input_file, line );
            check_sanity_line(line, file_face) ;
            std::istringstream is (line);
            is >> trash ;
            is >> i ;
            is >> j ;
            is >> k ;
            IOCellType cell({i, j, k}) ;
            add_cell(cell) ;
        }
        input_file.close() ;
        cout << "--- " << f_nfaces << " faces" << endl ;
    }
    
    void add_tets()
    {
        const std::string file_ele = ftets_from_prefix(_prefix) ;
        std::ifstream input_file ( file_ele );
        size_t f_ntets ;
        // Open the input file
        if ( not input_file . good () )
        {
            std::cerr << "File " << file_ele << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }
        // First line is the number of edges
        std::size_t line_number = 0;
        if (not input_file.eof())
        {
            ++line_number ;
            std::string line;
            getline( input_file, line );
            check_sanity_line(line, file_ele) ;
            // First number is the number of nodes
            std::istringstream is (line);
            is >> f_ntets ;
        }
        while ( !(input_file.eof()) && (f_ntets>0))
        {
            size_t trash, i, j, k, l ;
            ++ line_number;
            --f_ntets ;
            std::string line;
            getline( input_file, line );
            check_sanity_line(line, file_ele) ;
            std::istringstream is (line);
            is >> trash ;
            is >> i ;
            is >> j ;
            is >> k ;
            is >> l ;
            IOCellType cell({i, j, k, l}) ;
            add_cell(cell) ;
        }
        input_file.close() ;
        cout << "--- " << f_ntets << " tets" << endl ;
    }
private:
    std::string _prefix ;
    std::string fnodes_from_prefix(const std::string &prefix) {return prefix+".1.node"; }
    std::string fedges_from_prefix(const std::string &prefix) {return prefix+".1.edge"; }
    std::string ffaces_from_prefix(const std::string &prefix) {return prefix+".1.face"; }
    std::string ftets_from_prefix(const std::string &prefix) {return prefix+".1.ele"; }
} ;



// ----- mesh BB

//Mesh_object mesh_BB(const IONodeType &BBmin, const IONodeType &BBmax) ;

class IcosphereObject : public Mesh_object
{
public:
    using Index=size_t ;
    using Lookup=std::map<std::pair<Index, Index>, Index>;
    
    IcosphereObject(size_t subdivisions, const IONodeType &c = IONodeType({0, 0, 0}), double r=1.) : Mesh_object(2, vertices_ico, triangles_ico)
    {
        for (size_t i=0; i<subdivisions; ++i)
        {
            subdivide();
        }
        rigid_transformation(c, r) ;
    }
    
    
    // Icosahedron
    inline static const float X=.525731112119133606f;
    inline static const float Z=.850650808352039932f;
    inline static const float N=0.f;
    
    inline static const std::vector<IONodeType> vertices_ico =
    {
        IONodeType({-X,N,Z}), IONodeType({X,N,Z}), IONodeType({-X,N,-Z}), IONodeType({X,N,-Z}),
        IONodeType({N,Z,X}), IONodeType({N,Z,-X}), IONodeType({N,-Z,X}), IONodeType({N,-Z,-X}),
        IONodeType({Z,X,N}), IONodeType({-Z,X, N}), IONodeType({Z,-X,N}), IONodeType({-Z,-X, N})
    };
    
    inline static const std::vector<IOCellType> triangles_ico =
    {
        {0,4,1},{0,9,4},{9,5,4},{4,5,8},{4,8,1},
        {8,10,1},{8,3,10},{5,3,8},{5,2,3},{2,7,3},
        {7,10,3},{7,6,10},{7,11,6},{11,0,6},{0,1,6},
        {6,1,10},{9,0,11},{9,11,2},{9,2,5},{7,2,11}
    };
    
    // Methods
    Index vertex_for_edge(Lookup& lookup, Index first, Index second)
    {
        Lookup::key_type key(first, second);
        if (key.first>key.second)
            std::swap(key.first, key.second);
        
        auto inserted=lookup.insert({key, nodes.size()});
        if (inserted.second)
        {
            IONodeType& edge0(nodes[first]);
            IONodeType& edge1(nodes[second]);
            IONodeType& point(edge0+edge1);
            normalize(point) ;
            add_node(point);
        }
        
        return inserted.first->second;
    }
    
    void subdivide()
    {
        Lookup lookup;
        std::vector<IOCellType> result ;
        
        for (IOCellType & each:cells)
        {
            std::array<Index, 3> mid;
            std::vector<size_t> each_vertex ;
            for (size_t c : each)
                each_vertex.push_back(c) ;
            for (size_t edge=0; edge<3; ++edge)
            {
                mid[edge]=vertex_for_edge(lookup,
                                          each_vertex[edge], each_vertex[(edge+1)%3]);
            }
            result.push_back({each_vertex[0], mid[0], mid[2]});
            result.push_back({each_vertex[1], mid[1], mid[0]});
            result.push_back({each_vertex[2], mid[2], mid[1]});
            result.push_back({mid[0], mid[1], mid[2]});
        }
        clear_cells() ;
        for (IOCellType c : result)
            add_cell(c) ;
    }
    
    void rigid_transformation(const IONodeType &c, double r)
    {
        for (size_t i=0; i<nvertices; ++i)
        {
            nodes[i] *= r ;
            nodes[i] += c ;
        }
    }
};

// ///////////////////

// Generic Cub_object class - for 2D/3D cubical meshes
class Cub_object
{
public:
    int dim = 0 ; // Dimension of the complex
    vector<size_t> ncubs ; // Number of cubs in each dimension
    vector<size_t> N ; // Size of BB along each dimension
    std::vector<IOCubCellType> cubs ;
    bool khalimsky ;
    
    Cub_object() : dim(0), ncubs(vector<size_t>(4)), N(vector<size_t>(3)), khalimsky(false) {}
    // TODO: check !!! 4/3 ...
    Cub_object(int d, const std::vector<IOCubCellType> &vcubs, bool khal = false) : dim(d), ncubs(vector<size_t>(4)), N(vector<size_t>(3)), cubs(vcubs), khalimsky(khal)
    { check_dimension() ;}
    Cub_object(const Cub_object &m) : dim(m.dim), ncubs(m.ncubs), N(m.N), cubs(m.cubs), khalimsky(m.khalimsky) {}
    
    // Mesh operations
    void clear_cubs() { cubs.clear() ; for (size_t i=0; i<4; ++i) ncubs[i] = 0 ; }
    void add_cub(const IOCubCellType &c) {cubs.push_back(c); ++ncubs[cub_dim(c)] ;}
    void frame() // Enlarge the bouding box to add 1 voxel around
    {
        // Enlarge the BB along each dimension
        for (size_t i=0; i<N.size(); ++i)
        {
            if (khalimsky)
                N.at(i) += 4 ;
            else
                N.at(i) += 2 ;
        }
        // Shift cells coordinates
        for (size_t i = 0; i<cubs.size(); ++i)
        {
            for (int k = 0; k<dim; ++k)
            {
                if (khalimsky)
                    cubs.at(i).at(k) += 2 ;
                else
                    ++cubs.at(i).at(k) ;
            }
        }
    }
    
    // PGM
    bool read_pgm(const std::string &filename, bool khal = false)
    {
        khalimsky = khal ;
        size_t row = 0, col = 0, numrows = 0, numcols = 0;
        ifstream infile(filename);
        if(!infile.is_open()) {
            cerr << "Warning: cannot open file " << filename << endl ;
            // failed to open the file
            return false;
        }
        
        stringstream ss;
        string inputLine = "";
        
        // First line : version
        getline(infile,inputLine);
        if(inputLine.compare("P2") != 0) cerr << "Version error" << endl;
        else cout << "Version : " << inputLine << endl;
        
        
        // Second line: dimensions
        getline(infile,inputLine);
        vector<size_t> sizes ;
        
        size_t tmp ;
        stringstream sseizes (inputLine);
        while (sseizes >> tmp)
            sizes.push_back(tmp) ;
        
        cout << "dimensions : " ;
        for (size_t i=0; i<sizes.size(); ++i)
            cout << sizes.at(i) << " " ;
        cout << endl ;
        dim = sizes.size() ;
        N = vector<size_t>(dim) ;
        for (size_t i=0; i<dim; ++i)
            N.at(i) = sizes.at(i) ;
        
        // Throw away next data (255)
        
        getline(infile,inputLine);
        
        // Remainder: data
        // Read by decreasing dimension
        // In order to get : read row, then column, then depth...
        
        // Continue with a stringstream
        ss << infile.rdbuf();
        
        size_t NN = 1 ;
        for (int i=0; i<dim; ++i)
            NN = NN * N.at(i) ;
        
        // Following lines : data
        for(size_t i = 0; i < NN; ++i)
        {
            ss >> tmp ;
            if (tmp) // pixel on
            {
                cubs.push_back(index_to_coords(i, khal)) ;
                ++ncubs[dim] ;
            }
        }
        
        if (khal)
        {
            // Change N to Khalimsky maxima
            for(int i=0; i<dim; ++i)
                N.at(i) = 2*N.at(i)+1 ;
        }
        
        infile.close();
    }
    
    bool write_pgm(const std::string &filename) ;
    // CUB
    bool read_cub(const std::string &filename, bool khalimsky = false)
    {
        // 0 - open input file
        std::ifstream infile(filename);
        if(!infile.is_open()) {
            // failed to open the file
            cerr << "Warning: file " << filename << " does not exist" << endl ;
            return false;
        }
        std::size_t line_number = 0;
        while ( !(infile.eof()) )
        {
            ++line_number;
            std::string line;
            getline( infile, line );
            // Check that line is sanitized. If not, throw.
            for ( size_t i = 0; i < line.size(); ++ i ) {
                if ( ! ( std::isspace(line[i]) || std::isdigit(line[i]) ) ) {
                    std::cerr << "Fatal Error:\n  Cannot parse line #" << line_number << " of " << filename << "\n";
                    std::cerr << " --> " << line << "\n";
                    throw std::runtime_error("File Parsing Error: Invalid file");
                }
            }
            std::istringstream is( line );
            
            if (line_number == 1) // Header 1
                is >> dim ;
            else if (line_number == 2) // Header 2
            {
                for (int i=0; i<dim; ++i)
                    is >> N[i] ;
            }
            else // Cub line
            {
                IOCubCellType cub ;
                
                size_t v;
                while ( is >> v )
                    cub.push_back(v);
                if (cub.size() != dim)
                    cout << "Discard line " << line_number << " cub of invalid dimension" << endl ;
                else
                {
                    // Add this cub to cubs
                    ++ncubs[cub_dim(cub)] ;
                    
                    if (khalimsky)
                        cubs.push_back(cub) ;
                    else
                    {
                        size_t ind(khal_to_index(cub)) ;
                        cubs.push_back(index_to_coords(ind, false)) ;
                    }
                }
            }
        }
        infile.close() ;
        return true ;
    }
    
    bool write_cub(const std::string &filename) ;
    
    void print_infos (size_t level = 0) const
    {
        cout << "Cub_object infos - dim : " << dim << ", cubs : " << cubs.size() << endl ;
        for (int q=0; q < dim; ++q)
            cout << "\tSize along dim " << q << " : " << N[q] << endl ;
        for (int q = 0; q <= dim; ++q)
            cout << "Cubs of dim " << q << " : " << ncubs[q] << endl ;
        if (level == 1) // Print coordinates of cubes
        {
            for (size_t i=0; i<cubs.size(); ++i)
            {
                for (size_t j : cubs.at(i))
                    cout << j << " " ;
                cout << endl ;
            }
        }
    }
    
private:
    void check_dimension()
    {
        if (khalimsky)
        {
            for (IOCubCellType c : cubs)
                ++ncubs[cub_dim(c)] ;
        }
        else
        {
            ncubs[dim] = cubs.size() ;
        }
    }
    
    inline int cub_dim (IOCubCellType c)
    {
        int q = 0 ;
        for (size_t i : c)
        {
            if (i%2 == 1)
                ++q ;
        }
        return q ;
    }
    
    inline IOCubCellType index_to_coords(size_t i, bool khalimsky = false)
    {
        IOCubCellType coords ;
        // Convert index in binary object to size_t coordinates
        for (int q = 0; q < dim; ++q)
        {
            coords.push_back(i % N[q]) ;
            i = i / (N[q]) ;
        }
        if (khalimsky)
        {
            // Convert int coordinates to Khalimsky
            for (int q = 0; q < dim; ++q)
            {
                coords.at(q) = 2*coords.at(q)+1 ;
            }
        }
        return coords ;
    }
    
    inline size_t khal_to_index (const IOCubCellType& coords)
    {
        size_t cell_index(0);
        for (int i = 0; i < dim; ++i) {
            cell_index += coords[i] * N[i];
        }
    }
} ;

} /* end namespace HDVF */
} /* end namespace CGAL */


#endif /* tools_io_hpp */
