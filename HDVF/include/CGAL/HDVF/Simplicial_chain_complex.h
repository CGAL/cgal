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

#ifndef CGAL_HDVF_SIMPLICIAL_CHAIN_COMPLEX_H
#define CGAL_HDVF_SIMPLICIAL_CHAIN_COMPLEX_H

//#include <CGAL/license/HDVF.h>

#include <vector>
#include <map>
#include <CGAL/HDVF/tools_io.h>
#include <CGAL/HDVF/Abstract_simplicial_chain_complex.h>
#include <CGAL/OSM/OSM.h>

namespace CGAL {
namespace HDVF {

// Forward declaration of SimpComplexTools
template<typename CoefficientType> class Duality_simplicial_complex_tools ;

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Simplicial_chain_complex` refines the `Abstract_simplicial_chain_complex` class by assigning coordinates to vertices (i.e.\ 0-simplices). Hence, vtk output is available.

 \cgalModels{GeometricChainComplex}

 \tparam CoefficientType a model of the `Ring` concept.
 */

template<typename CoefficientType>
class Simplicial_chain_complex : public Abstract_simplicial_chain_complex<CoefficientType> {
public:
    /** \brief Type of vertices coordinates */
    typedef std::vector<double> Point ;

protected:
    /** \brief Vector of vertices coordinates */
    std::vector<Point> _coords ;

private:
    /** \brief Vector of VTK types associated to cells in each dimension
     e.g. {1, 3, 5, 10} */
    static const std::vector<int> VTK_simptypes ;

public:

    /**
     * \brief Default constructor: builds an empty  simplicial complex.
     */
    Simplicial_chain_complex() {} ;
    // Constructor from a simplicial complex stored in a .simp file and read into a MeshObject
    /**
     * \brief Constructor from a Mesh_object.
     *
     * Builds a simplicial complex from a Mesh_object and saves the vector of vertices coordinates.
     */
    Simplicial_chain_complex(const Mesh_object& mesh, std::vector<Point> coords) : Abstract_simplicial_chain_complex<CoefficientType>(mesh), _coords(coords) {} ;

    /**
     * \brief Affectation operator for simplicial complexes.
     *
     * Stores a copy of a simplicial complex in *this.
     *
     * \param[in] complex The simplicial complex which will be copied.
     */
    Simplicial_chain_complex& operator= (const Simplicial_chain_complex& complex)
    {
        this->Abstract_simplicial_chain_complex<CoefficientType>::operator=(complex) ;
        _coords = complex._coords ;
        return *this ;
    }

    /** \brief Friend class `Duality_simplicial_complex_tools` computes the complementary simplicial complex for Alexander duality */
    friend Duality_simplicial_complex_tools<CoefficientType> ;

    /** \brief Get the vector of vertices coordinates  */
    const std::vector<Point>& get_vertices_coords() const
    {
        return _coords ;
    }

    /** \brief Get the coordinates of the ith dimension-0 simplex

     * \warning This does not come to return vertices indices, as dimension 0 simplices enumerate vertices in any order. For instance, if an abstract simplicial complex is build from 3 vertices {1,2,3} such that the enumeration of dimension 0 simplicies is:
     *  id0: 3, id1 : 2, id2: 1
     * then the bottom_faces of the 1-simplex {1,2} are two 0-simplices with id 2 and 1.
     */
    Point get_vertex_coords (size_t i) const
    {
        const Simplex simpl(this->_ind2simp.at(0).at(i)) ;
        const std::set<size_t> verts(simpl.get_vertices()) ;
        const size_t id(*(verts.cbegin())) ;
        return _coords.at(id);
    }

    // VTK export

    /**
     * \brief Method exporting a simplicial complex (plus, optionally, labels) to a VTK file.
     *
     * The method generates legacy text VTK files. Labels are exported as such in a VTK property, together with CellID property, containing the index of each cell.
     *
     * \tparam LabelType Type of labels provided (default: int).
     *
     * \param[in] K Simplicial complex exported.
     * \param[in] filename Output file root (output filenames will be built from this root).
     * \param[in] labels Pointer to a vector of labels in each dimension. (*labels).at(q) is the set of integer labels of cells of dimension q. If labels is NULL, only CellID property is exported.
     * \param[in] label_type_name Typename used in vtk export (e.g. "int" or "unsigned_long", see <a href = "https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html">VTK manual </a>).
     */
    template <typename LabelType = int>
    static void chain_complex_to_vtk(const Simplicial_chain_complex &K, const std::string &filename, const std::vector<std::vector<LabelType> > *labels=NULL, string label_type_name = "int")
    {
        typedef Simplicial_chain_complex<CoefficientType> ComplexType;
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
        const std::vector<ComplexType::Point>& coords(K.get_vertices_coords()) ;
        for (size_t n = 0; n < nnodes; ++n)
        {
            vector<double> p(coords.at(n)) ;
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
            // Size : size of a cell of dimension q : q+1
            for (int q=0; q<=K._dim; ++q)
            {
                ncells_tot += K.nb_cells(q) ;
                const size_t size_cell = q+1 ;
                size_cells_tot += (size_cell+1)*K.nb_cells(q) ;
            }
            out << "CELLS " << ncells_tot << " " << size_cells_tot << endl ;
            // Output cells by increasing dimension
            for (int q=0; q<=K._dim; ++q)
            {
                const size_t size_cell = q+1 ;
                for (size_t id =0; id < K.nb_cells(q); ++id)
                {
                    Simplex verts(K._ind2simp.at(q).at(id)) ;
                    out << size_cell << " " ;
                    for (typename Simplex::const_iterator it = verts.cbegin(); it != verts.cend(); ++it)
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
     * \brief Method exporting a chain over a simplicial complex to a VTK file.
     *
     * The method generates legacy text VTK files. All the cells of the chain with non zero coefficient are exported. If a cellId is provided, labels are exported in a VTK property (2 for all cells, 0 for cell of index cellId).  The index of each cell is exported in a CellID property.
     *
     * \param[in] K Simplicial complex exported.
     * \param[in] filename Output file root (output filenames will be built from this root).
     * \param[in] chain Sparse_chain exported (all the cells with non-zero coefficients in the chain are exported to vtk).
     * \param[in] q Dimension of the cells of the chain.
     * \param[in] cellId If cellID is not -1 (that is MAX_SIZE_T), labels are exported to distinguish cells of the chain (label 2) from cellId cell (label 0).
     */
    static void chain_complex_chain_to_vtk(const Simplicial_chain_complex &K, const std::string &filename, const OSM::Sparse_chain<CoefficientType, OSM::COLUMN>& chain, int q, size_t cellId = -1) ;
};

// Initialization of static VTK_simptypes
template <typename CoefficientType> const
std::vector<int> Simplicial_chain_complex<CoefficientType>::VTK_simptypes({1, 3, 5, 10});


// chain_complex_chain_to_vtk
template <typename CoefficientType>
void Simplicial_chain_complex<CoefficientType>::chain_complex_chain_to_vtk(const Simplicial_chain_complex &K, const std::string &filename, const OSM::Sparse_chain<CoefficientType, OSM::COLUMN>& chain, int q, size_t cellId)
{
    typedef Simplicial_chain_complex<CoefficientType> ComplexType ;
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
    size_t nnodes = K._coords.size() ;
    out << "POINTS " << nnodes << " double" << endl ;
    const std::vector<ComplexType::Point>& coords(K.get_vertices_coords()) ;
    for (size_t n = 0; n < nnodes; ++n)
    {
        vector<double> p(coords.at(n)) ;
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
            const size_t size_cell = q+1 ;
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
            const size_t size_cell = q+1 ;
            for (size_t id =0; id < K.nb_cells(q); ++id)
            {
                if (!chain.is_null(id))
                {
                    Simplex verts(K._ind2simp.at(q).at(id)) ;
                    out << size_cell << " " ;
                    for (typename Simplex::const_iterator it = verts.cbegin(); it != verts.cend(); ++it)
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
        for (size_t i : ids)
            out << i << " " ;
        out << endl ;
    }
    out.close() ;
}

} /* end namespace HDVF */
} /* end namespace CGAL */

#endif // CGAL_HDVF_SIMPLICIAL_CHAIN_COMPLEX_H
