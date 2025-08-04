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

#ifndef CGAL_HDVF_GEOMETRIC_CHAIN_COMPLEX_TOOLS_H
#define CGAL_HDVF_GEOMETRIC_CHAIN_COMPLEX_TOOLS_H

#include <CGAL/license/HDVF.h>

#include <vector>
#include <map>
#include <stdexcept>
#include <unordered_set>
#include <CGAL/HDVF/Simplicial_chain_complex.h>
#include <CGAL/HDVF/Cubical_chain_complex.h>
#include <CGAL/HDVF/Hdvf_core.h>
#include <CGAL/HDVF/Hdvf_persistence.h>
#include <CGAL/HDVF/Hdvf_duality.h>
#include <CGAL/HDVF/Sub_chain_complex_mask.h>
#include <CGAL/HDVF/Mesh_object_io.h>
#include <CGAL/HDVF/Cub_object_io.h>
#include <CGAL/HDVF/Tet_object_io.h>
#include <CGAL/HDVF/Icosphere_object_io.h>
#include <CGAL/OSM/OSM.h>

/**
 * \class SimpComplexTools
 * \brief Provides tools for SimpComplexes ((co)homology, duality, persistent homology).
 * \brief Friend class of SimpComplex
 */

namespace CGAL {
namespace HDVF {

// GeometricChainComplex IO

// Hdvf vtk export

/*!
 * \brief Export all the `Hdvf` information of a geometric chain complex to vtk files.
 *
 * Export PSC labels and homology/cohomology generators (depending on HDVF options) associated to each critical cell to vtk files.
 *
 * \param[in] hdvf Reference to the HDVF exported.
 * \param[in] complex Underlying geometric chain complex.
 * \param[in] filename Prefix of all generated files.
 * \param[in] co_faces Export the cohomology generator or its co-faces (sometimes more convenient for visualisation).
 *
 * Below, a sample mesh with, (left) homology generators, (right) two examples of cohomology generators (corresponding generators/co-generators bear similar colours):
 *
 * <img src="HDVF_dtorus_homs.png" align="center" width=25%/>
 * <img src="HDVF_dtorus_cohom1.png" align="center" width=25%/>
 * <img src="HDVF_dtorus_cohom2.png" align="center" width=25%/>
 *
 * The same generators displayed through their co-faces:
 *
 * <img src="HDVF_dtorus_cohom1_co.png" align="center" width=25%/>
 * <img src="HDVF_dtorus_cohom2_co.png" align="center" width=25%/>
 *
 * All homology / cohomology generators:
 *
 *<img src="HDVF_dtorus_all.png" align="center" width=30%/>
 */

template <typename CoefType, typename ComplexType, template <typename, int> typename _ChainType = OSM::Sparse_chain, template <typename, int> typename _SparseMatrixType = OSM::Sparse_matrix, typename VertexIdType = size_t>
void hdvf_geometric_chain_complex_output_vtk (Hdvf_core<CoefType, ComplexType, _ChainType, _SparseMatrixType> &hdvf, ComplexType &complex, string filename = "test", bool co_faces = false)
{
    typedef Hdvf_core<CoefType, ComplexType, _ChainType, _SparseMatrixType> HDVF_type;
    // Export PSC labelling
    string outfile(filename+"_PSC.vtk") ;
    vector<vector<int> > labels = hdvf.get_psc_labels() ;
    ComplexType::chain_complex_to_vtk(complex, outfile, &labels) ;

    if (hdvf.get_hdvf_opts() != OPT_BND)
    {
        // Export generators of all critical cells
        vector<vector<size_t> > criticals(hdvf.get_flag(CRITICAL)) ;
        for (int q = 0; q <= complex.dim(); ++q)
        {
            for (size_t c : criticals.at(q))
            {
                // Homology generators
                if (hdvf.get_hdvf_opts() & (OPT_FULL | OPT_G))
                {
                    string outfile_g(filename+"_G_"+to_string(c)+"_dim_"+to_string(q)+".vtk") ;
                    //                    vector<vector<size_t> > labels = hdvf.export_label(G,c,q) ;
                    OSM::Sparse_chain<CoefType,OSM::COLUMN> chain(hdvf.get_homology_chain(c,q)) ;
                    ComplexType::chain_complex_chain_to_vtk(complex, outfile_g, chain, q, c) ;
                }
                // Cohomology generators
                if (hdvf.get_hdvf_opts() & (OPT_FULL | OPT_F))
                {
                    string outfile_f(filename+"_FSTAR_"+to_string(c)+"_dim_"+to_string(q)+".vtk") ;
                    OSM::Sparse_chain<CoefType,OSM::COLUMN> chain(hdvf.get_cohomology_chain(c, q)) ;
                    if (!co_faces)
                    {
                        ComplexType::chain_complex_chain_to_vtk(complex, outfile_f, chain, q, c) ;
                    }
                    else
                    {
                        if (q < complex.dim())
                        {
                            ComplexType::chain_complex_chain_to_vtk(complex, outfile_f, complex.cofaces_chain(chain, q), q+1, c) ;
                        }
                    }
                }
            }
        }
    }
}

// Hdvf_persistence vtk export

/** \brief Export all the `HDVF_persistence` information of a geometric chain complex to vtk files.
 *
 * Export PSC labels and homology/cohomology generators (depending on HDVF options) associated to each persistent intervals to vtk files.
 *
 * \param[in] per_hdvf Reference to the persistent HDVF exported.
 * \param[in] complex Underlying geometric chain complex.
 * \param[in] filename Prefix of all generated files.
 * \param[in] co_faces Export the cohomology generator or its co-faces (sometimes more convenient for visualisation).
 */

template <typename CoefType, typename ComplexType, typename DegType, typename FiltrationType>
void hdvf_persistence_geometric_chain_complex_output_vtk (Hdvf_persistence<CoefType, ComplexType, DegType, FiltrationType> &per_hdvf, ComplexType &complex, string filename = "per", bool co_faces = false)
{
    if (!per_hdvf.with_export())
        throw("Cannot export persistent generators to vtk: with_export is off!") ;

    using perHDVFType = Hdvf_persistence<CoefType, ComplexType, DegType, FiltrationType>;
    using PerHole = PerHoleT<DegType> ;

    // Export the filtration
    string out_file_filtration = filename+"_filtration.vtk" ;
    vector<vector<size_t> > filtr_labels = per_hdvf.get_filtration().export_filtration();
    ComplexType::template chain_complex_to_vtk<size_t>(complex, out_file_filtration, &filtr_labels, "unsigned_long") ;

    // Iterate over persistence diagram (iterator over non zero intervals)
    // Batch informations are stored in file filename_infos.txt
    std::ofstream info_file(filename+"_infos.txt") ;
    size_t i = 0 ;
    for (typename perHDVFType::iterator it = per_hdvf.begin(); it != per_hdvf.end(); ++it)
    {
        typename perHDVFType::PerIntervalInformation hole_data(*it) ;
        const PerHole hole(hole_data.hole) ;
        // Export informations of this hole
        info_file << i << " -> " << " --- duration : " << per_hdvf.hole_duration(hole) << " -- " << hole << std::endl ;

        if (per_hdvf.hole_duration(hole)>0) // Export vtk of "finite" holes
        {
            // Build name associated to the ith hole : filename_i
            string out_file = filename+"_"+to_string(i) ;
            // Export PSC labels to vtk
            ComplexType::chain_complex_to_vtk(complex, out_file+"_PSC.vtk", &hole_data.labelsPSC) ;
            const CellsPerInterval per_int_cells(std::get<1>(hole_data.hole)) ;
            // Export homology generators (g)
            if (per_hdvf.get_hdvf_opts()  & (OPT_FULL | OPT_G))
            {
                // First generator : filename_i_g_sigma_q.vtk
                {
                    const size_t id(per_int_cells.first.first) ;
                    const int dim(per_int_cells.first.second) ;
                    string tmp(out_file+"_g_"+to_string(id)+"_"+to_string(dim)+".vtk") ;
                    ComplexType::chain_complex_chain_to_vtk(complex, tmp, hole_data.g_chain_sigma, dim);
                }
                // Second generator : filename_i_g_tau_q+1.vtk
                {
                    const size_t id(per_int_cells.second.first) ;
                    const int dim(per_int_cells.second.second) ;
                    string tmp(out_file+"_g_"+to_string(id)+"_"+to_string(dim)+".vtk") ;
                    ComplexType::chain_complex_chain_to_vtk(complex, tmp, hole_data.g_chain_tau, dim);
                }
            }
            // Export cohomology generators (fstar)
            if ((per_hdvf.get_hdvf_opts() == OPT_FULL) || (per_hdvf.get_hdvf_opts() == OPT_F))
            {
                // First generator : filename_i_fstar_sigma_q.vtk
                {
                    const size_t id(per_int_cells.first.first) ;
                    const int dim(per_int_cells.first.second) ;
                    string tmp(out_file+"_fstar_"+to_string(id)+"_"+to_string(dim)+".vtk") ;
                    if (co_faces)
                    {
                        ComplexType::chain_complex_chain_to_vtk(complex, tmp, complex.cofaces_chain(hole_data.fstar_chain_sigma, dim), dim+1);
                    }
                    else
                        ComplexType::chain_complex_chain_to_vtk(complex, tmp, hole_data.fstar_chain_sigma, dim);
                }
                // Second generator : filename_i_fstar_tau_q+1.vtk
                {
                    const size_t id(per_int_cells.second.first) ;
                    const int dim(per_int_cells.second.second) ;
                    string tmp(out_file+"_fstar_"+to_string(id)+"_"+to_string(dim)+".vtk") ;
                    if (co_faces)
                    {
                        ComplexType::chain_complex_chain_to_vtk(complex, tmp, complex.cofaces_chain(hole_data.fstar_chain_tau, dim), dim+1);
                    }
                    else
                        ComplexType::chain_complex_chain_to_vtk(complex, tmp, hole_data.fstar_chain_tau, dim);
                }
            }
        }
        ++i ;
    }
    info_file.close() ;
    hdvf_geometric_chain_complex_output_vtk<CoefType, ComplexType>(per_hdvf, complex, filename+"_inf", co_faces) ;
}

// Hdvf_duality vtk export

/** \brief Export all the `HDVF_duality` information of a geometric chain complex to vtk files.

 Export PSC labels and homology/cohomology generators (depending on HDVF options) associated to each persistent intervals to vtk files.

 \param[in] hdvf Reference to the HDVF exported.
 \param[in] complex Underlying geometric chain complex.
 \param[in] filename Prefix of all generated files.
 \param[in] co_faces Export the cohomology generator or its co-faces (sometimes more convenient for visualisation).
 */

template <typename CoefType, typename ComplexType, typename VertexIdType = size_t>
void hdvf_duality_geometric_chain_complex_output_vtk (Hdvf_duality<CoefType, ComplexType> &hdvf, ComplexType &complex, string filename = "test", bool co_faces = false)
{
    typedef Hdvf_duality<CoefType, ComplexType> HDVF_type;
    // Export PSC labelling
    string outfile(filename+"_PSC.vtk") ;
    vector<vector<int> > labels = hdvf.get_psc_labels() ;
    ComplexType::chain_complex_to_vtk(complex, outfile, &labels) ;

    if (hdvf.get_hdvf_opts() != OPT_BND)
    {
        // Export generators of all critical cells
        vector<vector<size_t> > criticals(hdvf.get_flag(CRITICAL)) ;
        for (int q = 0; q <= complex.dim(); ++q)
        {
            for (size_t c : criticals.at(q))
            {
                // Homology generators
                if (hdvf.get_hdvf_opts() & (OPT_FULL | OPT_G))
                {
                    string outfile_g(filename+"_G_"+to_string(c)+"_dim_"+to_string(q)+".vtk") ;
                    //                    vector<vector<size_t> > labels = hdvf.export_label(G,c,q) ;
                    OSM::Sparse_chain<CoefType,OSM::COLUMN> chain(hdvf.get_homology_chain(c,q)) ;
                    ComplexType::chain_complex_chain_to_vtk(complex, outfile_g, chain, q, c) ;
                }
                // Cohomology generators
                if (hdvf.get_hdvf_opts() & (OPT_FULL | OPT_F))
                {
                    string outfile_f(filename+"_FSTAR_"+to_string(c)+"_dim_"+to_string(q)+".vtk") ;
                    OSM::Sparse_chain<CoefType,OSM::COLUMN> chain(hdvf.get_cohomology_chain(c, q)) ;
                    if (!co_faces)
                    {
                        ComplexType::chain_complex_chain_to_vtk(complex, outfile_f, chain, q, c) ;
                    }
                    else // Compute co-faces
                    {
                        if (q < complex.dim())
                        {
                            // Restrict the cofaces of the cohomology generator to the current sub chain complex
                            OSM::Sparse_chain<CoefType,OSM::COLUMN> cofaces_chain(complex.cofaces_chain(chain, q)) ;
                            Sub_chain_complex_mask<CoefType,ComplexType> sub(hdvf.get_current_mask());
                            sub.screen_chain(cofaces_chain, q+1);
                            // Display
                            ComplexType::chain_complex_chain_to_vtk(complex, outfile_f, cofaces_chain, q+1, c) ;
                        }
                    }
                }
            }
        }
    }
}

// Duality tools

// =========== Simplicial chain complex tools
// Build K sub-chain complex of L (homeomorphic to B^n)

/*!
 \ingroup PkgHDVFAlgorithmClasses

The class `Duality_simplicial_complex_tools` is dedicated to Alexander duality for 3D surface meshes. Starting from a simplicial chain complex (encoding a 3D surface mesh), it provides methods to embed the complex into a larger icosphere and generate a 3D constrained Delaunay triangulation.

Technically, starting from a Simplicial_chain_complex `_K`, the method `simp_complex_bb` builds a Simplicial_chain_complex L and a Sub_chain_complex_mask K.
- L : complex built out of _K together with a closing icosphere, meshed by tetgen (constrained Delaunay triangulation)
- K (Sub_chain_complex_mask) : Sub_chain_complex_mask identifying _K inside L

\tparam CoefficientType a model of the `Ring` concept providing the ring used to compute homology.
*/

template<typename CoefficientType>
class Duality_simplicial_complex_tools {
public:
    /** \brief Type of simplicial chain complex encoding `L`. */
    typedef Simplicial_chain_complex<CoefficientType> ComplexType ;
    /** \brief Type of sub chain complex mask encoding the sub complex `K`. */
    typedef Sub_chain_complex_mask<CoefficientType, ComplexType> SubCCType ;
    /** \brief Default constructor. */
    Duality_simplicial_complex_tools() {}

    /** \brief Type returned by `simp_complex_bb`.
     *
     * The structure contains a triple:
     * - A simplicial chain complex `L` (homeomorphic to \f$\mathbb B^3\f$).
     * - A sub chain complex mask `K` (encoding the initial mesh)
     * - The vector of vertices coordinates of `L`
     */
    typedef struct {
        ComplexType& L ;
        SubCCType& K ;
        std::vector<IONodeType> nodes ;
    } TripleRes ;

    /** \brief Generate a subcomplex \f$K\f$K and a complex \f$L\f$ with \f$K\subseteq L\f$ from a simplicial complex `_K`.
     *
     * `_K` is embedded into a larger icosphere and a 3D constrained Delaunay triangulation is generated. Then \f$K\f$, \f$L\f$ and vertices coordinates are extracted and stored in a `TripleRes` structure.
     *
     * The following figures shows the resulting complex with an initial `Twirl` mesh (right - sectional view):
     * <img src="HDVF_twirl_view1.png" align="center" width=35%/>
     * <img src="HDVF_twirl_view2.png" align="center" width=30%/>
     *
     * \param[in] _K Simplicial chain complex (working mesh).
     * \param[in] BB_ratio Ratio of the "closing" icosphere diameter with respect to the diameter of the object's bounding box.
     * \param[in] out_file_prefix Prefix of tetgen intermediate files (default: "file_K_closed.off").
     */
    static TripleRes simplicial_chain_complex_bb (const ComplexType& _K, double BB_ratio=1.5, const string& out_file_prefix = "file_K_closed.off")
    {
        // Export _K to a MeshObject to add the icosphere and mesh with tetGen
        Mesh_object_io mesh_L = Duality_simplicial_complex_tools::export_meshObject(_K) ;

        // Closing K by adding the icosphere
        //  Compute a bounding icosphere
        IONodeType bary = mesh_L.barycenter() ;
        double r = mesh_L.radius(bary) ;
        Icosphere_object_io ico(2,bary, BB_ratio*r) ;
        ico.print_infos() ;

        // Add it to the mesh
        mesh_L.push_back(ico) ;
        mesh_L.print_infos() ;

        // Write this mesh to an off file for Tetgen
        mesh_L.write_off(out_file_prefix) ;
        //        const std::string tetgen_path(CMAKE_TETGEN_PATH) ;
        // WARNING: use CGAL 3D triangulation to get rid of tetgen...
        const std::string tetgen_path("/Users/umenohana/Dropbox/G-Mod/TopAlg/code/tetgen1.6.0/build/") ;
        std::string tetgen_command = tetgen_path+"tetgen -pqkcYfe "+out_file_prefix+".off" ;
        system(tetgen_command.c_str()) ;

        // Read the mesh built by tetgen for L
        Tet_object_io tetL(out_file_prefix) ;

        // Build the associated SimpComplex
        ComplexType& L = *new ComplexType(tetL, tetL.get_nodes()) ;

        // Build the Sub_chain_complex_mask encoding _K inside L
        SubCCType& K(*new SubCCType(L, false)) ;
        // Visit all cells of _K and activate the corresponding bit in K
        for (int q=0; q<=_K.dim(); ++q)
        {
            for (size_t i=0; i<_K.nb_cells(q); ++i)
            {
                const std::set<size_t>& simplex(_K._ind2simp.at(q).at(i).get_vertices()) ;
                const size_t id(L._simp2ind.at(q)[simplex]) ;
                K.set_bit_on(q, id) ;
            }
        }
        TripleRes t = {L,K,tetL.nodes} ;
        return t ;
    }

    /** \brief Export a SimpComplex to a MeshObject  */
    static Mesh_object_io& export_meshObject(const ComplexType& _CC)
    {
        std::vector<IOCellType> vcells ;
        for (int q = 0; q <= _CC.dim(); ++q)
            for (size_t i = 0; i<_CC.nb_cells(q); ++i)
            {
                const Simplex& s(_CC._ind2simp.at(q).at(i)) ;
                vcells.push_back(s.get_vertices()) ;
            }

        Mesh_object_io &m = *(new Mesh_object_io(-3, _CC.get_vertices_coords(), vcells)) ;
        return m ;
    }
} ;

// =========== Cubical chain complex tools
// Build K sub-chain complex of L (homeomorphic to B^n)
// Adjust L size

/*!
 \ingroup PkgHDVFAlgorithmClasses

The class `Duality_cubical_complex_tools` is dedicated to Alexander duality for 3D binary volumes.

Starting from a Cubical_chain_complex `_K`, the method `cubical_chain_complex_BB` builds a Cubical_chain_complex L and Sub_chain_complex_mask K.
- L : complex built of the "full" bounding box of _K
- K (Sub_chain_complex_mask) : Sub_chain_complex_mask identifying _K inside L

Use the `frame` method from the `Cub_object_io` class to enlarge the bounding box (via a 1 pixel dilatation) if necessary.

\tparam CoefficientType a model of the `Ring` concept providing the ring used to compute homology.
*/

template<typename CoefficientType>
class Duality_cubical_complex_tools {
public:
    /** \brief Type of cubical complexes used for the initial complex and \f$L\f$. */
    typedef Cubical_chain_complex<CoefficientType> ComplexType ;
    /** \brief Type of sub chain complex mask used to encode the sub-complex  \f$K\f$. */
    typedef Sub_chain_complex_mask<CoefficientType, ComplexType> SubCCType ;
    // Constructor
    Duality_cubical_complex_tools() {}


    /** \brief Generate a subcomplex \f$K\f$K and a complex \f$L\f$ with \f$K\subseteq L\f$ from a cubical complex `_K`.
     *
     * `L` is the bounding box of `_K` (homeomorphic to a ball) and \f$K\f$ is a sub chain complex mask encoding `_K`.
     *
     * The following figures shows the resulting complex with an initial simple cubical complex (right - sectional view):
     * <img src="HDVF_eight_view1.png" align="center" width=35%/>
     * <img src="HDVF_eight_view2.png" align="center" width=30%/>
     *
     * \param[in] _K Initial cubical chain complex (working mesh).
     */
    static std::pair<ComplexType&, SubCCType&> cubical_chain_complex_bb (const ComplexType& _K)
    {
        Cub_object_io tmp ;
        tmp.dim = _K.dim() ;
        tmp.N = _K._size_bb ;
        tmp.ncubs.resize(tmp.dim+1) ;
        // Visit all boolean indices in the BB of _CC and insert corresponding Cells
        for (size_t i=0; i<_K._P.at(_K.dim()); ++i)
        {
            const std::vector<size_t> tmpkhal(_K.ind2khal(i)) ;
            const int dtmp(_K.calculate_dimension(tmpkhal)) ;
            (tmp.ncubs)[dtmp] += 1 ;
            tmp.cubs.push_back(tmpkhal) ;
        }
        ComplexType& L(*new ComplexType(tmp, ComplexType::PRIMAL)) ;

        // Build the Sub_chain_complex_mask corresponding to _CC
        SubCCType& K(*new SubCCType(L, false)) ;
        // Visit all cells of _CC and activate the corresponding bit in K
        for (int q=0; q<=_K.dim(); ++q)
        {
            for (size_t i=0; i<_K.nb_cells(q); ++i)
            {
                const std::vector<size_t> khal(_K.ind2khal(_K._base2bool.at(q).at(i))) ;
                const size_t j = L._bool2base.at(q).at(L.khal2ind(khal)) ;
                K.set_bit_on(q,j) ;
            }
        }
        return std::pair<ComplexType&, SubCCType&>(L,K) ;
    }
} ;

} /* end namespace HDVF */
} /* end namespace CGAL */

#endif // CGAL_HDVF_GEOMETRIC_CHAIN_COMPLEX_TOOLS_H
