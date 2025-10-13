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
namespace IO {

// GeometricChainComplex IO

// Hdvf vtk export

/*!
 * \brief Exports all the `Hdvf` information of a geometric chain complex to vtk files.
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

template <typename ChainComplex, template <typename, int> typename ChainType = OSM::Sparse_chain, template <typename, int> typename SparseMatrixType = OSM::Sparse_matrix, typename VertexIdType = size_t>
void write_VTK (Homological_discrete_vector_field::Hdvf_core<ChainComplex, ChainType, SparseMatrixType>& hdvf, ChainComplex &complex, std::string filename = "test", bool co_faces = false)
{
    typedef typename ChainComplex::Coefficient_ring Coefficient_ring;
    typedef Homological_discrete_vector_field::Hdvf_core<ChainComplex, ChainType, SparseMatrixType> HDVF_parent;
    // Export PSC labelling
    std::string outfile(filename+"_PSC.vtk") ;
    std::vector<std::vector<int> > labels = hdvf.psc_labels() ;
    ChainComplex::chain_complex_to_vtk(complex, outfile, &labels) ;

    if (hdvf.hdvf_opts() != Homological_discrete_vector_field::OPT_BND)
    {
        // Export generators of all critical cells
        std::vector<std::vector<size_t> > criticals(hdvf.psc_flags(Homological_discrete_vector_field::CRITICAL)) ;
        for (int q = 0; q <= complex.dimension(); ++q)
        {
            for (size_t c : criticals.at(q))
            {
                // Homology generators
                if (hdvf.hdvf_opts() & (Homological_discrete_vector_field::OPT_FULL | Homological_discrete_vector_field::OPT_G))
                {
                    std::string outfile_g(filename+"_hom_"+std::to_string(c)+"_dim_"+std::to_string(q)+".vtk") ;
                    //                    std::vector<std::vector<size_t> > labels = hdvf.export_label(G,c,q) ;
                    OSM::Sparse_chain<Coefficient_ring,OSM::COLUMN> chain(hdvf.homology_chain(c,q)) ;
                    ChainComplex::chain_to_vtk(complex, outfile_g, chain, q, c) ;
                }
                // Cohomology generators
                if (hdvf.hdvf_opts() & (Homological_discrete_vector_field::OPT_FULL | Homological_discrete_vector_field::OPT_F))
                {
                    std::string outfile_f(filename+"_cohom_"+std::to_string(c)+"_dim_"+std::to_string(q)+".vtk") ;
                    OSM::Sparse_chain<Coefficient_ring,OSM::COLUMN> chain(hdvf.cohomology_chain(c, q)) ;
                    if (!co_faces)
                    {
                        ChainComplex::chain_to_vtk(complex, outfile_f, chain, q, c) ;
                    }
                    else
                    {
                        if (q < complex.dimension())
                        {
                            ChainComplex::chain_to_vtk(complex, outfile_f, complex.cofaces_chain(chain, q), q+1, c) ;
                        }
                    }
                }
            }
        }
    }
}

// Hdvf_persistence vtk export

/** \brief Exports all the `HDVF_persistence` information of a geometric chain complex to vtk files.
 *
 * Export PSC labels and homology/cohomology generators (depending on HDVF options) associated to each persistent intervals to vtk files.
 *
 * \param[in] per_hdvf Reference to the persistent HDVF exported.
 * \param[in] complex Underlying geometric chain complex.
 * \param[in] filename Prefix of all generated files.
 * \param[in] co_faces Export the cohomology generator or its co-faces (sometimes more convenient for visualisation).
 */

template <typename ChainComplex, typename Degree, typename FiltrationType>
void write_VTK (Homological_discrete_vector_field::Hdvf_persistence<ChainComplex, Degree, FiltrationType> &per_hdvf, ChainComplex &complex, std::string filename = "per", bool co_faces = false)
{
    typedef typename ChainComplex::Coefficient_ring Coefficient_ring;
    if (!per_hdvf.with_export())
        throw("Cannot export persistent generators to vtk: with_export is off!") ;

    using HDVF_type = Homological_discrete_vector_field::Hdvf_persistence<ChainComplex, Degree, FiltrationType>;
    using Persistent_interval = HDVF_type::Persistent_interval;

    // Export the filtration
    std::string out_file_filtration = filename+"_filtration.vtk" ;
    std::vector<std::vector<size_t> > filtr_labels = per_hdvf.get_filtration().export_filtration();
    ChainComplex::template chain_complex_to_vtk<size_t>(complex, out_file_filtration, &filtr_labels, "unsigned_long") ;

    // Iterate over persistence diagram (iterator over non zero intervals)
    // Batch informations are stored in file filename_infos.txt
    std::ofstream info_file(filename+"_infos.txt") ;
    size_t i = 0 ;
    for (typename HDVF_type::iterator it = per_hdvf.begin(); it != per_hdvf.end(); ++it)
    {
        typename HDVF_type::Persistent_hole hole_data(*it) ;
        const Persistent_interval hole(hole_data.hole) ;
        // Export informations of this hole
        /* V1 */
        info_file << i << " -> " << " --- duration : " << hole.duration() << " -- " ;
        hole.insert(info_file);
        info_file << std::endl ;
        /* V2 */
        //        info_file << i << " -> " << " --- duration : " << hole.duration() << " -- " << hole << std::endl ;

        if (hole.duration()>0) // Export vtk of "finite" holes
        {
            // Build name associated to the ith hole : filename_i
            std::string out_file = filename+"_"+std::to_string(i) ;
            // Export PSC labels to vtk
            ChainComplex::chain_complex_to_vtk(complex, out_file+"_PSC.vtk", &hole_data.labelsPSC) ;

            // Export homology generators (g)
            if (per_hdvf.hdvf_opts()  & (Homological_discrete_vector_field::OPT_FULL | Homological_discrete_vector_field::OPT_G))
            {
                // First generator : filename_i_g_sigma_q.vtk
                {
                    const size_t id(hole_data.hole.cell_birth.first) ;
                    const int dim(hole_data.hole.cell_birth.second) ;
                    std::string tmp(out_file+"_hom_"+std::to_string(id)+"_"+std::to_string(dim)+".vtk") ;
                    ChainComplex::chain_to_vtk(complex, tmp, hole_data.homology_chain_birth, dim);
                }
                // Second generator : filename_i_g_tau_q+1.vtk
                {
                    const size_t id(hole_data.hole.cell_death.first) ;
                    const int dim(hole_data.hole.cell_death.second) ;
                    std::string tmp(out_file+"_hom_"+std::to_string(id)+"_"+std::to_string(dim)+".vtk") ;
                    ChainComplex::chain_to_vtk(complex, tmp, hole_data.homology_chain_death, dim);
                }
            }
            // Export cohomology generators (fstar)
            if ((per_hdvf.hdvf_opts() == Homological_discrete_vector_field::OPT_FULL) || (per_hdvf.hdvf_opts() == Homological_discrete_vector_field::OPT_F))
            {
                // First generator : filename_i_fstar_sigma_q.vtk
                {
                    const size_t id(hole_data.hole.cell_birth.first) ;
                    const int dim(hole_data.hole.cell_birth.second) ;
                    std::string tmp(out_file+"_cohom_"+std::to_string(id)+"_"+std::to_string(dim)+".vtk") ;
                    if (co_faces)
                    {
                        ChainComplex::chain_to_vtk(complex, tmp, complex.cofaces_chain(hole_data.cohomology_chain_birth, dim), dim+1);
                    }
                    else
                        ChainComplex::chain_to_vtk(complex, tmp, hole_data.cohomology_chain_birth, dim);
                }
                // Second generator : filename_i_fstar_tau_q+1.vtk
                {
                    const size_t id(hole_data.hole.cell_death.first) ;
                    const int dim(hole_data.hole.cell_death.second) ;
                    std::string tmp(out_file+"_cohom_"+std::to_string(id)+"_"+std::to_string(dim)+".vtk") ;
                    if (co_faces)
                    {
                        ChainComplex::chain_to_vtk(complex, tmp, complex.cofaces_chain(hole_data.cohomology_chain_death, dim), dim+1);
                    }
                    else
                        ChainComplex::chain_to_vtk(complex, tmp, hole_data.cohomology_chain_death, dim);
                }
            }
        }
        ++i ;
    }
    info_file.close() ;
    // Export infinite holes by calling Hdvf_core write_VTK
    (static_cast<void (*) (Homological_discrete_vector_field::Hdvf_core<ChainComplex, CGAL::OSM::Sparse_chain,CGAL::OSM::Sub_sparse_matrix>&, ChainComplex&, std::string, bool)>(CGAL::IO::write_VTK<ChainComplex,CGAL::OSM::Sparse_chain,CGAL::OSM::Sub_sparse_matrix, size_t>))(per_hdvf, complex, filename+"_inf", co_faces);
}

// Hdvf_duality vtk export

/** \brief Exports all the `HDVF_duality` information of a geometric chain complex to vtk files.

 Export PSC labels and homology/cohomology generators (depending on HDVF options) associated to each persistent intervals to vtk files.

 \param[in] hdvf Reference to the HDVF exported.
 \param[in] complex Underlying geometric chain complex.
 \param[in] filename Prefix of all generated files.
 \param[in] co_faces Export the cohomology generator or its co-faces (sometimes more convenient for visualisation).
 */

template <typename ChainComplex, typename VertexIdType = size_t>
void write_VTK (Homological_discrete_vector_field::Hdvf_duality<ChainComplex> &hdvf, ChainComplex &complex, std::string filename = "test", bool co_faces = false)
{
    typedef typename ChainComplex::Coefficient_ring Coefficient_ring;
    typedef Homological_discrete_vector_field::Hdvf_duality<ChainComplex> HDVF_parent;
    // Export PSC labelling
    std::string outfile(filename+"_PSC.vtk") ;
    std::vector<std::vector<int> > labels = hdvf.psc_labels() ;
    ChainComplex::chain_complex_to_vtk(complex, outfile, &labels) ;

    if (hdvf.hdvf_opts() != Homological_discrete_vector_field::OPT_BND)
    {
        // Export generators of all critical cells
        std::vector<std::vector<size_t> > criticals(hdvf.psc_flags(Homological_discrete_vector_field::CRITICAL)) ;
        for (int q = 0; q <= complex.dimension(); ++q)
        {
            for (size_t c : criticals.at(q))
            {
                // Homology generators
                if (hdvf.hdvf_opts() & (Homological_discrete_vector_field::OPT_FULL | Homological_discrete_vector_field::OPT_G))
                {
                    std::string outfile_g(filename+"_hom_"+std::to_string(c)+"_dim_"+std::to_string(q)+".vtk") ;
                    //                    std::vector<std::vector<size_t> > labels = hdvf.export_label(G,c,q) ;
                    OSM::Sparse_chain<Coefficient_ring,OSM::COLUMN> chain(hdvf.homology_chain(c,q)) ;
                    ChainComplex::chain_to_vtk(complex, outfile_g, chain, q, c) ;
                }
                // Cohomology generators
                if (hdvf.hdvf_opts() & (Homological_discrete_vector_field::OPT_FULL | Homological_discrete_vector_field::OPT_F))
                {
                    std::string outfile_f(filename+"_cohom_"+std::to_string(c)+"_dim_"+std::to_string(q)+".vtk") ;
                    OSM::Sparse_chain<Coefficient_ring,OSM::COLUMN> chain(hdvf.cohomology_chain(c, q)) ;
                    if (!co_faces)
                    {
                        ChainComplex::chain_to_vtk(complex, outfile_f, chain, q, c) ;
                    }
                    else // Compute co-faces
                    {
                        if (q < complex.dimension())
                        {
                            // Restrict the cofaces of the cohomology generator to the current sub chain complex
                            OSM::Sparse_chain<Coefficient_ring,OSM::COLUMN> cofaces_chain(complex.cofaces_chain(chain, q)) ;
                            Homological_discrete_vector_field::Sub_chain_complex_mask<ChainComplex> sub(hdvf.get_current_mask());
                            sub.screen_chain(cofaces_chain, q+1);
                            // Display
                            ChainComplex::chain_to_vtk(complex, outfile_f, cofaces_chain, q+1, c) ;
                        }
                    }
                }
            }
        }
    }
}
} /* end namespace IO */

namespace Homological_discrete_vector_field {

// Duality tools

// =========== Simplicial chain complex tools
// Build K sub-chain complex of L (homeomorphic to B^n)

/*!
 \ingroup PkgHDVFAlgorithmClasses

The class `Duality_simplicial_complex_tools` is dedicated to Alexander duality for 3D surface meshes. Starting from a simplicial chain complex (encoding a 3D surface mesh), it provides methods to embed the complex into a larger icosphere and generate a 3D constrained Delaunay triangulation.

Technically, starting from a Simplicial_chain_complex `_K`, the method `simplicial_chain_complex_bb()` builds a `Simplicial_chain_complex`  `L` and a `Sub_chain_complex_mask`  `K`.
- L : complex built out of _K together with a closing icosphere, meshed by tetgen (constrained Delaunay triangulation)
- K (Sub_chain_complex_mask) : Sub_chain_complex_mask identifying _K inside L

\tparam CoefficientRing a model of the `IntegralDomainWithoutDivision` concept providing the ring used to compute homology.

\tparam Traits a geometric traits class model of the `HDVFTraits` concept.
*/

template<typename CoefficientRing, typename Traits>
class Duality_simplicial_complex_tools {
public:
    typedef typename Traits::Point Point ;
    /** \brief Type of simplicial chain complex encoding `L`. */
    typedef Simplicial_chain_complex<CoefficientRing,Traits> Chain_complex ;
    /** \brief Type of sub chain complex mask encoding the sub complex `K`. */
    typedef Sub_chain_complex_mask<Chain_complex> Sub_chain_complex ;
    /** \brief Default constructor. */
    Duality_simplicial_complex_tools() {}

    /** \brief Type returned by `simplicial_chain_complex_bb()`.
     *
     * The structure contains a triple:
     * - A simplicial chain complex `L` (homeomorphic to \f$\mathbb B^3\f$).
     * - A sub chain complex mask `K` (encoding the initial mesh)
     * - The vector of vertex coordinates of `L`
     */
    typedef struct {
        Chain_complex& L ;
        Sub_chain_complex& K ;
        std::vector<Point> nodes ;
    } Complex_duality_data ;

    /** \brief Generates a subcomplex \f$K\f$K and a complex \f$L\f$ with \f$K\subseteq L\f$ from a simplicial complex `_K`.
     *
     * `_K` is embedded into a larger icosphere and a 3D constrained Delaunay triangulation is generated. Then \f$K\f$, \f$L\f$ and vertex coordinates are extracted and stored in a `Complex_duality_data` structure.
     *
     * The following figures shows the resulting complex with an initial `Twirl` mesh (right - sectional view):
     * <img src="HDVF_twirl_view1.png" align="center" width=35%/>
     * <img src="HDVF_twirl_view2.png" align="center" width=30%/>
     *
     * \param[in] _K Simplicial chain complex (working mesh).
     * \param[in] BB_ratio Ratio of the "closing" icosphere diameter with respect to the diameter of the object's bounding box.
     * \param[in] out_file_prefix Prefix of tetgen intermediate files (default: "file_K_closed.off").
     */
    static Complex_duality_data simplicial_chain_complex_bb (const Chain_complex& _K, double BB_ratio=1.5, const std::string& out_file_prefix = "file_K_closed.off")
    {

        std::cerr << "-- Starting simplicial_chain_complex_bb" << std::endl;
        std::cerr << "Imported mesh" << std::endl;
        std::cout << _K;
        // Export _K to a MeshObject to add the icosphere and mesh with tetGen
        Mesh_object_io mesh_L = Duality_simplicial_complex_tools::export_meshObject(_K) ;
        std::cerr << "Mesh_object_io from mesh" << std::endl;
        mesh_L.print_infos();

        // Closing K by adding the icosphere
        //  Compute a bounding icosphere
        Point center = mesh_L.centroid() ;
        double r = mesh_L.radius(center) ;
        Icosphere_object_io<Traits> ico(2,center, BB_ratio*r) ;
        std::cerr << "Icosphere generated" << std::endl;
        ico.print_infos() ;

        // Add it to the mesh
        mesh_L.push_back(ico) ;
        std::cerr << "Mesh concatenation" << std::endl;
        mesh_L.print_infos() ;

        // Write this mesh to an off file for Tetgen
        mesh_L.write_off(out_file_prefix) ;
        //        const std::string tetgen_path(CMAKE_TETGEN_PATH) ;
        // WARNING: use CGAL 3D triangulation to get rid of tetgen...
        const std::string tetgen_path("/Users/umenohana/Dropbox/G-Mod/TopAlg/code/tetgen1.6.0/build/") ;
        std::string tetgen_command = tetgen_path+"tetgen -pqkcYfe "+out_file_prefix+".off" ;
        system(tetgen_command.c_str()) ;

        // Read the mesh built by tetgen for L
        Tet_object_io<Traits> tetL(out_file_prefix) ;

        // Build the associated SimpComplex
        Chain_complex& L = *new Chain_complex(tetL) ;
        std::cout << "------ L:" << L;

        // Build the Sub_chain_complex_mask encoding _K inside L
        Sub_chain_complex& K(*new Sub_chain_complex(L, false)) ;
        // Visit all cells of _K and activate the corresponding bit in K
        for (int q=0; q<=_K.dimension(); ++q)
        {
            for (size_t i=0; i<_K.number_of_cells(q); ++i)
            {
//                const std::vector<size_t>& simplex(_K._ind2simp.at(q).at(i).get_vertices())
                const Simplex&  simplex(_K.index_to_cell(i,q)) ;
//                const size_t id(L._simp2ind.at(q)[simplex]) ;
                const size_t id(L.cell_to_index(simplex));
                K.set_bit_on(q, id) ;
            }
        }
        Complex_duality_data t = {L,K,tetL.nodes} ;
        return t ;
    }

    /** \brief Exports a SimpComplex to a MeshObject  */
    static Mesh_object_io<Traits>& export_meshObject(const Chain_complex& _CC)
    {
        std::vector<Io_cell_type> vcells ;
        for (int q = 0; q <= _CC.dimension(); ++q)
            for (size_t i = 0; i<_CC.number_of_cells(q); ++i)
            {
                const Simplex& s(_CC.index_to_cell(i,q)) ;
                vcells.push_back(s.get_vertices()) ;
            }

        std::vector<Point> coords;
        for (auto it = _CC.get_vertices_coords().begin(); it != _CC.get_vertices_coords().end(); ++it)
            coords.push_back(Point(*it));
        Mesh_object_io<Traits> &m = *(new Mesh_object_io<Traits>(-3, coords, vcells)) ;
        return m ;
    }
} ;

// =========== Cubical chain complex tools
// Build K sub-chain complex of L (homeomorphic to B^n)
// Adjust L size

/*!
 \ingroup PkgHDVFAlgorithmClasses

The class `Duality_cubical_complex_tools` is dedicated to Alexander duality for 3D binary volumes.

Starting from a Cubical_chain_complex `_K`, the method `cubical_chain_complex_bb()` builds a `Cubical_chain_complex` `L` and `Sub_chain_complex_mask` `K`.
- L : complex built of the "full" bounding box of _K
- K (Sub_chain_complex_mask) : Sub_chain_complex_mask identifying _K inside L

Use the `frame` method from the `Cub_object_io` class to enlarge the bounding box (via a 1 pixel dilatation) if necessary.

\tparam CoefficientRing a model of the `IntegralDomainWithoutDivision` concept providing the ring used to compute homology.
*/

template<typename CoefficientRing, typename Traits>
class Duality_cubical_complex_tools {
public:
    /** \brief Type of cubical complexes used for the initial complex and \f$L\f$. */
    typedef Cubical_chain_complex<CoefficientRing, Traits> Chain_complex ;
    /** \brief Type of sub chain complex mask used to encode the sub-complex  \f$K\f$. */
    typedef Sub_chain_complex_mask<Chain_complex> Sub_chain_complex ;
    // Constructor
    Duality_cubical_complex_tools() {}


    /** \brief Generates a subcomplex \f$K\f$K and a complex \f$L\f$ with \f$K\subseteq L\f$ from a cubical complex `_K`.
     *
     * `L` is the bounding box of `_K` (homeomorphic to a ball) and \f$K\f$ is a sub chain complex mask encoding `_K`.
     *
     * The following figures shows the resulting complex with an initial simple cubical complex (right - sectional view):
     * <img src="HDVF_eight_view1.png" align="center" width=35%/>
     * <img src="HDVF_eight_view2.png" align="center" width=30%/>
     *
     * \param[in] _K Initial cubical chain complex (working mesh).
     */
    static std::pair<Chain_complex&, Sub_chain_complex&> cubical_chain_complex_bb (const Chain_complex& _K)
    {
        Cub_object_io tmp_L ;
        tmp_L.dim = _K.dimension() ;
        tmp_L.N = _K.size_bb() ;
        tmp_L.ncubs.resize(tmp_L.dim+1) ;
        // Visit all boolean indices in the BB of _CC and insert corresponding cells in the Cub_object_io of L
        for (size_t i=0; i<_K.size(); ++i)
        {
            const std::vector<size_t> tmpkhal(_K.bindex_to_cell(i)) ;
            const int dtmp(_K.dimension(tmpkhal)) ;
            (tmp_L.ncubs)[dtmp] += 1 ;
            tmp_L.cubs.push_back(tmpkhal) ;
        }
        Chain_complex& L(*new Chain_complex(tmp_L, Chain_complex::PRIMAL)) ;

        // Build the Sub_chain_complex_mask corresponding to _CC
        Sub_chain_complex& K(*new Sub_chain_complex(L, false)) ;
        // Visit all cells of _CC and activate the corresponding bit in K
        for (int q=0; q<=_K.dimension(); ++q)
        {
            for (size_t i=0; i<_K.number_of_cells(q); ++i)
            {
                const std::vector<size_t> khal(_K.index_to_cell(i, q)) ;
                const size_t j = L.cell_to_index(khal);
                K.set_bit_on(q,j) ;
            }
        }
        return std::pair<Chain_complex&, Sub_chain_complex&>(L,K) ;
    }
} ;

} /* end namespace Homological_discrete_vector_field */
} /* end namespace CGAL */

#endif // CGAL_HDVF_GEOMETRIC_CHAIN_COMPLEX_TOOLS_H
