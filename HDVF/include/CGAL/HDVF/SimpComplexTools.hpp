#ifndef SIMPCOMPLEXTOOLS_HPP
#define SIMPCOMPLEXTOOLS_HPP

#include <vector>
#include <map>
#include <stdexcept>
#include <unordered_set>
//#include "config.h"
#include "Simplex.hpp"
#include "Abstract_simplicial_chain_complex.hpp"
#include "Hdvf_core.h"
#include "Hdvf_persistence.h"
#include "Hdvf_duality.h"
#include "Sub_chain_complex_mask.h"
#include "tools_io.hpp"
#include "CGAL/OSM/OSM.hpp"

/**
 * \class SimpComplexTools
 * \brief Provides tools for SimpComplexes ((co)homology, duality, persistent homology).
 * \brief Friend class of SimpComplex
 *
 * \author Bac A.
 * \version 0.2.0
 * \date 06/11/2024
 */

namespace CGAL {
namespace HDVF {

// SimpComplex

/** \brief Export HDVF for a Simplicial_complex to a vtk file.
 *
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

template <typename CoefType, template <typename, int> typename _ChainType = OSM::Sparse_chain, template <typename, int> typename _SparseMatrixType = OSM::Sparse_matrix, typename VertexIdType = size_t>
void Simp_output_vtk (Hdvf_core<CoefType, Simplicial_chain_complex<CoefType>, _ChainType, _SparseMatrixType> &hdvf, Simplicial_chain_complex<CoefType> &complex, string filename = "test", bool co_faces = false)
{
    typedef Simplicial_chain_complex<CoefType> ComplexType;
    typedef Hdvf_core<CoefType, Simplicial_chain_complex<CoefType>, _ChainType, _SparseMatrixType> HDVF_type;
    // Export PSC labelling
    string outfile(filename+"_PSC.vtk") ;
    vector<vector<int> > labels = hdvf.export_psc_labels() ;
    ComplexType::Simplicial_chain_complex_to_vtk(complex, outfile, &labels) ;
    
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
                    OSM::Sparse_chain<CoefType,OSM::COLUMN> chain(hdvf.export_homology_chain(c,q)) ;
                    ComplexType::Simplicial_chain_complex_chain_to_vtk(complex, outfile_g, chain, q, c) ;
                }
                // Cohomology generators
                if (hdvf.get_hdvf_opts() & (OPT_FULL | OPT_F))
                {
                    string outfile_f(filename+"_FSTAR_"+to_string(c)+"_dim_"+to_string(q)+".vtk") ;
                    OSM::Sparse_chain<CoefType,OSM::COLUMN> chain(hdvf.export_cohomology_chain(c, q)) ;
                    if (!co_faces)
                    {
                        ComplexType::Simplicial_chain_complex_chain_to_vtk(complex, outfile_f, chain, q, c) ;
                    }
                    else
                    {
                        if (q < complex.dim())
                        {
                            ComplexType::Simplicial_chain_complex_chain_to_vtk(complex, outfile_f, complex.cofaces_chain(chain, q), q+1, c) ;
                        }
                    }
                }
            }
        }
    }
}

// Duality

template<typename CoefficientType>
class Duality_simplicial_complex_tools {
public:
    typedef Simplicial_chain_complex<CoefficientType> _ComplexType ;
    typedef Sub_chain_complex_mask<CoefficientType, _ComplexType> _SubCCType ;
    // Constructor
    Duality_simplicial_complex_tools() {}
    
    /** \brief Build a SimpComplex L and Sub_chain_complex_mask K from a SimpComplex _K.
     * L : complex built out of _K and a closing icosphere meshed by tetgen
     * K (Sub_chain_complex_mask) : Sub_chain_complex_mask identifying _K inside L
     */
    
    /** \brief Build the bounding box CubComplex of _CC  */
    typedef struct ztriple {
        _ComplexType& L ;
        _SubCCType& K ;
        std::vector<IONodeType> nodes ;
    } TripleRes ;
    
    static TripleRes SimpComplexBB (const _ComplexType& _K, double BB_ratio=1.5, const string& out_file_prefix = "file_K_closed.off")
    {
        // Export _CC to a MeshObject to add the icosphere and mesh with tetGen
        Mesh_object mesh_L = Duality_simplicial_complex_tools::export_meshObject(_K) ;
        
        // Closing K by adding the icosphere
        //  Compute a bounding icosphere
        IONodeType bary = mesh_L.barycenter() ;
        double r = mesh_L.radius(bary) ;
        IcosphereObject ico(2,bary, BB_ratio*r) ;
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
        TetObject tetL(out_file_prefix) ;
        
        // Build the associated SimpComplex
        _ComplexType& L = *new _ComplexType(tetL, tetL.nodes) ;
        
        // Build the Sub_chain_complex_mask encoding _K inside L
        _SubCCType& K(*new _SubCCType(L, false)) ;
        // Visit all cells of _K and activate the corresponding bit in K
        for (int q=0; q<=_K.dim(); ++q)
        {
            for (size_t i=0; i<_K.nb_cells(q); ++i)
            {
                const std::set<size_t>& simplex(_K._ind2simp.at(q).at(i).getVertices()) ;
                const size_t id(L._simp2ind.at(q)[simplex]) ;
                K.set_bit_on(q, id) ;
            }
        }
        TripleRes t = {L,K,tetL.nodes} ;
        return t ;
    }
    
    /** \brief Export a SimpComplex to a MeshObject  */
    static Mesh_object& export_meshObject(const _ComplexType& _CC)
    {
        std::vector<IOCellType> vcells ;
        for (int q = 0; q <= _CC.dim(); ++q)
            for (size_t i = 0; i<_CC.nb_cells(q); ++i)
            {
                const Simplex& s(_CC._ind2simp.at(q).at(i)) ;
                vcells.push_back(s.getVertices()) ;
            }
        
        Mesh_object &m = *(new Mesh_object(-3, _CC.get_vertices_coords(), vcells)) ;
        return m ;
    }
} ;

template <typename CoefType, typename VertexIdType = size_t>
void Dual_simp_output_vtk (Hdvf_duality<CoefType, Simplicial_chain_complex<CoefType> > &hdvf, Simplicial_chain_complex<CoefType> &complex, string filename = "test", bool co_faces = false)
{
    typedef Simplicial_chain_complex<CoefType> ComplexType;
    typedef Hdvf_duality<CoefType, Simplicial_chain_complex<CoefType> > HDVF_type;
    // Export PSC labelling
    string outfile(filename+"_PSC.vtk") ;
    vector<vector<int> > labels = hdvf.export_psc_labels() ;
    ComplexType::Simplicial_chain_complex_to_vtk(complex, outfile, &labels) ;
    
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
                    OSM::Sparse_chain<CoefType,OSM::COLUMN> chain(hdvf.export_homology_chain(c,q)) ;
                    ComplexType::Simplicial_chain_complex_chain_to_vtk(complex, outfile_g, chain, q, c) ;
                }
                // Cohomology generators
                if (hdvf.get_hdvf_opts() & (OPT_FULL | OPT_F))
                {
                    string outfile_f(filename+"_FSTAR_"+to_string(c)+"_dim_"+to_string(q)+".vtk") ;
                    OSM::Sparse_chain<CoefType,OSM::COLUMN> chain(hdvf.export_cohomology_chain(c, q)) ;
                    if (!co_faces)
                    {
                        ComplexType::Simplicial_chain_complex_chain_to_vtk(complex, outfile_f, chain, q, c) ;
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
                            ComplexType::Simplicial_chain_complex_chain_to_vtk(complex, outfile_f, cofaces_chain, q+1, c) ;
                        }
                    }
                }
            }
        }
    }
}

// Persistent homology


/** \brief Export persistent information to vtk files */
template <typename CoefType, typename DegType, typename FiltrationType>
void Per_Simp_output_vtk (Hdvf_persistence<CoefType, Simplicial_chain_complex<CoefType>, DegType, FiltrationType> &per_hdvf, Simplicial_chain_complex<CoefType> &complex, string filename = "per", bool co_faces = false)
{
    if (!per_hdvf.with_export())
        throw("Cannot export persistent generators to vtk: with_export is off!") ;
    
    using perHDVFType = Hdvf_persistence<CoefType, Simplicial_chain_complex<CoefType>, DegType, FiltrationType> ;
    using ComplexType = Simplicial_chain_complex<CoefType> ;
    using PerHole = PerHoleT<DegType> ;
    
    // Export the filtration
    string out_file_filtration = filename+"_filtration.vtk" ;
    vector<vector<size_t> > filtr_labels = per_hdvf.get_filtration().export_filtration();
    ComplexType::template Simplicial_chain_complex_to_vtk<size_t>(complex, out_file_filtration, &filtr_labels, "unsigned_long") ;
    
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
            ComplexType::Simplicial_chain_complex_to_vtk(complex, out_file+"_PSC.vtk", &hole_data.labelsPSC) ;
            const CellsPerInterval per_int_cells(std::get<1>(hole_data.hole)) ;
            // Export homology generators (g)
            if (per_hdvf.get_hdvf_opts()  & (OPT_FULL | OPT_G))
            {
                // First generator : filename_i_g_sigma_q.vtk
                {
                    const size_t id(per_int_cells.first.first) ;
                    const int dim(per_int_cells.first.second) ;
                    string tmp(out_file+"_g_"+to_string(id)+"_"+to_string(dim)+".vtk") ;
                    ComplexType::Simplicial_chain_complex_chain_to_vtk(complex, tmp, hole_data.g_chain_sigma, dim);
                }
                // Second generator : filename_i_g_tau_q+1.vtk
                {
                    const size_t id(per_int_cells.second.first) ;
                    const int dim(per_int_cells.second.second) ;
                    string tmp(out_file+"_g_"+to_string(id)+"_"+to_string(dim)+".vtk") ;
                    ComplexType::Simplicial_chain_complex_chain_to_vtk(complex, tmp, hole_data.g_chain_tau, dim);
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
                        ComplexType::Simplicial_chain_complex_chain_to_vtk(complex, tmp, complex.cofaces_chain(hole_data.fstar_chain_sigma, dim), dim+1);
                    }
                    else
                        ComplexType::Simplicial_chain_complex_chain_to_vtk(complex, tmp, hole_data.fstar_chain_sigma, dim);
                }
                // Second generator : filename_i_fstar_tau_q+1.vtk
                {
                    const size_t id(per_int_cells.second.first) ;
                    const int dim(per_int_cells.second.second) ;
                    string tmp(out_file+"_fstar_"+to_string(id)+"_"+to_string(dim)+".vtk") ;
                    if (co_faces)
                    {
                        ComplexType::Simplicial_chain_complex_chain_to_vtk(complex, tmp, complex.cofaces_chain(hole_data.fstar_chain_tau, dim), dim+1);
                    }
                    else
                        ComplexType::Simplicial_chain_complex_chain_to_vtk(complex, tmp, hole_data.fstar_chain_tau, dim);
                }
            }
        }
        ++i ;
    }
    info_file.close() ;
    Simp_output_vtk<CoefType>(per_hdvf, complex, filename+"_inf", co_faces) ;
}

} /* end namespace HDVF */
} /* end namespace CGAL */

#endif // SimpComplexTools_HPP
