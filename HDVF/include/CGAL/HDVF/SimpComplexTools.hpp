#ifndef SIMPCOMPLEXTOOLS_HPP
#define SIMPCOMPLEXTOOLS_HPP

#include <vector>
#include <map>
#include <stdexcept>
#include <unordered_set>
//#include "config.h"
#include "Simplex.hpp"
#include "Abstract_simplicial_chain_complex.hpp"
#include "hdvf.hpp"
#include "per_hdvf.hpp"
#include "sub_chain_complex.hpp"
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

// SimpComplex

/** \brief vtk export for SimpComplex */

template <typename CoefType, template <typename, int> typename _ChainType = OSM::Chain, template <typename, int> typename _SparseMatrixType = OSM::SparseMatrix, typename VertexIdType = int>
void Simp_output_vtk (HDVF<CoefType, Simplicial_chain_complex<CoefType>, _ChainType, _SparseMatrixType> &hdvf, Simplicial_chain_complex<CoefType> &complex, string filename = "test")
{
    typedef Simplicial_chain_complex<CoefType> ComplexType;
    // Export PSC labelling
    string outfile(filename+"_PSC.vtk") ;
    vector<vector<int> > labels = hdvf.export_label(PSC) ;
    ComplexType::Simplicial_chain_complex_to_vtk(complex, outfile, &labels) ;
    
    if (hdvf.get_hdvf_opts() != OPT_BND)
    {
        // Export generators of all critical cells
        vector<vector<int> > criticals(hdvf.get_flag(CRITICAL)) ;
        for (int q = 0; q <= complex.dim(); ++q)
        {
            for (int c : criticals.at(q))
            {
                if (hdvf.get_hdvf_opts() & (OPT_FULL | OPT_G))
                {
                    string outfile_g(filename+"_G_"+to_string(c)+"_dim_"+to_string(q)+".vtk") ;
//                    vector<vector<int> > labels = hdvf.export_label(G,c,q) ;
                    OSM::Chain<CoefType,OSM::COLUMN> chain(hdvf.export_GChain(c,q)) ;
                    ComplexType::Simplicial_chain_complex_chain_to_vtk(complex, outfile_g, chain, q, c) ;
                }
                if (q < complex.dim())
                {
                    if (hdvf.get_hdvf_opts() & (OPT_FULL | OPT_F))
                    {
                        string outfile_f(filename+"_FSTAR_"+to_string(c)+"_dim_"+to_string(q)+".vtk") ;
                        OSM::Chain<CoefType,OSM::COLUMN> chain(hdvf.export_FSTARChain(c,q)) ;
                        ComplexType::Simplicial_chain_complex_chain_to_vtk(complex, outfile_f, chain, q+1, c) ;
                    }
                }
            }
        }
    }
}

// Duality

template<typename CoefficientType>
class Duality_SimpComplexTools {
public:
    typedef Simplicial_chain_complex<CoefficientType> _ComplexType ;
    typedef SubChainComplex<CoefficientType, _ComplexType> _SubCCType ;
    // Constructor
    Duality_SimpComplexTools() {}
    
    /** \brief Build a SimpComplex L and SubChainComplex K from a SimpComplex _K.
     * L : complex built out of _K and a closing icosphere meshed by tetgen
     * K (SubChainComplex) : SubChainComplex identifying _K inside L
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
        Mesh_object mesh_L = Duality_SimpComplexTools::export_meshObject(_K) ;
        
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
        
        // Build the SubChainComplex encoding _K inside L
        _SubCCType& K(*new _SubCCType(L, false)) ;
        // Visit all cells of _K and activate the corresponding bit in K
        for (int q=0; q<=_K.dim(); ++q)
        {
            for (int i=0; i<_K.nb_cells(q); ++i)
            {
                const std::set<int>& simplex(_K._ind2simp.at(q).at(i).getVertices()) ;
                const int id(L._simp2ind.at(q)[simplex]) ;
                K.set_BitOn(q, id) ;
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
            for (int i = 0; i<_CC.nb_cells(q); ++i)
            {
                const Simplex& s(_CC._ind2simp.at(q).at(i)) ;
                vcells.push_back(s.getVertices()) ;
            }
                
        Mesh_object &m = *(new Mesh_object(-3, _CC.get_nodes_coords(), vcells)) ;
        return m ;
    }
} ;

// Persistent homology


/** \brief Export persistent information to vtk files */
template <typename CoefType, typename DegType>
void Per_Simp_output_vtk (PersistentHDVF<CoefType, Simplicial_chain_complex<CoefType>, DegType> &per_hdvf, Simplicial_chain_complex<CoefType> &complex, string filename = "per")
{
    if (!per_hdvf.with_export())
        throw("Cannot export persistent generators to vtk: with_export is off!") ;
    
    using perHDVFType = PersistentHDVF<CoefType, Simplicial_chain_complex<CoefType>, DegType> ;
    using ComplexType = Simplicial_chain_complex<CoefType> ;
    using PerHole = PerHoleT<DegType> ;
    // Iterate over persistence diagram (iterator over non zero intervals)
    // Batch informations are stored in file filename_infos.txt
    std::ofstream info_file(filename+"_infos.txt") ;
    int i = 0 ;
    for (typename perHDVFType::iterator it = per_hdvf.begin(); it != per_hdvf.end(); ++it)
    {
        typename perHDVFType::Exp_infos hole_data(*it) ;
        const PerHole hole(hole_data.hole) ;
        // Export informations of this hole
        info_file << i << " -> " << " --- duration : " << per_hdvf.hole_duration(hole) << " -- " << hole << std::endl ;
        
        if (per_hdvf.hole_duration(hole)>0) // Export vtk of "finite" holes
        {
            // Build name associated to the ith hole : filename_i
            string out_file = filename+"_"+to_string(i) ;
            // Export PSC labels to vtk
            ComplexType::Simplicial_chain_complex_to_vtk(complex, out_file+"_PSC.vtk", &hole_data.labelsPSC) ;
            const PerIntervalCells per_int_cells(std::get<1>(hole_data.hole)) ;
            // Export homology generators (g)
            if (per_hdvf.get_hdvf_opts()  & (OPT_FULL | OPT_G))
            {
                // First generator : filename_i_g_sigma_q.vtk
                {
                    const int id(per_int_cells.first.first), dim(per_int_cells.first.second) ;
                    string tmp(out_file+"_g_"+to_string(id)+"_"+to_string(dim)+".vtk") ;
                    ComplexType::Simplicial_chain_complex_chain_to_vtk(complex, tmp, hole_data.g_chain_sigma, dim);
                }
                // Second generator : filename_i_g_tau_q+1.vtk
                {
                    const int id(per_int_cells.second.first), dim(per_int_cells.second.second) ;
                    string tmp(out_file+"_g_"+to_string(id)+"_"+to_string(dim)+".vtk") ;
                    ComplexType::Simplicial_chain_complex_chain_to_vtk(complex, tmp, hole_data.g_chain_tau, dim);
                }
            }
            // Export cohomology generators (fstar)
            if ((per_hdvf.get_hdvf_opts() == OPT_FULL) || (per_hdvf.get_hdvf_opts() == OPT_F))
            {
                // First generator : filename_i_fstar_sigma_q.vtk
                {
                    const int id(per_int_cells.first.first), dim(per_int_cells.first.second) ;
                    string tmp(out_file+"_fstar_"+to_string(id)+"_"+to_string(dim)+".vtk") ;
                    ComplexType::Simplicial_chain_complex_chain_to_vtk(complex, tmp, hole_data.fstar_chain_sigma, dim+1);
                }
                // Second generator : filename_i_fstar_tau_q+1.vtk
                {
                    const int id(per_int_cells.second.first), dim(per_int_cells.second.second) ;
                    string tmp(out_file+"_fstar_"+to_string(id)+"_"+to_string(dim)+".vtk") ;
                    ComplexType::Simplicial_chain_complex_chain_to_vtk(complex, tmp, hole_data.fstar_chain_tau, dim+1);
                }
            }
        }
        ++i ;
    }
    info_file.close() ;
    Simp_output_vtk<CoefType>(per_hdvf, complex, filename+"_inf") ;
}

#endif // SimpComplexTools_HPP
