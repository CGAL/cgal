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

#ifndef Duality_cubical_complex_tools_HPP
#define Duality_cubical_complex_tools_HPP

#include <vector>
#include <map>
#include <stdexcept>
#include <unordered_set>
#include "Cubical_chain_complex.hpp"
#include "Hdvf_core.h"
#include "Hdvf_persistence.h"
#include "Sub_chain_complex_mask.h"
#include "tools_io.hpp"
#include "CGAL/OSM/OSM.hpp"

/**
 * \class Duality_cubical_complex_tools
 * \brief Provides tools for CubComplexes ((co)homology, duality, persistent homology).
 * \brief Friend class of CubComplex
 *
 * \author Bac A.
 * \version 0.2.0
 * \date 06/11/2024
 */

// CubComplex

namespace CGAL {
namespace HDVF {

/** \brief Export HDVF for a Cubical_complex to a vtk file.
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

template <typename CoefType, template <typename, int> typename _ChainType = OSM::Sparse_chain, template <typename, int> typename _SparseMatrixType = OSM::Sparse_matrix>
void Cub_output_vtk (Hdvf_core<CoefType, Cubical_chain_complex<CoefType>, _ChainType, _SparseMatrixType> &hdvf, Cubical_chain_complex<CoefType> &complex, string filename = "test", bool co_faces = false)
{
    typedef Hdvf_core<CoefType, Cubical_chain_complex<CoefType>, _ChainType, _SparseMatrixType> HDVF_type;
    typedef Cubical_chain_complex<CoefType> ComplexType;
    
    // Export PSC labelling
    vector<vector<int> > labels = hdvf.export_psc_labels() ;
    string outfilePSC(filename+"_PSC.vtk") ;
    ComplexType::Cubical_chain_complex_to_vtk(complex, outfilePSC, &labels) ;
    
    // Export generators of all critical cells if available
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
                    OSM::Sparse_chain<CoefType,OSM::COLUMN> chain(hdvf.export_homology_chain(c,q)) ;
                    ComplexType::Cubical_chain_complex_chain_to_vtk(complex, outfile_g, chain, q, c) ;
                }
                // Cohomology generators
                if (hdvf.get_hdvf_opts() & (OPT_FULL | OPT_F))
                {
                    string outfile_f(filename+"_FSTAR_"+to_string(c)+"_dim_"+to_string(q)+".vtk") ;
                    OSM::Sparse_chain<CoefType,OSM::COLUMN> chain(hdvf.export_cohomology_chain(c, q)) ;
                    if (!co_faces)
                    {
                        ComplexType::Cubical_chain_complex_chain_to_vtk(complex, outfile_f, chain, q, c) ;
                    }
                    else
                    {
                        if (q < complex.dim())
                        {
                            ComplexType::Cubical_chain_complex_chain_to_vtk(complex, outfile_f, complex.cofaces_chain(chain, q), q+1, c) ;
                        }
                    }
                }
            }
        }
    }
}

// Duality

template<typename T>
class Duality_cubical_complex_tools {
public:
    typedef Cubical_chain_complex<T> _ComplexType ;
    typedef Sub_chain_complex_mask<T, _ComplexType> _SubCCType ;
    // Constructor
    Duality_cubical_complex_tools() {}
    
    
    /** \brief Build a CubComplex L and Sub_chain_complex_mask K from a CubComplex.
     * L : full bounding box
     * K (Sub_chain_complex_mask) : CubComplex
     */
    
    // Tools
    
    /** \brief Build the bounding box CubComplex of _CC  */
    static std::pair<_ComplexType&, _SubCCType&> CubComplexBB (const _ComplexType& _CC)
    {
        Cub_object tmp ;
        tmp.dim = _CC.dim() ;
        tmp.N = _CC._size_bb ;
        tmp.ncubs.resize(tmp.dim+1) ;
        //        tmp.ncubs = tmp.N ;
        // Visit all boolean indices in the BB of _CC and insert corresponding Cells
        for (size_t i=0; i<_CC._P.at(_CC.dim()); ++i)
        {
            const std::vector<size_t> tmpkhal(_CC.ind2khal(i)) ;
            const int dtmp(_CC.calculate_dimension(tmpkhal)) ;
            (tmp.ncubs)[dtmp] += 1 ;
            tmp.cubs.push_back(tmpkhal) ;
        }
        _ComplexType& L(*new _ComplexType(tmp, _ComplexType::PRIMAL)) ;
        
        // Build the Sub_chain_complex_mask corresponding to _CC
        _SubCCType& K(*new _SubCCType(L, false)) ;
        // Visit all cells of _CC and activate the corresponding bit in K
        for (int q=0; q<=_CC.dim(); ++q)
        {
            for (size_t i=0; i<_CC.nb_cells(q); ++i)
            {
                const std::vector<size_t> khal(_CC.ind2khal(_CC._base2bool.at(q).at(i))) ;
                const size_t j = L._bool2base.at(q).at(L.khal2ind(khal)) ;
                K.set_bit_on(q,j) ;
            }
        }
        return std::pair<_ComplexType&, _SubCCType&>(L,K) ;
    }
} ;

template <typename CoefType, typename VertexIdType = size_t>
void Dual_cub_output_vtk (Hdvf_duality<CoefType, Cubical_chain_complex<CoefType> > &hdvf, Cubical_chain_complex<CoefType> &complex, string filename = "test", bool co_faces = false)
{
    typedef Hdvf_duality<CoefType, Simplicial_chain_complex<CoefType> > HDVF_type;
    typedef Cubical_chain_complex<CoefType> ComplexType;
    
    // Export PSC labelling
    vector<vector<int> > labels = hdvf.export_psc_labels() ;
    string outfilePSC(filename+"_PSC.vtk") ;
    ComplexType::Cubical_chain_complex_to_vtk(complex, outfilePSC, &labels) ;
    
    // Export generators of all critical cells if available
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
                    OSM::Sparse_chain<CoefType,OSM::COLUMN> chain(hdvf.export_homology_chain(c,q)) ;
                    ComplexType::Cubical_chain_complex_chain_to_vtk(complex, outfile_g, chain, q, c) ;
                }
                // Cohomology generators
                if (hdvf.get_hdvf_opts() & (OPT_FULL | OPT_F))
                {
                    string outfile_f(filename+"_FSTAR_"+to_string(c)+"_dim_"+to_string(q)+".vtk") ;
                    OSM::Sparse_chain<CoefType,OSM::COLUMN> chain(hdvf.export_cohomology_chain(c, q)) ;
                    if (!co_faces)
                    {
                        ComplexType::Cubical_chain_complex_chain_to_vtk(complex, outfile_f, chain, q, c) ;
                    }
                    else
                    {
                        if (q < complex.dim())
                        {
                            OSM::Sparse_chain<CoefType,OSM::COLUMN> cofaces_chain(complex.cofaces_chain(chain, q)) ;
                            Sub_chain_complex_mask<CoefType,ComplexType> sub(hdvf.get_current_mask());
                            sub.screen_chain(cofaces_chain, q+1);
                            ComplexType::Cubical_chain_complex_chain_to_vtk(complex, outfile_f, cofaces_chain, q+1, c) ;
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
void Per_Cub_output_vtk (Hdvf_persistence<CoefType, Cubical_chain_complex<CoefType>, DegType, FiltrationType> &per_hdvf, Cubical_chain_complex<CoefType> &complex, string filename = "per", bool co_faces = false)
{
    if (!per_hdvf.with_export())
        throw("Cannot export persistent generators to vtk: with_export is off!") ;
    
    using perHDVFType = Hdvf_persistence<CoefType, Cubical_chain_complex<CoefType>, DegType, FiltrationType> ;
    using ComplexType = Cubical_chain_complex<CoefType> ;
    using PerHole = PerHoleT<DegType> ;
    
    // Export the filtration
    string out_file_filtration = filename+"_filtration.vtk" ;
    vector<vector<size_t> > filtr_labels = per_hdvf.get_filtration().export_filtration();
    ComplexType::template Cubical_chain_complex_to_vtk<size_t>(complex, out_file_filtration, &filtr_labels, "unsigned_long") ;
    
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
        
        if (per_hdvf.hole_duration(hole)>=0)
        {
            // Build name associated to the ith hole : filename_i
            string out_file = filename+"_"+to_string(i) ;
            // Export PSC labels to vtk
            ComplexType::Cubical_chain_complex_to_vtk(complex, out_file+"_PSC.vtk", &hole_data.labelsPSC) ;
            const CellsPerInterval per_int_cells(std::get<1>(hole_data.hole)) ;
            // Export homology generators (g)
            if (per_hdvf.get_hdvf_opts()  & (OPT_FULL | OPT_G))
            {
                // First generator : filename_i_g_sigma_q.vtk
                {
                    const size_t id(per_int_cells.first.first);
                    const int dim(per_int_cells.first.second) ;
                    string tmp(out_file+"_g_"+to_string(id)+"_"+to_string(dim)+".vtk") ;
                    ComplexType::Cubical_chain_complex_chain_to_vtk(complex, tmp, hole_data.g_chain_sigma, dim);
                }
                // Second generator : filename_i_g_tau_q+1.vtk
                {
                    const size_t id(per_int_cells.second.first) ;
                    const int dim(per_int_cells.second.second) ;
                    string tmp(out_file+"_g_"+to_string(id)+"_"+to_string(dim)+".vtk") ;
                    ComplexType::Cubical_chain_complex_chain_to_vtk(complex, tmp, hole_data.g_chain_tau, dim);
                }
            }
            // Export cohomology generators (fstar)
            if ((per_hdvf.get_hdvf_opts() == OPT_FULL) || (per_hdvf.get_hdvf_opts() == OPT_F))
            {
                
                // First generator : filename_i_fstar_sigma_q.vtk
                {
                    const size_t id(per_int_cells.first.first);
                    const int dim(per_int_cells.first.second) ;
                    string tmp(out_file+"_fstar_"+to_string(id)+"_"+to_string(dim)+".vtk") ;
                    if (co_faces)
                    {
                        ComplexType::Cubical_chain_complex_chain_to_vtk(complex, tmp, complex.cofaces_chain(hole_data.fstar_chain_sigma, dim), dim+1);
                    }
                    else
                        ComplexType::Cubical_chain_complex_chain_to_vtk(complex, tmp, hole_data.fstar_chain_sigma, dim);
                }
                // Second generator : filename_i_fstar_tau_q+1.vtk
                {
                    const size_t id(per_int_cells.second.first) ;
                    const int dim(per_int_cells.second.second) ;
                    string tmp(out_file+"_fstar_"+to_string(id)+"_"+to_string(dim)+".vtk") ;
                    if (co_faces)
                    {
                        ComplexType::Cubical_chain_complex_chain_to_vtk(complex, tmp, complex.cofaces_chain(hole_data.fstar_chain_tau, dim), dim+1);
                    }
                    else
                        ComplexType::Cubical_chain_complex_chain_to_vtk(complex, tmp, hole_data.fstar_chain_tau, dim);
                }
            }
        }
        ++i ;
    }
    info_file.close() ;
    Cub_output_vtk<CoefType>(per_hdvf, complex, filename+"_inf", co_faces) ;
}

} /* end namespace HDVF */
} /* end namespace CGAL */

#endif // Duality_cubical_complex_tools_HPP
