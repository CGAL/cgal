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

template <typename CoefType, template <typename, int> typename _ChainType = OSM::Sparse_chain, template <typename, int> typename _SparseMatrixType = OSM::Sparse_matrix>
void Cub_output_vtk (Hdvf_core<CoefType, Cubical_chain_complex<CoefType>, _ChainType, _SparseMatrixType> &hdvf, Cubical_chain_complex<CoefType> &complex, string filename = "test")
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
        vector<vector<int> > criticals(hdvf.get_flag(CRITICAL)) ;
        for (int q = 0; q <= complex.dim(); ++q)
        {
            for (int c : criticals.at(q))
            {
                if (hdvf.get_hdvf_opts() & (OPT_FULL | OPT_G))
                {
                    string outfile_g(filename+"_G_"+to_string(c)+"_dim_"+to_string(q)+".vtk") ;
                    //                    vector<vector<int> > labels = hdvf.export_label(G,c,q) ;
                    OSM::Sparse_chain<CoefType,OSM::COLUMN> chain(hdvf.export_homology_chain(c,q)) ;
                    ComplexType::Cubical_chain_complex_chain_to_vtk(complex, outfile_g, chain, q, c) ;
                }
                if (q < complex.dim())
                {
                    if (hdvf.get_hdvf_opts() & (OPT_FULL | OPT_F))
                    {
                        string outfile_f(filename+"_FSTAR_"+to_string(c)+"_dim_"+to_string(q)+".vtk") ;
                        OSM::Sparse_chain<CoefType,OSM::COLUMN> chain(hdvf.export_cohomology_chain(c,q)) ;
                        ComplexType::Cubical_chain_complex_chain_to_vtk(complex, outfile_f, chain, q+1, c) ;
                    }
                }
            }
        }
    }
}
//        vector<vector<int> > criticals(hdvf.get_flag(HDVF_coreT::CRITICAL)) ;
//        for (int q = 0; q <= complex.dim(); ++q)
//        {
//            for (int c : criticals.at(q))
//            {
//                if (hdvf.get_hdvf_opts() & (OPT_FULL | OPT_G))
//                {
//                    string outfile(filename+"_G_"+to_string(c)+"_dim_"+to_string(q)+".vtk") ;
//                    vector<vector<int> > labels = hdvf.export_label(G,c,q) ;
//                    Cubical_chain_complex_write_to_vtk<CoefType>(complex, outfile, &labels, BOOLEAN_LABELS, "bool") ;
//                }
//                if (q < complex.dim())
//                {
//                    if (hdvf.get_hdvf_opts() & (OPT_FULL | OPT_F))
//                    {
//                        string outfile_f(filename+"_FSTAR_"+to_string(c)+"_dim_"+to_string(q)+".vtk") ;
//                        vector<vector<int> > labels_fstar = hdvf.export_label(FSTAR,c,q) ;
//                        Cubical_chain_complex_write_to_vtk<CoefType>(complex, outfile_f, &labels_fstar, BOOLEAN_LABELS, "bool") ;
//                    }
//                }
//            }
//        }
//    }
//}

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
        for (int i=0; i<_CC._P.at(_CC.dim()); ++i)
        {
            const std::vector<int> tmpkhal(_CC.ind2khal(i)) ;
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
            for (int i=0; i<_CC.nb_cells(q); ++i)
            {
                const std::vector<int> khal(_CC.ind2khal(_CC._base2bool.at(q).at(i))) ;
                const int j = L._bool2base.at(q).at(L.khal2ind(khal)) ;
                K.set_bit_on(q,j) ;
            }
        }
        return std::pair<_ComplexType&, _SubCCType&>(L,K) ;
    }
} ;

// Persistent homology

/** \brief Export persistent information to vtk files */
template <typename CoefType, typename DegType, typename FiltrationType>
void Per_Cub_output_vtk (Hdvf_persistence<CoefType, Cubical_chain_complex<CoefType>, DegType, FiltrationType> &per_hdvf, Cubical_chain_complex<CoefType> &complex, string filename = "per")
{
    if (!per_hdvf.with_export())
        throw("Cannot export persistent generators to vtk: with_export is off!") ;
    
    using perHDVFType = Hdvf_persistence<CoefType, Cubical_chain_complex<CoefType>, DegType, FiltrationType> ;
    using ComplexType = Cubical_chain_complex<CoefType> ;
    using PerHole = PerHoleT<DegType> ;
    
    // Export the filtration
    string out_file_filtration = filename+"_filtration.vtk" ;
    vector<vector<int> > filtr_labels = per_hdvf.get_filtration().export_filtration();
    ComplexType::Cubical_chain_complex_to_vtk(complex, out_file_filtration, &filtr_labels) ;
    
    // Iterate over persistence diagram (iterator over non zero intervals)
    // Batch informations are stored in file filename_infos.txt
    std::ofstream info_file(filename+"_infos.txt") ;
    int i = 0 ;
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
                    const int id(per_int_cells.first.first), dim(per_int_cells.first.second) ;
                    string tmp(out_file+"_g_"+to_string(id)+"_"+to_string(dim)+".vtk") ;
                    ComplexType::Cubical_chain_complex_chain_to_vtk(complex, tmp, hole_data.g_chain_sigma, dim);
                }
                // Second generator : filename_i_g_tau_q+1.vtk
                {
                    const int id(per_int_cells.second.first), dim(per_int_cells.second.second) ;
                    string tmp(out_file+"_g_"+to_string(id)+"_"+to_string(dim)+".vtk") ;
                    ComplexType::Cubical_chain_complex_chain_to_vtk(complex, tmp, hole_data.g_chain_tau, dim);
                }
            }
            // Export cohomology generators (fstar)
            if ((per_hdvf.get_hdvf_opts() == OPT_FULL) || (per_hdvf.get_hdvf_opts() == OPT_F))
            {
                // First generator : filename_i_fstar_sigma_q.vtk
                {
                    const int id(per_int_cells.first.first), dim(per_int_cells.first.second) ;
                    string tmp(out_file+"_fstar_"+to_string(id)+"_"+to_string(dim)+".vtk") ;
                    ComplexType::Cubical_chain_complex_chain_to_vtk(complex, tmp, hole_data.fstar_chain_sigma, dim);
                }
                // Second generator : filename_i_fstar_tau_q+1.vtk
                {
                    const int id(per_int_cells.second.first), dim(per_int_cells.second.second) ;
                    string tmp(out_file+"_fstar_"+to_string(id)+"_"+to_string(dim)+".vtk") ;
                    ComplexType::Cubical_chain_complex_chain_to_vtk(complex, tmp, hole_data.fstar_chain_tau, dim);
                }
            }
        }
        ++i ;
    }
    info_file.close() ;
    Cub_output_vtk<CoefType>(per_hdvf, complex, filename+"_inf") ;
}

} /* end namespace HDVF */
} /* end namespace CGAL */

#endif // Duality_cubical_complex_tools_HPP
