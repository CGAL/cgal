/**
 * \file per_hom.hpp
 * \brief Persistent homology with HDVF.
 * \author Bac A.
 * \version 0.1.0
 * \date 22/08/2024
 *
 * Persistent homology with HDVF.
 */

#ifndef PER_HOM_HPP
#define PER_HOM_HPP

#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include <functional>
#include "CGAL/OSM/OSM.hpp"
#include "CGAL/HDVF/SubSparseMatrix.hpp"
#include "CGAL/HDVF/hdvf_core.hpp"

/**
 * \class Filtration class
 * \brief Implementation of the filtration.
 *
 * \tparam _CoefficientType The chain's coefficient types (default is OSM::ZCoefficient)
 * \tparam _ComplexType The type of complex  (default is SimpComplex)
 *
 * \author Bac A.
 * \version 0.1.0
 * \date 22/08/2024
 */

namespace CGAL {

template <typename _CoefType, typename _ComplexType, typename _DegType>
class Filtration
{
public:
    typedef std::pair<int, int> CellDim ;
    typedef struct zFiltration_iter_type {
        CellDim cell_dim ;
        _DegType degree ;
    } Filtration_iter_type ;
    
protected:
    const _ComplexType& _K ;
    std::vector<CellDim> _filtration ;
    std::vector<_DegType> _deg ;
    
    // Cell -> (filtration index, basis index according to the filtration)
    std::map<CellDim,int> _cell_to_t ;
    
    typedef OSM::SparseMatrix<_CoefType,OSM::COLUMN> CMatrix ;
    typedef OSM::SparseMatrix<_CoefType,OSM::ROW> RMatrix ;
    typedef OSM::Chain<_CoefType,OSM::COLUMN> CChain ;
    typedef OSM::Chain<_CoefType,OSM::ROW> RChain ;
public:
    // Constructor from a filtration on all cells
    Filtration(const _ComplexType& K, const std::vector<CellDim>& filtration, const std::vector<_DegType>& deg) :
    _K(K), _filtration(filtration), _deg(deg)
    {
        check_filtration() ;
        build_filtration_structure() ;
    }
    
    // Constructor of an empty filtration
    Filtration(const _ComplexType& K) : _K(K)
    {
    }
    
    // Copy constructor
    Filtration(const Filtration& f) : _K(f._K), _filtration(f._filtration), _cell_to_t(f._cell_to_t) {}
    
    // Build lower (resp. upper)-star filtration
    // From the degrees of vertices
    void star_filtration(const std::vector<_DegType>& deg, bool lower = true) ;
    // From a function over vertices
    void star_filtration(std::function<_DegType(int)>& deg_fun, bool lower = true) ;
    
    // iterator (iterate a filtration)
    struct iterator
    {
        // Iterator tags
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = Filtration_iter_type;
        
        // Iterator constructors
        iterator(const Filtration& f, int i=0) : _i(i), _f(f) {}
        
        // Operators
        value_type operator*() const
        {
            Filtration_iter_type res ;
            res.cell_dim = _f._filtration.at(_i) ;
            res.degree = _f._deg.at(_i) ;
            return res ;
        }
        
        int time () { return _i ; }
        CellDim cell_dim () { return _f._filtration.at(_i); }
        _DegType degree () { return _f._deg.at(_i); }
        
        // Prefix increment
        // Find next hole of positive duration
        iterator& operator++()
        {
            ++_i;
            return *this;
        }
        
        // Postfix increment
        iterator operator++(int) { iterator tmp = *this; ++(*this); return tmp; }
        
        friend bool operator== (const iterator& a, const iterator& b) { return a._i == b._i; };
        friend bool operator!= (const iterator& a, const iterator& b) { return a._i != b._i; };
        
    private:
        int _i ; // Index along _persist
        const Filtration& _f ; // Filtration iterated
    };
    iterator begin() { return iterator(*this, 0) ; }
    iterator end() { return iterator(*this, _filtration.size()) ; }
    
    // getters
    CellDim get_cell_dim (int i) { return _filtration.at(i); }
    _DegType get_degree (int i) { return _deg.at(i); }
    
    // Output filtration
    friend ostream & operator<<(ostream & out, const Filtration &f)
    {
        const int N(f._filtration.size()) ;
        for (int i=0; i<N; ++i)
            // Filtration
        {
            out << i << " -> (" << f._filtration.at(i).first << "- dim " << f._filtration.at(i).second << " , " << f._deg.at(i) << ") " << std::endl ;
        }
        return out ;
    }
    
protected:
    void build_filtration_structure() ;
    bool check_filtration() ;
    
    template <typename CoefT, typename ComplexT, typename DegT>
    friend class PersistentHDVF ;
};

template <typename _CoefType, typename _ComplexType, typename _DegType>
void Filtration<_CoefType, _ComplexType, _DegType>::build_filtration_structure()
{
    for (int i = 0; i<_filtration.size(); ++i)
    {
        const CellDim c(_filtration.at(i)) ;
        // c : filtration index i, index in the basis reordered by filtration : j
        _cell_to_t[c] = i ;
    }
}

template <typename _CoefType, typename _ComplexType, typename _DegType>
bool Filtration<_CoefType, _ComplexType, _DegType>::check_filtration()
{
    bool valid = true ;
    for (int i=0; i<_filtration.size() && valid; ++i)
    {
        if (i>0)
            valid = valid & (_deg.at(i) >= _deg.at(i-1)) ;
        CellDim c(_filtration.at(i)) ;
        cout << i << " -> " << c.first << " dim "<< c.second << endl ;
        const int q = c.second ;
        if (q>0)
        {
            CChain dc = _K.d(c.first, q) ;
            cout << "bnd : " << dc << endl ;
            for (typename CChain::iterator it = dc.begin(); it != dc.end() && valid; ++it)
            {
                // Faces of c
                const CellDim face(it->first,q-1) ;
                // Check if the face c belongs to the filtration
                auto it_face(_cell_to_t.find(face)) ;
                valid = valid & (it_face != _cell_to_t.end()) ;
                if (!valid)
                    cout << "face not found" << endl ;
                if (it_face != _cell_to_t.end())
                    valid = valid & (_cell_to_t[face] < i) ;
                if (!valid)
                    cout << "face " << it->first << " at time : " << _cell_to_t[face] << " with i : " << i << endl ;
            }
        }
        if (!valid)
            cout << "check failed : " << i << endl ;
    }
    return valid ;
}

template <typename _CoefType, typename _ComplexType, typename _DegType>
void Filtration<_CoefType, _ComplexType, _DegType>::star_filtration(const std::vector<_DegType> &deg, bool lower)
{
    if (deg.size() != _K.nb_cells(0))
        throw "Star filtration error : deg should provide one value by vertex" ;
    
    // Create filtration and degrees for all cells according to deg
    // -> lower star: maximum degree of vertices
    // -> upper star: minimum degree of vertices
    std::vector<CellDim> tmp_filtration ;
    std::vector<_DegType> tmp_deg ;
    std::vector<int> tmp_perm ;
    // Init tmp_perm (for vertices)
    //    cout << "vertices degrees : " << endl ;
    for (int i=0; i<deg.size(); ++i)
    {
        tmp_perm.push_back(i) ;
        tmp_filtration.push_back(CellDim(i,0)) ;
        tmp_deg.push_back(deg.at(i)) ;
        //        cout << i << " -> " << deg.at(i) << endl ;
    }
    // For all other cells
    for (int q=1; q<=_K.dim(); ++q)
    {
        for (int i=0; i<_K.nb_cells(q); ++i)
        {
            tmp_filtration.push_back(CellDim(i,q)) ;
            
            // Compute corresponding degree
            // Vertices of the cell
            // TODO : bottom_faces should produce a set
            std::vector<int> verts(_K.bottom_faces(i,q)) ;
            //            cout << "cell " << i << " dim " << q << " : " << verts ;
            // Compute the degree of the cell
            _DegType d = deg.at(verts.at(0)) ;
            for (int j=1; j<verts.size(); ++j)
            {
                const _DegType tmp_d(deg.at(verts.at(j))) ;
                if (lower) // lower star : maximum of degrees
                {
                    if (tmp_d > d)
                        d = tmp_d ;
                }
                else // upper
                {
                    if (tmp_d < d)
                        d = tmp_d ;
                }
            }
            tmp_deg.push_back(d) ;
            tmp_perm.push_back(tmp_perm.size()) ;
            //            cout << " -> " << d << endl ;
        }
    }
    // Sort filtration
    
    // Create sorting function : lexicographic order over (deg, dim)
    // Test if cell i < cell j
    auto f_sort = [&tmp_filtration, &tmp_deg] (int i, int j)
    {
        // degree of i is lower than degree of j or (they are equal and the dimension of i is lower than the dimension of j)
        return ((tmp_deg[i] < tmp_deg[j]) || ((tmp_deg[i] == tmp_deg[j]) && (tmp_filtration[i].second < tmp_filtration[j].second))) ;
    } ;
    // Sort -> create the right permutation
    std::sort(tmp_perm.begin(), tmp_perm.end(), f_sort) ;
    // Insert cells in the filtration and degrees accordingly
    for (int i=0; i<tmp_perm.size(); ++i)
    {
        _filtration.push_back(tmp_filtration.at(tmp_perm.at(i))) ;
        _deg.push_back(tmp_deg.at(tmp_perm.at(i))) ;
        _cell_to_t[tmp_filtration.at(tmp_perm.at(i))] = i ;
    }
    
    // Check filtration
    //    cout << "check_filtration : " << check_filtration() << endl ;
}

template <typename _CoefType, typename _ComplexType, typename _DegType>
void Filtration<_CoefType, _ComplexType, _DegType>::star_filtration(std::function<_DegType(int)>& deg_fun, bool lower)
{
    std::vector<_DegType> deg ;
    for (int i=0; i<_K.nb_cells(0); ++i)
    {
        deg.push_back(deg_fun(i)) ;
    }
    star_filtration(deg) ;
}

/* Standard functions for lower star filtration */
/* DegType = double */

//std::function<double(int)> deg_fun_x = [&complex](int i)
//{
//    std::vector<double> Xi(complex.get_vertex_coords(i)) ;
//    return (Xi.at(0)) ;
//} ;
//
//std::function<double(int)> deg_fun_z = [&complex](int i)
//{
//    std::vector<double> Xi(complex.get_vertex_coords(i)) ;
//    return (Xi.at(2)) ;
//} ;

std::function<double(const std::vector<double>&)> f_x = [](const std::vector<double>& v)
{
    return (v.at(0)) ;
} ;

std::function<double(const std::vector<double>&)> f_y = [](const std::vector<double>& v)
{
    return (v.at(1)) ;
} ;

std::function<double(const std::vector<double>&)> f_z = [](const std::vector<double>& v)
{
    return (v.at(2)) ;
} ;

template<typename ComplexType>
std::function<double(int)>  deg_fun (const ComplexType& complex, std::function<double(const vector<double>&)>& f)
{
    std::function<double(int)> deg_fun_f = [&complex, &f](int i)
    {
        const std::vector<double> Xi(complex.get_vertex_coords(i)) ;
        //        cout << "deg_fun_f, i: " << i ;
        //        cout << " Xi : " ;
        //        for (double c : Xi)
        //            cout << c << " " ;
        //        cout << " - degree: " << f(Xi) << endl ;
        return f(Xi) ;
    } ;
    return deg_fun_f ;
}


/* ************************************************************************ */

/**
 * \class Persistent HDVF
 * \brief Implementation of persistent HDVF (inherits HDVF).
 *
 * The PersistentHDVF class contains all functions to build persistent HDVF and create perfect persistent HDVF.
 *
 * \tparam _CoefficientType The chain's coefficient types (default is OSM::ZCoefficient)
 * \tparam _ComplexType The type of complex  (default is SimpComplex)
 *
 * \author Bac A.
 * \version 0.1.0
 * \date 22/08/2024
 */

// Persistent interval types
typedef std::pair<int, int> PerInterval ;
typedef std::pair<int, int> CellDim ;
typedef std::pair<CellDim, CellDim> PerIntervalCells ;

template <typename _DegType>
using PerDegIntervalT =  std::pair<_DegType,_DegType> ;

template <typename _DegType>
using PerHoleT = std::tuple<PerInterval, PerIntervalCells, PerDegIntervalT<_DegType> > ;

template <typename _DegType>
ostream& operator<< (ostream& out, const PerHoleT<_DegType>& hole)
{
    // time (cell, dim) -> time (cell, dim) / time to live
    const PerInterval per_int(std::get<0>(hole)) ;
    const PerIntervalCells per_int_cells(std::get<1>(hole)) ;
    const PerDegIntervalT<_DegType> per_int_deg(std::get<2>(hole)) ;
    
    //    const _DegType duration((per_int_deg.second>=0)?per_int_deg.second-per_int_deg.first:per_int_deg.second) ;
    const _DegType duration(per_int_deg.second-per_int_deg.first) ;
    if (duration >= 0) // finite interval
    {
        out << "[" << per_int.first << " (" << per_int_cells.first.first << ", " << per_int_cells.first.second << ") -> " ;
        out << per_int.second << " (" << per_int_cells.second.first << ", " << per_int_cells.second.second << ") / duration: " << duration << "]" << endl ;
    }
    else
    {
        out << "[" << per_int.first << " (" << per_int_cells.first.first << ", " << per_int_cells.first.second << ") -> " ;
        out << "inf]" << endl ;
    }
    return out ;
}


template<typename _CoefficientType, typename _ComplexType, typename _DegType>
class PersistentHDVF : public HDVF_core<_CoefficientType,_ComplexType, OSM::Chain, OSM::SubSparseMatrix>
{
    // Matrices types
    typedef OSM::Chain<_CoefficientType, OSM::COLUMN> CChain;
    typedef OSM::Chain<_CoefficientType, OSM::ROW> RChain;
    typedef OSM::SubSparseMatrix<_CoefficientType, OSM::COLUMN> CMatrix;
    typedef OSM::SubSparseMatrix<_CoefficientType, OSM::ROW> RMatrix;
    
    // HDVF type
    typedef HDVF_core<_CoefficientType, _ComplexType, OSM::Chain, OSM::SubSparseMatrix> HDVF_type ;
    
    // Persistence diagram types
    typedef PerDegIntervalT<_DegType> PerDegInterval;
    typedef PerHoleT<_DegType> PerHole;
    
public:
    // Export types
    typedef std::vector<std::vector<int> > ExpLabels ;
    typedef CChain ExpChain ;
    
    // For persistent diagram iterator
    typedef struct zexp_infos {
        PerHole hole ;
        ExpLabels labelsPSC ;
        ExpChain g_chain_sigma, g_chain_tau, fstar_chain_sigma, fstar_chain_tau ;
    } Exp_infos ;
    
protected:
    const Filtration<_CoefficientType, _ComplexType, _DegType> &_f ;
    // Permutation between K indices and persistent indices (order given by the filtration)
    std::vector<std::vector<int> > _K_to_per, _per_to_K ;
    // Persistent pairs
    std::vector<PerHole> _persist ;
    // Export data
    bool _with_export ;
    std::vector<ExpLabels> _export_labels ;
    std::vector<std::pair<ExpChain, ExpChain> > _export_g, _export_fstar ;
    
    // Current "time" in the filtration
    int _t ;
    // Current "times" along each dimension
    std::vector<int> _t_dim ;
    // Corresponding masks (for SubSparseMatrices)
    std::vector<OSM::Bitboard> _masks ;
    // Persistence computed
    bool _computation_over ;
    
public:
    /** \brief PersistentHDVF Constructor */
    PersistentHDVF(const _ComplexType& K, const Filtration<_CoefficientType, _ComplexType, _DegType>& f, int hdvf_opt = OPT_BND, bool with_export = false) ;
    
    // Persistence step
    PairCell findNextPair(bool &found) ;
    
    void stepPersist() ;
    
    void computePersistentHomology()
    {
        for (int i=0; i < _f._filtration.size(); ++i)
        {
            this->progress_bar(i, _f._filtration.size()) ;
            stepPersist() ;
        }
        
        // Compute "infinite" holes
        vector<vector<int> > criticals(this->get_flag(CRITICAL)) ;
        for (int q=0; q < criticals.size(); ++q)
        {
            for (int i : criticals.at(q))
            {
                // i : persistence index
                const PairCell p = {i, -1, q} ;
                const int ki(_per_to_K.at(q).at(i)) ; // K index
                const CellDim c(ki,q) ;
                const int ti(_f._cell_to_t.at(c)) ;
                const _DegType di(_f._deg.at(i)) ;
                PerInterval per_int(ti,ti-1) ;
                CellDim inf(-1,q+1) ;
                PerIntervalCells per_int_cell(c,inf) ;
                PerDegInterval per_deg_int(di,di-1) ;
                PerHole hole(per_int, per_int_cell, per_deg_int) ;
                _persist.push_back(hole) ;
                
                // If export is on, store export data
                if (_with_export)
                    export_perHDVF(p) ;
            }
        }
        _computation_over = true ;
    }
    
    bool with_export () { return _with_export ; }
    
    /* ********** */
    
    _DegType hole_duration (const PerHole hole) const
    {
        const PerDegInterval per_int_deg(std::get<2>(hole)) ;
        return per_int_deg.second - per_int_deg.first ;
    }
    
    friend ostream& operator<< (ostream& out, const PersistentHDVF& per_hdvf)
    {
        int i = 0 ;
        for (PerHole hole : per_hdvf._persist)
        {
            if (abs(per_hdvf.hole_duration(hole)) > 0)
                out << i << " --- duration : " << per_hdvf.hole_duration(hole) << " -- " << hole << endl ;
            ++i ;
            //            // time (cell, dim) -> time (cell, dim) / time to live
            //            const PerInterval per_int(std::get<0>(hole)) ;
            //            const PerIntervalCells per_int_cells(std::get<1>(hole)) ;
            //            const PerDegInterval per_int_deg(std::get<2>(hole)) ;
            //
            //            const _DegType duration((per_int_deg.second>=0)?per_int_deg.second-per_int_deg.first:-per_int_deg.second) ;
            //            if (duration >= 0) // finite interval
            //            {
            //                    out << "[" << per_int.first << " (" << per_int_cells.first.first << ", " << per_int_cells.first.second << ") -> " ;
            //                    out << per_int.second << " (" << per_int_cells.second.first << ", " << per_int_cells.second.second << ") / duration: " << duration << "]" << endl ;
            //            }
            //            else
            //            {
            //                out << "[" << per_int.first << " (" << per_int_cells.first.first << ", " << per_int_cells.first.second << ") -> " ;
            //                out << "inf]" << endl ;
            //            }
        }
        return out ;
    }
    
    ostream& print_perHDVFInfo (ostream& out)
    {
        out << "Filtration: " << _f << endl ;
        out << "_K_to_per and _per_to_K" << endl ;
        for (int q=0; q<=this->_K.dim(); ++q)
        {
            out << "-> dim " << q << endl ;
            out << "index_per -(_per_to_K)-> index_K -(_K_to_per)-> index_per" << endl ;
            for (int i=0; i<this->_K.nb_cells(q); ++i)
            {
                const int id_K(_per_to_K.at(q).at(i)) ;
                out << i << " -> " << id_K << " -> " << _K_to_per.at(q).at(id_K) << endl ;
            }
        }
        return out ;
    }
    
    // Method to generate labels for visualisation
    // For PSC, no additional arguments, for FSTAR and G, specify the critical cell concerned (with persistent index)
    vector<vector<int> > export_label (ExportType type, int cell=0, int dim=0) const
    {
        cout << "per export label" << endl ;
        vector<vector<int> > labels(this->_K.dim()+1) ;
        if (type == PSC)
        {
            for (int q=0; q<=this->_K.dim(); ++q)
            {
                for (int i = 0; i<this->_K.nb_cells(q); ++i) // _K indices
                {
                    // if index of i in the filtration is lower than time in dimension q -> display
                    if (_K_to_per.at(q).at(i) <= _t_dim.at(q))
                    {
                        if (this->_flag.at(q).at(i) == PRIMARY)
                            labels.at(q).push_back(-1) ;
                        else if (this->_flag.at(q).at(i) == SECONDARY)
                            labels.at(q).push_back(1) ;
                        else
                            labels.at(q).push_back(0) ;
                    }
                    else // i has not yet been reached
                        labels.at(q).push_back(2) ;
                }
            }
        }
        else
        {
            if ((type == FSTAR) && (this->_hdvf_opt & (OPT_FULL | OPT_F)))
            {
                for (int q=0; q<=this->_K.dim(); ++q)
                {
                    labels.at(q).resize(this->_K.nb_cells(q)) ;
                }
                
                if (dim < this->_K.dim())
                {
                    int cellK = _per_to_K.at(dim).at(cell) ; // index of the cell in K
                    
                    labels.at(dim).at(cellK) = 1 ;
                    // Cross f*
                    const typename HDVF_type::RChain& fstar_cell = OSM::cgetRow(this->_F_row.at(dim), cell) ;
                    for (typename HDVF_type::RChain::const_iterator it = fstar_cell.cbegin(); it != fstar_cell.cend(); ++it)
                    {
                        int fstarcellK = _per_to_K.at(dim).at(it->first) ; // index of the cell of f* in K
                        // Set the cofaces of idcellK in dimension dim+1
                        typename HDVF_type::RChain cofaces(this->_K.cod(fstarcellK,dim)) ;
                        for (typename HDVF_type::RChain::const_iterator it2 =  cofaces.cbegin(); it2 != cofaces.cend(); ++it2)
                        {
                            labels.at(dim+1).at(it2->first) = 1 ; // set the flag to 1 only for FSTAR(cell) cells
                        }
                    }
                }
                else
                    throw "Error: cannot export f* for a cell of maximal dimension" ;
            }
            else if (this->_hdvf_opt & (OPT_FULL | OPT_G)) // G
            {
                //                cout << "Export : _subCC : " << _subCC << endl ;
                for (int q=0; q<=this->_K.dim(); ++q)
                {
                    labels.at(q).resize(this->_K.nb_cells(q)) ;
                }
                
                int cellK = _per_to_K.at(dim).at(cell) ; // index of the cell in K
                labels.at(dim).at(cellK) = 1 ;
                // Cross g
                const typename HDVF_type::CChain& g_cell = OSM::cgetColumn(this->_G_col.at(dim), cell) ;
                for (typename HDVF_type::CChain::const_iterator it = g_cell.cbegin(); it != g_cell.cend(); ++it)
                {
                    // Index of the cell in K
                    int gcellK = _per_to_K.at(dim).at(it->first) ;
                    labels.at(dim).at(gcellK) = 1 ; // set the flag to 1 only for G(cell) cells
                }
            }
        }
        return labels ;
    }
    
    // Method to generate PSC labels for visualisation
    virtual vector<vector<int> > export_labelsPSC () const
    {
        vector<vector<int> > labels(this->_K.dim()+1) ;
        for (int q=0; q<=this->_K.dim(); ++q)
        {
            for (int i = 0; i<this->_K.nb_cells(q); ++i)
            {
                const int id_per(_K_to_per.at(q).at(i)) ;
                if (id_per <=_t_dim.at(q))
                {
                    if (this->_flag.at(q).at(id_per) == PRIMARY)
                        labels.at(q).push_back(-1) ;
                    else if (this->_flag.at(q).at(id_per) == SECONDARY)
                        labels.at(q).push_back(1) ;
                    else if (this->_flag.at(q).at(id_per) == CRITICAL)
                        labels.at(q).push_back(0) ;
                    else // NONE
                        labels.at(q).push_back(2) ;
                }
                else
                    labels.at(q).push_back(2) ;
            }
        }
        return labels ;
    }
    
    // Method to export a G chain for visualisation
    // with cell indices in initial _K
    virtual CChain export_GChain (int cell, int dim) const
    {
        if ((dim<0) || (dim>this->_K.dim()))
            throw "Error : export_GChain with dim out of range" ;
        //        if (_K_to_per.at(dim).at(cell) > _t_dim.at(dim))
        //            throw "Error : export_GChain with 'future' cell wrt persistence" ;
        
        if (this->_hdvf_opt & (OPT_FULL | OPT_G))
        {
            // Get g(cell, dim) with per indices
            CChain g_cell(OSM::getColumn(this->_G_col.at(dim), cell)) ;
            // Add 1 to the cell
            g_cell[cell] = 1 ;
            // Compute the chain with _K indices
            CChain g_cell_K(g_cell.dimension()) ;
            for (typename CChain::const_iterator it = g_cell.begin(); it != g_cell.end(); ++it)
            {
                const int i(_per_to_K.at(dim).at(it->first)) ;
                g_cell_K[i] = it->second ;
            }
            return g_cell_K ;
        }
        else
            throw "Error : trying to export g_chain without proper HDVF option" ;
    }
    
    // Method to export a FSTAR chain for visualisation
    // with cell indices in initial _K
    virtual CChain export_FSTARChain (int cell, int dim) const
    {
        if ((dim<0) || (dim>this->_K.dim()))
            throw "Error : export_GChain with dim out of range" ;
        //        if (_K_to_per.at(dim).at(cell) > _t_dim.at(dim))
        //            throw "Error : export_GChain with 'future' cell wrt persistence" ;
        
        if (this->_hdvf_opt & (OPT_FULL | OPT_F))
        {
            // Get fstar(cell, dim) with per indices
            RChain fstar_cell(OSM::getRow(this->_F_row.at(dim), cell)) ;
            // Add 1 to the cell
            fstar_cell[cell] = 1 ;
            // Compute the cofaces of the chain with _K indices
            // Compute the cofaces
            if (dim < this->_K.dim())
            {
                CChain fstar_cofaces(this->_K.nb_cells(dim+1)) ;
                for (typename RChain::const_iterator it = fstar_cell.begin(); it != fstar_cell.end(); ++it)
                {
                    // Set the cofaces of indices_K(it->first) in dimension dim+1
                    // belonging to _K(_t)
                    const int i(_per_to_K.at(dim).at(it->first)) ;
                    RChain cofaces(this->_K.cod(i,dim)) ;
                    for (typename RChain::const_iterator it2 =  cofaces.cbegin(); it2 != cofaces.cend(); ++it2)
                    {
                        const int id(it2->first) ;
                        if (_K_to_per.at(dim+1).at(id) <=_t_dim.at(dim+1))
                            fstar_cofaces[id] = 1 ;
                    }
                }
                return fstar_cofaces ;
            }
            else
                return CChain(0) ;
        }
        else
            throw "Error : trying to export g_chain without proper HDVF option" ;
    }
    
    // Export method (at a given time, for a pair of critical cells - before pairing)
    void export_perHDVF(PairCell p)
    {
        // Export labels
        ExpLabels labels(this->export_labelsPSC()) ;
        _export_labels.push_back(labels) ;
        // Export g (according to options)
        if (this->_hdvf_opt & (OPT_FULL | OPT_G))
        {
            CChain chain_sigma(export_GChain(p.sigma, p.dim)) ;
            CChain chain_tau ;
            if (p.tau >= 0)
                chain_tau = export_GChain(p.tau, p.dim+1) ;
            _export_g.push_back(std::pair<CChain,CChain>(chain_sigma, chain_tau)) ;
        }
        // Export fstar (according to options)
        if (this->_hdvf_opt & (OPT_FULL | OPT_F))
        {
            CChain chain_sigma(export_FSTARChain(p.sigma, p.dim)) ;
            CChain chain_tau ;
            if (p.tau >= 0)
                chain_tau = export_FSTARChain(p.tau, p.dim+1) ;
            _export_fstar.push_back(std::pair<CChain,CChain>(chain_sigma, chain_tau)) ;
        }
    }
    // Empty export method (for persistent intervals below the threshold)
    void export_perHDVF()
    {
        // Export labels
        ExpLabels labels ;
        _export_labels.push_back(labels) ;
        // Export g (according to options)
        if (this->_hdvf_opt & (OPT_FULL | OPT_G))
        {
            CChain chain_sigma, chain_tau ;
            _export_g.push_back(std::pair<CChain,CChain>(chain_sigma, chain_tau)) ;
        }
        // Export fstar (according to options)
        if (this->_hdvf_opt & (OPT_FULL | OPT_F))
        {
            CChain chain_sigma, chain_tau ;
            _export_fstar.push_back(std::pair<CChain,CChain>(chain_sigma, chain_tau)) ;
        }
    }
    
    struct iterator
    {
        // Iterator tags
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = Exp_infos;
        using pointer           = Exp_infos*;  // or also value_type*
        
        // Iterator constructors
        iterator(const PersistentHDVF& per_hdvf, int i=0) : _i(i), _per_hdvf(per_hdvf)
        {
            if (_i == 0)
            {
                // Iterate over holes of duration > 0
                while ((_i<_per_hdvf._persist.size()) && (_per_hdvf.hole_duration(_per_hdvf._persist.at(_i)) == 0))
                {
                    ++_i ;
                }
            }
        } ;
        
        // Operators
        value_type operator*() const
        {
            Exp_infos res ;
            res.hole = _per_hdvf._persist.at(_i) ;
            res.labelsPSC = _per_hdvf._export_labels.at(_i) ;
            if (_per_hdvf._hdvf_opt & (OPT_G | OPT_FULL))
            {
                res.g_chain_sigma = _per_hdvf._export_g.at(_i).first ;
                res.g_chain_tau = _per_hdvf._export_g.at(_i).second ;
            }
            if (_per_hdvf._hdvf_opt & (OPT_F | OPT_FULL))
            {
                res.fstar_chain_sigma = _per_hdvf._export_fstar.at(_i).first ;
                res.fstar_chain_tau = _per_hdvf._export_fstar.at(_i).second ;
            }
            return res ;
        }
        pointer operator->()
        {
            Exp_infos res ;
            res.hole = _per_hdvf._persist.at(_i) ;
            res.labelsPSC = _per_hdvf._export_labels.at(_i) ;
            if (_per_hdvf._hdvf_opt & (OPT_G | OPT_FULL))
            {
                res.g_chain_sigma = _per_hdvf._export_g.at(_i).first ;
                res.g_chain_tau = _per_hdvf._export_g.at(_i).second ;
            }
            if (_per_hdvf._hdvf_opt & (OPT_F | OPT_FULL))
            {
                res.fstar_chain_sigma = _per_hdvf._export_fstar.at(_i).first ;
                res.fstar_chain_tau = _per_hdvf._export_fstar.at(_i).second ;
            }
            return res ;
        }
        
        // Prefix increment
        // Find next hole of positive duration
        iterator& operator++()
        {
            ++_i;
            while ((_i<_per_hdvf._persist.size()) && (_per_hdvf.hole_duration(_per_hdvf._persist.at(_i)) == 0))
            {
                ++_i ;
            }
            return *this;
        }
        
        // Postfix increment
        iterator operator++(int) { iterator tmp = *this; ++(*this); return tmp; }
        
        friend bool operator== (const iterator& a, const iterator& b) { return a._i == b._i; };
        friend bool operator!= (const iterator& a, const iterator& b) { return a._i != b._i; };
        
    private:
        int _i ; // Index along _persist
        const PersistentHDVF& _per_hdvf ; // per_hdvf iterated
    };
    iterator begin() { return iterator(*this, 0) ; }
    iterator end() { return iterator(*this, _persist.size()) ; }    
} ;


template<typename _CoefficientType, typename _ComplexType, typename _DegType>
PersistentHDVF<_CoefficientType, _ComplexType, _DegType>::PersistentHDVF(const _ComplexType& K, const Filtration<_CoefficientType, _ComplexType, _DegType>& f, int hdvf_opt, bool with_export) : HDVF_core<_CoefficientType,_ComplexType, OSM::Chain, OSM::SubSparseMatrix>(K,hdvf_opt), _f(f), _with_export(with_export), _t(-1), _computation_over(false)
{
    // Initialisation of _t_dim, _K_to_per and _per_to_K
    _t_dim.resize(this->_K.dim()+1, 0) ;
    _K_to_per.resize(this->_K.dim()+1) ;
    _per_to_K.resize(this->_K.dim()+1) ;
    for (int q=0; q<=this->_K.dim(); ++q)
    {
        _K_to_per.at(q).resize(this->_K.nb_cells(q)) ;
        _t_dim.at(q) = -1 ;
    }
    
    for (int i = 0; i<_f._filtration.size(); ++i)
    {
        const CellDim c(_f._filtration.at(i));
        const int q(c.second) ;
        const int ind_K_i(c.first) ;
        const int ind_per_i(_per_to_K.at(q).size()) ;
        _per_to_K.at(q).push_back(ind_K_i) ;
        _K_to_per.at(q).at(ind_K_i) = ind_per_i ;
    }
    
    // Init _masks
    _masks.resize(this->_K.dim()+1) ;
    for (int q=0; q<this->_K.dim()+1; ++q)
    {
        _masks.at(q) = OSM::Bitboard(this->_K.nb_cells(q)) ;
    }
    
    // Init boundary matrices
    vector<CMatrix> _DD_per(this->_K.dim()+1) ;
    
    // Copy _DD_col with filtration order (for dimensions q>0)
    for (int q = 0 ; q <= this->_K.dim(); ++q)
    {
        const std::pair<int, int> s(this->_DD_col.at(q).dimensions()) ;
        _DD_per.at(q) = CMatrix(s.first, s.second) ;
    }
    // Set empty mask for _DD_per[0]
    _DD_per.at(0).complement();
    
    for (int q = 1 ; q <= this->_K.dim(); ++q)
    {
        // Cross _DD_col.at(q) and set _DD_per.at(q) coefficients on the fly
        for (OSM::Bitboard::iterator it_col = this->_DD_col.at(q).begin(); it_col != this->_DD_col.at(q).end(); ++it_col)
        {
            const int j(*it_col) ;
            const CChain& col(OSM::cgetColumn(this->_DD_col.at(q), j)) ;
            for (typename CChain::const_iterator it = col.begin(); it != col.end(); ++it)
            {
                const int i(it->first) ;
                const _CoefficientType v(it->second) ;
                // Cells in the _K basis : i(dim q-1) / j(dim q)
                // Convert to indices in the persistent order
                const int pi(_K_to_per.at(q-1).at(i)), pj(_K_to_per.at(q).at(j)) ;
                _DD_per.at(q).set_coef(pi, pj, v) ;
            }
        }
    }
    this->_DD_col = _DD_per ;
    
    // Init _DD_col mask (empty for all cells)
    for (int q = 1 ; q <= this->_K.dim(); ++q)
        this->_DD_col.at(q).set_sub(_masks.at(q)) ;
}

template<typename _CoefficientType, typename _ComplexType, typename _DegType>
PairCell PersistentHDVF<_CoefficientType, _ComplexType, _DegType>::findNextPair(bool &found)
{
    PairCell p ;
    // Get current cell (in the basis K)
    CellDim c(_f._filtration.at(_t)) ;
    const int q(c.second), sigma(_K_to_per.at(q).at(c.first))  ;
    // Search for pairing
    found = false;
    
    if (q >= 1)
    {
        // Compute bounded max while iterating over the Chain
        const CChain& tmp2(OSM::cgetColumn(this->_DD_col.at(q), sigma)) ;
        std::size_t tmax = _t_dim.at(q-1) ;
        std::size_t i ;
        for (typename CChain::const_iterator it = tmp2.cbegin(); it != tmp2.cend(); ++it)
        {
            if ((it->first <= tmax) && (abs(it->second) == 1)) // possible pairing
            {
                if (!found) // for first cell met
                {
                    found = true ;
                    i = it->first ;
                }
                else // update the i with the maximum cell index (persistance)
                {
                    if (it->first > i)
                        i = it->first ;
                }
            }
        }
        if (found)
        {
            p.sigma = i ;
            p.tau = sigma ;
            p.dim = q-1 ;
        }
    }
    return p;
}

template<typename _CoefficientType, typename _ComplexType, typename _DegType>
void PersistentHDVF<_CoefficientType, _ComplexType, _DegType>::stepPersist()
{
    // Compute next persistent pair
    
    ++_t ; // Step forward in the filtration
    const int q_current(_f._filtration.at(_t).second) ; // Get the dimension of the new current cell
    ++_t_dim.at(q_current) ; // Update time in the dimension of the current cell
    _masks.at(q_current).setOn(_t_dim.at(q_current)) ; // Update mask accordingly
    this->_DD_col.at(q_current).set_bitOn(_t_dim.at(q_current)) ; // Update _DD_col mask
    
    // Search for pairing
    bool found ;
    PairCell p(findNextPair(found)) ;
    if (found)
    {
        // Corresponding persistent interval
        const int q(p.dim) ;
        // indices of both cells in the _K basis
        const int ki(_per_to_K.at(q).at(p.sigma)), kj(_per_to_K.at(q+1).at(p.tau)) ;
        CellDim ci(ki, q), cj(kj, q+1) ; // cells of the interval - in the K basis
        int ti(_f._cell_to_t.at(ci)), tj(_f._cell_to_t.at(cj)) ; // times of the interval
        PerInterval interval(ti, tj) ;
        PerIntervalCells interval_cells(ci, cj) ;
        PerDegInterval interval_deg(_f._deg.at(ti), _f._deg.at(tj)) ;
        PerHole hole(interval, interval_cells, interval_deg) ;
        // Add this interval
        _persist.push_back(hole) ;
        
        // If export is on, store export data for significant persistant intervals
        if (_with_export)
        {
            if ((interval_deg.second-interval_deg.first)>0)
                export_perHDVF(p) ;
            else
                export_perHDVF() ;
        }
        
        //        if (interval_deg.second-interval_deg.first>0)
        //        {
        //            cout << "========> time : " << _t << endl ;
        //            cout << "pair found : " << p.sigma << ", " << p.tau << " - dim " << p.dim << endl ;
        //            this->print_reduction() ;
        //            
        //        }
        
        // Prepare for next step
        this->A(p.sigma, p.tau, p.dim) ; // Update the reduction
        //        if (interval_deg.second-interval_deg.first>0)
        //        {
        //            cout << "====> after A : " << endl ;
        //            this->print_reduction() ;
        //            
        //        }
    }
}

}

#endif // PER_HOM_HPP
