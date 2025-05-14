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

#ifndef CGAL_FILTRATION_H
#define CGAL_FILTRATION_H

#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include <functional>

/**
 * \class Filtration_lower_star class
 * \brief Implementation of the filtration.
 *
 * \tparam CoefficientType The chain's coefficient types (default is OSM::ZCoefficient)
 * \tparam ComplexType The type of complex
 */

namespace CGAL {
namespace HDVF {

template <typename CoefficientType, typename ComplexType, typename DegreeType>
class Filtration_lower_star
{
public:
    typedef std::pair<int, int> CellDim ;
    typedef struct zFiltration_iter_type {
        CellDim cell_dim ;
        DegreeType degree ;
    } FiltrationIterValue ;
    
protected:
    const ComplexType& _K ;
    std::vector<CellDim> _filtration ;
    std::vector<DegreeType> _deg ;
    
    // Cell -> (filtration index, basis index according to the filtration)
    std::map<CellDim,int> _cell_to_t ;
    
    typedef OSM::SparseMatrix<CoefficientType,OSM::COLUMN> CMatrix ;
    typedef OSM::SparseMatrix<CoefficientType,OSM::ROW> RMatrix ;
    typedef OSM::Chain<CoefficientType,OSM::COLUMN> CChain ;
    typedef OSM::Chain<CoefficientType,OSM::ROW> RChain ;
public:
    // Constructor of an empty filtration
    Filtration_lower_star(const ComplexType& K) : _K(K)
    {
    }
    
    // Copy constructor
    Filtration_lower_star(const Filtration_lower_star& f) : _K(f._K), _filtration(f._filtration), _cell_to_t(f._cell_to_t) {}
    
    // Constructor from a filtration on all cells
    Filtration_lower_star(const ComplexType& K, const std::vector<CellDim>& filtration, const std::vector<DegreeType>& deg) :
    _K(K), _filtration(filtration), _deg(deg)
    {
        check_filtration() ;
        build_filtration_structure() ;
    }
    
    // Constructor from the vector of degrees of vertices
    Filtration_lower_star(const ComplexType& K, const std::vector<DegreeType>& deg) : _K(K)
    {
        star_filtration(deg);
    }
    
    // Constructor from a function over vertices
    Filtration_lower_star(const ComplexType& K, std::function<DegreeType(int)>& deg_fun) : _K(K)
    {
        star_filtration(deg_fun);
    }
    
    // Build lower-star filtration
    // From the degrees of vertices
    void star_filtration(const std::vector<DegreeType>& deg) ;
    // From a function over vertices
    void star_filtration(std::function<DegreeType(int)>& deg_fun) ;
    
    // iterator (iterate a filtration)
    struct iterator
    {
        // Iterator tags
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = FiltrationIterValue;
        
        // Iterator constructors
        iterator(const Filtration_lower_star& f, int i=0) : _i(i), _f(f) {}
        
        // Operators
        value_type operator*() const
        {
            FiltrationIterValue res ;
            res.cell_dim = _f._filtration.at(_i) ;
            res.degree = _f._deg.at(_i) ;
            return res ;
        }
        
        int time () const { return _i ; }
        CellDim cell_dim () const { return _f._filtration.at(_i); }
        DegreeType degree () const { return _f._deg.at(_i); }
        
        // Prefix increment
        // Find next hole of positive duration
        iterator& operator++()
        {
            ++_i;
            return *this;
        }
        
        // Postfix increment
        iterator operator++(int) { iterator tmp = *this; ++(*this); return tmp; }
        
        friend bool operator== (const iterator& a, const iterator& b)  { return a._i == b._i; };
        friend bool operator!= (const iterator& a, const iterator& b)  { return a._i != b._i; };
        
    private:
        int _i ; // Index along _persist
        const Filtration_lower_star& _f ; // Filtration_lower_star iterated
    };
    iterator begin() { return iterator(*this, 0) ; }
    iterator end() { return iterator(*this, _filtration.size()) ; }
    
    // getters
    int get_filtration_size () const { return _filtration.size();}
    CellDim get_cell_dim (int i) const { return _filtration.at(i); }
    DegreeType get_degree (int i) const { return _deg.at(i); }
    
    // Output filtration
    friend ostream & operator<<(ostream & out, const Filtration_lower_star &f)
    {
        const int N(f._filtration.size()) ;
        for (int i=0; i<N; ++i)
            // Filtration_lower_star
        {
            out << i << " -> (" << f._filtration.at(i).first << "- dim " << f._filtration.at(i).second << " , " << f._deg.at(i) << ") " << std::endl ;
        }
        return out ;
    }
    
    bool check_filtration() const ;
    
protected:
    void build_filtration_structure() ;
    
    template <typename CoefT, typename ComplexT, typename DegT, typename FiltrT>
    friend class Hdvf_persistence ;
};

template <typename CoefficientType, typename ComplexType, typename DegreeType>
void Filtration_lower_star<CoefficientType, ComplexType, DegreeType>::build_filtration_structure()
{
    for (int i = 0; i<_filtration.size(); ++i)
    {
        const CellDim c(_filtration.at(i)) ;
        // c : filtration index i, index in the basis reordered by filtration : j
        _cell_to_t[c] = i ;
    }
}

template <typename CoefficientType, typename ComplexType, typename DegreeType>
bool Filtration_lower_star<CoefficientType, ComplexType, DegreeType>::check_filtration() const
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

template <typename CoefficientType, typename ComplexType, typename DegreeType>
void Filtration_lower_star<CoefficientType, ComplexType, DegreeType>::star_filtration(const std::vector<DegreeType> &deg)
{
    if (deg.size() != _K.nb_cells(0))
        throw "Star filtration error : deg should provide one value by vertex" ;
    
    // Create filtration and degrees for all cells according to deg
    // -> lower star: maximum degree of vertices
    // -> upper star: minimum degree of vertices
    std::vector<CellDim> tmp_filtration ;
    std::vector<DegreeType> tmp_deg ;
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
            DegreeType d = deg.at(verts.at(0)) ;
            for (int j=1; j<verts.size(); ++j)
            {
                const DegreeType tmp_d(deg.at(verts.at(j))) ;
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

template <typename CoefficientType, typename ComplexType, typename DegreeType>
void Filtration_lower_star<CoefficientType, ComplexType, DegreeType>::star_filtration(std::function<DegreeType(int)>& deg_fun)
{
    std::vector<DegreeType> deg ;
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

} /* end namespace HDVF */
} /* end namespace CGAL */

#endif // CGAL_FILTRATION_H
