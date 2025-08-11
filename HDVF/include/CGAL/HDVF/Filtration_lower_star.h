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

#ifndef CGAL_HDVF_FILTRATION_LOWER_STAR_H
#define CGAL_HDVF_FILTRATION_LOWER_STAR_H

#include <CGAL/license/HDVF.h>

#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include <functional>
#include <CGAL/HDVF/Filtration_core.h>

namespace CGAL {
namespace HDVF {

/* Standard functions for lower star filtration */
/* DegType = double */

//std::function<double(size_t)> deg_fun_x = [&complex](size_t i)
//{
//    std::vector<double> Xi(complex.get_vertex_coords(i)) ;
//    return (Xi.at(0)) ;
//} ;
//
//std::function<double(size_t)> deg_fun_z = [&complex](size_t i)
//{
//    std::vector<double> Xi(complex.get_vertex_coords(i)) ;
//    return (Xi.at(2)) ;
//} ;

/** \brief For lower star filtration along x: function mapping coordinates to x */
std::function<double(const std::vector<double>&)> f_x = [](const std::vector<double>& v)
{
    return (v.at(0)) ;
} ;

/** \brief For lower star filtration along y: function mapping coordinates to y */
std::function<double(const std::vector<double>&)> f_y = [](const std::vector<double>& v)
{
    return (v.at(1)) ;
} ;

/** \brief For lower star filtration along z: function mapping coordinates to z */
std::function<double(const std::vector<double>&)> f_z = [](const std::vector<double>& v)
{
    return (v.at(2)) ;
} ;

/** \brief Degree function from a coordinates to scalar map. */
template<typename ComplexType>
std::function<double(size_t)>  deg_fun (const ComplexType& complex, std::function<double(const std::vector<double>&)>& f)
{
    std::function<double(size_t)> deg_fun_f = [&complex, &f](size_t i)
    {
        const std::vector<double> Xi(complex.get_vertex_coords(i)) ;
        return f(Xi) ;
    } ;
    return deg_fun_f ;
}

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Filtration_lower_star` implements the lower star filtration on a given complex implementing the concept `AbstractChainComplex`.

 A filtration associated to a chain complex `K` associates to each cell of `K` a scalar value (called degree) such that the degree of a cell is larger than the degrees of its faces.

 Let `Degree_type` be a scalar type. The lower star filtration is a filtration obtained from a map \f$\mathrm{deg}\,:\, K_0 \to \mathrm{Degree_type}\f$ (where \f$K_0\f$ denotes the set of vertices of \f$K\f$) associating a degree to each vertex of the complex \f$K\f$.

 The map is extended to cells of any dimension by setting, for a cell \f$\sigma\f$:
 \f[\mathrm{deg}(\sigma) = \max_{\substack{v\in K_0\\v\text{ face of }\sigma}} \mathrm{deg}(v)
 \f]

 For geometric complexes, standard lower star filtrations are obtained by taking as a degree function the \f$x\f$, \f$y\f$ or \f$z\f$ coordinate of vertices. The image below illustrates such a filtration ((left) lower star filtration with a \f$z\f$ degree map on vertices, (right) lower star filtration with a \f$y\f$ degree map on vertices).

 <img src="lower_star_filtration_z.png" align="center" width=15%/>
 <img src="lower_star_filtration_y.png" align="center" width=15%/>

 The `Filtration_lower_star` class provides constructors taking as input:
 - either the vector of vertex degrees
 - or a function mapping each vertex to its degree.

 \cgalModels{Filtration}

 \tparam CoefficientType a model of the `Ring` concept (ring used for homology computation).
 \tparam ComplexType a model of the `AbstractChainComplex` concept (type of the underlying chain complex).
 \tparam Degree_type the scalar type of degrees.
 */

template <typename CoefficientType, typename ComplexType, typename Degree_type>
class Filtration_lower_star : public Filtration_core<CoefficientType, ComplexType, Degree_type>
{
private:
    /*! Type of parent filtration instance. */
    typedef Filtration_core<CoefficientType, ComplexType, Degree_type> FiltrationCoreT;
    /*! Type of cell identifier (cell index and dimension). */
    typedef FiltrationCoreT::Cell_index_dimension Cell_index_dimension;
    /*! Type of value returned by the iterator. */
    typedef FiltrationCoreT::Filtration_iter_value Filtration_iter_value;
public:
    /*!
     * \brief Copy constructor.
     *
     * Builds a filtration by copy from another.
     *
     * \param[in] f An initial lower star filtration.
     */
    Filtration_lower_star(const Filtration_lower_star& f) : FiltrationCoreT(f) {}

    /*! \brief Constructor from vertex degrees.
     *
     * The constructor computes all cells degrees as the minimum of the degrees of their vertices and sorts all the cells of the complex to fulfill the filtration ordering constraints.
     *
     * \param[in] K Constant reference to the underlying complex.
     * \param[in] deg Vector of vertex degrees.
     */
    Filtration_lower_star(const ComplexType& K, const std::vector<Degree_type>& deg) : FiltrationCoreT(K)
    {
        star_filtration(deg);
    }

    /*! \brief Constructor from a function mapping vertices to degrees.
     *
     * The constructor computes all cells degrees as the minimum of the degrees of their vertices (obtained through `deg_fun`) and sorts all the cells of the complex to fulfill the filtration ordering constraints.
     *
     * \param[in] K Constant reference to the underlying complex.
     * \param[in] deg_fun Function mapping vertices of `K` to their degree.
     */
    Filtration_lower_star(const ComplexType& K, std::function<Degree_type(size_t)>& deg_fun) : FiltrationCoreT(K)
    {
        star_filtration(deg_fun);
    }


protected:
    // Build lower-star filtration
    /*! \brief Function building the filtration from the vector of vertices degrees.
     */
    void star_filtration(const std::vector<Degree_type>& deg) ;

    /*! \brief Function building the filtration from a function mapping vertices to their degree.
     */
    void star_filtration(std::function<Degree_type(size_t)>& deg_fun) ;
};



template <typename CoefficientType, typename ComplexType, typename Degree_type>
void Filtration_lower_star<CoefficientType, ComplexType, Degree_type>::star_filtration(const std::vector<Degree_type> &deg)
{
    if (deg.size() != this->_K.number_of_cells(0))
        throw "Star filtration error : deg should provide one value by vertex" ;

    // Create filtration and degrees for all cells according to deg
    // -> lower star: maximum degree of vertices
    // -> upper star: minimum degree of vertices
    std::vector<Cell_index_dimension> tmp_filtration ;
    std::vector<Degree_type> tmp_deg ;
    std::vector<size_t> tmp_perm ;
    for (size_t i=0; i<deg.size(); ++i)
    {
        tmp_perm.push_back(i) ;
        tmp_filtration.push_back(Cell_index_dimension(i,0)) ;
        tmp_deg.push_back(deg.at(i)) ;
    }
    // For all other cells
    for (int q=1; q<=this->_K.dimension(); ++q)
    {
        for (size_t i=0; i<this->_K.number_of_cells(q); ++i)
        {
            tmp_filtration.push_back(Cell_index_dimension(i,q)) ;

            // Compute corresponding degree
            // Vertices of the cell
            std::vector<size_t> verts(this->_K.bottom_faces(i,q)) ;
            // Compute the degree of the cell
            Degree_type d = deg.at(verts.at(0)) ;
            for (size_t j=1; j<verts.size(); ++j)
            {
                const Degree_type tmp_d(deg.at(verts.at(j))) ;
                if (tmp_d > d)
                    d = tmp_d ;
                // If upper star filtration (for cohomology)
                //                    if (tmp_d < d)
                //                        d = tmp_d ;

            }
            tmp_deg.push_back(d) ;
            tmp_perm.push_back(tmp_perm.size()) ;
        }
    }
    // Sort filtration

    // Create sorting function : lexicographic order over (deg, dim)
    // Test if cell i < cell j
    auto f_sort = [&tmp_filtration, &tmp_deg] (size_t i, size_t j)
    {
        // degree of i is lower than degree of j or (they are equal and the dimension of i is lower than the dimension of j)
        return ((tmp_deg[i] < tmp_deg[j]) || ((tmp_deg[i] == tmp_deg[j]) && (tmp_filtration[i].second < tmp_filtration[j].second))) ;
    } ;
    // Sort -> create the right permutation
    std::sort(tmp_perm.begin(), tmp_perm.end(), f_sort) ;
    // Insert cells in the filtration and degrees accordingly
    for (size_t i=0; i<tmp_perm.size(); ++i)
    {
        this->_filtration.push_back(tmp_filtration.at(tmp_perm.at(i))) ;
        this->_deg.push_back(tmp_deg.at(tmp_perm.at(i))) ;
        this->_cell_to_t[tmp_filtration.at(tmp_perm.at(i))] = i ;
    }
}

template <typename CoefficientType, typename ComplexType, typename Degree_type>
void Filtration_lower_star<CoefficientType, ComplexType, Degree_type>::star_filtration(std::function<Degree_type(size_t)>& deg_fun)
{
    std::vector<Degree_type> deg ;
    for (size_t i=0; i<this->_K.number_of_cells(0); ++i)
    {
        deg.push_back(deg_fun(i)) ;
    }
    star_filtration(deg) ;
}

} /* end namespace HDVF */
} /* end namespace CGAL */

#endif // CGAL_HDVF_FILTRATION_LOWER_STAR_H
