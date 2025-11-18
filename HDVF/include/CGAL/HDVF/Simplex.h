
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


#ifndef CGAL_HDVF_SIMPLEX_H
#define CGAL_HDVF_SIMPLEX_H

#include <CGAL/license/HDVF.h>

#include <set>
#include <vector>
#include <utility>
#include <algorithm>
#include <iostream>

namespace CGAL {
namespace Homological_discrete_vector_field {

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Simplex` is used by the class `Abstract_simplicial_chain_complex` to implement the structure de simplex (i.e.\ cells of a simplicial complex).

 Simplices are described by the *ordered vector* of the indices of their vertices (see the documentation of `Abstract_simplicial_chain_complex` for examples).

 */

class Simplex {
    // Friend class: Abstract_simplicial_chain_complex
    template<typename _CoefficientType>
    friend class Abstract_simplicial_chain_complex ;

public:
    /** \brief Type of simplices representation.
     *
     * Simplices are stored as a sorted vector of vertex indices (a simplex of dimension \f$q\f$ has \f$q+1\f$ vertices and the vector must be **sorted**).
     **/
    typedef std::vector<size_t> Simplex_data;
private:
    // Indices of the simplex vertices - sorted vector
    Simplex_data _vertices;

    public:
    /** \brief Constructor from a vector of vertex indices.
     *
     * \param vertices Vector of the simplex vertex indices (a simplex of dimension \f$q\f$ has \f$q+1\f$ vertices and the vector must be **sorted**).
     * \param sort_data If `true` vectors of vertex indices are sorted, if `false` indices are assumed to be sorted (faster).
     */
    Simplex(const Simplex_data& vertices, bool sort_data = false) : _vertices(vertices) {
        if (sort_data)
            std::sort(_vertices.begin(), _vertices.end());
    }

    /** \brief Constant iterator over the vertices of a simplex. */
    typedef Simplex_data::const_iterator const_iterator ;

    /** \brief Returning a constant iterator on the first vertex of the simplex. */
    const_iterator cbegin ()
    {
        return _vertices.cbegin() ;
    }
    /** \brief Returning a constant past-the-end iterator of the simplex. */
    const_iterator cend () { return _vertices.cend() ; }

    /** \brief Returns the dimension of a simplex.
     *
     * A simplex of dimension \f$q\f$ has \f$q+1\f$ vertices.
     */
    int dimension() const
    {
        return _vertices.size() - 1;
    }

    /** \brief Returns the boundary of a simplex.
     *
     * The method returns the vector of the simplices in the boundary of the object.
     */

    std::vector<Simplex> boundary() const
    {
        std::vector<Simplex> result;
        result.reserve(_vertices.size());

        auto it = _vertices.begin();
        for (size_t i = 0; i < _vertices.size(); ++i) {
            Simplex_data simplex_vertices;
            std::copy_if(_vertices.begin(), _vertices.end(), std::inserter(simplex_vertices, simplex_vertices.begin()),
                         [it](const size_t& vertex) { return vertex != *it; });
            result.emplace_back(simplex_vertices);
            ++it;
        }

        return result;
    }

    /** \brief Gets the set of vertex indices of the simplex. */
    const Simplex_data& vertices() const
    {
        return _vertices ;
    }

    /** \brief Comparison operator.
     *
     * Compare the object with another simplex according to the lexicographical order on vertex indices sets.
     *
     * \param other Compare `this` with `other` (returns  `this < other`).
     */
    bool operator<(const Simplex& other) const {
        return _vertices < other._vertices;
    }

    /** \brief Equality operator.
     *
     * As vertex indices must be store in increasing order, comparison just comes to compare the ordered vector of indices.
     */
    bool operator==(const Simplex& other) const {
        return _vertices == other._vertices;
    }


    /** \brief Outputs a simplex.
     */
    friend std::ostream& operator<<(std::ostream& out, const Simplex& simplex)
    {
        out << "<";
        bool first = true;
        for (size_t vertex : simplex._vertices) {
            if (!first) {

                out << " ";
            }
            out << vertex;
            first = false;
        }
        out << ">";
        return out;
    }
};

} /* end namespace Homological_discrete_vector_field */
} /* end namespace CGAL */

#endif // CGAL_HDVF_SIMPLEX_H
