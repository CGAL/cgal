
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

#include <set>
#include <vector>
#include <utility>
#include <algorithm>
#include <iostream>

namespace CGAL {
namespace HDVF {

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Simplex` is used by the class `Abstract_simplicial_chain_complex` to implement the structure de simplex (i.e.\ cells of a simplicial complex).

 Simplices are described by the set of the indices of their vertices (see the documentation of `Abstract_simplicial_chain_complex` for examples).

 */

class Simplex {
    // Friend class: Abstract_simplicial_chain_complex
    template<typename _CoefficientType>
    friend class Abstract_simplicial_chain_complex ;

    private:
    // Indices of the simplex vertices
    std::set<size_t> _vertices;

    public:
    /** \brief Constant iterator over the vertices of a simplex. */
    typedef std::set<size_t>::const_iterator const_iterator ;

    /** \brief Method returning a constant iterator on the first vertex of the simplex. */
    const_iterator cbegin ()
    {
        return _vertices.cbegin() ;
    }
    /** \brief Method returning a constant iterator after the last vertex of the simplex. */
    const_iterator cend () { return _vertices.cend() ; }

    /** \brief Constructor from a set of vertices indices.
     *
     * \param[in] vertices Set of the simplex vertices indices (a simplex of dimension \f$q\f$ has \f$q+1\f$ vertices.
     */
    Simplex(const std::set<size_t>& vertices) : _vertices(vertices) {}

    /** \brief Dimension of a simplex.
     *
     * A simplex of dimension \f$q\f$ has \f$q+1\f$ vertices.
     */
    int dimension() const
    {
        return _vertices.size() - 1;
    }

    /** \brief Boundary of a simplex.
     *
     * The method returns the vector of the simplices in the boundary of the object.
     */
    std::vector<Simplex> boundary() const
    {
        std::vector<Simplex> result;
        result.reserve(_vertices.size());

        auto it = _vertices.begin();
        for (size_t i = 0; i < _vertices.size(); ++i) {
            std::set<size_t> simplex_vertices;
            std::copy_if(_vertices.begin(), _vertices.end(), std::inserter(simplex_vertices, simplex_vertices.begin()),
                         [it](const size_t& vertex) { return vertex != *it; });
            result.emplace_back(simplex_vertices);
            ++it;
        }

        return result;
    }

    /** \brief Get the set of vertices indices of the simplex. */
    const std::set<size_t>& get_vertices() const
    {
        return _vertices ;
    }

    /** \brief Comparison operator.
     *
     * Compare the object with another simplex according to the lexicographical order on vertices indices sets.
     *
     * \param[in] other Compare `this` with `other` (returns  `this < other`).
     */
    bool operator<(const Simplex& other) const {
        return _vertices < other._vertices;
    }

    /** \brief Output a simplex.
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

} /* end namespace HDVF */
} /* end namespace CGAL */

#endif // CGAL_HDVF_SIMPLEX_H
