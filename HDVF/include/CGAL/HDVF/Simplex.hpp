#ifndef SIMPLEX_HPP
#define SIMPLEX_HPP

#include <set>
#include <vector>
#include <utility> 
#include <algorithm> 
#include <iostream>

class Simplex {
    template<typename _CoefficientType>
    friend class Abstract_simplicial_complex ;
private:
    std::set<int> _vertices;

public:
    typedef std::set<int>::const_iterator const_iterator ;
    const_iterator cbegin () { return _vertices.cbegin() ; }
    const_iterator cend () { return _vertices.cend() ; }
    
    // Constructor
    Simplex(const std::set<int>& vertices) : _vertices(vertices) {}

    int dimension() const
    {
        return _vertices.size() - 1;
    }

    // Boundary operator
    std::vector<Simplex> boundary() const
    {
        std::vector<Simplex> result;
        result.reserve(_vertices.size());

        auto it = _vertices.begin();
        for (size_t i = 0; i < _vertices.size(); ++i) {
            std::set<int> simplex_vertices;
            std::copy_if(_vertices.begin(), _vertices.end(), std::inserter(simplex_vertices, simplex_vertices.begin()),
                         [it](const int& vertex) { return vertex != *it; });
            result.emplace_back(simplex_vertices);
            ++it;
        }

        return result;
    }

    // Method to get the vertices of the simplex (for testing)
    const std::set<int>& getVertices() const
    {
        return _vertices ;
    }

    // Comparison operator used for sorting simplices (especially in the Complex class)
    bool operator<(const Simplex& other) const {
        return _vertices < other._vertices;
    }

    // Operator <<
    friend std::ostream& operator<<(std::ostream& out, const Simplex& simplex)
    {
            out << "<";
            bool first = true;
            for (int vertex : simplex._vertices) {
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


#endif // SIMPLEX_HPP
