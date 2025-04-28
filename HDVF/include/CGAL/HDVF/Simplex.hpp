#ifndef SIMPLEX_HPP
#define SIMPLEX_HPP

#include <set>
#include <vector>
#include <utility> 
#include <algorithm> 
#include <iostream>

template <typename VertexIdType = int>
class Simplex {
    template<typename _CoefficientType, typename _VertexIdType>
    friend class Abstract_simplicial_complex ;
private:
    std::set<int> vertices; 
    bool inside;

public:
    typedef std::set<int>::const_iterator const_iterator ;
    const_iterator cbegin () { return vertices.cbegin() ; }
    const_iterator cend () { return vertices.cend() ; }
    
    // Constructor
    Simplex(const std::set<int>& vertices, bool inside = true);

    int dimension() const;

    // Boundary operator
    std::vector<Simplex> boundary() const;

    // Method to get the vertices of the simplex (for testing)
    const std::set<int>& getVertices() const;

    // Comparison operator used for sorting simplices (especially in the Complex class)
    bool operator<(const Simplex& other) const {
        return vertices < other.vertices;
    }

    // Operator <<
    friend std::ostream& operator<<(std::ostream& out, const Simplex& simplex); 



};

// Implementations

template <typename VertexIdType>
Simplex<VertexIdType>::Simplex(const std::set<int>& vertices, bool inside)
    : vertices(vertices), inside(inside) {}

template <typename VertexIdType>
int Simplex<VertexIdType>::dimension() const {
    return vertices.size() - 1;
}

template <typename VertexIdType>
std::vector<Simplex<VertexIdType> > Simplex<VertexIdType>::boundary() const {
    std::vector<Simplex> result;
    result.reserve(vertices.size());

    auto it = vertices.begin();
    for (size_t i = 0; i < vertices.size(); ++i) {
        std::set<int> simplex_vertices;
        std::copy_if(vertices.begin(), vertices.end(), std::inserter(simplex_vertices, simplex_vertices.begin()),
                     [it](const int& vertex) { return vertex != *it; });
        result.emplace_back(simplex_vertices, inside);
        ++it;
    }

    return result;
}

template <typename VertexIdType>
const std::set<int>& Simplex<VertexIdType>::getVertices() const {
//    return std::vector<int>(vertices.begin(), vertices.end());
    return vertices ;
}

template <typename VertexIdType>
std::ostream& operator<<(std::ostream& out, const Simplex<VertexIdType>& simplex)
{
        out << "<";
        bool first = true;
        for (int vertex : simplex.vertices) {
            if (!first) {
                
                out << " ";
            }
            out << vertex;
            first = false;
        }
        out << ">";
        return out;
}

#endif // SIMPLEX_HPP
