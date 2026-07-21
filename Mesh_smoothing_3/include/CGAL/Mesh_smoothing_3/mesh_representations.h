#pragma once

#include <vector>
#include <array>
#include <Eigen/Eigen>
#include "default_shapes.h"

#include <iostream>

namespace Mesh_smoothing_3 {

namespace default_structures {

    // Empty minimal mesh representation
    class Empty_mesh {
    public:
        using Cell_descriptor = std::size_t;
        using Vertex_descriptor = std::size_t;
        using Point_3 = Eigen::Vector3d;

        std::size_t nb_cells() const { return 0; }
        std::size_t nb_vertices() const { return 0; }

        Point_3 vertex_coordinates(Vertex_descriptor vertex) const { return {0.,0.,0.}; }
        void set_new_vertex_coordinates(Vertex_descriptor vertex, Point_3 coord) {}   // only non const

        std::vector<Cell_descriptor> cell_range() const { return {}; } // should return a range of Cell_descriptor
        std::array<Vertex_descriptor, 4> cell_vertices(Cell_descriptor cell) const { return {0,0,0,0}; } // can return anything of size 4 with [int] operator
        std::array<Point_3, 4> cell_reference_shape(Cell_descriptor cell) const {
            return Shapes::VTK_TETRAHEDRON<Point_3>();
        }
    };


    // Empty minimal boundary representation
    template<typename Vertex_descriptor>
    class Empty_boundary {
        public:
            using Face_descriptor = std::size_t;
            using Normal_3 = Eigen::Vector3d;
            using Surface_patch_index = unsigned;
            std::size_t nb_faces() const { return 0; }
            std::vector<Face_descriptor> face_range() const { return {}; }
            std::size_t nb_face_vertices(Face_descriptor face) const { return 0; }
            Surface_patch_index patch_id(Face_descriptor) const { return 0; }
            std::vector<Vertex_descriptor> face_vertices(Face_descriptor face) const { return {}; }
    };


    // empty minimal edge network representation
    template<typename Vertex_descriptor>
    class Empty_edge_network {
        public:
            using Edge_descriptor = std::size_t;
            using Curve_index = unsigned;
            std::size_t nb_edges() const { return 0; }
            std::vector<Edge_descriptor> edge_range() const { return {}; }
            Curve_index curve_id(Edge_descriptor) const { return 0; }
            Vertex_descriptor edge_vertex(Edge_descriptor edge, unsigned i) const { return Vertex_descriptor(); }
    };

}

namespace utils {

    // Allows iterating over a range of unsigned integers with a simple syntax (for (unsigned i : Contiguous_unsigned_range{0, n}) { ... })
    struct Contiguous_unsigned_range {
        std::size_t i, n;
        void operator++() { ++i; }
        bool operator!=(Contiguous_unsigned_range const & rhs) const { return i != rhs.i; }
        std::size_t operator*() const { return i; }
        auto begin() { return Contiguous_unsigned_range{0,n}; }
        auto end()   { return Contiguous_unsigned_range{n,n}; }
    };
}

namespace basic_structures {
    class Tetrahedral_mesh {
    public:
        using Cell_descriptor = std::size_t;
        using Vertex_descriptor = std::size_t;
        using Point_3 = Eigen::Vector3d;

        std::size_t nb_cells() const { return _tetrahedra.size(); }
        std::size_t nb_vertices() const { return _points.size(); }

        Point_3 vertex_coordinates(Vertex_descriptor vertex) const { return _points[vertex]; }
        void set_new_vertex_coordinates(Vertex_descriptor vertex, Point_3 coord) { _points[vertex] = coord; }   // only non const

        utils::Contiguous_unsigned_range cell_range() const { return utils::Contiguous_unsigned_range{0, nb_cells()}; }
        std::array<Vertex_descriptor, 4> cell_vertices(Cell_descriptor cell) const { return _tetrahedra[cell]; }
        std::array<Point_3, 4> cell_reference_shape(Cell_descriptor cell) const {
            return Shapes::VTK_TETRAHEDRON<Point_3>();
        }
    public:
        std::vector<Point_3> _points;
        std::vector<std::array<std::size_t, 4>> _tetrahedra;
    };

    class Triangle_boundary {
    public:
        using Face_descriptor = std::size_t;
        using Normal_3 = Eigen::Vector3d;
        using Surface_patch_index = unsigned;
        std::size_t nb_faces() const { return _triangles.size(); }
        utils::Contiguous_unsigned_range face_range() const { return utils::Contiguous_unsigned_range{0, nb_faces()}; }
        std::size_t nb_face_vertices(Face_descriptor) const { return 3; }
        Surface_patch_index patch_id(Face_descriptor) const { return 0; }
        auto face_vertices(Face_descriptor face) const { return _triangles[face]; }

    public:
        std::vector<std::array<std::size_t, 3>> _triangles;
    };

    class Simple_edge_network {
    public:
        using Edge_descriptor = std::size_t;
        using Curve_index = unsigned;
        std::size_t nb_edges() const { return _edge_vertices.size(); }
        utils::Contiguous_unsigned_range edge_range() const { return utils::Contiguous_unsigned_range{0, nb_edges()}; }
        Curve_index curve_id(Edge_descriptor edge) const { return _id[edge]; }
        std::size_t edge_vertex(Edge_descriptor edge, unsigned i) const { return _edge_vertices[edge][i]; }
    public:
        void add_edge(std::size_t v0, std::size_t v1, unsigned id = 0) {
            _edge_vertices.push_back({v0, v1});
            _id.push_back(id);
        }
        std::vector<std::array<std::size_t, 2>> _edge_vertices;
        std::vector<unsigned> _id;
    };
}


namespace helper_structures {

    // Templated structure for representing mixed-element meshes
    template<
        typename InputCellDescriptor,
        typename VertexDescriptor,
        typename Point3,
        typename InputCellRangeType
    >
    class Mixed_element_mesh {
    public:
        // you redefine these in your derived class and then you need to call assemble() before using the structure
        using Shape = Mesh_smoothing_3::Shapes::Base_element_shape_reference<Point3>;

        virtual std::size_t nb_vertices() const = 0;

        virtual Point3 vertex_coordinates(VertexDescriptor vertex) const = 0;
        virtual void set_new_vertex_coordinates(VertexDescriptor vertex, Point3 coord) = 0;

        virtual InputCellRangeType input_cell_range() const = 0;
        virtual Shape const * get_element_shape(InputCellDescriptor cell) const = 0; // you can return nullptr if you want to ignore the cell
        virtual VertexDescriptor get_cell_vertex(InputCellDescriptor cell, unsigned local_Vertex_descriptor) const = 0;

        bool has_reference_mesh = false;
        virtual Point3 get_ref_vertex_coordinates(VertexDescriptor vertex) const { return Point3(); } // redefine if has_reference_mesh == true


    public:
        using Cell_descriptor = std::size_t;
        using Vertex_descriptor = VertexDescriptor;
        using Point_3 = Point3;

        std::size_t nb_cells() const { return optimization_tet_2_input_element.size(); };

        utils::Contiguous_unsigned_range cell_range() const { return utils::Contiguous_unsigned_range{0, nb_cells()}; }
        std::array<Vertex_descriptor, 4> cell_vertices(Cell_descriptor cell) const {
            std::array<Vertex_descriptor, 4> sub_decomposition;
            auto input_element = optimization_tet_2_input_element[cell].first;
            std::size_t tet_number = optimization_tet_2_input_element[cell].second;
            for (std::size_t i = 0; i < 4; ++i) {
                sub_decomposition[i] = get_cell_vertex(input_element, get_element_local_vert(input_element, tet_number, i));
            }
            return sub_decomposition;
        }
        std::array<Point_3, 4> cell_reference_shape(Cell_descriptor cell) const {
            return get_element_ref_shape(optimization_tet_2_input_element[cell].first, optimization_tet_2_input_element[cell].second);
        }

    public:
        unsigned get_nb_inner_tetrahedra(InputCellDescriptor cell) const {
            if (_element_shape.at(cell) == nullptr) return 0;
            return _element_shape.at(cell)->nb_inner_tetrahedra();
        };

        unsigned get_element_local_vert(InputCellDescriptor cell, unsigned tet, unsigned tet_vert) const {
            if (_element_shape.at(cell) == nullptr) return 0;
            return _element_shape.at(cell)->inner_tetrahedra_local_vert(tet, tet_vert);
        };

        std::array<Point_3, 4> get_element_ref_shape(InputCellDescriptor cell, unsigned tet) const {
            if (_element_shape.at(cell) == nullptr) return {Point_3(), Point_3(), Point_3(), Point_3()};
            if (!has_reference_mesh) {
                return _element_shape.at(cell)->inner_tetrahedra_reference_shape(tet);
            }
            else {
                return {
                    get_ref_vertex_coordinates(get_cell_vertex(cell, get_element_local_vert(cell, tet, 0))),
                    get_ref_vertex_coordinates(get_cell_vertex(cell, get_element_local_vert(cell, tet, 1))),
                    get_ref_vertex_coordinates(get_cell_vertex(cell, get_element_local_vert(cell, tet, 2))),
                    get_ref_vertex_coordinates(get_cell_vertex(cell, get_element_local_vert(cell, tet, 3)))
                };
            }
        };

        std::vector<std::pair<InputCellDescriptor, unsigned>> optimization_tet_2_input_element;
        std::unordered_map<InputCellDescriptor, Shape const *> _element_shape;

        void assemble() {
            for (auto c : input_cell_range()) {
                _element_shape.emplace(c, get_element_shape(c));
            }
            optimization_tet_2_input_element.clear();
            for (auto cell_descriptor : input_cell_range()) {
                for (unsigned i = 0; i < get_nb_inner_tetrahedra(cell_descriptor); ++i) {
                    optimization_tet_2_input_element.push_back({cell_descriptor, i});
                }
            }
        }
    };

    // Templated structure for polygonal boundary representation
    template<
        typename VertexDescriptor = std::size_t,
        typename FaceDescriptor = std::size_t,
        typename NormalType = Eigen::Vector3d
    >
    class Polygonal_boundary {
    public:
        using Face_descriptor = FaceDescriptor;
        using Normal_3 = NormalType;
        using Vertex_descriptor = VertexDescriptor;
        using Surface_patch_index = unsigned;
        std::size_t nb_faces() const { return _face_vertices.size(); }
        utils::Contiguous_unsigned_range face_range() const { return utils::Contiguous_unsigned_range{0, nb_faces()}; }
        std::size_t nb_face_vertices(Face_descriptor face) const { return _face_vertices[face].size(); }
        Surface_patch_index patch_id(Face_descriptor face) const { return _id[face]; }
        auto face_vertices(Face_descriptor face) const { return _face_vertices[face]; }

    public:
        void add_polygon(std::vector<Vertex_descriptor> const &polygon, unsigned id = 0) {
            _face_vertices.push_back(polygon);
            _id.push_back(id);
        }
        std::vector<std::vector<Vertex_descriptor>> _face_vertices;
        std::vector<unsigned> _id;
    };

}


// WORK IN PROGRESS
namespace cgal_types {
template <typename Triangulation_3>
class Triangulation_3_wrapper {
    using Cell_descriptor = Triangulation_3::Cell_handle;
    using Vertex_descriptor = Triangulation_3::Vertex_handle;
    using Point_3 = Triangulation_3::Geom_traits::Point_3;

    std::size_t nb_cells() const { return tr.number_of_finite_cells(); }
    std::size_t nb_vertices() const { return tr.number_of_vertices(); }

    decltype(auto) vertex_coordinates(Vertex_descriptor vertex) const { return tr.point(vertex); }
    void set_new_vertex_coordinates(Vertex_descriptor vertex, const Point_3& coord) {
        vertex->set_point(coord);
    }

    auto cell_range() const { return tr.finite_cell_handles(); }
    std::array<Vertex_descriptor, 4> cell_vertices(Cell_descriptor cell) const {
        std::array<Vertex_descriptor, 4> vertices;
        for (int i = 0; i < 4; ++i) {
            vertices[static_cast<unsigned>(i)] = cell->vertex(i);
        }
        return vertices;
    }
    std::array<Point_3, 4> cell_reference_shape(Cell_descriptor cell) const {
        return Shapes::VTK_TETRAHEDRON<Point_3>();
    }

    Triangulation_3 &tr;
};

template <typename C3t3>
class C3t3_wrapper {
public:
    using Cell_descriptor = C3t3::Cell_handle;
    using Vertex_descriptor = C3t3::Vertex_handle;
    using Face_descriptor = C3t3::Facet;
    using Edge_descriptor = C3t3::Edge;
    using Normal_3 = C3t3::Triangulation::Geom_traits::Vector_3;
    using Point_3 = C3t3::Triangulation::Geom_traits::Point_3;
    using Weighted_point_3 = C3t3::Triangulation::Geom_traits::Weighted_point_3;

    using Surface_patch_index = C3t3::Surface_patch_index;
    using Curve_index = C3t3::Curve_index;

    std::size_t nb_cells() const { return c3t3.number_of_cells(); }
    std::size_t nb_faces() const { return c3t3.number_of_edges(); }
    std::size_t nb_edges() const { return c3t3.number_of_facets(); }
    std::size_t nb_vertices() const { return c3t3.triangulation().number_of_vertices(); }

    decltype(auto) vertex_coordinates(Vertex_descriptor vertex) const {
        return c3t3.triangulation().point(vertex).point(); // c3t3 holds weighted points
    }
    void set_new_vertex_coordinates(Vertex_descriptor vertex, const Point_3& coord) { vertex->set_point(Weighted_point_3{coord}); }
    auto cell_range() const { return c3t3.cells_in_complex(); }
    std::array<Vertex_descriptor, 4> cell_vertices(Cell_descriptor cell) const {
        std::array<Vertex_descriptor, 4> vertices;
        for (int i = 0; i < 4; ++i) {
            vertices[static_cast<unsigned>(i)] = cell->vertex(i);
        }
        return vertices;
    }
    std::array<Point_3, 4> cell_reference_shape(Cell_descriptor cell) const {
        return Shapes::VTK_TETRAHEDRON<Point_3>();
    }

    auto face_range() const { return c3t3.facets_in_complex(); }
    std::size_t nb_face_vertices(Face_descriptor face) const { return 3; }
    Surface_patch_index patch_id(Face_descriptor face) const { return c3t3.surface_patch_index(face); }
    std::vector<Vertex_descriptor> face_vertices(Face_descriptor face) const {
        std::vector<Vertex_descriptor> vertices(3);
        for (int i = 1; i < 4; ++i) {
            vertices[static_cast<unsigned>(i-1)] = face.first->vertex((face.second + i)%4);
        }
        return vertices;
    }

    auto edge_range() const { return c3t3.edges_in_complex(); }
    Curve_index curve_id(Edge_descriptor edge) const { return c3t3.curve_index(edge); } 
    Vertex_descriptor edge_vertex(Edge_descriptor edge, unsigned i) const {
        assert(i < 2);
        return edge.first->vertex(i == 0 ? edge.second : edge.third);
    }

    C3t3 &c3t3;
};

}

}

