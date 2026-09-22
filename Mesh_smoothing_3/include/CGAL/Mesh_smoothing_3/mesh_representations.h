// Copyright (c) 2026  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : François Protais

#ifndef CGAL_MESH_SMOOTHING_3_MESH_REPRESENTATIONS_H
#define CGAL_MESH_SMOOTHING_3_MESH_REPRESENTATIONS_H

#include <CGAL/license/Mesh_smoothing_3.h>

#include <CGAL/Mesh_smoothing_3/default_shapes.h>
#include <CGAL/Mesh_smoothing_3/internal/type_definitions.h>

#include <CGAL/Element_topo.h>

#include <Eigen/Eigen>

#include <vector>
#include <array>
#include <iostream>

namespace CGAL {

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

        Point_3 vertex_coordinates(Vertex_descriptor) const { return {0.,0.,0.}; }
        void  set_vertex_coordinates(Vertex_descriptor, Point_3) {}   // only non const

        std::vector<Cell_descriptor> cell_range() const { return {}; } // should return a range of Cell_descriptor
        std::array<Vertex_descriptor, 4> cell_vertices(Cell_descriptor) const { return {0,0,0,0}; } // can return anything of size 4 with [int] operator
        std::array<Point_3, 4> cell_reference_shape(Cell_descriptor) const {
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
            std::size_t nb_face_vertices(Face_descriptor) const { return 0; }
            Surface_patch_index patch_id(Face_descriptor) const { return 0; }
            std::vector<Vertex_descriptor> face_vertices(Face_descriptor) const { return {}; }
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
            Vertex_descriptor edge_vertex(Edge_descriptor, unsigned) const { return Vertex_descriptor(); }
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
        void  set_vertex_coordinates(Vertex_descriptor vertex, Point_3 coord) { _points[vertex] = coord; }   // only non const

        utils::Contiguous_unsigned_range cell_range() const { return utils::Contiguous_unsigned_range{0, nb_cells()}; }
        std::array<Vertex_descriptor, 4> cell_vertices(Cell_descriptor cell) const { return _tetrahedra[cell]; }
        std::array<Point_3, 4> cell_reference_shape(Cell_descriptor) const {
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

// concept without being a concept
struct Example_mixed_mesh {
    using Cell_descriptor = std::size_t;
    using Vertex_descriptor = std::size_t;
    using Point_3 = Eigen::Vector3d;
    using Input_cell_descriptor = std::size_t;
    using Element_shape_type = unsigned;

    using Shape = Mesh_smoothing_3::Shapes::Base_element_shape_reference<Point_3>;

    std::size_t nb_vertices() const { return 0; }
    std::size_t nb_input_cells() const { return 0; }

    Point_3 vertex_coordinates(Vertex_descriptor /*vertex*/) const { return Point_3(); }
    void  set_vertex_coordinates(Vertex_descriptor /*vertex*/, Point_3 /*coord*/) {}

    std::vector<Input_cell_descriptor> input_cell_range() const { return {}; }

    // pointers are for multi-typing of shapes
    Shape const * get_shape(Element_shape_type /*index*/) const { return nullptr; } // you can return nullptr if you want to ignore the cell
    Element_shape_type get_element_shape_id(Input_cell_descriptor /*cell*/) const { return 0; }
    Vertex_descriptor get_cell_vertex(Input_cell_descriptor /*cell*/, unsigned /*local_vertex_index*/) const { return 0; }

    bool has_reference_mesh = false;
    Point_3 get_ref_vertex_coordinates(Vertex_descriptor) const { return Point_3(); } // redefine if has_reference_mesh == true
};

template<typename MixedMesh>
class Mixed_mesh_wrapper {
public:
    Mixed_mesh_wrapper(MixedMesh &mesh)
    : _mesh(mesh)
    {
        optimization_tet_2_input_element.reserve(mesh.nb_input_cells());
        optimization_tet_2_input_element.clear();
        for (Input_cell_descriptor const &cell_descriptor : _mesh.input_cell_range()) {
            for (unsigned i = 0; i < get_nb_inner_tetrahedra(cell_descriptor); ++i) {
                optimization_tet_2_input_element.push_back({cell_descriptor, i});
            }
        }
    }

    using Cell_descriptor = typename MixedMesh::Cell_descriptor;
    using Vertex_descriptor = typename MixedMesh::Vertex_descriptor;
    using Point_3 = typename MixedMesh::Point_3;

    using Input_cell_descriptor = typename MixedMesh::Input_cell_descriptor;
    using Element_shape_type = typename MixedMesh::Element_shape_type;
    using Shape = typename MixedMesh::Shape;

    std::size_t nb_cells() const { return optimization_tet_2_input_element.size(); };
    std::size_t nb_vertices() const { return _mesh.nb_vertices(); }

    Point_3 vertex_coordinates(Vertex_descriptor vertex) const { return _mesh.vertex_coordinates(vertex); }
    void  set_vertex_coordinates(Vertex_descriptor vertex, Point_3 coord) { _mesh.set_vertex_coordinates(vertex, coord); }   // only non const


    utils::Contiguous_unsigned_range cell_range() const { return utils::Contiguous_unsigned_range{0, nb_cells()}; }
    std::array<Vertex_descriptor, 4> cell_vertices(Cell_descriptor cell) const {
        std::array<Vertex_descriptor, 4> sub_decomposition;
        auto input_element = optimization_tet_2_input_element[cell].first;
        std::size_t tet_number = optimization_tet_2_input_element[cell].second;
        for (std::size_t i = 0; i < 4; ++i) {
            sub_decomposition[i] = _mesh.get_cell_vertex(input_element, get_element_local_vert(input_element, tet_number, i));
        }
        return sub_decomposition;
    }
    std::array<Point_3, 4> cell_reference_shape(Cell_descriptor cell) const {
        return get_element_ref_shape(optimization_tet_2_input_element[cell].first, optimization_tet_2_input_element[cell].second);
    }

private:

    unsigned get_nb_inner_tetrahedra(Input_cell_descriptor cell) const {
        Shape const * shape = _mesh.get_shape(_mesh.get_element_shape_id(cell));
        if (shape == nullptr) return 0;
        return shape->nb_inner_tetrahedra();
    };

    unsigned get_element_local_vert(Input_cell_descriptor cell, unsigned tet, unsigned tet_vert) const {
        Shape const * shape = _mesh.get_shape(_mesh.get_element_shape_id(cell));
        assert(shape != nullptr);
        return shape->inner_tetrahedra_local_vert(tet, tet_vert);
    };

    std::array<Point_3, 4> get_element_ref_shape(Input_cell_descriptor cell, unsigned tet) const {
        Shape const * shape = _mesh.get_shape(_mesh.get_element_shape_id(cell));
        assert(shape != nullptr);
        if (!_mesh.has_reference_mesh) {
            return shape->inner_tetrahedra_reference_shape(tet);
        }
        else {
            return {
                _mesh.get_ref_vertex_coordinates(_mesh.get_cell_vertex(cell, get_element_local_vert(cell, tet, 0))),
                _mesh.get_ref_vertex_coordinates(_mesh.get_cell_vertex(cell, get_element_local_vert(cell, tet, 1))),
                _mesh.get_ref_vertex_coordinates(_mesh.get_cell_vertex(cell, get_element_local_vert(cell, tet, 2))),
                _mesh.get_ref_vertex_coordinates(_mesh.get_cell_vertex(cell, get_element_local_vert(cell, tet, 3)))
            };
        }
    };

    MixedMesh &_mesh;
    std::vector<std::pair<Input_cell_descriptor, unsigned>> optimization_tet_2_input_element;

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


namespace cgal_types {


template <typename C3t3>
class C3t3_wrapper {
public:
    C3t3_wrapper(C3t3 &c3t3)
      :c3t3(c3t3)
    {}

    using Cell_descriptor = typename C3t3::Cell_handle;
    using Vertex_descriptor = typename C3t3::Vertex_handle;
    using Face_descriptor = typename C3t3::Facet;
    using Edge_descriptor = typename C3t3::Edge;
    using Normal_3 = typename C3t3::Triangulation::Geom_traits::Vector_3;
    using Point_3 = typename C3t3::Triangulation::Geom_traits::Point_3;
    using Weighted_point_3 = typename C3t3::Triangulation::Geom_traits::Weighted_point_3;
    using Surface_patch_index = std::pair<typename C3t3::Surface_patch_index, Face_descriptor>;
    using Curve_index = std::pair<typename C3t3::Curve_index, Edge_descriptor>;
    using Construct_point_3 = typename C3t3::Triangulation::Geom_traits::Construct_point_3;

    std::size_t nb_cells() const { return c3t3.number_of_cells(); }
    std::size_t nb_faces() const { return c3t3.number_of_facets(); }
    std::size_t nb_edges() const { return c3t3.number_of_edges(); }
    std::size_t nb_vertices() const { return c3t3.triangulation().number_of_vertices(); }

    decltype(auto) vertex_coordinates(Vertex_descriptor vertex) const {
        return Mesh_smoothing_3_internal::get_point<C3t3>(vertex); // c3t3 holds weighted points
    }
    void set_vertex_coordinates(Vertex_descriptor vertex, const Point_3& coord)
    {
        const auto old_point = vertex->point();
        const auto new_point = Construct_point_3()(coord);
        if constexpr (std::is_same_v<std::decay_t<decltype(old_point)>, Weighted_point_3>)
            vertex->set_point(Weighted_point_3(new_point, old_point.weight()));
        else
            vertex->set_point(new_point);
    }
    auto cell_range() const { return c3t3.cells_in_complex(); }
    std::array<Vertex_descriptor, 4> cell_vertices(Cell_descriptor cell) const {
        std::array<Vertex_descriptor, 4> vertices;
        for (int i = 0; i < 4; ++i) {
            vertices[static_cast<unsigned>(i)] = cell->vertex(i);
        }
        return vertices;
    }
    std::array<Point_3, 4> cell_reference_shape(Cell_descriptor) const {
        return Shapes::VTK_TETRAHEDRON<Point_3>();
    }

    auto face_range() const { return c3t3.facets_in_complex(); }
    std::size_t nb_face_vertices(Face_descriptor) const { return 3; }
    Surface_patch_index patch_id(Face_descriptor face) const { return {c3t3.surface_patch_index(face), face}; }
    std::vector<Vertex_descriptor> face_vertices(Face_descriptor face) const {
        std::vector<Vertex_descriptor> vertices(3);
        for (int i = 1; i < 4; ++i) {
            vertices[static_cast<unsigned>(i-1)] = face.first->vertex((face.second + i)%4);
        }
        return vertices;
    }

    auto edge_range() const { return c3t3.edges_in_complex(); }
    Curve_index curve_id(Edge_descriptor edge) const { return {c3t3.curve_index(edge), edge}; }
    Vertex_descriptor edge_vertex(Edge_descriptor edge, unsigned i) const {
        assert(i < 2);
        return edge.first->vertex(i == 0 ? edge.second : edge.third);
    }

    C3t3 &c3t3;
};


template <typename LCC>
class LCC_mixed_mesh {
public:
    using Cell_descriptor = std::size_t;
    using Vertex_descriptor = typename LCC::Vertex_attribute_descriptor;
    using Point_3 = typename LCC::Point;
    using Input_cell_descriptor = typename LCC::Dart_descriptor;
    using Element_shape_type = typename CGAL::CMap::Element_topo::cell_topo;

    using Shape = Mesh_smoothing_3::Shapes::Base_element_shape_reference<Point_3>;

    using Input_cell_range = std::vector<Input_cell_descriptor>;

    std::size_t nb_vertices() const { return lcc.number_of_vertex_attributes(); }
    std::size_t nb_input_cells() const { return input_cells.size(); }

    Point_3 vertex_coordinates(Vertex_descriptor vertex) const { return vertex->point(); }
    void  set_vertex_coordinates(Vertex_descriptor vertex, Point_3 coord) { vertex->point() = coord; }

    std::vector<Input_cell_descriptor> input_cell_range() const { return input_cells; }

    Shape const * get_shape(Element_shape_type type) const {
        using namespace CGAL::CMap::Element_topo;
        switch(type)
        {
        case TETRAHEDRON:
            return &tet_ref;
        case HEXAHEDRON:
            return &hex_ref;
        case PYRAMID:
            return &pyr_ref;
        case PRISM:
            return &wedge_ref;
        default:
            // Unsupported LCC volume.
            return nullptr;
        }
    }
    Element_shape_type get_element_shape_id(Input_cell_descriptor cell) const {
        Input_cell_descriptor d;
        return CGAL::CMap::Element_topo::get_cell_topo<3>(lcc, cell, d);
    }

    Vertex_descriptor get_cell_vertex(Input_cell_descriptor cell, unsigned local_vertex) const {
        std::vector<Vertex_descriptor> const vertices = vtk_cell_vertices(cell);
        CGAL_assertion(local_vertex < vertices.size());
        return vertices[local_vertex];
     }



    bool has_reference_mesh = false;
    Point_3 get_ref_vertex_coordinates(Vertex_descriptor) const { return Point_3(); } // redefine if has_reference_mesh == true

public:
    LCC_mixed_mesh(LCC &lcc)
      :lcc(lcc)
    {
        static_assert(LCC::dimension == 3);
        static_assert(LCC::ambient_dimension == 3);

        wedge_ref.inverse = true;

        // Important: the iterator itself is the Dart_descriptor.
        auto cells = lcc.template one_dart_per_cell<3>();
        input_cells.reserve(std::distance(cells.begin(), cells.end()));
        for(auto cell = cells.begin(); cell != cells.end(); ++cell) input_cells.push_back(cell);
    }
private:
    LCC &lcc;
    std::vector<Input_cell_descriptor> input_cells;

    Shapes::VTK_TETRAHEDRON<Point_3> tet_ref;
    Shapes::VTK_HEXAHEDRON<Point_3> hex_ref;
    Shapes::VTK_PYRAMID<Point_3> pyr_ref;
    Shapes::VTK_WEDGE<Point_3> wedge_ref;

private:

    // Recover the local vertices using exactly the same convention as CGAL::IO::write_VTK().
    std::vector<Vertex_descriptor>
    vtk_cell_vertices(Input_cell_descriptor cell) const {
        using namespace CGAL::CMap::Element_topo;
        Input_cell_descriptor sd;
        Element_shape_type type = get_cell_topo<3>(lcc, cell, sd);

        std::vector<Vertex_descriptor> vertices;
        switch(type)
        {
        case TETRAHEDRON:
            vertices = {
                lcc.vertex_attribute(sd),
                lcc.vertex_attribute(lcc.template beta<1>(sd)),
                lcc.vertex_attribute(lcc.template beta<0>(sd)),
                lcc.vertex_attribute(lcc.template beta<2, 0>(sd))
            };
            break;

        case PYRAMID:
            vertices = {
                lcc.vertex_attribute(sd),
                lcc.vertex_attribute(lcc.template beta<1>(sd)),
                lcc.vertex_attribute(lcc.template beta<1, 1>(sd)),
                lcc.vertex_attribute(lcc.template beta<0>(sd)),
                lcc.vertex_attribute(lcc.template beta<2, 0>(sd))
            };
            break;

        case PRISM:
        {
            vertices = {
                lcc.vertex_attribute(sd),
                lcc.vertex_attribute(lcc.template beta<1>(sd)),
                lcc.vertex_attribute(lcc.template beta<0>(sd))
            };

            const auto sd2 = lcc.template beta<2, 1, 1, 2>(sd);
            vertices.push_back(lcc.vertex_attribute(lcc.template beta<1>(sd2)));
            vertices.push_back(lcc.vertex_attribute(sd2));
            vertices.push_back(lcc.vertex_attribute(lcc.template beta<0>(sd2)));
            break;
        }
        case HEXAHEDRON:
        {
            auto current = sd;
            // VTK vertices 0..3.
            for(unsigned i = 0; i < 4; ++i) {
                vertices.push_back(lcc.vertex_attribute(current));
                current = lcc.template beta<1>(current);
            }

            // VTK vertices 4..7.
            auto opposite = lcc.template beta<2, 1, 1, 2, 1>(current);
            for(unsigned i = 0; i < 4; ++i) {
                vertices.push_back(lcc.vertex_attribute(opposite));
                opposite = lcc.template beta<0>(opposite);
            }
            break;
        }
        default:
            break;
        }

        return vertices;
    }
};

template <typename LCC>
class LCC_surface_mesh {
public:
    using Vertex_descriptor = typename LCC::Vertex_attribute_descriptor;
    using Normal_3 = typename LCC::Vector;
    using Face_descriptor = std::size_t;
    using Surface_patch_index = unsigned;


    std::size_t nb_faces() const {
        return boundary_faces.size();
    }

    utils::Contiguous_unsigned_range face_range() const {
        return {0, nb_faces()};
    }

    std::size_t nb_face_vertices(Face_descriptor face) const {
        const auto first = boundary_faces[face];
        auto dart = first;
        std::size_t n = 0;
        do {
            ++n;
            dart = lcc.template beta<1>(dart);
        } while(dart != first);
        return n;
    }

    Surface_patch_index patch_id(Face_descriptor face) const {
        // For now, each boundary polygon is its own patch.
        return face;
    }

    std::vector<Vertex_descriptor>
    face_vertices(Face_descriptor face) const {
        std::vector<Vertex_descriptor> vertices;
        const auto first = boundary_faces[face];
        auto dart = first;
        do {
            vertices.push_back(lcc.vertex_attribute(dart));
            dart = lcc.template beta<1>(dart);
        } while(dart != first);
        return vertices;
    }

    LCC_surface_mesh(LCC &lcc): lcc(lcc) {
        auto faces = lcc.template one_dart_per_cell<2>();

        for(auto face = faces.begin(); face != faces.end(); ++face) {
            if(lcc.template is_free<3>(face))
                boundary_faces.push_back(face);
        }
    }

private:
    LCC &lcc;
    std::vector<typename LCC::Dart_descriptor> boundary_faces;
};

} } } // end of CGAL::Mesh_smoothing_3::default_structures namespace

#endif // CGAL_MESH_SMOOTHING_3_MESH_REPRESENTATIONS_H
