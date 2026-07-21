#pragma once

#include <vector>
#include <array>
#include <Eigen/Eigen>

namespace Mesh_smoothing_3 {

namespace Shapes {
    
    template < 
        typename Point3  = Eigen::Vector3d
    >
    class Base_element_shape_reference {
    public:
        virtual unsigned nb_inner_tetrahedra() const { return 0; }

        virtual unsigned nb_vertices() const { return 0; }

        virtual unsigned inner_tetrahedra_local_vert(unsigned tet, unsigned tet_vert) const { return 0; }

        virtual Point3 vertex_reference_coordinates(unsigned local_vert) const { return Point3(); }

        std::array<Point3, 4> inner_tetrahedra_reference_shape(unsigned tet) const {
            return {
                vertex_reference_coordinates(inner_tetrahedra_local_vert(tet, 0)),
                vertex_reference_coordinates(inner_tetrahedra_local_vert(tet, 1)),
                vertex_reference_coordinates(inner_tetrahedra_local_vert(tet, inverse?3:2)),
                vertex_reference_coordinates(inner_tetrahedra_local_vert(tet, inverse?2:3))
            };
        }

        bool inverse = false;
    };

    // This heavy class set-up is to obtain something generic that does not store data while still working with Eigen that does not allow constexpr allocations 
    template <
        typename Point3  = Eigen::Vector3d,
        unsigned NbVertices = 0,
        unsigned NbTetrahedra = 0
    >
    class Array_based_shape_reference : public Base_element_shape_reference<Point3>  {
    public:
        std::array<Point3, NbVertices> const * reference_points_ptr = nullptr;
        std::array<std::array<unsigned, 4>, NbTetrahedra> const * inner_tetrahedra_ptr = nullptr;

        unsigned nb_inner_tetrahedra() const override { return inner_tetrahedra_ptr->size(); }

        unsigned nb_vertices() const override { return reference_points_ptr->size(); }

        unsigned inner_tetrahedra_local_vert(unsigned tet, unsigned tet_vert) const override { return (*inner_tetrahedra_ptr)[tet][tet_vert]; }

        Point3 vertex_reference_coordinates(unsigned local_vert) const override { return (*reference_points_ptr)[local_vert]; }

};

    // following: https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html#legacy-file-examples
    template <
        typename Point3  = Eigen::Vector3d
    >
    class VTK_TETRAHEDRON : public Array_based_shape_reference<Point3, 4, 1>  {
    public:
        // cannot be constexpr because of Eigen initialization constraints
        static const std::array<Point3, 4> reference_points;
        static const std::array<std::array<unsigned, 4>, 1> inner_tetrahedra;
        VTK_TETRAHEDRON() {
            this->reference_points_ptr = &reference_points;
            this->inner_tetrahedra_ptr = &inner_tetrahedra;
        }
        operator std::array<Point3, 4> () const { return this->inner_tetrahedra_reference_shape(0);  }
    };
    
    template <typename Point3> const std::array<Point3, 4> VTK_TETRAHEDRON<Point3>::reference_points = {{
        Point3{std::sqrt(8./9),0,-1./3},
        Point3{-std::sqrt(2./9),std::sqrt(2./3),-1./3},
        Point3{-std::sqrt(2./9),-std::sqrt(2./3),-1./3},
        Point3{0.,0.,1.}
    }};
    template <typename Point3> const std::array<std::array<unsigned, 4>, 1> VTK_TETRAHEDRON<Point3>::inner_tetrahedra = {{
        {0,1,2,3}
    }};

    template <
        typename Point3  = Eigen::Vector3d
    >
    class VTK_PYRAMID : public Array_based_shape_reference<Point3, 5, 4>  {
    public:
        // cannot be constexpr because of Eigen initialization constraints
        static const std::array<std::array<unsigned, 4>, 4> inner_tetrahedra;
        static const std::array<Point3, 5> reference_points;
        VTK_PYRAMID() {
            this->reference_points_ptr = &reference_points;
            this->inner_tetrahedra_ptr = &inner_tetrahedra;
        }
    };
    template <typename Point3> const std::array<Point3, 5> VTK_PYRAMID<Point3>::reference_points = {{
        Point3{0.,0.,0.},
        Point3{1.,0.,0.},
        Point3{1.,1.,0.},
        Point3{0.,1.,0.},
        Point3{0.5,0.5,0.5}
    }};
    template <typename Point3> const std::array<std::array<unsigned, 4>, 4> VTK_PYRAMID<Point3>::inner_tetrahedra = {{
        {0,1,3,4},
        {1,2,0,4},
        {2,3,1,4},
        {3,0,2,4}
    }};

    template <
        typename Point3  = Eigen::Vector3d
    >
    class VTK_WEDGE : public Array_based_shape_reference<Point3, 6, 6>  {
    public:
        // cannot be constexpr because of Eigen initialization constraints
        static const std::array<Point3, 6> reference_points;
        static const std::array<std::array<unsigned, 4>, 6> inner_tetrahedra;
        VTK_WEDGE() {
            this->reference_points_ptr = &reference_points;
            this->inner_tetrahedra_ptr = &inner_tetrahedra;
        }
    };
    template <typename Point3> const std::array<Point3, 6> VTK_WEDGE<Point3>::reference_points = {{
        Point3{0.,0.,0.},
        Point3{1.,0.,0.},
        Point3{0.5,0.,std::sqrt(3)/2},
        Point3{0.,1.,0.},
        Point3{1.,1.,0.},
        Point3{0.5,1.,std::sqrt(3)/2},
    }};
    template <typename Point3> const std::array<std::array<unsigned, 4>, 6> VTK_WEDGE<Point3>::inner_tetrahedra = {{
        {0,1,3,2},
        {1,4,0,2},
        {2,1,0,5},
        {3,4,5,0},
        {4,5,3,1},
        {5,3,4,2},
    }};

    template <
        typename Point3  = Eigen::Vector3d
    >
    class VTK_HEXAHEDRON : public Array_based_shape_reference<Point3, 8, 8>  {
    public:
        // cannot be constexpr because of Eigen initialization constraints
        static const std::array<Point3, 8> reference_points;
        static const std::array<std::array<unsigned, 4>, 8> inner_tetrahedra;
        VTK_HEXAHEDRON() {
            this->reference_points_ptr = &reference_points;
            this->inner_tetrahedra_ptr = &inner_tetrahedra;
        }
    };
    template <typename Point3> const std::array<Point3, 8> VTK_HEXAHEDRON<Point3>::reference_points = {{
        Point3{0.,0.,0.},
        Point3{1.,0.,0.},
        Point3{1.,1.,0.},
        Point3{0.,1.,0.},
        Point3{0.,0.,1.},
        Point3{1.,0.,1.},
        Point3{1.,1.,1.},
        Point3{0.,1.,1.}
    }};
    template <typename Point3> const std::array<std::array<unsigned, 4>, 8> VTK_HEXAHEDRON<Point3>::inner_tetrahedra = {{
        {0,1,3,4},
        {1,2,0,5},
        {2,3,1,6},
        {3,2,0,7},
        {4,7,5,0},
        {5,4,7,1},
        {6,5,4,2},
        {7,6,4,3},
        // {0,2,7,5},
        // {1,3,4,6},
    }};

    // is equivalent to VTK voxels
    template <
        typename Point3  = Eigen::Vector3d
    >
    class GEOGRAM_HEXAHEDRON : public Array_based_shape_reference<Point3, 8, 8>  {
    public:
        static const std::array<Point3, 8> reference_points;
        static const std::array<std::array<unsigned, 4>, 8> inner_tetrahedra;
        GEOGRAM_HEXAHEDRON() {
            this->reference_points_ptr = &reference_points;
            this->inner_tetrahedra_ptr = &inner_tetrahedra;
        }
    };
    template <typename Point3> const std::array<Point3, 8> GEOGRAM_HEXAHEDRON<Point3>::reference_points = {{
        Point3{0.,0.,0.},
        Point3{1.,0.,0.},
        Point3{0.,1.,0.},
        Point3{1.,1.,0.},
        Point3{0.,0.,1.},
        Point3{1.,0.,1.},
        Point3{0.,1.,1.},
        Point3{1.,1.,1.}
    }};
    template <typename Point3> const std::array<std::array<unsigned, 4>, 8> GEOGRAM_HEXAHEDRON<Point3>::inner_tetrahedra = {{
        {0,1,2,4},
        {1,3,0,5},
        {3,2,1,7},
        {2,0,3,6},
        {4,6,5,0},
        {5,4,7,1},
        {7,5,6,3},
        {6,7,4,2},
        // {0,3,6,5},
        // {1,2,4,7},
    }};
}



}

