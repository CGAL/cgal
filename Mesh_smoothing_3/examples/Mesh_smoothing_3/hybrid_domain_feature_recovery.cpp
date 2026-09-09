#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/Mesh_smoothing_3/boundary_aware_mesh_smoothing.h>
#include <CGAL/Mesh_smoothing_3/projectors.h>

#include <CGAL/centroid.h>
#include <CGAL/Kernel/global_functions_3.h>

#include <CGAL/IO/File_medit.h>

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>

// Kernel
using K = CGAL::Exact_predicates_inexact_constructions_kernel;

// Implicit domain
using Implicit_domain = CGAL::Labeled_mesh_domain_3<K>;

// Polyhedral domain
using Polyhedron = CGAL::Polyhedron_3<K>;
using Polyhedron_domain = CGAL::Polyhedral_mesh_domain_3<Polyhedron, K>;
using Polyhedron_projector = CGAL::Mesh_smoothing_3::Polyhedral_mesh_domain_projector<Polyhedron_domain>;

class Hybrid_domain
{
    const Implicit_domain& implicit_domain;
    const Polyhedron_domain& polyhedron_domain;
    Polyhedron_projector polyhedron_projector;

public:
    Hybrid_domain(const Implicit_domain& implicit_domain,
                  const Polyhedron_domain& polyhedron_domain)
        : implicit_domain(implicit_domain)
        , polyhedron_domain(polyhedron_domain)
        , polyhedron_projector(polyhedron_domain)
    {}

    // Types required by MeshDomain_3
    using Surface_patch_index = int;
    using Subdomain_index = int;
    using Index = int;

    using R = K;
    using Point_3 = K::Point_3;
    using Vector_3 = K::Vector_3;
    using FT = K::FT;
    using Intersection = std::tuple<Point_3, Index, int>;

    // Type required by ConstructTangentSpace
    using Tangent_space = CGAL::Mesh_smoothing_3::Tangent_space<K>;

    CGAL::Bbox_3 bbox() const
    {
        return implicit_domain.bbox() + polyhedron_domain.bbox();
    }

    struct Construct_initial_points
    {
        Construct_initial_points(const Hybrid_domain& domain)
            : r_domain_(domain)
        {}

        template<class OutputIterator>
        OutputIterator operator()(OutputIterator pts, const int n = 20) const
        {
            using Implicit_Index = Implicit_domain::Index;
            std::vector<std::pair<Point_3, Implicit_Index> > implicit_points_vector;

            Implicit_domain::Construct_initial_points cstr_implicit_initial_points =
                r_domain_.implicit_domain.construct_initial_points_object();

            cstr_implicit_initial_points(
                std::back_inserter(implicit_points_vector), n / 2);

            for(const auto& p : implicit_points_vector)
                *pts++ = std::make_pair(p.first, 2);

            using Polyhedron_Index = Polyhedron_domain::Index;
            std::vector<std::pair<Point_3, Polyhedron_Index> > polyhedron_points_vector;

            Polyhedron_domain::Construct_initial_points cstr_polyhedron_initial_points =
                r_domain_.polyhedron_domain.construct_initial_points_object();

            cstr_polyhedron_initial_points(
                std::back_inserter(polyhedron_points_vector), n / 2);

            for(const auto& p : polyhedron_points_vector)
                *pts++ = std::make_pair(p.first, 1);

            return pts;
        }

    private:
        const Hybrid_domain& r_domain_;
    };

    Construct_initial_points construct_initial_points_object() const
    {
        return Construct_initial_points(*this);
    }

    struct Is_in_domain
    {
        Is_in_domain(const Hybrid_domain& domain)
            : r_domain_(domain)
        {}

        std::optional<Subdomain_index> operator()(const Point_3& p) const
        {
            const std::optional<Subdomain_index> subdomain_index =
                r_domain_.implicit_domain.is_in_domain_object()(p);

            if(subdomain_index)
                return 2;

            return r_domain_.polyhedron_domain.is_in_domain_object()(p);
        }

    private:
        const Hybrid_domain& r_domain_;
    };

    Is_in_domain is_in_domain_object() const
    {
        return Is_in_domain(*this);
    }

    struct Construct_intersection
    {
        Construct_intersection(const Hybrid_domain& domain)
            : r_domain_(domain)
        {}

        template<typename Query>
        Intersection operator()(const Query& query) const
        {
            using boost::get;

            const Implicit_domain::Intersection implicit_inter =
                r_domain_.implicit_domain.construct_intersection_object()(query);

            if(get<2>(implicit_inter) != 0)
                return Intersection(get<0>(implicit_inter), 2, get<2>(implicit_inter));

            const Polyhedron_domain::Intersection polyhedron_inter =
                r_domain_.polyhedron_domain.construct_intersection_object()(query);

            if(get<2>(polyhedron_inter) != 0)
            {
                const Point_3 inter_point = get<0>(polyhedron_inter);

                if(!r_domain_.implicit_domain.is_in_domain_object()(inter_point))
                    return Intersection(inter_point, 1, get<2>(polyhedron_inter));
            }

            return Intersection();
        }

    private:
        const Hybrid_domain& r_domain_;
    };

    Construct_intersection construct_intersection_object() const
    {
        return Construct_intersection(*this);
    }

    Index index_from_surface_patch_index(const Surface_patch_index& index) const
    {
        return index;
    }

    Index index_from_subdomain_index(const Subdomain_index& index) const
    {
        return index;
    }

    Surface_patch_index surface_patch_index(const Index& index) const
    {
        return index;
    }

    Subdomain_index subdomain_index(const Index& index) const
    {
        return index;
    }

    template<typename Patch_face>
    Tangent_space
    patch_face_projection_plane(const Patch_face& patch_face,
                                const std::vector<Point_3>& face_points) const
    {
        const Point_3 face_center =
            CGAL::centroid(face_points.begin(), face_points.end());

        // Patch 1 is the polyhedral domain.
        if(patch_face.first == 1)
            return polyhedron_projector.patch_face_projection_plane(patch_face, face_points);

        // Patch 2 is the implicit sphere:
        // ||p - (1,1,1)||^2 - 1 = 0.
        CGAL_assertion(patch_face.first == 2);

        const Point_3 sphere_center(1., 1., 1.);
        const Vector_3 radial = face_center - sphere_center;

        CGAL_assertion(radial.squared_length() != FT(0));

        const FT length = CGAL::sqrt(radial.squared_length());
        const Vector_3 normal = radial / length;
        const Point_3 projected = sphere_center + normal;

        return Tangent_space{projected, normal};
    }

    template<typename Curve_edge>
    Tangent_space
    curve_edge_projection_line(const Curve_edge& curve_edge,
                               const std::array<Point_3, 2>& edge_points) const
    {
        return polyhedron_projector.curve_edge_projection_line(curve_edge, edge_points);
    }
};

using Domain = Hybrid_domain;

// Triangulation
using Tr = CGAL::Mesh_triangulation_3<Domain, K>::type;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;

// Criteria
using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;
using Facet_criteria = Mesh_criteria::Facet_criteria;
using Cell_criteria = Mesh_criteria::Cell_criteria;

using Point = K::Point_3;
using FT = K::FT;

FT sphere_centered_at_111(const Point& p)
{
    const FT dx = p.x() - 1;
    const FT dy = p.y() - 1;
    const FT dz = p.z() - 1;

    return dx * dx + dy * dy + dz * dz - 1;
}

namespace params = CGAL::parameters;

int main()
{
    const std::string fname = CGAL::data_file_path("meshes/cube.off");

    Polyhedron polyhedron;
    std::ifstream input(fname);
    input >> polyhedron;

    if(input.bad())
    {
        std::cerr << "Error: Cannot read file " << fname << std::endl;
        return EXIT_FAILURE;
    }

    Polyhedron_domain polyhedron_domain(polyhedron);

    Implicit_domain sphere_domain =
        Implicit_domain::create_implicit_mesh_domain(
            sphere_centered_at_111,
            K::Sphere_3(K::Point_3(1, 1, 1), K::FT(2)));

    Domain domain(sphere_domain, polyhedron_domain);

    Facet_criteria facet_criteria(30, 0.08, 0.025);
    Cell_criteria cell_criteria(2, 0.1);
    Mesh_criteria criteria(facet_criteria, cell_criteria);

    C3t3 c3t3 =
        CGAL::make_mesh_3<C3t3>(
            domain,
            criteria,
            params::no_perturb().no_exude());

    dump_c3t3(c3t3, "hybrid_initial");


    const auto result =
        CGAL::boundary_aware_mesh_smoothing(
            c3t3,
            domain,
            CGAL::parameters::verbose(true));

    std::cout << "Number of inverted elements: "
              << result.nb_invalid_elements << '\n';
    std::cout << "Number of vertex updates: "
              << result.nb_vertex_updates << '\n';
    std::cout << "Number of metric evaluations: "
              << result.nb_metric_evaluations << '\n';
    std::cout << "Smoothing time: "
              << result.total_time << " s." << '\n';

    dump_c3t3(c3t3, "hybrid_smoothed");

    return EXIT_SUCCESS;
}