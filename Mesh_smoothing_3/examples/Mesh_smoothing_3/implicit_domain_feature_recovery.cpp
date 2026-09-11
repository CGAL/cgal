#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Implicit_to_labeling_function_wrapper.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Mesh_3/Dump_c3t3.h>

#include <CGAL/Mesh_smoothing_3/boundary_aware_mesh_smoothing.h>
#include <CGAL/Mesh_smoothing_3/Projectors.h>

#include <utility>
#include <vector>
#include <string>

// Kernel
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = K::FT;
using Point_3 = K::Point_3;
using Vector_3 = K::Vector_3;

using Value_and_gradient = std::pair<FT, Vector_3>;

// -----------------------------------------------------------------------------
// Implicit functions
// -----------------------------------------------------------------------------

Value_and_gradient torus_distance_and_gradient(const Point_3& p)
{
    const FT major_radius = FT(1.5);
    const FT minor_radius = FT(0.5);

    const FT x = p.x();
    const FT y = p.y();
    const FT z = p.z();

    const FT radial = CGAL::sqrt(x * x + z * z);
    const FT dr = radial - major_radius;
    const FT tube_distance = CGAL::sqrt(dr * dr + y * y);

    const FT value = tube_distance - minor_radius;

    if (radial < FT(1e-8) || tube_distance < FT(1e-8))
        return {value, Vector_3(FT(1), FT(0), FT(0))};

    const Vector_3 gradient(
        dr * x / (tube_distance * radial),
        y / tube_distance,
        dr * z / (tube_distance * radial));

    return {value, gradient};
}

Value_and_gradient sphere_distance_and_gradient(const Point_3& p)
{
    const FT x = p.x();
    const FT y = p.y();
    const FT z = p.z();

    const FT distance_to_center =
        CGAL::sqrt(x * x + y * y + z * z);

    const FT radius = CGAL::sqrt(FT(3));

    if (distance_to_center < FT(1e-8))
        return {-radius, Vector_3(FT(1), FT(0), FT(0))};

    return {
        distance_to_center - radius,
        Vector_3(
            x / distance_to_center,
            y / distance_to_center,
            z / distance_to_center)
    };
}

// Scalar versions used by Mesh_3.
FT torus_function(const Point_3& p)
{
    return torus_distance_and_gradient(p).first;
}

FT sphere_function(const Point_3& p)
{
    return sphere_distance_and_gradient(p).first;
}

// -----------------------------------------------------------------------------
// Mesh domain
// -----------------------------------------------------------------------------

using Function = FT (*)(const Point_3&);

using Function_wrapper =
    CGAL::Implicit_multi_domain_to_labeling_function_wrapper<Function>;

using Function_vector = Function_wrapper::Function_vector;
using Mesh_domain = CGAL::Labeled_mesh_domain_3<K>;

// Triangulation
using Tr = CGAL::Mesh_triangulation_3<Mesh_domain>::type;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;

// Criteria
using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;
using Facet_criteria = Mesh_criteria::Facet_criteria;
using Cell_criteria = Mesh_criteria::Cell_criteria;

int main()
{
    namespace params = CGAL::parameters;

    // The same domain as mesh_implicit_domains_2.cpp:
    //
    //     torus_function > 0
    //     sphere_function < 0
    //
    Function_vector functions;
    functions.push_back(&torus_function);
    functions.push_back(&sphere_function);

    std::vector<std::string> signs;
    signs.push_back("+-");

    Mesh_domain domain(
        Function_wrapper(functions, signs),
        K::Sphere_3(
            CGAL::ORIGIN,
            CGAL::square(FT(5))),
        params::relative_error_bound(1e-6));

    // Same mesh criteria as the original example.
    Facet_criteria facet_criteria(
        30,     // angle
        0.2,    // size
        0.02);  // approximation

    Cell_criteria cell_criteria(
        2.,     // radius-edge ratio
        0.4);   // size

    Mesh_criteria criteria(
        facet_criteria,
        cell_criteria);

    C3t3 c3t3 =
        CGAL::make_mesh_3<C3t3>(
            domain,
            criteria,
            params::no_exude().no_perturb());

    CGAL::dump_c3t3(c3t3, "implicit_initial");

    // -------------------------------------------------------------------------
    // Projection target
    // -------------------------------------------------------------------------
    //
    // The domain corresponds to:
    //
    //     torus > 0  &&  sphere < 0
    //
    // Reorient both functions so that the desired domain is negative:
    //
    //     -torus < 0 && sphere < 0.
    //
    // max(-torus, sphere) therefore represents the entire boundary as one
    // implicit function. No patch identifier is used by the projector.
    //
    auto projection_function = [](const Point_3& p) -> Value_and_gradient {
        const auto [torus_distance, torus_gradient] =
            torus_distance_and_gradient(p);

        const auto [sphere_distance, sphere_gradient] =
            sphere_distance_and_gradient(p);

        // Desired domain:
        //
        //   outside torus  -> -torus_distance < 0
        //   inside sphere  ->  sphere_distance < 0
        //
        const FT outside_torus_distance = -torus_distance;

        if(outside_torus_distance > sphere_distance)
        {
            return {
                outside_torus_distance,
                -torus_gradient
            };
        }

        return {
            sphere_distance,
            sphere_gradient
        };
    };

    using Projector =
        CGAL::Mesh_smoothing_3::Signed_distance_function_projector<
            K,
            decltype(projection_function)>;

    // tolerance is expressed in the units of the input geometry.
    Projector projector(
        projection_function,
        10,     // maximum projection iterations
        1e-8);  // positional tolerance

    const auto result =
        CGAL::boundary_aware_mesh_smoothing(
            c3t3,
            projector,
            CGAL::parameters::verbose(true));

    std::cout << "Number of inverted elements: "
              << result.nb_invalid_elements << '\n';
    std::cout << "Number of vertex updates: "
              << result.nb_vertex_updates << '\n';
    std::cout << "Number of metric evaluations: "
              << result.nb_metric_evaluations << '\n';
    std::cout << "Smoothing time: "
              << result.total_time << " s." << '\n';

    CGAL::dump_c3t3(c3t3, "implicit_smoothed");

    return EXIT_SUCCESS;
}