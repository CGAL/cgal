#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERPOLATED_CURVATURE_MEASURES_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERPOLATED_CURVATURE_MEASURES_H
#endif

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <vector>
#include <numeric>
#include <CGAL/assertions.h>

namespace CGAL {

namespace Polygon_mesh_processing {

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;

// enum to specify which measure is computed
enum MEASURE_INDEX {
    MU0_AREA_MEASURE,
    MU1_MEAN_CURVATURE_MEASURE,
    MU2_GAUSSIAN_CURVATURE_MEASURE
};

Epic::FT interpolated_mu_i_triangle(const Epic::Vector_3 x0, const Epic::Vector_3 x1, const Epic::Vector_3 x2,
    const Epic::Vector_3 u0, const Epic::Vector_3 u1, const Epic::Vector_3 u2, MEASURE_INDEX mu_i)
{
    Epic::Vector_3 um;
    switch (mu_i)
    {
    case MU0_AREA_MEASURE:

        um = (u0 + u1 + u2) / 3.0;
        
        return 0.5 * um * CGAL::cross_product(x1 - x0, x2 - x0);

    case MU1_MEAN_CURVATURE_MEASURE:

        um = (u0 + u1 + u2) / 3.0;
        
        return 0.5 * um * (CGAL::cross_product(u2 - u1, x0)
                         + CGAL::cross_product(u0 - u2, x1)
                         + CGAL::cross_product(u1 - u0, x2));

    case MU2_GAUSSIAN_CURVATURE_MEASURE:

        return 0.5 * u0 * CGAL::cross_product(u1, u2);

    default: return 0;
    }
}

Epic::FT interpolated_mu_i_face(std::vector<Epic::Vector_3>& x, std::vector<Epic::Vector_3>& u, MEASURE_INDEX mu_i)
{
    std::size_t n = x.size();
    CGAL_precondition(u.size() == n);
    CGAL_precondition(n >= 3);

    if (n == 3)
        return interpolated_mu_i_triangle(x[0], x[1], x[2],
            u[0], u[1], u[2], mu_i);

    /// If Quad measure formulas (Bilinear Interpolation) proved to be better,
    /// they will be implemented and called here
    //else if(n == 4)
    //    return interpolated_mu0_quad(x[0], x[1], x[2], x[3]
    //                                 u[0], u[1], u[2], u[3]);
    else {
        Epic::FT mu0 = 0;

        // getting barycenter of points
        Epic::Vector_3 xm = std::accumulate(x.begin(), x.end(), Epic::Vector_3(0, 0, 0));
        xm /= n;

        // getting unit average normal of points
        Epic::Vector_3 um = std::accumulate(u.begin(), u.end(), Epic::Vector_3(0,0,0));
        um /= sqrt(um * um);

        // summing each triangle's measure after triangulation by barycenter split.
        for (std::size_t i = 0; i < n; i++)
        {
            mu0 += interpolated_mu_i_triangle(x[i], x[(i + 1) % n], xm,
                u[i], u[(i + 1) % n], um, mu_i);
        }
        return mu0;
    }
}

}
}
