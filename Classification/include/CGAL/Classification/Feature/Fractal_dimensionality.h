#ifndef CGAL_CLASSIFICATION_FEATURE_FRACTAL_DIMENSION_H
#define CGAL_CLASSIFICATION_FEATURE_FRACTAL_DIMENSION_H
#include <CGAL/Classification.h>
#include <CGAL/license/Classification.h>

#include <iostream>
#include <math.h>

#include <Eigen/Dense>
#include <boost/multi_array.hpp>

#include <CGAL/Kd_tree.h>
#include <CGAL/Splitters.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix2d;
typedef Eigen::Matrix<double, 1, Eigen::Dynamic> VectorXd;
typedef Eigen::Matrix<int, 1, Eigen::Dynamic> VectorXi;
typedef boost::multi_array<int, 3> Matrix3i;
typedef Matrix3i::array_view<3>::type Matrix3iV;
typedef Matrix3i::index_range range;


namespace CGAL {
namespace   Classification {
namespace     Feature {

// given f() takes extra arguments
template<class T, class F, class... Args>
typename std::enable_if<(T::dimensionality == 1), void>::type IterateArrayView(T& array, F f, Args& ...args) {
    for (auto& element : array) {
        f(element, args...);
    }
}

template<class T, class F, class... Args>
typename std::enable_if<(T::dimensionality > 1), void>::type IterateArrayView(T& array, F f, Args& ...args) {
    for (auto element : array) {
        IterateArrayView<decltype(element), F, Args...>(element, f, args...);
    }
}

// User-defined feature
template <typename GeomTraits, typename PointRange, typename PointMap, typename ConcurrencyTag = CGAL::Parallel_if_available_tag>
class Fractal_dimensionality : public CGAL::Classification::Feature_base
{
    // Feature description
    const PointRange& input;
    PointMap point_map;
    
    std::vector<float> feature_values; // fractal dimensionality

public:
    Fractal_dimensionality(const PointRange& input, PointMap point_map, float feature_scale, float r1 = -1, float r2 = -1, std::size_t n = 30) : input(input), point_map(point_map)
    {
        assert((r1 >= 0) && (r2 > r1));

        this->set_name("Fractal_dimensionality");

        typedef GeomTraits::Point_3 Point;
        typedef GeomTraits::Iso_cuboid_3 Iso_cuboid_3;

        typedef CGAL::Search_traits_3<GeomTraits> SearchTraits_3;
        typedef CGAL::Sliding_midpoint<SearchTraits_3> Splitter;
        typedef CGAL::Kd_tree<SearchTraits_3, Splitter, CGAL::Tag_true> Tree;
        typedef CGAL::Fuzzy_iso_box<SearchTraits_3> Iso_box;

        double l = feature_scale; // the length of neighborhoods (default: 6)
        if (r1 == -1)
            r1 = feature_scale * 0.3f; // the smallest ratio (default: 0.3)
        if (r2 == -1)
            r2 = feature_scale * 0.8f; // the biggest ratio (default: 0.8)
        double step = (r2 - r1) / n;

        // Build K-D Tree
        Tree m_tree;
        for (auto pt = point_map.begin(); pt != point_map.end(); pt++)
            m_tree.insert(*pt);

        double gb_xmin = m_tree.bounding_box().min_coord(0), gb_xmax = m_tree.bounding_box().max_coord(0);
        double gb_ymin = m_tree.bounding_box().min_coord(1), gb_ymax = m_tree.bounding_box().max_coord(1);
        double gb_zmin = m_tree.bounding_box().min_coord(2), gb_zmax = m_tree.bounding_box().max_coord(2);

        // Initialize containers to save results
        std::vector<Matrix2d> results;

        // Initialize containers to save bounding boxes
        Matrix2d bboxes(input.size() + 1, 6);
        bboxes.row(0) << gb_xmin, gb_ymin, gb_zmin, gb_xmin, gb_ymin, gb_zmin;


        std::vector<Point> m_nbs;
        int global_index = 1;
        for (auto pt = point_map.begin(); pt != point_map.end(); pt++) {
            Matrix2d tmp(n, 4);
            results.push_back(tmp);

            // Find neighbourhoods
            Point pt_left(pt->x() - l / 2, pt->y() - l / 2, pt->z() - l / 2);
            Point pt_right(pt->x() + l / 2, pt->y() + l / 2, pt->z() + l / 2);
            Iso_box m_box(pt_left, pt_right);
            m_tree.search(std::back_inserter(m_nbs), m_box);

            // Find bounding box
            Iso_cuboid_3 my_box = CGAL::bounding_box(m_nbs.begin(), m_nbs.end());
            double my_xmin = my_box.xmin(), my_xmax = my_box.xmax();
            double my_ymin = my_box.ymin(), my_ymax = my_box.ymax();
            double my_zmin = my_box.zmin(), my_zmax = my_box.zmax();

            bboxes.row(global_index++) << my_xmin, my_ymin, my_zmin, my_xmax, my_ymax, my_zmax;
            m_nbs.clear();
        }

        CGAL::for_each<ConcurrencyTag>
            (CGAL::make_counting_range<std::size_t>(0, n),
                [&](const std::size_t& i) -> bool {
                    double cell_size = r1 + i * step;
                    size_t gb_xsize = static_cast<size_t> (max(ceil((gb_xmax - gb_xmin) / cell_size), 1.));
                    size_t gb_ysize = static_cast<size_t> (max(ceil((gb_ymax - gb_ymin) / cell_size), 1.));
                    size_t gb_zsize = static_cast<size_t> (max(ceil((gb_zmax - gb_zmin) / cell_size), 1.));

                    Matrix3i my_voxel = Matrix3i(boost::extents[gb_xsize][gb_ysize][gb_zsize]);
                    std::fill_n(my_voxel.data(), my_voxel.num_elements(), 0);

                    for (auto pt = point_map.begin(); pt != point_map.end(); pt++) {
                        size_t my_xid = min((size_t)((pt->x() - gb_xmin) / cell_size), gb_xsize - 1);
                        size_t my_yid = min((size_t)((pt->y() - gb_ymin) / cell_size), gb_ysize - 1);
                        size_t my_zid = min((size_t)((pt->z() - gb_zmin) / cell_size), gb_zsize - 1);
                        my_voxel[my_xid][my_yid][my_zid] += 1;
                    }

                    std::size_t global_index = 0;
                    Matrix2d limits(1, 6);
                    limits.row(0) << gb_xsize - 1, gb_ysize - 1, gb_zsize - 1, gb_xsize - 1, gb_ysize - 1, gb_zsize - 1;

                    for (auto pt = point_map.begin(); pt != point_map.end(); pt++) {
                        // Calculate range
                        VectorXd ranges_d(6);
                        ranges_d << ((bboxes.row(global_index + 1) - bboxes.row(0)) / cell_size).cwiseMin(limits);
                        VectorXi ranges_i(6);
                        ranges_i = ranges_d.cast <int>();

                        // Create view from voxel matrix
                        Matrix3iV my_view = my_voxel[boost::indices[range(ranges_i[0], ranges_i[3] + 1)][range(ranges_i[1], ranges_i[4] + 1)][range(ranges_i[2], ranges_i[5] + 1)]];
                        double total_pts = 0.0;
                        IterateArrayView(my_view, [](int& elem, double& total_pts) {total_pts += elem; }, total_pts);

                        // Calculate Fractal Dimension
                        double d0 = 0., d1 = 0., d2 = 0.;
                        IterateArrayView(my_view, [](int& elem, double& d0, double& d1, double& d2, double& total_pts) {
                            if (elem > 0) {
                                d0 += 1;
                                d1 += elem / total_pts * log(elem / total_pts);
                                d2 += pow(elem / total_pts, 2);
                            }
                            }, d0, d1, d2, total_pts);

                        d0 = log(d0);
                        d2 = log(d2);

                        results[global_index].row(i) << cell_size, d0, d1, d2;
                        global_index += 1;
                    }
                });

        // Pair Counting
        std::vector<double> d0_map, d1_map, d2_map;
        d0_map.reserve(input.size());
        d1_map.reserve(input.size());
        d2_map.reserve(input.size());
        feature_values.resize(input.size(), 0.f);

        CGAL::for_each<ConcurrencyTag>
            (CGAL::make_counting_range<std::size_t>(0, input.size()),
                [&](const std::size_t& i) -> bool {
                    // Perform line fitting
                    Matrix2d A(n, 2);
                    for (int j = 0; j < n; ++j) {
                        A.row(j) << log(results[i](j, 0)), -results[i](j, 1); // d0
                        //A.row(j) << log(results[i](j, 0)), results[i](j, 2); // d1
                        //A.row(j) << log(results[i](j, 0)), results[i](j, 3); // d2
                    }
                    A.col(0).array() -= A.col(0).mean();
                    A.col(1).array() -= A.col(1).mean();
                    Eigen::JacobiSVD<Matrix2d> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
                    Matrix2d V = svd.matrixV();
                    d0_map[i] = V(1, 0) / V(0, 0); // slope
                    d0_map[i] = (d0_map[i] < 0) ? 0 : ((d0_map[i] > 3) ? 3 : d0_map[i]); // clamp
                    feature_values[i] = d0_map[i];
                });
       
    }

    float value(std::size_t pt_index) {
        return feature_values[pt_index];
    }
};

} // namespace Feature
} // namespace Classification
} // namespace CGAL
#endif CGAL_CLASSIFICATION_FEATURE_FRACTAL_DIMENSION_H