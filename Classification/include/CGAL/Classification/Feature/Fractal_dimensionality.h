#include <CGAL/Classification.h>
#include <CGAL/license/Classification.h>

/*
#include <omp.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classification.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Real_timer.h>
// Sphere generation
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Aff_transformation_3.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/alpha_wrap_3.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>
*/


#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <math.h>
#include <Eigen/Dense>
#include <boost/multi_array.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classification.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/bounding_box.h>

#include <CGAL/Kd_tree.h>
#include <CGAL/Splitters.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>

#include <CGAL/Real_timer.h>



namespace Classification = CGAL::Classification;

typedef Classification::Label_handle                                            Label_handle;
typedef Classification::Feature_handle                                          Feature_handle;
typedef Classification::Label_set                                               Label_set;
typedef Classification::Feature_set                                             Feature_set;
typedef Classification::Point_set_feature_generator<Kernel, Point_set, Pmap>    Feature_generator;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix2d;
typedef Eigen::Matrix<double, 1, Eigen::Dynamic> VectorXd;
typedef Eigen::Matrix<int, 1, Eigen::Dynamic> VectorXi;
typedef boost::multi_array<int, 3> Matrix3i;
typedef Matrix3i::array_view<3>::type Matrix3iV;
typedef Matrix3i::index_range range;



namespace CGAL {
namespace   Classification {
namespace     Feature {


class Fractal_0 : public CGAL::Classification::Feature_base {
    std::vector<float> d0;
public:
    Fractal_0(std::vector<float>& d0) : d0(d0) {
        this->set_name("fractal0");
    }

    float value(std::size_t pt_index) {
        return d0[pt_index];
    }
};

class Fractal_1 : public CGAL::Classification::Feature_base {
    std::vector<float> d1;
public:
    Fractal_1(std::vector<float>& d1) : d1(d1) {
        this->set_name("fractal1");
    }

    float value(std::size_t pt_index) {
        return d1[pt_index];
    }
};

class Fractal_2 : public CGAL::Classification::Feature_base {
    std::vector<float> d2;
public:
    Fractal_2(std::vector<float>& d2) : d2(d2) {
        this->set_name("fractal2");
    }

    float value(std::size_t pt_index) {
        return d2[pt_index];
    }
};

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

template <typename GeomTraits, typename PointRange, typename PointMap>
std::tuple< std::vector<float>, std::vector<float>, std::vector<float> > compute_fractal_dimensionality(const PointRange& input, PointMap point_map, float feature_scale) {
    typedef GeomTraits::Point_3 Point;
    typedef PointRange Point_set;
    typedef GeomTraits::Iso_cuboid_3 Iso_cuboid_3;

    typedef PointRange::Point_map Pmap;
    typedef PointRange::Property_map<int> Imap;
    typedef PointRange::Property_map<unsigned char> UCmap;
    typedef PointRange::Property_map<double> Dmap;

    typedef CGAL::Real_timer Timer;

    typedef CGAL::Search_traits_3<GeomTraits> SearchTraits_3;
    typedef CGAL::Sliding_midpoint<SearchTraits_3> Splitter;
    typedef CGAL::Kd_tree<SearchTraits_3, Splitter, CGAL::Tag_true> Tree;
    typedef CGAL::Fuzzy_iso_box<SearchTraits_3> Iso_box;
    typedef CGAL::Fuzzy_sphere<SearchTraits_3> Sphere;

    using FloatMap = typename PointRange::template Property_map<float>;
    //std::vector<typename FloatMap::value_type> d0, d1, d2;
    //std::vector<float> d0, d1, d2;
    
    double l = 6; // the length of neighborhoods (default: 6) // RMD: Should this be scaled by the "feature_scale" parameter?
    float r1 = feature_scale * 1.0f; // the smallest ratio (default: 0.3)
    float r2 = feature_scale * 10.0f; // the biggest ratio (default: 0.8)
    int n = 10; // the number of ratios (default: 30)
    double step = (r2 - r1) / n;

    std::cout << "feature scale: " << feature_scale << std::endl;
    std::cerr << "l : " << l << std::endl;
    std::cerr << "r1: " << r1 << std::endl;
    std::cerr << "r2: " << r2 << std::endl;
    std::cerr << "n: " << n << std::endl;
    std::cout << "step: " << step << std::endl;

    Timer t, t_total;
    t_total.reset();
    t_total.start();

    // Build K-D Tree
    Tree m_tree;
    std::cerr << "Building Tree..." << std::endl;
    for (auto pt = point_map.begin(); pt != point_map.end(); pt++)
        m_tree.insert(*pt);
    std::cerr << "Tree is built!" << std::endl;

    double gb_xmin = m_tree.bounding_box().min_coord(0), gb_xmax = m_tree.bounding_box().max_coord(0);
    double gb_ymin = m_tree.bounding_box().min_coord(1), gb_ymax = m_tree.bounding_box().max_coord(1);
    double gb_zmin = m_tree.bounding_box().min_coord(2), gb_zmax = m_tree.bounding_box().max_coord(2);

    t.reset();
    t.start();
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
    t.stop();
    std::cerr << "Found neighbors in " << t.time() << "s" << std::endl;

    for (int i = 0; i < n; i++) {
        std::cerr << i << std::endl;

        t.reset();
        t.start();
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
        t.stop();
        std::cerr << "Build graph in " << t.time() << "s" << std::endl;

        global_index = 0;
        Matrix2d limits(1, 6);
        limits.row(0) << gb_xsize - 1, gb_ysize - 1, gb_zsize - 1, gb_xsize - 1, gb_ysize - 1, gb_zsize - 1;


        t.reset();
        t.start();
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
        t.stop();
        std::cerr << "Calculate Fractal Dimension in " << t.time() << "s" << std::endl;
    }

    // Pair Counting
    std::vector<float> D0, D1, D2;
    D0.resize(input.size(), 0.f);
    D1.resize(input.size(), 0.f);
    D2.resize(input.size(), 0.f);
    t.reset();
    t.start();
    Matrix2d combinations(n * (n - 1) / 2, 4);
    for (int i = 0; i < input.size(); i++) {
        int idx = 0;
        for (int j = 0; j < n - 1; j++)
            for (int k = j + 1; k < n; k++)
                combinations.row(idx++) = results[i].row(k) - results[i].row(j);

        double thresh = 0.1;
        combinations.col(1) = combinations.col(1).array() / combinations.col(0).array();
        combinations.col(2) = combinations.col(2).array() / combinations.col(0).array();
        combinations.col(3) = combinations.col(3).array() / combinations.col(0).array();
        combinations = (combinations.array().abs() < thresh).select(Matrix2d::Zero(n * (n - 1) / 2, 4), combinations);

        D0[i] = -10 * combinations.col(1).sum() / (combinations.col(1).array() != 0).count();
        D1[i] = 10 * combinations.col(2).sum() / (combinations.col(2).array() != 0).count();
        D2[i] = 10 * combinations.col(3).sum() / (combinations.col(3).array() != 0).count();
    }
    t.stop();
    t_total.stop();
    std::cerr << "Pair Counting in " << t.time() << "s" << std::endl;


    std::cerr << "All done in " << t_total.time() << "s" << std::endl;
    return std::make_tuple(D0, D1, D2);
}


}}}