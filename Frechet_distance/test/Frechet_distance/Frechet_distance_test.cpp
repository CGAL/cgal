#include <CGAL/Frechet_distance.h>
#include <CGAL/Frechet_distance_traits_2.h>
#include <CGAL/Frechet_distance/Neighbor_search.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Epeck_d.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Dimension.h>

#include <algorithm>
#include <cassert>
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

template <class TestKernel, class TestTraits, class TestPoint>
struct Test_struct
{

//
// helpers
//
using Test_distance_t = double;
using TestCurve = std::vector<TestPoint>;
using TestCurves = std::vector<TestCurve>;

struct FrechetDistanceQuery {
    std::size_t id1, id2;
    Test_distance_t distance;
    bool decision;
};
using FrechetDistanceQueries = std::vector<FrechetDistanceQuery>;

struct FrechetDistanceNearNeighborsDSQuery {
        std::size_t id;
        Test_distance_t distance;
        std::vector<std::size_t> expected_result; // TODO: should be curve ids
};
using FrechetDistanceNearNeighborsDSQueries =
    std::vector<FrechetDistanceNearNeighborsDSQuery>;


static void readCurve(std::ifstream& curve_file, TestCurve& curve)
{
    // Read everything into a stringstream.
    std::stringstream ss;
    ss << curve_file.rdbuf();
    CGAL::set_ascii_mode(ss);

    TestPoint p;
    auto ignore_count = (std::numeric_limits<std::streamsize>::max)();
    while (ss >> p) {
        ss.ignore(ignore_count, '\n');

        if ((!curve.empty()) && (p == curve.back())) {
            continue;
        }
        curve.push_back(p);
    }
}

static TestCurves readCurves(std::string const& curve_directory)
{
    TestCurves curves;
    std::vector<std::string> curve_filenames;

    // read filenames of curve files
    std::ifstream file(curve_directory + "dataset.txt");
    assert(file);

    std::string line;
    while (std::getline(file, line)) {
        curve_filenames.push_back(line);
    }

    // read curves
    curves.reserve(curve_filenames.size());
    for (auto const& curve_filename : curve_filenames) {
        std::ifstream curve_file(curve_directory + curve_filename);
        assert(curve_file);

        curves.emplace_back();
        readCurve(curve_file, curves.back());

        if (curves.back().empty()) {
            curves.pop_back();
        }
    }

    return curves;
}


static FrechetDistanceQueries readFrechetDistanceQueries(std::string const& query_file)
{
    FrechetDistanceQueries queries;

    std::ifstream file(query_file);
    assert(file);

    std::string line;
    while (std::getline(file, line)) {
        queries.emplace_back();
        auto& query = queries.back();

        std::stringstream ss(line);
        ss >> query.id1 >> query.id2 >> query.distance >> query.decision;
    }

    return queries;
}

static FrechetDistanceNearNeighborsDSQueries readFrechetDistanceNearNeighborsDSQueries(
    std::string const& query_file)
{
        FrechetDistanceNearNeighborsDSQueries queries;

        std::ifstream file(query_file);
        assert(file);

        std::string line;
        while (std::getline(file, line)) {
                queries.emplace_back();
                auto& query = queries.back();

                std::stringstream ss(line);
                ss >> query.id >> query.distance;

                CGAL::Frechet_distance::internal::CurveID result_id;
                while (ss >> result_id) {
                        query.expected_result.push_back(result_id);
                }
                std::sort(query.expected_result.begin(),
                query.expected_result.end());
        }

        return queries;
}


//
// tests
//
template<bool force_filtering=false>
static double testFrechetDistance()
{
    namespace params = CGAL::parameters;
    std::string curve_directory = "./data/curves/";
    std::string query_directory = "./data/queries/";
    std::vector<std::string> datasets;

    auto const dimension = TestTraits::Dimension::value;
    if (dimension == 2) {
        // datasets = {"sigspatial", "OV"};
        datasets = { "sigspatial" };
    }
    else if (dimension == 3) {
        datasets = { "generated_3d" };
    }
    else if (dimension == 100) {
        datasets = { "generated_100d" };
    }

    CGAL::Real_timer timer;
    for (auto const& dataset : datasets) {
        auto curves = readCurves(curve_directory + dataset + "/");
        auto queries =
            readFrechetDistanceQueries(query_directory + dataset + ".txt");

        for (auto const& query : queries) {
          /*
            std::cout << CGAL::bounded_error_Frechet_distance(curves[query.id1], curves[query.id2], 0.001)
                      << std::endl;
          */
            timer.start();
            auto decision =
                ! CGAL::is_Frechet_distance_larger(
                    curves[query.id1], curves[query.id2], query.distance,
                    params::force_filtering(std::bool_constant<force_filtering>())
                           .geom_traits(TestTraits()));
            timer.stop();

            if (decision != query.decision) {
                std::cout << "Wrong decision on query." << std::endl;
                exit(- 1);
            }
        }
    }
    return timer.time();
}

static double testFrechetDistanceNearNeighborsDS()
{
        std::string curve_directory = "./data/curves/";
        std::vector<std::string> datasets = { "sigspatial" };
        std::string query_directory = "./data/ds_queries/";

        CGAL::Real_timer timer;
        for (auto const& dataset: datasets) {
                auto curves = readCurves(curve_directory + dataset + "/");
                auto queries =
                readFrechetDistanceNearNeighborsDSQueries(query_directory + dataset +
                ".txt");

                CGAL::Frechet_distance::Neighbor_search<TestCurve, TestTraits> ds;
                ds.insert(curves);

                for (auto const& query: queries) {
                        auto result = ds.get_close_curves(curves[query.id],
                                                        query.distance);
                         std::sort(result.begin(), result.end());

                        timer.start();
                        if (!std::equal(result.begin(), result.end(), query.expected_result.begin(), query.expected_result.end())) {
                                std::cout << "Wrong result on query." << std::endl;
                                exit(- 1);
                        }
                        timer.stop();
                }
        }

        return timer.time();
}

};

int main(int argc, char** argv)
{
  std::set<int> test_set;
  for (int i=1;i<argc;++i)
    test_set.insert(atoi(argv[i]));

  if (test_set.empty() || test_set.count(0))
  {
    using Kernel = CGAL::Simple_cartesian<double>;
    using Traits = CGAL::Frechet_distance_traits_2<Kernel>;
    using Point = Kernel::Point_2;
    std::cout <<"Simple_cartesian<double>\n";
    double t1=Test_struct<Kernel, Traits, Point>::testFrechetDistanceNearNeighborsDS();
    double t2=Test_struct<Kernel, Traits, Point>::testFrechetDistance();
    std::cout << t1 << " " << t2 << "\n";
  }

  if (test_set.empty() || test_set.count(0))
  {
    using Kernel = CGAL::Simple_cartesian<double>;
    using Traits = CGAL::Frechet_distance_traits_2<Kernel>;
    using Point = Kernel::Point_2;
    std::cout <<"Simple_cartesian<double> (force_filtering)\n";
    double t1=Test_struct<Kernel, Traits, Point>::testFrechetDistanceNearNeighborsDS();
    double t2=Test_struct<Kernel, Traits, Point>::testFrechetDistance<true>();
    std::cout << t1 << " " << t2 << "\n";
  }

  if (test_set.empty() || test_set.count(1))
  {
    using Kernel = CGAL::Epick;
    using Traits = CGAL::Frechet_distance_traits_2<Kernel>;
    using Point = Kernel::Point_2;
    std::cout <<"CGAL::Epick\n";
    double t1=Test_struct<Kernel, Traits, Point>::testFrechetDistanceNearNeighborsDS();
    double t2=Test_struct<Kernel, Traits, Point>::testFrechetDistance();
    std::cout << t1 << " " << t2 << "\n";
  }

  if (test_set.empty() || test_set.count(2))
  {
    using Kernel = CGAL::Epeck;
    using Traits = CGAL::Frechet_distance_traits_2<Kernel>;
    using Point = Kernel::Point_2;
    std::cout <<"CGAL::Epeck\n";
    double t1=Test_struct<Kernel, Traits, Point>::testFrechetDistanceNearNeighborsDS();
    double t2=Test_struct<Kernel, Traits, Point>::testFrechetDistance();
    std::cout << t1 << " " << t2 << "\n";
  }

  if (test_set.empty() || test_set.count(3))
  {
    using Kernel = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;
    using Traits = CGAL::Frechet_distance_traits_2<Kernel>;
    using Point = Kernel::Point_2;
    std::cout <<"Exact_predicates_exact_constructions_kernel_with_sqrt (force filtering)\n";
    double t1=Test_struct<Kernel, Traits, Point>::testFrechetDistanceNearNeighborsDS();
    double t2=Test_struct<Kernel, Traits, Point>::testFrechetDistance<true>();
    std::cout << t1 << " " << t2 << "\n";
  }

  if (test_set.empty() || test_set.count(4))
  {
    using Kernel = CGAL::Simple_cartesian<CGAL::Exact_rational>;
    using Traits = CGAL::Frechet_distance_traits_2<Kernel>;
    using Point = Kernel::Point_2;
    std::cout <<"Simple_cartesian<Exact_rational> (force filtering)\n";
    double t1=Test_struct<Kernel, Traits, Point>::testFrechetDistanceNearNeighborsDS();
    double t2=Test_struct<Kernel, Traits, Point>::testFrechetDistance<true>();
    std::cout << t1 << " " << t2 << "\n";
  }

  if (test_set.empty() || test_set.count(5))
  {
    using Kernel = CGAL::Simple_cartesian<double>;
    using Traits = CGAL::Frechet_distance_traits_3<Kernel>;
    using Point = Kernel::Point_3;
    std::cout <<"Simple_cartesian<double> in 3D\n";
    double t1=Test_struct<Kernel, Traits, Point>::testFrechetDistance();
    std::cout << t1 << "\n";
  }

  if (test_set.empty() || test_set.count(6))
  {
    using Kernel = CGAL::Simple_cartesian<double>;
    using Traits = CGAL::Frechet_distance_traits_3<Kernel>;
    using Point = Kernel::Point_3;
    std::cout <<"Simple_cartesian<double> in 3D (force filtering)\n";
    double t1=Test_struct<Kernel, Traits, Point>::testFrechetDistance<true>();
    std::cout << t1 << "\n";
  }

  if (test_set.empty() || test_set.count(7))
  {
    using Kernel = CGAL::Epick_d<CGAL::Dimension_tag<100>>;
    using Traits = CGAL::Frechet_distance_traits_d<Kernel>;
    using Point = Kernel::Point_d;
    std::cout <<"CGAL::Epick_d\n";
    double t1=Test_struct<Kernel, Traits, Point>::testFrechetDistance();
    std::cout << t1 << "\n";
  }

  if (test_set.empty() || test_set.count(8))
  {
    using Kernel = CGAL::Epeck_d<CGAL::Dimension_tag<100>>;
    using Traits = CGAL::Frechet_distance_traits_d<Kernel>;
    using Point = Kernel::Point_d;
    std::cout <<"CGAL::Epeck_d\n";
    double t1=Test_struct<Kernel, Traits, Point>::testFrechetDistance();
    std::cout << t1 << "\n";
  }

  return 0;
}
