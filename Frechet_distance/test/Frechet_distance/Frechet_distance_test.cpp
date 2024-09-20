#include <CGAL/Frechet_distance.h>
#include <CGAL/Frechet_distance_traits_2.h>
#include <CGAL/Frechet_distance_near_neighbors_ds.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Timer.h>

#include <algorithm>
#include <cassert>
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{

//
// helpers
//
using Test_distance_t = double;
//using TestKernel = CGAL::Simple_cartesian<double>;
//using TestKernel = CGAL::Simple_cartesian<CGAL::Exact_rational>;
using TestKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
//using TestKernel = CGAL::Exact_predicates_exact_constructions_kernel;
//using TestKernel = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;

using TestTraits = CGAL::Frechet_distance_traits_2<TestKernel>;
using TestPoint = TestKernel::Point_2;
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


void readCurve(std::ifstream& curve_file, TestCurve& curve)
{
    // Read everything into a stringstream.
    std::stringstream ss;
    ss << curve_file.rdbuf();
    CGAL::set_ascii_mode(ss);

    TestPoint p;
    auto ignore_count = std::numeric_limits<std::streamsize>::max();
    while (ss >> p) {
        ss.ignore(ignore_count, '\n');

        if ((!curve.empty()) && (p == curve.back())) {
            continue;
        }
        curve.push_back(p);
    }
}

TestCurves readCurves(std::string const& curve_directory)
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


FrechetDistanceQueries readFrechetDistanceQueries(std::string const& query_file)
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

FrechetDistanceNearNeighborsDSQueries readFrechetDistanceNearNeighborsDSQueries(
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

                CGAL::Frechet_distance_::internal::CurveID result_id;
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

void testFrechetDistance()
{
    std::string curve_directory = "./data/curves/";
    // std::vector<std::string> datasets = {"sigspatial", "OV"};
    std::vector<std::string> datasets = { "sigspatial" };
    std::string query_directory = "./data/queries/";
    CGAL::Timer timer;
    for (auto const& dataset : datasets) {
        auto curves = readCurves(curve_directory + dataset + "/");
        auto queries =
            readFrechetDistanceQueries(query_directory + dataset + ".txt");

        for (auto const& query : queries) {
          /*
                        std::cout
                            << CGAL::approximate_Frechet_distance<TestCurve,
                                                              TestTraits>(
                                   curves[query.id1], curves[query.id2], 0.001)
                            << std::endl;
          */
            timer.start();
            auto decision =
                ! CGAL::is_Frechet_distance_larger<TestTraits>(
                    curves[query.id1], curves[query.id2], query.distance);
            timer.stop();
            if (decision != query.decision) {
                std::cout << "Wrong decision on query." << std::endl;
                exit(- 1);
            }
        }
    }
    std::cout << timer.time() << "sec." << std::endl;
}

void testFrechetDistanceNearNeighborsDS()
{
        std::string curve_directory = "./data/curves/";
        std::vector<std::string> datasets = { "sigspatial" };
        std::string query_directory = "./data/ds_queries/";

        for (auto const& dataset: datasets) {
                auto curves = readCurves(curve_directory + dataset + "/");
                auto queries =
                readFrechetDistanceNearNeighborsDSQueries(query_directory + dataset +
                ".txt");

                CGAL::FrechetDistanceNearNeighborsDS<TestCurve, TestTraits> ds;
                ds.insert(curves);

                for (auto const& query: queries) {
                        auto result = ds.get_close_curves(curves[query.id],
                        query.distance); std::sort(result.begin(), result.end());

                        if (!std::equal(result.begin(), result.end(), query.expected_result.begin(), query.expected_result.end())) {
                                std::cout << "Wrong result on query." << std::endl;
                                exit(- 1);
                        }
                }
        }

}

}  // end anonymous namespace


int main()
{
    // TODO: add actualy query data for DS
      //  std::cout << "testFrechetDistanceNearNeighborsDS start" << std::endl;
      testFrechetDistanceNearNeighborsDS();
      //  std::cout << "testFrechetDistanceNearNeighborsDS done" << std::endl;
    std::cout << "testFrechetDistance start" << std::endl;
    testFrechetDistance();
    std::cout << "testFrechetDistance done" << std::endl;
    return 0;
}
