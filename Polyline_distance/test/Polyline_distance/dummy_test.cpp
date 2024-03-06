// Copyright (c) 2019 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : André Nusser <anusser@mpi-inf.mpg.de>
//                 Marvin Künnemann <marvin@mpi-inf.mpg.de>
//                 Karl Bringmann <kbringma@mpi-inf.mpg.de>
//
// =============================================================================

#include <CGAL/Frechet_distance.h>
#include <CGAL/Frechet_distance_near_neighbors_ds.h>
#include <CGAL/Simple_cartesian.h>
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
using distance_t = double;
using Kernel = CGAL::Simple_cartesian<double>;
using Traits = CGAL::Polyline_traits_2<Kernel, double>;
// using NT = Kernel::FT;
using Point = Kernel::Point_2;
using TestCurve = std::vector<Point>;
using TestCurves = std::vector<TestCurve>;

struct FrechetDistanceQuery {
    std::size_t id1, id2;
    distance_t distance;
    bool decision;
};
using FrechetDistanceQueries = std::vector<FrechetDistanceQuery>;

struct FrechetDistanceNearNeighborsDSQuery {
    std::size_t id;
    distance_t distance;
    std::vector<std::size_t> expected_result;  // TODO: should be curve ids
};
using FrechetDistanceNearNeighborsDSQueries =
    std::vector<FrechetDistanceNearNeighborsDSQuery>;

void readCurve(std::ifstream& curve_file, TestCurve& curve)
{
    // Read everything into a stringstream.
    std::stringstream ss;
    ss << curve_file.rdbuf();
    CGAL::set_ascii_mode(ss);

    Point p;
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

        CGAL::Polyline_distance::internal::CurveID result_id;
        while (ss >> result_id) {
            query.expected_result.push_back(result_id);
        }
        std::sort(query.expected_result.begin(), query.expected_result.end());
    }

    return queries;
}

//
// tests
//

void testFrechetDistance()
{
    std::string curve_directory = "../data/curves/";
    // std::vector<std::string> datasets = {"sigspatial", "OV"};
    std::vector<std::string> datasets = {"sigspatial"};
    std::string query_directory = "../data/queries/";
    CGAL::Timer timer;
    for (auto const& dataset : datasets) {
        auto curves = readCurves(curve_directory + dataset + "/");
        auto queries =
            readFrechetDistanceQueries(query_directory + dataset + ".txt");

        for (auto const& query : queries) {
            timer.start();
            auto decision =
                CGAL::continuous_Frechet_distance_less_than<TestCurve, Traits>(
                    curves[query.id1], curves[query.id2], query.distance);
            timer.stop();
            if (decision != query.decision) {
                assert(false);
                ERROR("Wrong decision on query.");
            }
        }
    }
    std::cout << timer.time() << "sec." << std::endl;
}

void testFrechetDistanceNearNeighborsDS()
{
    std::string curve_directory = "../data/curves/";
    std::vector<std::string> datasets = {"sigspatial"};
    std::string query_directory = "../data/ds_queries/";

    for (auto const& dataset : datasets) {
        auto curves = readCurves(curve_directory + dataset + "/");
        auto queries = readFrechetDistanceNearNeighborsDSQueries(
            query_directory + dataset + ".txt");

        CGAL::FrechetDistanceNearNeighborsDS<TestCurve> ds;
        ds.insert(curves);

        for (auto const& query : queries) {
            auto result = ds.get_close_curves(curves[query.id], query.distance);
            std::sort(result.begin(), result.end());

            if (!std::equal(result.begin(), result.end(),
                            query.expected_result.begin(),
                            query.expected_result.end())) {
                assert(false);
                ERROR("Wrong result on query.");
            }
        }
    }
}

}  // end anonymous namespace

int main()
{
    // TODO: add actualy query data for DS
    std::cout << "testFrechetDistanceNearNeighborsDS start" << std::endl;
    testFrechetDistanceNearNeighborsDS();
    std::cout << "testFrechetDistanceNearNeighborsDS done" << std::endl;
    std::cout << "testFrechetDistance start" << std::endl;
    testFrechetDistance();
    std::cout << "testFrechetDistance done" << std::endl;
}
