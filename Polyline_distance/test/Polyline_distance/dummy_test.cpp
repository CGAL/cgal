#include <CGAL/Frechet_distance.h>
#include <CGAL/Cartesian.h>

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

//
// helpers
//

using Kernel = CGAL::Cartesian<double>;
// using NT = Kernel::FT;
using Point = Kernel::Point_2;
using Curve = std::vector<Point>;
using Curves = std::vector<Curve>;

struct Query {
	std::size_t id1, id2;
	distance_t distance;
	bool decision;
};
using Queries = std::vector<Query>;

void readCurve(std::ifstream& curve_file, Curve& curve)
{
	// Read everything into a stringstream.
	std::stringstream ss;
	ss << curve_file.rdbuf();
	CGAL::set_ascii_mode(ss);

	Point p;
	auto ignore_count = std::numeric_limits<std::streamsize>::max();
	while (ss >> p) {
		ss.ignore(ignore_count, '\n');

		if (!curve.empty() && p == curve.back()) {
			continue;
		}
		curve.push_back(p);
	}
}

Curves readCurves(std::string const& curve_directory)
{
	Curves curves;
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
	for (auto const& curve_filename: curve_filenames) {
		std::ifstream curve_file(curve_directory + curve_filename);
		assert(curve_file);

		curves.emplace_back();
		readCurve(curve_file, curves.back());

		if (curves.back().empty()) { curves.pop_back(); }
	}

	return curves;
}

Queries readQueries(std::string const& query_file)
{
	Queries queries;

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

//
// tests
//

int testCorrectness()
{
	std::string curve_directory = "data/curves/";
	std::vector<std::string> datasets = { "sigspatial", "OV" };
	std::string query_directory = "data/queries/";

	for (auto const& dataset: datasets) {
		auto curves = readCurves(curve_directory + dataset + "/");
		auto queries = readQueries(query_directory + dataset + ".txt");

		for (auto const& query: queries) {
			auto decision = continuous_Frechet_distance_less_than(curves[query.id1], curves[query.id2], query.distance);
			assert(decision == query.decision);
		}
	}
}

} // end anonymous namespace

int main()
{
	testCorrectness();
}
