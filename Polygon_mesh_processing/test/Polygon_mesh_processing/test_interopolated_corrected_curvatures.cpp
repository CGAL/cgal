#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/Curvatures/interpolated_corrected_curvature_measures.h>
#include <CGAL/Polygon_mesh_processing/random_perturbation.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>

#include <boost/graph/graph_traits.hpp>

#include <iostream>
#include <unordered_map>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel EpicKernel;
typedef CGAL::Surface_mesh<EpicKernel::Point_3> SMesh;
typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<SMesh>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;

void test(std::string mesh_path, EpicKernel::FT rel_expansion_radius, EpicKernel::FT rel_noise_magnitude) {
	SMesh pmesh;
	const std::string filename = CGAL::data_file_path(mesh_path);

	if (!CGAL::IO::read_polygon_mesh(filename, pmesh))
	{
		std::cerr << "Invalid input file." << std::endl;
	}

	bool created = false;

	SMesh::Property_map<vertex_descriptor, EpicKernel::FT> mean_curvature_map, gaussian_curvature_map;
	boost::tie(mean_curvature_map, created) = pmesh.add_property_map<vertex_descriptor, EpicKernel::FT>("v:mean_curvature_map", 0);
	assert(created);

	boost::tie(gaussian_curvature_map, created) = pmesh.add_property_map<vertex_descriptor, EpicKernel::FT>("v:gaussian_curvature_map", 0);
	assert(created);

	// getting the max and min edge lengthes
	const auto edge_range = CGAL::edges(pmesh);

	const auto edge_length_comparator = [&, pmesh](auto l, auto r) {
		return PMP::edge_length(l, pmesh) < PMP::edge_length(r, pmesh);
	};

	const edge_descriptor longest_edge = *std::max_element(edge_range.begin(), edge_range.end(), edge_length_comparator);
	const EpicKernel::FT max_edge_length = PMP::edge_length(longest_edge, pmesh);

	const edge_descriptor shortest_edge = *std::min_element(edge_range.begin(), edge_range.end(), edge_length_comparator);
	const EpicKernel::FT min_edge_length = PMP::edge_length(shortest_edge, pmesh);


	if (rel_noise_magnitude > 0)
	{
		if (!CGAL::is_triangle_mesh(pmesh))
			return;

		SMesh::Property_map<vertex_descriptor, EpicKernel::Vector_3> vnm;
		boost::tie(vnm, created) = pmesh.add_property_map<vertex_descriptor, EpicKernel::Vector_3>("v:vnm", { 0 , 0 , 0 });
		assert(created);
		
		CGAL::Polygon_mesh_processing::compute_vertex_normals(pmesh, vnm);

		PMP::random_perturbation(pmesh, rel_noise_magnitude * min_edge_length, CGAL::parameters::random_seed(0));
		PMP::interpolated_corrected_mean_curvature(
			pmesh,
			mean_curvature_map,
			CGAL::parameters::ball_radius(rel_expansion_radius * max_edge_length).vertex_normal_map(vnm)
		);
		PMP::interpolated_corrected_gaussian_curvature(
			pmesh,
			gaussian_curvature_map,
			CGAL::parameters::ball_radius(rel_expansion_radius * max_edge_length).vertex_normal_map(vnm)
		);
	}
	else {
		PMP::interpolated_corrected_mean_curvature(
			pmesh,
			mean_curvature_map,
			CGAL::parameters::ball_radius(rel_expansion_radius * max_edge_length)
		);
		PMP::interpolated_corrected_gaussian_curvature(
			pmesh,
			gaussian_curvature_map,
			CGAL::parameters::ball_radius(rel_expansion_radius * max_edge_length)
		);

	}

	

	//PMP::interpolated_corrected_mean_curvature(
	//	pmesh,
	//	mean_curvature_map,
	//	CGAL::parameters::ball_radius(rel_expansion_radius * max_edge_length)
	//);
	//PMP::interpolated_corrected_gaussian_curvature(
	//	pmesh,
	//	gaussian_curvature_map,
	//	CGAL::parameters::ball_radius(rel_expansion_radius * max_edge_length)
	//);


	const EpicKernel::FT max_mean_curvature = *std::max_element(mean_curvature_map.begin(), mean_curvature_map.end());
	const EpicKernel::FT min_mean_curvature = *std::min_element(mean_curvature_map.begin(), mean_curvature_map.end());
	const EpicKernel::FT max_gaussian_curvature = *std::max_element(gaussian_curvature_map.begin(), gaussian_curvature_map.end());
	const EpicKernel::FT min_gaussian_curvature = *std::min_element(gaussian_curvature_map.begin(), gaussian_curvature_map.end());

	std::cout << "# " << mesh_path << ":\n" 
		<< "expansion radius ratio to max length / expansion radius = " << rel_expansion_radius << " / " << rel_expansion_radius * max_edge_length << ",\n" 
		<< "max perturbation ratio to minlength / max perturbation = " << rel_noise_magnitude << " / " << rel_noise_magnitude * min_edge_length << "\n"
		<< "mean curvature: min = " << min_mean_curvature << " <-> " << max_mean_curvature << " = max" << "\n"
		<< "gaussian curvature: min = " << min_gaussian_curvature << " <-> " << max_gaussian_curvature << " = max" << "\n\n\n";


}

int main()
{
	const std::vector<std::string> mesh_pathes_to_test = {
		"meshes/icc_test/Sphere Quads + Tris.obj",
		"meshes/icc_test/Sphere Quads + Tris 100352.obj",
		"meshes/icc_test/Sphere Tris Ico.obj",
		"meshes/icc_test/Sphere Tris Tet.obj",
		"meshes/icc_test/Sphere Tris Oct.obj",
		"meshes/icc_test/Sphere Quads.obj",
		"meshes/icc_test/Sphere Quads Remesh.obj",
		"meshes/icc_test/Sphere Ngons + Quads + Tris.obj",
		"meshes/icc_test/Cube with fillet Quads.obj",
		"meshes/cylinder.off",
		"meshes/icc_test/Lantern Tris.obj",
		"meshes/icc_test/Lantern Quads.obj"
	};

	const std::vector<EpicKernel::FT> rel_expansion_radii = { 0, 0.1, 0.5, 1 };
	const std::vector<EpicKernel::FT> rel_noise_magnitudes = { 0, 0.5, 0.9 };

	for (auto mesh_path : mesh_pathes_to_test) {
		for (EpicKernel::FT rel_expansion_radius : rel_expansion_radii)
			for (EpicKernel::FT rel_noise_magnitude : rel_noise_magnitudes)
			{
				test(mesh_path, rel_expansion_radius, rel_noise_magnitude);
			}

		std::cout << "_________________________________________________________________________________\n\n";
	}

}
