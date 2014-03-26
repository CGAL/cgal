//A demo for the scale space.
//Copyright (C) 2013  INRIA - Sophia Antipolis
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Author(s):      Thijs van Lankveld


#include <algorithm>
#include <vector>
#include <string>

//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "readFile.h"
#include "writeFile.h"
#include "quad.h"

#include "Scale_space_surface_constructer.h"


int main(int argc, char** argv) {
	if (argc > 1 && (!std::strcmp(argv[1], "-help") || !std::strcmp(argv[1], "--help"))) {
		std::cout << "usage:\nDigne.exe [input].xyzn [outfile] [#neighbors] [#iterations] [#samples for r-estimate]" << std::endl;
		exit(0);
	}

	char* input = "kitten.xyzn";
	char* outFile = "kitten_o";
	unsigned int nn = 30;
	unsigned int iter = 4;
	unsigned int samples = 400;
	if (argc > 1)
		input = argv[1];
	if (argc > 2)
		outFile = argv[2];
	if (argc > 3)
		sscanf(argv[3], "%u", &nn);
	if (argc > 4)
		sscanf(argv[4], "%u", &iter);
	if (argc > 5)
		sscanf(argv[5], "%u", &samples);

	//typedef CGAL::Exact_predicates_exact_constructions_kernel			Kernel;
	typedef CGAL::Exact_predicates_inexact_constructions_kernel			Kernel;

	typedef Kernel::Point_3												Point_3;
	typedef Kernel::Vector_3											Vector_3;
	typedef Kernel::Triangle_3											Triangle_3;

	typedef std::vector<Point_3>										PointCollection;
	typedef std::vector<Vector_3>										VectorCollection;
	typedef float4														Color;
	typedef std::vector<Color>											ColorCollection;

	// Read the input file.
	PointCollection points;
	VectorCollection normals;
	ColorCollection colors;
	std::string is(input);
	std::transform(is.begin(), is.end(), is.begin(), std::tolower);
	bool loaded = false;
	if (is.rfind(".xyzn") == is.length()-5)
		loaded = loadXYZN<Kernel>(input, std::back_inserter(points), std::back_inserter(normals));
/*	else if (is.rfind(".ive") == is.length()-4 || is.rfind(".osg") == is.length()-4)
		loaded = loadOSG<Kernel, Color>(input, std::back_inserter(points), std::back_inserter(normals), std::back_inserter(colors));*/
	if (!loaded) {
		std::cout << "Error loading input." << std::endl;
		exit(1);
	}
	std::cout << "Input loaded:" << std::endl;
	std::cout << points.size() << " points, " << normals.size() << " normals, " << colors.size() << " colors" << std::endl;
	if (points.size() != normals.size())
		normals.resize(points.size());

	// Construct the mesh in a scale space.
	typedef Scale_space_surface_constructer<Kernel>		Scale_space_surface_constructer;
	typedef Scale_space_surface_constructer::Triangle	Triangle;
	std::list<Triangle> triangles;
	Scale_space_surface_constructer sssc(iter);
	sssc(points.begin(), points.end(), std::back_inserter(triangles));

	std::cout << "Save to: " << outFile << std::endl;
	std::string fs(outFile);

	// Save the points.
	std::cout << "Save as xyzn" << std::endl;
	std::ofstream fout;
	fout.open((fs+".xyzn").c_str());
	if (fout) {
		VectorCollection::const_iterator nit = normals.begin();
		for (Scale_space_surface_constructer::PointCollection::const_iterator pit = sssc.moved().begin(); pit != sssc.moved().end(); ++pit, ++nit) {
			fout << *pit << " " << *nit << std::endl;
		}
		fout.close();
	}
	else
		std::cout << "Error writing to .xyzn output." << std::endl;
	
	
	// Construct the normals of the triangles.
	Kernel::Construct_cross_product_vector_3 cross;
	VectorCollection tri_norm;
	tri_norm.reserve(triangles.size());
	for (std::list<Triangle>::const_iterator tit = triangles.begin(); tit != triangles.end(); ++tit) {
		const Point_3& p1 = tit->vertex(0);
		const Point_3& p2 = tit->vertex(1);
		const Point_3& p3 = tit->vertex(1);

		Vector_3 n = cross((p2 - p1), (p3 - p1));
		tri_norm.push_back(n / CGAL::sqrt(CGAL::to_double(n.squared_length())));
	}

	// Link the colors of the points.
	std::map<Point_3, Color> toColor;
	PointCollection::const_iterator pit = points.begin();
	for (ColorCollection::const_iterator cit = colors.begin(); cit != colors.end(); ++cit)
		toColor[*pit++] = *cit;
	ColorCollection tri_colors;
	tri_colors.reserve(3*triangles.size());
	if (!toColor.empty())
		for (std::list<Triangle>::const_iterator tit = triangles.begin(); tit != triangles.end(); ++tit)
			for (unsigned int i = 0; i < 3; ++i)
				tri_colors.push_back(toColor[tit->vertex(i)]);
	
	// Save the shapes (OSG style).
/*	bool saveIve = false;
	std::cout << "Save as IVE" << std::endl;
	osg::ref_ptr<osg::Switch> shape = new osg::Switch;
	if (savePointsIVE(*shape, sssc.moved().begin(), sssc.moved().end(), normals.begin(), normals.end(), colors.begin(), colors.end()))
		if (saveTrianglesIVE(*shape, facets.begin(), facets.end(), tri_norm.begin(), tri_norm.end(), tri_colors.begin(), tri_colors.end()))
			if (saveIVE(shape.get(), (fs+".ive").c_str()))
				saveIve = true;
	if (!saveIve)
		std::cout << "Error writing to .ive output." << std::endl;*/
	
	// Save the shapes (OFF style).
	std::cout << "Save as OFF" << std::endl;
	if (!saveOFF((fs+".off").c_str(), triangles.begin(), triangles.end(), tri_norm.begin(), tri_norm.end(), colors.begin(), colors.end()))
		std::cout << "Error writing to .off output." << std::endl;
}