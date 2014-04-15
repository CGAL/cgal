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
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <CGAL/Scale_space_surface_reconstructer_3.h>

//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

//typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef CGAL::Scale_space_surface_reconstructer_3< Kernel > Reconstructer;

typedef Reconstructer::Point                                Point;
typedef std::vector< Point >                                Pointset;

typedef Reconstructer::Tripleset	                        Tripleset;

const unsigned int LINE_SIZE = 1024;


void readLine( std::ifstream& fin, char* line ) {
	fin.getline( line, LINE_SIZE );
	while( fin && line[0] == '#' ) {
		// Comment line.
		std::cout << "Comment: " << line << std::endl;
		fin.getline( line, LINE_SIZE );
	}
}

bool readOFF( std::string input, Pointset& points ) {
	// Read the input file.
    std::ifstream fin( input.c_str() );

	// Any .off file should start with the line "OFF"
    char line[LINE_SIZE];
	readLine( fin, line );
	if( std::strncmp(line, "OFF", 3) != 0 ) return false;

	// The next line contains the object counts.
	unsigned int n_points, n_polygons, n_segments;
	fin >> n_points >> n_polygons >> n_segments;

	// Collect the points.
    points.reserve( n_points );
	for( unsigned int n = 0; fin && n < n_points; ++n ) {
        Point p;
		fin >> p;
        points.push_back( p );
	}
	if( !fin ) return false;

    // Ignore the rest
    fin.close();
    return true;
}

bool writeOFF( std::string output, const Pointset& points, const Tripleset& triples ) {
    // Write the output file.
	std::ofstream fout( output );
	fout << "OFF" << std::endl;
    fout << points.size() << " " << triples.size() << " 0" << std::endl;

	// Write the points.
	for( Pointset::const_iterator pit = points.begin(); pit != points.end(); ++pit)
		fout << *pit << std::endl;

	// Write the triples.
    for( Tripleset::const_iterator tit = triples.begin(); tit != triples.end(); ++tit )
		fout << "3  " << *tit << std::endl;

	return true;
}

int main(int argc, char** argv) {
	std::string base = "kitten";
	unsigned int neighbors = 30;
	unsigned int iterations = 4;
	unsigned int samples = 200;
	if (argc > 1)
		base = argv[1];
	if (argc > 2)
		sscanf(argv[2], "%u", &neighbors);
	if (argc > 3)
		sscanf(argv[3], "%u", &iterations);
	if (argc > 4)
		sscanf(argv[4], "%u", &samples);
    
    std::string input = base + ".off";
	std::string output_ss = base + "_ss.off";
	std::string output_sm = base + "_sm.off";
    
    // Read the data.
    std::cout << "Input: " << input << std::flush;
	Pointset points;
    if( !readOFF( input, points ) ) {
        std::cerr << std::endl << "Error reading " << input << std::endl;
        exit(-1);
    }
	std::cout << " loaded: " << points.size() << " points." << std::endl;

	// Construct the mesh in a scale space.
	Reconstructer reconstruct;

	reconstruct.compute_surface( points.begin(), points.end(), neighbors, iterations, samples );
    std::cout << "Reconstruction done." << std::endl;

    // Write the reconstruction.
    std::cout << "Output: " << output_ss << std::flush;
    if( !writeOFF( output_ss, points, reconstruct.surface() ) ) {
        std::cerr << std::endl << "Error writing " << output_ss << std::endl;
        exit(-1);
    }
	std::cout << " written." << std::endl;

    // Write the reconstruction.
    std::cout << "Output: " << output_sm << std::flush;
    Pointset smoothed( reconstruct.scale_space_begin(), reconstruct.scale_space_end() );
    if( !writeOFF( output_sm, smoothed, reconstruct.surface() ) ) {
        std::cerr << std::endl << "Error writing " << output_ss << std::endl;
        exit(-1);
    }
	std::cout << " written." << std::endl;
}