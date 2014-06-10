#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <CGAL/Scale_space_surface_reconstruction_3.h>

//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

//typedef CGAL::Exact_predicates_exact_constructions_kernel     Kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;

typedef CGAL::Scale_space_surface_reconstruction_3< Kernel >    Reconstruction;

typedef Reconstruction::Point                                   Point;
typedef std::vector< Point >                                    Pointset;

typedef Reconstruction::Const_triple_iterator	                TripleIterator;

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

bool writeOFF( std::string output, const Pointset& points, TripleIterator triples_begin, TripleIterator triples_end, size_t n_triples ) {
    // Write the output file.
	std::ofstream fout( output.c_str() );
	fout << "OFF" << std::endl;
    fout << points.size() << " " << n_triples << " 0" << std::endl;

	// Write the points.
	for( Pointset::const_iterator pit = points.begin(); pit != points.end(); ++pit)
		fout << *pit << std::endl;

	// Write the triples.
    for( TripleIterator tit = triples_begin; tit != triples_end; ++tit ) {
		fout << "3  " << *tit << std::endl;
    }

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
	Reconstruction reconstruct( neighbors, samples );
	reconstruct.reconstruct_surface( points.begin(), points.end(), iterations );
    std::cout << "Reconstruction done." << std::endl;

    // Write the reconstruction.
    std::cout << "Output: " << output_ss << std::flush;
    if( !writeOFF( output_ss, points, reconstruct.surface_begin(), reconstruct.surface_end(), reconstruct.number_of_triangles() ) ) {
        std::cerr << std::endl << "Error writing " << output_ss << std::endl;
        exit(-1);
    }
	std::cout << " written." << std::endl;

    // Write the reconstruction.
    std::cout << "Output: " << output_sm << std::flush;
    Pointset smoothed( reconstruct.scale_space_begin(), reconstruct.scale_space_end() );
    if( !writeOFF( output_sm, smoothed, reconstruct.surface_begin(), reconstruct.surface_end(), reconstruct.number_of_triangles() ) ) {
        std::cerr << std::endl << "Error writing " << output_ss << std::endl;
        exit(-1);
    }
	std::cout << " written." << std::endl;
}