#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include <boost/xpressive/xpressive_static.hpp>
#include <boost/tuple/tuple.hpp>

using namespace boost::xpressive;


typedef boost::tuple<std::string, std::string, std::string> Vertex;
typedef boost::tuple<int, int, int> Face;
typedef std::vector<Vertex> Vertices;
typedef std::vector<Face> Faces;
Vertices vertices;
Faces faces;

int usage(const char* argv0, bool failure = true)
{
  std::cerr << "Usage:\n"
	    << "  " << argv0 << " INPUT.xyz OUTPUT.off\n";
  if(failure)
    return EXIT_FAILURE;
  else
    return EXIT_SUCCESS;
}

int incorrect_input(std::string error = "")
{
  std::cerr << "Incorrect input file.\n"
	    << "  " << error << "\n";
  std::cerr << "So far: " << vertices.size() << " vertices, "
	    << faces.size() << " faces\n";
  return EXIT_FAILURE;
}

int main(int argc, char** argv)
{
  if(argc != 3)
    return usage(argv[0]);
  mark_tag vertex_index(1), vertex_x(2), vertex_y(3), vertex_z(4);

  sregex tface_re = bos >> *space >> "TFACE" >> *space >> eos;

  sregex vertex_re = bos >> *space >> "VRTX" 
			 >> +space >> (vertex_index=+_d)
			 >> +space >> (vertex_x=+(digit|'-'|'+'|'.'))
			 >> +space >> (vertex_y=+(digit|'-'|'+'|'.'))
			 >> +space >> (vertex_z=+(digit|'-'|'+'|'.'))
			 >> eos;

  sregex triangle_re = bos >> *space >> "TRGL"
			   >> +space >> (s1=+digit)
			   >> +space >> (s2=+digit)
			   >> +space >> (s3=+digit)
			   >> eos;
  sregex end_re = bos >> *space >> "END" >> *space >> eos;

  std::ifstream input(argv[1]);
  std::ofstream output(argv[2]);

  if(!input) {
    std::cerr << "Cannot read \"" << argv[1] << "\"!\n";
    return EXIT_FAILURE;
  }
  if(!output) {
    std::cerr << "Cannot write to \"" << argv[2] << "\"!\n";
    return EXIT_FAILURE;
  }

  std::string line;
  std::getline(input, line);
  smatch results;
  while(input && ! regex_match(line, tface_re)) // search line "TFACE"
  {
    std::getline(input, line);
  }
  std::getline(input, line);
  while(input && regex_match(line, results, vertex_re)) {
    vertices.push_back(boost::make_tuple(results[vertex_x],
					 results[vertex_y],
					 results[vertex_z]));
    std::getline(input, line);
  }
  while(input && regex_match(line, results, triangle_re)) {
    std::stringstream s;
    int i, j, k;
    s << results[1] << " " << results[2] << " " << results[3];
    s >> i >> j >> k;
    faces.push_back(boost::make_tuple(i, j, k));
    std::getline(input, line);
  }
  if(!input || !regex_match(line, end_re))
    return incorrect_input("premature end of file!");

  output << "OFF " << vertices.size() << " " << faces.size() << " " << "0\n";
  for(Vertices::const_iterator vit = vertices.begin(), vend = vertices.end();
      vit != vend; ++vit)
    output << boost::get<0>(*vit) << " "
	   << boost::get<1>(*vit) << " "
	   << boost::get<2>(*vit) << "\n";
  for(Faces::const_iterator fit = faces.begin(), fend = faces.end();
      fit != fend; ++fit)
    output << "3 " << boost::get<0>(*fit)
	   << " " << boost::get<1>(*fit)
	   << " " << boost::get<2>(*fit) << "\n";
  if(output)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
};
