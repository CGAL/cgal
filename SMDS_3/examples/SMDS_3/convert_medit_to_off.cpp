#include <ios>
#include <string>
#include <iostream>
#include <fstream>

int main(int argc, char** argv)
{
  if(argc != 3)
  {
    std::cerr << "Usage: " << argv[0] << " input.medit output.off" << std::endl;
    return 1;
  }
  std::ifstream input(argv[1]);
  if(!input)
  {
    std::cerr << "Cannot open file " << argv[1] << std::endl;
    return 1;
  }
  std::ofstream output(argv[2]);
  if(!output)
  {
    std::cerr << "Cannot open file " << argv[2] << std::endl;
    return 1;
  }
  bool verbose = true;

  auto& is = input;

  int dim;
  int nv, nf;
  std::string word;
  std::ifstream::pos_type vertices_pos;

  is >> word >> dim; // MeshVersionFormatted 1
  is >> word >> dim; // Dimension 3

  if(verbose)
  {
    std::cout << "Reading .mesh file..." << std::endl;
  }

  std::string line;
  while(std::getline(is, line) && line != "End")
  {
    // remove trailing whitespace, in particular a possible '\r' from Windows
    // end-of-line encoding
    if(!line.empty() && std::isspace(line.back())) {
      line.pop_back();
    }
    if (line.size() > 0 && line.at(0) == '#') {
      continue;
    }

    if(line.find("Vertices") != std::string::npos)
    {
      is >> nv;
      if(verbose)
        std::cerr << "Reading "<< nv << " vertices" << std::endl;

      vertices_pos = is.tellg();
      continue;
    }

    if(line.find("Triangles") != std::string::npos)
    {
      is >> nf;
      if(verbose)
        std::cerr << "Reading "<< nf << " triangles" << std::endl;
      break;
    }
  }
  if(!input) {
    std::cerr << "Issue after reading the number of vertices and triangles\n";
    return 1;
  }
  input.seekg(vertices_pos);
  if(!input) {
    std::cerr << "Issue after seekg\n";
    return 1;
  }
  output << "OFF\n";
  output << nv << " " << nf << " 0\n";
  for(int i = 0; i < nv; ++i) {
    std::string x, y, z, index;
    input >> x >> y >> z >> index;
    output << x << ' ' << y << ' ' << z << '\n';
    if(i < 10) std::cerr << x << ' ' << y << ' ' << z << '\n';
  }
  while(std::getline(is, line) && line != "End")
  {
    if(line.find("Triangles") != std::string::npos)
    {
      is >> nf;
      if(verbose)
        std::cerr << "Reading "<< nf << " triangles" << std::endl;
      for(int i=0; i<nf; ++i)
      {
        int n[3];
        int surface_patch_id;
        if(!(is >> n[0] >> n[1] >> n[2] >> surface_patch_id))
        {
          if(verbose)
            std::cerr << "Issue while reading triangles" << std::endl;
          return 1;
        }
        output << "3 " << n[0] - 1 << ' ' << n[1] - 1 << ' ' << n[2] - 1 << '\n';
        if(i < 10) std::cerr << n[0] - 1 << ' ' << n[1] - 1 << ' ' << n[2] - 1 << '\n';
      }
    }
  }
}
