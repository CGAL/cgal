#include <iostream>
#include <fstream>
#include <string>

void usage(char* name)
{
  std::cerr << "Usage:\n"
	    << name << " [FILE_BASENAME]\n"
	    << "    will convert FILE_BASENAME.off to FILE_BASENAME.mesh.\n";
  exit(1);
}

int main(int argc, char** argv)
{
  if(argc < 2)
    usage(argv[0]);

  std::string file_base = argv[1];
  
  std::ifstream off((file_base + ".off").c_str());
  std::ofstream mesh((file_base + ".mesh").c_str());
  
  std::string dummy_s;

  int number_of_points;
  int number_of_faces;
  int dummy_i;

  mesh << "MeshVersionFormatted 1\n"
       << "Dimension 3\n"
       << "Vertices\n";
  
  off >> dummy_s; // OFF
  off >> number_of_points >> number_of_faces >> dummy_i;
  mesh << number_of_points << "\n";

  while( number_of_points > 0 && ! off.eof() )
  {
    std::string x, y, z;
    off >> x >> y >> z;
    if(off)
    {
      mesh << x << " " << y << " " << z << " 1\n";
      --number_of_points;
    }
  }

  if( number_of_points != 0)
  {
    std::cerr << "Error: bad number of points in OFF file "
              << file_base << ".off\n";
    exit(2);
  }

  mesh << "Triangles" << "\n"
       << number_of_faces << "\n";

  while( ! off.eof() )
  {
    int i0, i1, i2, i3;
    off >> i0 >> i1 >> i2 >> i3;
    if(off)
    {
      mesh << (i1 + 1) << " "
           << (i2 + 1) << " "
           << (i3 + 1) << " 0\n";
      --number_of_faces;
    }
  }

  mesh << "Tetrahedra\n0\nEnd\n";

  mesh.close();
  if( number_of_faces != 0)
  {
    std::cerr << "Error: bad number of faces in OFF file "
              << file_base << ".off\n";
    exit(2);
  }
}
