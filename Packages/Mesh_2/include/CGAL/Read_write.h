#include <iostream>
#include <fstream>

#include <CGAL/IO/File_header_extended_OFF.h> 
// to skip comments and EOF in the function read_poly



//   // IO

//   // write and read the constrained edges in the format:
//   //   number_of_edges
//   //   segment1
//   //   segment2
//   //   ...
//   void write(std::ostream &f) const;
//   void read(std::istream &f, bool dont_refine = false);

//   // write and read a mesh in the Triangle .poly format 
//   // (see http://www-2.cs.cmu.edu/~quake/triangle.poly.html)
//   void write_poly(std::ostream &f) const;
//   void read_poly(std::istream &f, bool dont_refine = false);

CGAL_BEGIN_NAMESPACE

// IO
//the function that writes a file
template <class Tr>
void
write(const Conform_triangulation_2<Tr>& mesh, std::ostream &f)
{
  typedef typename Conform_triangulation_2<Tr>::Finite_edges_iterator
    Finite_edges_iterator;

  f << mesh.number_of_constrained_edges() << std::endl;
  for(Finite_edges_iterator eit = mesh.finite_edges_begin();
      eit!=mesh.finite_edges_end();
      ++eit)
    if((*eit).first->is_constrained((*eit).second)) 
      {
	f << (*eit).first->vertex(mesh.cw((*eit).second))->point() << " "
	  << (*eit).first->vertex(mesh.ccw((*eit).second))->point() <<std::endl;
      }
}

//the function that reads a file
template <class Tr>
void
read(Conform_triangulation_2<Tr>& mesh, std::istream &f)
{
  typedef typename Conform_triangulation_2<Tr>::Point Point;
  int nedges = 0;
  mesh.clear();
  f>>nedges;
  for(int n = 0; n<nedges; n++) {
    Point p1, p2;
    f >> p1 >> p2;
    mesh.insert_constraint(p1, p2);
  }
}

//the function that write a Shewchuk Triangle .poly file
template <class Tr>
void
write_poly(const Conform_triangulation_2<Tr>& mesh, std::ostream &f)
{
  typedef Conform_triangulation_2<Tr> Triangulation;
  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Finite_vertices_iterator
    Finite_vertices_iterator;
  typedef typename Triangulation::Finite_edges_iterator
    Finite_edges_iterator;


  std::map<Vertex_handle, unsigned int> index;

  // write vertices
  f << "# Shewchuk Triangle .poly file, produced by the CGAL::Mesh_2 package"
    << std::endl
    << "# Neither attributes nor boundary markers are used." << std::endl
    << mesh.number_of_vertices() << " " << 2 << " " 
    << 0 << " " << 0 << std::endl;

  f << std::endl;

  unsigned int vertices_counter = 0;
  for(Finite_vertices_iterator vit = mesh.finite_vertices_begin();
      vit != mesh.finite_vertices_end();
      ++vit)
    {
      f << ++vertices_counter << " " << vit->point() << std::endl;
      index[vit] = vertices_counter;
    }

  f << std::endl;

  // write constrained edges

  f << mesh.number_of_constrained_edges() << " " << 0 << std::endl;
  unsigned int edges_counter = 0;
  for(Finite_edges_iterator eit = mesh.finite_edges_begin();
      eit != mesh.finite_edges_end();
      ++eit)
    if((*eit).first->is_constrained((*eit).second)) 
      f << ++edges_counter << " "
	<< index[(*eit).first->vertex(mesh.cw((*eit).second))] << " "
	<< index[(*eit).first->vertex(mesh.ccw((*eit).second))] 
	<< std::endl;

  f << std::endl;

//   // write seeds, assuming that the seeds unmarks faces
//   unsigned int seeds_counter = 0;
//   f << mesh.seeds.size() << std::endl;
//   for(typename Seeds::const_iterator sit = seeds.begin();
//       sit!=seeds.end(); ++sit)
//     f << ++seeds_counter << " " << *sit << std::endl;
}

//the function that reads a Shewchuk Triangle .poly file
template <class Tr>
void
read_poly(Conform_triangulation_2<Tr>& mesh, std::istream &f)
{
  typedef Conform_triangulation_2<Tr> Triangulation;
  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Point Point;

  mesh.clear();

  unsigned int number_of_points;
  skip_comment_OFF(f);
  f >> number_of_points;
  skip_until_EOL(f);
  skip_comment_OFF(f);
  
  // read vertices
  std::vector<Vertex_handle> vertices(number_of_points);
  for(unsigned int i = 0; i < number_of_points; ++i)
    {
      unsigned int j;
      Point p;
      f >> j >> p;
      skip_until_EOL(f); skip_comment_OFF(f);
      vertices[--j] = mesh.insert(p);
    }

  // read segments
  unsigned int number_of_segments;
  f >> number_of_segments;
  skip_until_EOL(f); skip_comment_OFF(f);
  for(unsigned int k = 0; k < number_of_segments; ++k)
    {
      unsigned int l, v1, v2;
      f >> l >> v1 >> v2;
      skip_until_EOL(f); skip_comment_OFF(f);
      mesh.insert_constraint(vertices[--v1], vertices[--v2]);
    }

  // read holes
  unsigned int number_of_holes;
  f >> number_of_holes;
  for(unsigned int m = 0; m < number_of_holes; ++m)
    {
      unsigned int n;
      Point p;
      f >> n >> p;
      skip_until_EOL(f); skip_comment_OFF(f);
      //seeds.push_back(p);
    }
}

CGAL_END_NAMESPACE

