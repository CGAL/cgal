#define CGAL_RECONSTRUCTION_FROM_PARALLEL_SLICES_3_DEBUG
//#define CGAL_ALLOW_NON_MANIFOLD_INPUT

#include <CGAL/Reconstruction_from_parallel_slices_3.h>
#include <CGAL/Reconstruction_from_parallel_slices_3/contour_providers.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <iostream>
#include <fstream>
#include <CGAL/self_intersect.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

//to test the reader from polygon range provided to SWIG
void file_polygon_list(char* fname,std::list<CGAL::Polygon_2<Kernel> >& polygons,std::list<double>& zcoords,const int constant_coordinate){
  typedef Kernel::Point_2 Point_2;
  std::ifstream input(fname);
  int nbpt;
  double coords[3];
  CGAL_assertion_code( double first_coords[3] );
  
  while (input >> nbpt && input){
    input >> coords[0] >> coords[1] >> coords[2];
    CGAL_assertion_code( first_coords[0]=coords[0] );
    CGAL_assertion_code( first_coords[1]=coords[1] );
    CGAL_assertion_code( first_coords[2]=coords[2] );

    polygons.push_back(CGAL::Polygon_2<Kernel>());
    zcoords.push_back(coords[constant_coordinate]);
    CGAL::Polygon_2<Kernel>& polygon=polygons.back();

    for (int i=1;i<nbpt;++i){
      input >> coords[0] >> coords[1] >> coords[2];
      polygon.push_back(Point_2(coords[(constant_coordinate+1)%3],coords[(constant_coordinate+2)%3]));
    }
    CGAL_assertion(first_coords[0]==coords[0]);
    CGAL_assertion(first_coords[1]==coords[1]);
    CGAL_assertion(first_coords[2]==coords[2]);    
  }
}

void file_vector_vector_Point_3(char* fname,std::vector< std::vector<Kernel::Point_3> >&  contours){
  typedef Kernel::Point_3 Point_3;
  std::ifstream input(fname);
  int nbpt;
  double x,y,z;
  
  while (input >> nbpt && input){
    contours.push_back(std::vector<Point_3>());
    contours.back().reserve(nbpt);
    
    input >> x >> y >> z;
    contours.back().push_back(Point_3(x,y,z));
    
    for (int i=1;i<nbpt;++i){
      input >> x >> y >> z;
      contours.back().push_back(Point_3(x,y,z));
    }
  }
}

int main(int argc, char* argv[])
{
  if(argc<2){
    std::cerr << "Usage: " << argv[0] << "all_in_one.cgal [012]" << std::endl;
    std::cerr << "The second optional argument indicates the axis perpendicular to the contour planes." <<std::endl;
    return 0;
  }
  
  unsigned int coord=2;
  if (argc==3)
    coord=atoi(argv[2]);
  
  #ifdef CGAL_ALLOW_NON_MANIFOLD_INPUT
  typedef CGAL::Slice_writer_into_file<Kernel::Point_3> Slice_writer;
  Slice_writer output("graph.off");
  #else
  typedef CGAL::Incremental_slice_writer_into_polyhedron<Polyhedron,Kernel> Slice_writer;
  Polyhedron polyhedron;
  Slice_writer output(polyhedron);
  #endif

  CGAL::Reconstruction_from_parallel_slices_3<Slice_writer> reconstruction;
  
  switch(coord){
    case 0:
    {
      CGAL::All_polygons_in_one_file_axis_aligned_planes<Kernel::Point_2,0> reader(argv[1]);
      
      //using vector of vector of Point_3
      //std::vector< std::vector<Kernel::Point_3> >  contours;
      //file_vector_vector_Point_3(argv[1],contours);
      //CGAL::Polygon_as_vector_of_Point_3_in_axis_aligned_planes<Kernel> reader(contours,0);
      
      //using vector of Polygon
      //std::list<CGAL::Polygon_2<Kernel> > polygons;
      //std::list<double> zcoords;
      //file_polygon_list(argv[1],polygons,zcoords,0);
      //CGAL::Polygon_range_in_axis_aligned_planes<std::list<CGAL::Polygon_2<Kernel> >::iterator,std::list<double>::iterator> 
      //  reader(polygons.begin(),polygons.end(),zcoords.begin(),zcoords.end());
      
      reconstruction.run(reader,output,0);
      break;
    }
    case 1:
    {
      CGAL::All_polygons_in_one_file_axis_aligned_planes<Kernel::Point_2,1> reader(argv[1]);
      reconstruction.run(reader,output,1);
      break;      
    }
    case 2:
    {
      CGAL::All_polygons_in_one_file_axis_aligned_planes<Kernel::Point_2,2> reader(argv[1]);
      reconstruction.run(reader,output,2);
      break;      
    }
    default:
      std::cerr << "Error the specified dimension is incorrect! This should be 0, 1 or 2"<< std::endl; 
      exit(EXIT_FAILURE);
  }
  
  #ifndef CGAL_ALLOW_NON_MANIFOLD_INPUT
  std::ofstream fout("graph.off");
  if ( !polyhedron.is_closed() ){
    std::cout << "The polyhedron is not closed" << std::endl;
    polyhedron.normalize_border();
    for (Polyhedron::Edge_iterator  eit=polyhedron.border_edges_begin(),eit_end=polyhedron.edges_end();eit!=eit_end;++eit)
      std::cout << eit->vertex()->point() <<  " " << eit->opposite()->vertex()->point() << std::endl;
  }
  fout << polyhedron;
  fout.close();

  //check that there is not degenerate triangles
  int degen_triangle=0;
  for (Polyhedron::Facet_iterator it=polyhedron.facets_begin();it!=polyhedron.facets_end();++it)
  {
    std::set<Kernel::Point_3> points;
    points.insert( it->halfedge()->vertex()->point() );
    points.insert( it->halfedge()->next()->vertex()->point() );
    points.insert( it->halfedge()->opposite()->vertex()->point() );
    if (points.size()!=3) ++degen_triangle;
  }
  if (degen_triangle!=0) std::cout << degen_triangle << " degenerate triangles" << std::endl;
  
  //test for self-intersections
  std::size_t k=0;
  CGAL::Counting_output_iterator counter(&k);
  CGAL::self_intersect<Polyhedron,Kernel>(polyhedron,counter);
  if (k>0) std::cout << "The polyhedron is self-intersecting" << std::endl;
  
  //check all vertices are used
  std::set<Polyhedron::Vertex_handle> used_vertices;
  for (Polyhedron::Halfedge_iterator it=polyhedron.halfedges_begin();it!=polyhedron.halfedges_end();++it)
    used_vertices.insert(it->vertex());
  int unused_vert=0;
  for (Polyhedron::Vertex_iterator it=polyhedron.vertices_begin();it!=polyhedron.vertices_end();++it)
    if (used_vertices.find(it)==used_vertices.end()){
      ++unused_vert;
      std::cout << it->point() << std::endl;
    }
  if (unused_vert!=0) std::cout << unused_vert << " vertices not used" << std::endl;
    
  #endif
  
  return 0;
}
