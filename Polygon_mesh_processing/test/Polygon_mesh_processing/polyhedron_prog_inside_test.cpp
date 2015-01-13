#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_inside_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Timer.h>
#include <boost/lexical_cast.hpp>

//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron;

void random_points(char* fname, CGAL::Bbox_3 bbox, int n)
{
  std::cerr << "Write " << n << " points to file " << fname << std::endl;
  CGAL::Random rg;

  double grid_dx = bbox.xmax() - bbox.xmin();
  double grid_dy = bbox.ymax() - bbox.ymin();
  double grid_dz = bbox.zmax() - bbox.zmin();
  
  std::ofstream out(fname);
  out.precision(16);
  for(int i=0; i < n; i++){
    out << bbox.xmin() + rg.get_double()* grid_dx << " " 
        << bbox.ymin() + rg.get_double()* grid_dy << " "
        << bbox.zmin() + rg.get_double()* grid_dz << std::endl;
  }
  out.close();
}

int main(int argc, char* argv[])
{
  std::cerr.precision(17);

  std::ifstream points_file(argv[2]);

  std::vector<Point> points;
  std::copy(std::istream_iterator<Point>(points_file),
            std::istream_iterator<Point>(),
            std::back_inserter(points));


  int gridsize = 10;
  if(argc>3){
    gridsize = boost::lexical_cast<int>(argv[3]);
  }
  std::cerr << "gridsize = " << gridsize << std::endl;
  
  std::size_t nb_points=points.size();
  
  std::vector<bool> ray_res(nb_points);
  std::vector<bool> grid_res(nb_points);
  
  //using ray
  {
    Polyhedron polyhedron;
    std::ifstream polyhedron_file(argv[1]);
    polyhedron_file >> polyhedron;
    std::cerr << "|V| = " << polyhedron.size_of_vertices() << std::endl;
    
    CGAL::Timer timer;
    timer.start();
    CGAL::Point_inside_polyhedron_3<Polyhedron,K> inside_with_ray(polyhedron);
    timer.stop();
    std::cerr <<"Using ray"<< std::endl;
    std::cerr << "  Preprocessing took " << timer.time() << " sec." << std::endl;
    timer.reset();  
    int n_inside = 0;
    
    timer.start();
    for(std::size_t k=0; k<nb_points; ++k){
      ray_res[k] = (inside_with_ray(points[k]) == CGAL::ON_BOUNDED_SIDE);
      if(ray_res[k]){
        ++n_inside;
      }
    }
    timer.stop();
    std::cerr << "  " << n_inside << " points inside " << std::endl;
    std::cerr << "  " << points.size() - n_inside << " points outside "  << std::endl;
    std::cerr << " Queries took " << timer.time() << " sec." << std::endl; 

  }
  
  //using grid
  {
    Polyhedron polyhedron;
    std::ifstream polyhedron_file(argv[1]);
    polyhedron_file >> polyhedron;
    std::cerr << "|V| = " << polyhedron.size_of_vertices() << std::endl;
    
    CGAL::Timer timer;
    timer.start();
    CGAL::Point_inside_polyhedron_3<Polyhedron,K> inside_with_grid(polyhedron, gridsize);
    timer.stop();
    std::cerr <<"Using grid"<< std::endl;
    std::cerr << "  Preprocessing took " << timer.time() << " sec." << std::endl;
    timer.reset();

    if(argc>5){
      random_points(argv[4],inside_with_grid.bbox() ,boost::lexical_cast<int>(argv[5]) );
    }
    
    int n_inside = 0;
    timer.start();
    for(int k=0;k<nb_points;++k){
      grid_res[k]=inside_with_grid(points[k]);
      if(grid_res[k]){
        ++n_inside;
      }
    }
    timer.stop();
    std::cerr << "  " << n_inside << " points inside " << std::endl;
    std::cerr << "  " << points.size() - n_inside << " points outside "  << std::endl;
    std::cerr << "  Queries took " << timer.time() << " sec." << std::endl; 
  }

  for(int k=0;k<nb_points;++k){
    if(ray_res[k]!=grid_res[k]){
      std::cerr << "WARNING: Result is different for point " << k << std::endl;
    }
  }
  
  //using original code
  {
    Polyhedron polyhedron;
    std::ifstream polyhedron_file(argv[1]);
    polyhedron_file >> polyhedron;
    std::cerr << "|V| = " << polyhedron.size_of_vertices() << std::endl;

    
    std::cerr <<"Using ray (original code)"<< std::endl;
    CGAL::Point_inside_polyhedron_3<Polyhedron,K> inside_with_ray(polyhedron,0);
    
    CGAL::Timer timer;  
    int n_inside = 0;
    timer.start();
    for(int k=0;k<nb_points;++k)
      if(inside_with_ray(points[k],true)) ++n_inside;
    timer.stop();
    std::cerr << "  " << n_inside << " points inside " << std::endl;
    std::cerr << "  " << points.size() - n_inside << " points outside "  << std::endl;
    std::cerr << "  Queries took " << timer.time() << " sec." << std::endl; 
  }
  
  return 0;
}
