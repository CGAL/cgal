#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/IO/read_xyz_points.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Point_set_3<Kernel> Point_set;


struct  Point_push_pmap {
  Point_set& ps;
  Point_set::Index ind;

  Point_push_pmap(Point_set& ps, int ind=0)
    : ps(ps), ind(ind)
  {}

  inline friend void put(Point_push_pmap& pm, Point_set::Index& i, Point& p)
  {
    if(! pm.ps.surface_mesh().has_valid_index(pm.ind)){
      pm.ps.surface_mesh().add_vertex();
    }
    put(pm.ps.points(), pm.ind,p);
    i = pm.ind;
    ++pm.ind;
  }
};

struct  Normal_push_pmap {
  Point_set& ps;
   Point_set::Index ind;

  Normal_push_pmap(Point_set& ps, int ind=0)
    : ps(ps), ind(ind)
  {}

  inline friend void put(Normal_push_pmap& pm, Point_set::Index& i  , Vector& v)
  {
    if(! pm.ps.surface_mesh().has_valid_index(pm.ind)){
      pm.ps.surface_mesh().add_vertex();
    }
    put(pm.ps.normals(), pm.ind,v);
    i = pm.ind;
    ++pm.ind;
  }
};




int main (int argc, char** argv)
{
  Point_set point_set;

  if (point_set.has_normals())
    std::cerr << "Point set has normals" << std::endl;
  else
    std::cerr << "Point set doesn't have normals" << std::endl;

  point_set.add_normal_property();

  std::vector<Point_set::Index> indices;
  read_xyz_points_and_normals(std::ifstream("data.pwn"),
                              std::back_inserter(indices),
                              Point_push_pmap(point_set),
                              Normal_push_pmap(point_set),
                              Kernel());
  std::cerr << point_set.surface_mesh() << std::endl;

  for(int i =0; i < indices.size(); i++){
    std::cerr << indices[i] << std::endl;
    }
  Point_set::Index v0(0), v1(1);
  std::cerr << point_set.normal(v0) << "  " << point_set.normal(v1) << std::endl;
  point_set.surface_mesh().swap(v0,v1);
  std::cerr << point_set.normal(v0) << "  " << point_set.normal(v1) << std::endl;
  return 0;
    
  for (std::size_t i = 0; i < 10; ++ i)
    point_set.push_back (Point (rand () / (FT)RAND_MAX,
                                rand () / (FT)RAND_MAX,
                                rand () / (FT)RAND_MAX));


  for (std::size_t i = 0; i < point_set.size(); ++ i)
    std::cerr << "Point " << i << ": " << point_set[i] << std::endl;

  if (point_set.has_normals())
    std::cerr << "Point set has normals" << std::endl;
  else
    std::cerr << "Point set doesn't have normals" << std::endl;

  point_set.push_back (Point (rand() / (FT)RAND_MAX,
                              rand() / (FT)RAND_MAX,
                              rand() / (FT)RAND_MAX),
                       Vector (rand() / (FT)RAND_MAX,
                               rand() / (FT)RAND_MAX,
                               rand() / (FT)RAND_MAX));
  
  for (std::size_t i = 0; i < point_set.size(); ++ i)
    {
      point_set.normal(i) = Vector (rand () / (FT)RAND_MAX,
                                    rand () / (FT)RAND_MAX,
                                    rand () / (FT)RAND_MAX);
      std::cerr << "Normal " << i << ": " << point_set.normal(i) << std::endl;
    }

  point_set.remove_normal_property();
  
  if (point_set.has_normals())
    std::cerr << "Point set has normals" << std::endl;
  else
    std::cerr << "Point set doesn't have normals" << std::endl;


  Point_set::Point_pmap pm = point_set.points();
  Point_set::Vector_pmap nm = point_set.normals();
  
  
  std::cout << get(pm,Point_set::Index(0)) << std::endl;
  std::cout << get(nm,Point_set::Index(0)) << std::endl;
  
  
  return 0;
}
