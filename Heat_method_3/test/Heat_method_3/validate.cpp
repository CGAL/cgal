#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Heat_method_3/Heat_method_3.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <string>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef Kernel::Point_2                                      Point_2;
typedef CGAL::Surface_mesh<Point>                            Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;
typedef Surface_mesh::Property_map<vertex_descriptor,double> Vertex_distance_map;

typedef CGAL::Heat_method_3::Heat_method_3<Surface_mesh> Heat_method;
typedef CGAL::Heat_method_3::Heat_method_3<Surface_mesh, CGAL::Tag_true> Heat_method_idt;



int validate(char* fname)
{
  std::string s(fname);
  std::string base = s.substr(0,s.length()-4);
  Surface_mesh sm;

  std::ifstream in(fname);
  in >> sm;
  if(!in || num_vertices(sm) == 0) {
    std::cerr << "Problem loading the input data" << std::endl;
    return 1;
  }
  std::vector<int> sources;

  std::ifstream src((base+".src").c_str());
  int index;
  while(src >> index){
    std::cout << "source: " << index << std::endl;
    sources.push_back(index);
  }
  std::cout << std::endl;

  std::ifstream dists((base+".dists").c_str());
  std::vector<double> groundtruth;
  double d;
  while(dists >> d){
    groundtruth.push_back(d);
  }

  std::ifstream ref((base+".ref").c_str());
  std::vector<double> reference;
  double d2;
  while(ref >> d2){
    reference.push_back(d2);
  }


  Vertex_distance_map vdm = sm.add_property_map<vertex_descriptor,double>("v:dist",0).first;

  Heat_method hm(sm);
  for(int index : sources){
    hm.add_source(vertex_descriptor(index));
  }
  hm.fill_distance_map(vdm);

  Vertex_distance_map vdm_idt = sm.add_property_map<vertex_descriptor,double>("v:idt",0).first;

  Heat_method_idt hm_idt(sm);
  for(std::size_t i=0; i < sources.size(); i++){
    hm_idt.add_source(vertex_descriptor(index));
  }
  hm_idt.fill_distance_map(vdm_idt);

  std::cout.precision(17);
  int i = 0;
  BOOST_FOREACH(vertex_descriptor vd, vertices(sm)){
    std::cout << vd << " without idt: " << get(vdm, vd) << "  with idt: "  << get(vdm_idt, vd) << "  " << groundtruth[i] << " and exact polyhedral: " << reference[i] << std::endl;
    i++;
  }
  std::cout << "done" << std::endl;

  return 0;
}

int main(int argc, char*argv[])
{
  int res = 0;
  for(int i=1; i < argc; i++){
    std::cout << "validate("<< argv[i] << ")"<< std::endl;
    validate(argv[i]);
    res++;
  }
  return res;
}
