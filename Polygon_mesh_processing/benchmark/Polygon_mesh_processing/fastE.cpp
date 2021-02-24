#define TRACE
#include <fastenvelope/FastEnvelope.h>
#include <fastenvelope/Types.hpp>
#include <fstream>
#include <CGAL/Timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

int main(int argc, char* argv[]) {
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef typename K::Point_3 Point_3;

  std::vector<fastEnvelope::Vector3> env_vertices;
  std::vector<fastEnvelope::Vector3i> env_faces;

  std::ifstream in(argv[1]);

  double eps = std::stod(std::string(argv[2]));

  int ii, ij, ik;
  if(argc == 6){
    ii = std::stoi(std::string(argv[3]));
    ij = std::stoi(std::string(argv[4]));
    ik = std::stoi(std::string(argv[5]));
  }

  std::string off;
  int V, F, E;
  in >> off >> V >> F >> E;
  env_vertices.reserve(V);
  env_faces.reserve(F);

  double x,y,z;
  for(int i =0; i < V; i++){
    in >> x >> y >> z;
    fastEnvelope::Vector3 p = { x, y , z };
    env_vertices.push_back(p);
  }
  int three, vi, vj, vk;
  for(int i =0; i < F; ++i){
    in >> three >> vi >> vj >> vk;
    fastEnvelope::Vector3i f = {vi, vj, vk};
    env_faces.push_back(f);
  }

  CGAL::Timer t;
  t.start();
  fastEnvelope::FastEnvelope envelope(env_vertices, env_faces, eps);
  std::cout << t.time() << " sec." << std::endl;
  t.reset();

  if(argc == 6){
    fastEnvelope::Vector3 v0 = env_vertices[ii];
    fastEnvelope::Vector3 v1 = env_vertices[ij];
    fastEnvelope::Vector3 v2 = env_vertices[ik];
    std::array<fastEnvelope::Vector3, 3> tria = {v0, v1, v2};
    bool bbb = envelope.is_outside(tria);
    std::cout << t.time() << " sec." << std::endl;
  if(bbb){
      std::cout <<  "outside the envelope" << std::endl;
    }else{
      std::cout <<  "inside the envelope" << std::endl;
    }
    return 0;
  }

  int inside_count = 0;
  int outside_count = 0;

  std::ofstream inside("insideE.txt");
  std::ofstream outside("outsideE.txt");
  for(int i = 0; i <  env_vertices.size()  ; i+=10){
      for(int j = i+1; j < env_vertices.size(); j+=10){
        for(int k = j+1; k < env_vertices.size(); k+=10){
          if( ( i != j) && (i != k) && (j != k)){
            Point_3 p(env_vertices[i][0],env_vertices[i][1], env_vertices[i][2]);
            Point_3 q(env_vertices[j][0],env_vertices[j][1], env_vertices[j][2]);
            Point_3 r(env_vertices[k][0],env_vertices[k][1], env_vertices[k][2]);
            if(! CGAL::collinear(p,q,r)){
            std::array<fastEnvelope::Vector3, 3> f = { env_vertices[i],  env_vertices[j], env_vertices[k] };
              if(! envelope.is_outside(f)){
                inside_count++;
                //inside << i << " " << j << " "<< k <<std::endl;
              } else{
                outside_count++;
                //outside << i << " " << j << " "<< k <<std::endl;
              }
            }
          }
        }
      }

  }
  std::cout << inside_count << " " << outside_count << std::endl;
  std::cout << t.time() << " sec." << std::endl;


  return 0;
}
