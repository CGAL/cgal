#include <fastenvelope/FastEnvelope.h>
#include <fastenvelope/Types.hpp>
#include <fstream>
#include <CGAL/Timer.h>

 int main(int argc, char* argv[])
 {
 std::vector<fastEnvelope::Vector3> env_vertices;
  std::vector<fastEnvelope::Vector3i> env_faces;


  std::ifstream in(argv[1]);

  double eps = std::stod(std::string(argv[2]));

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

  fastEnvelope::Vector3 v0 = env_vertices[0];
  fastEnvelope::Vector3 v44 = env_vertices[44];
  bool b0_44 = envelope.is_outside(v0,v44);

  std::ofstream out("insideF.polylines.txt");
  int count = 0;
  for(int i = 0; i < env_vertices.size(); i++){
    fastEnvelope::Vector3 vi = env_vertices[i];
    if(! envelope.is_outside(v0,vi)){
      count++;
      //out << "2 " << v0 << " " << vi << std::endl;
    }
  }
  std::cout << count << " inside in " << t.time() << " sec." << std::endl;

  return 0;
}
