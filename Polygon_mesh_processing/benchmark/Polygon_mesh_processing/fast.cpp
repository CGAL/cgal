
// #define CGAL_PROFILE


//#define TRACE

// fast.cpp and fastE.cpp produce the same trace output to see where the executables take different paths
// As fastE operates on reordered faces it is important to use as input an off file
// where the faces are ordered in such a way that they do not get reordered.
// Therefore fastE.cpp writess the indices into a file named "sorted.off" which can be used
// to replace the faces in the original input.


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedral_envelope.h>
#include <CGAL/Timer.h>

#include <fstream>


// `fast` takes an off file and an offset as arguments
//  If called additionally with 3 more vertex indices it performs the envelope test with the triangle
//  Otherwise it tests for all vertex triples forming a non-degenerate trianges
//  and writes the triple in the file inside.txt or outside.txt

int main(int argc, char* argv[])
{
  typedef std::array<int, 3> Vector3i;
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef Kernel::Point_3 Point_3;
  typedef CGAL::Polyhedral_envelope<Kernel> Envelope;

  std::vector<Point_3> env_vertices;
  std::vector<Vector3i> env_faces;

  std::ifstream in(argv[1]);

  double eps = std::stod(std::string(argv[2]));

  int ii, ij, ik;

  int query_dim = -1;

  if(argc >= 4){
    query_dim++;
    ii = std::stoi(std::string(argv[3]));
  }

  if(argc >= 5){
    query_dim++;
    ij = std::stoi(std::string(argv[4]));
  }

  if(argc == 6){
    query_dim++;
    ik = std::stoi(std::string(argv[5]));
  }
  std::string off;
  int V, F, E;
  in >> off >> V >> F >> E;
  env_vertices.reserve(V);
  env_faces.reserve(F);

  Kernel::Point_3 p;
  for(int i =0; i < V; i++){
    in >> p;
    env_vertices.push_back(p);
  }
  int three, vi, vj, vk;
  for(int i =0; i < F; ++i){
    in >> three >> vi >> vj >> vk;
    Vector3i f = { vi, vj, vk };
    env_faces.push_back(f);
  }

  CGAL::Timer t;
  t.start();

  Envelope envelope(env_vertices, env_faces, eps);

  std::cout << t.time() << " sec." << std::endl;
  t.reset();


  if(query_dim != -1){

    bool bbb;

    if(query_dim == 0){
    Point_3 v0 = env_vertices[ii];
    bbb = envelope(v0);
    } else if (query_dim == 1){
      Point_3 v0 = env_vertices[ii];
      Point_3 v1 = env_vertices[ij];
      bbb = envelope(v0, v1);
    } else {
      Point_3 v0 = env_vertices[ii];
      Point_3 v1 = env_vertices[ij];
      Point_3 v2 = env_vertices[ik];
      /*
      {
        std::ofstream query("query.off");
        query << "OFF\n" << "3 1 0\n" << v0 << std::endl << v1 << std::endl << v2 << std::endl << "3 0 1 2" << std::endl;
      }
      */
      bbb = envelope(v0, v1, v2);

      std::cout << t.time() << " sec." << std::endl;
    }

    if(bbb){
      std::cout <<  "inside the envelope" << std::endl;
    }else{
      std::cout <<  "outside the envelope" << std::endl;
    }
    return 0;
  }


  int inside_count = 0;
  int outside_count = 0;

  std::ofstream inside("inside.txt");
  std::ofstream outside("outside.txt");

  for(int i = 0; i < env_vertices.size() ; i+=10){
      for(int j = i+1; j < env_vertices.size(); j+= 10){
        for(int k = j+1; k < env_vertices.size(); k+=10){
          if( ( i != j) && (i != k) && (j != k)){
            if(! CGAL::collinear(env_vertices[i], env_vertices[j],env_vertices[k])){
              if(envelope(env_vertices[i],  env_vertices[j], env_vertices[k])){
                inside_count++;
                // inside << i << " " << j << " "<< k <<std::endl;
              } else{
                outside_count++;
                // outside << i << " " << j << " "<< k <<std::endl;
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
