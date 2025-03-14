#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/lloyd_optimize_mesh_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_mesh_vertex_base_2<K>                Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K>                  Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>        Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds>  CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>            Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria>              Mesher;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;


using namespace std;

std::string result;

void CGAL_Remesh(double* vert_xyz_array, size_t vert_count, double criteria_a, double criteria_b, int iteration_number, double*& newVertices, int*& vCount, int*& newFaces, int*& fCount, int fileNum)
{
  CDT cdt;



  vector<Vertex_handle> cdt_Vh_Boundary;

  for (int i = 0; i < vert_count; ++i)
  {
    Vertex_handle vh = cdt.insert(Point(vert_xyz_array[3 * i + 0], vert_xyz_array[3 * i + 1]));
    cdt_Vh_Boundary.push_back(vh);
  }

  // insert Constrain
  for (int i = 0; i < cdt_Vh_Boundary.size() - 1; ++i)
  {
    cdt.insert_constraint(cdt_Vh_Boundary[i], cdt_Vh_Boundary[i + 1]);
  }
  cdt.insert_constraint(cdt_Vh_Boundary[cdt_Vh_Boundary.size() - 1], cdt_Vh_Boundary[0]);

  // refine and optimize mesh

  Mesher mesher(cdt);
  mesher.set_criteria(Criteria(criteria_a, criteria_b));
  mesher.refine_mesh();
  CGAL::lloyd_optimize_mesh_2(cdt, CGAL::parameters::max_iteration_number = iteration_number);


  // make index pair
  vector <CDT::Vertex_handle>  visitedVertices;   // collect visited vertices

  map <CDT::Vertex_handle, int> indexList;      // create a map to note the index

  int i = 0;
  for (CDT::Vertex_iterator v_it = cdt.vertices_begin(); v_it != cdt.vertices_end(); ++v_it)
  {

    CDT::Vertex_handle vh = v_it->handle();
    indexList[vh] = i;
    visitedVertices.push_back(vh);
    i++;
  }

  // Convert data into double array
  int vNum = cdt.number_of_vertices();

  newVertices = new double[vNum * 3];

  i = 0;
  for (CDT::Vertex_iterator vi = cdt.vertices_begin(); vi != cdt.vertices_end(); ++vi)
  {
    newVertices[i] = vi->point()[0];
    i += 1;
    newVertices[i] = vi->point()[1];
    i += 1;
    newVertices[i] = 0;
    i += 1;
  }


  int vertexCount = vNum;
  vCount = &vertexCount;

  int num_face_in_domain = 0;

  for (CDT::Face_iterator f_it = cdt.faces_begin(); f_it != cdt.faces_end(); ++f_it)
  {
    CDT::Face_handle face = f_it;
    if (face->is_in_domain())
    {
      num_face_in_domain++;
    }
  }


  newFaces = new int[num_face_in_domain * 3];


  i = 0;
  for (CDT::Face_iterator f_it = cdt.faces_begin(); f_it != cdt.faces_end(); ++f_it)
  {
    CDT::Face_handle face = f_it;

    if (face->is_in_domain())
    {
      newFaces[i] = int(indexList.find(face->vertex(0)->handle())->second);
      i += 1;
      newFaces[i] = int(indexList.find(face->vertex(1)->handle())->second);
      i += 1;
      newFaces[i] = int(indexList.find(face->vertex(2)->handle())->second);
      i += 1;
    }
  }
  int faceCount = num_face_in_domain;
  fCount = &faceCount;


  // print
  std::stringstream outputFile;

    if (fCount && newFaces) {
      outputFile << "\nRemeshed faces (" << *fCount << "):\n";
      for (int i = 0; i < (*fCount) * 3; i += 3) {
        outputFile << newFaces[i] << " " << newFaces[i + 1] << " " << newFaces[i + 2] << "\n";
      }
    }

    if(fileNum == 1)
    {
      result =  outputFile.str();
    }
    else
    {
      assert(result == outputFile.str());
    }

}

int main()
{
  //
  double vert_xyz_array[] = {
    3375.4981, 1935.35224056, 0.0,
    3350.77259333, 2066.38188194, 0.0,
    3210.83712383, 2054.53004190, 0.0,
    3068.88, 2060.98161842, 0.0,
    3034.24939361, 2066.55369658, 0.0,
    3025.6776, 2008.40156297, 0.0,
    3013.23241519, 1927.9864, 0.0,
    3033.36312437, 1924.87291062, 0.0,
    3078.68871988, 1917.86131994, 0.0,
    3124.01437021, 1910.85008365, 0.0,
    3167.66255113, 1908.86629111, 0.0,
    3169.83178711, 1908.76770020, 0.0,
    3215.64920401, 1906.68531674, 0.0,
    3260.96808047, 1913.74020504, 0.0,
    3306.03738982, 1922.24487242, 0.0,
    3351.10669914, 1930.74953992, 0.0
  };


  size_t vert_count = 16;
  double criteria_a = 0.125;
  double criteria_b = 36.691771392;
  int iteration_number = 20;

  double* newVertices = nullptr;
  int* vCount = nullptr;
  int* newFaces = nullptr;
  int* fCount = nullptr;

  for (int i = 1; i <= 2; ++i)
  {
    CGAL_Remesh(vert_xyz_array, vert_count, criteria_a, criteria_b, iteration_number, newVertices, vCount, newFaces, fCount, i);
  }


  std::cout << "\nDone";


  return 0;
}
