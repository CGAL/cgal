
#include <assert.h>
#include <iostream.h>
#include <fstream.h>
#include <CGAL/basic.h>
#include <CGAL/predicates_on_pointsC2.h> // for the orientation

#include "points.h"

#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>


typedef Euclidean_2 Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Triangulation_2<Gt,Tds>  Triangulation;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds>  Delaunay_triangulation;

int main(int argc, char* argv[])
{

  Triangulation T;

  PVector V;

  ifstream data(argv[1]);
  
  if(data.bad()){
    cout << "Problem with file " << argv[1] << endl;
  }
  data >> V;

  cout << V.size() << endl;

  int count = 0;
  cout << "Start insertion" << endl;
  for(int i = 0; i<V.size();i++){
    
    T.insert(V[i]);

    if(++count == 100){
      cout << ".";
      count = 0;
    }
  }
  cout << endl << "done" << endl;
  T.is_valid();

  return 0;
}
