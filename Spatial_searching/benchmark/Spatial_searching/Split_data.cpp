// $URL$
// $Id$

#include <CGAL/Cartesian_d.h>
#include <CGAL/Random.h>

#include <iostream>
#include <fstream>
#include <list>

typedef CGAL::Cartesian_d<double> K;
typedef K::Point_d Point_d;


int main() {

  CGAL::Random Rnd;
  int N;

  char filename[80];
  std::cout << "Enter input file name containing data points: \n" ;
  std::cin >> filename;

  char filename_out_query_points[80];
  std::cout << "Enter output file name containing query points: \n" ;
  std::cin >> filename_out_query_points;

  char filename_out_data_points[80];
  std::cout << "Enter output file name containing data points: \n" ;
  std::cin >> filename_out_data_points;

  std::ifstream in;
  std::ofstream out_query, out_data;

  int data_point_number, query_point_number;

  typedef std::list<Point_d> point_list;
  point_list all_points;

  in.open(filename);
  out_query.open(filename_out_query_points);
  out_data.open(filename_out_data_points);
  in >> N;
  in >> query_point_number;
  in >> data_point_number;
  int point_number=query_point_number+data_point_number;

  out_query << N; out_query << " "; out_query << query_point_number; out_query << std::endl;
  out_data << N; out_data << " "; out_data << data_point_number; out_data << std::endl;

  std::cout << "dimension = " << N << std::endl;
  std::cout << "data point number = " << data_point_number << std::endl;
  std::cout << "query point number = " << query_point_number << std::endl;

  std::vector<int> query_point(point_number);


  for (int ii = 0; ii < point_number; ii++) query_point[ii]=0;

  // random selection of data points
  for (int jj=0; jj < query_point_number; ) {
       int random_number=Rnd.get_int(0,point_number);
       if (query_point[random_number]==0) {
                query_point[random_number]=1;
                  jj++;
       }
 }


 for (int i = 0; i < point_number; i++) {
        std::vector<double> p(N);
        for (int j = 0; j < N; j++) {
          in >> p[j];
        }
        Point_d Pnt(N,p.begin(),p.end());
        all_points.push_back(Pnt);
  };


  // schrijf all points naar query points of data points
  int counter=0;
  for(point_list::iterator it = all_points.begin(); it != all_points.end(); ++it) {
     //
     if (query_point[counter]==1) {
             for (int j = 0; j < N; j++) {
                  out_query << (*it)[j]; out_query << " ";
             }
             out_query << std::endl;
     } else {
             for (int j = 0; j < N; j++) {
                  out_data << (*it)[j]; out_data << " ";
             }
             out_data << std::endl;
     }
     counter++;
  }


  out_query.close();
  out_data.close();
  return 0;
}
