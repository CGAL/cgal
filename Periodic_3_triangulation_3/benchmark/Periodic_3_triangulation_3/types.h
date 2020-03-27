#ifndef TYPES_H
#define TYPES_H

#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Random.h>
#include <CGAL/point_generators_3.h>

#include <CGAL/Timer.h>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

void print_usage(const std::string & bench_name) {
  std::cout<<bench_name
           <<" [insert | remove] [-i FILE | -c N M] "
           <<"[-f FILE] [-w FILE] [-r N] [-p N] [-s N] [--nospatial] "
           <<"[--noread] [--nowrite]" <<std::endl<<std::endl;
  std::cout<<"\t insert: \n\t\tbenchmark only the insert functionality"
      <<std::endl;
  std::cout<<"\t remove: \n\t\tbenchmark only the remove functionality"
      <<std::endl;
  std::cout<<"\t -i FILE: \n\t\tinserts test data points from FILE"
      <<std::endl;
  std::cout<<"\t -c N M: \n\t\tinserts N+M random points, N of which "
           <<"are created by in_cube an M of which created by on_sphere. "
           <<"Default: 5 10"<<std::endl;
  std::cout<<"\t -r N: \n\t\tremoves N vertices."<<std::endl;
  std::cout<<"\t -f FILE: \n\t\treads an initial triangulation from FILE, "
      <<"default: data/Triangulation.in"<<std::endl;
  std::cout<<"\t -w FILE: \n\t\twrites the computed triangulation to FILE, "
      <<"default: data/Triangulation.out"<<std::endl;
  std::cout<<"\t -p N: \n\t\tAdapt the floating point precision in the read "
      <<"and written files. Default: 18"<<std::endl;
  std::cout<<"\t -s N: \n\t\tSeed for the random point generation. Default: 7"
      <<std::endl;
  std::cout<<"\t --nospatial: \n\t\tDoes not do spatial sorting"<<std::endl;
  std::cout<<"\t --noread: \n\t\tDoes not read a triangulation from the disk "
      <<"but computes a completely new one"<<std::endl;
  std::cout<<"\t --nowrite: \n\t\tDoes not write the final triangulation to "
      <<"the disk"<<std::endl<<std::endl;
  exit(-1);
}

template <class Triangulation, class Stream, class InsertIterator>
Stream& read_point_set_from_skel(Stream &in, InsertIterator it) {
  std::string s;
  in >> s;
  CGAL_assertion( s == "SKEL" );
  int no_of_points, no_of_cells;
  in >> no_of_points >> no_of_cells ;
  typedef typename Triangulation::Geom_traits::FT Number_type;
  typedef typename Triangulation::Point Point;
  Number_type x,y,z;
  for (int i=0 ; i < no_of_points ; i++) {
    in>>x;
    in>>y;
    in>>z;
    *it++ = Point(x+0.5,y+0.5,z+0.5);
  }
  return in;
}

template <class Point, class InsertIterator>
void create_random_points(InsertIterator it, int c, int s, int seed = 7)
{
  CGAL::Random random(seed);
  typedef CGAL::Creator_uniform_3<double,Point>  Creator;

  CGAL::Random_points_in_cube_3<Point, Creator> in_cube(.5, random);
  for (int count=0; count < c; count++) {
    Point p = *in_cube; in_cube ++;
    p = Point(p.x()+.5,p.y()+.5,p.z()+.5);
    *it++ = p;
  }
  CGAL::Random_points_on_sphere_3<Point, Creator> on_sphere(.25, random);
  for (int count=0; count < s; count++) {
    Point p = *on_sphere; on_sphere ++;
    p = Point(p.x()+.5,p.y()+.5,p.z()+.5);
    *it++ = p;
  }
}

template <class Triangulation>
void test_remove(Triangulation &T, int n) {
  for (int i=0 ; i < n ; i++) {
    if (T.number_of_vertices() > 0)
      T.remove( ++T.all_vertices_begin() );
  }
}

template <class Triang>
int bench_triang(int argc, char* argv[], Triang T, bool periodic) {
  std::string bench_name,dfile, ifile, ofile;
  if (periodic) {
    bench_name = "Periodic_3_Delaunay_3_bench";
    ifile = "data/Periodic_3_triangulation.in";
    ofile = "data/Periodic_3_triangulation.out";
  } else {
    bench_name = "Delaunay_3_bench";
    ifile = "data/Triangulation.in";
    ofile = "data/Triangulation.out";
  }

  typedef typename Triang::Point Point;
  typedef CGAL::Timer Timer;

  bool do_write = true;
  bool do_read = true;
  bool do_insert = true;
  bool do_remove = true;
  bool do_spatial_sort = true;
  bool random_data = true;

  std::vector<Point> pts;
  unsigned int no_cube = 5;
  unsigned int no_sphe = 10;
  unsigned int no_rem = 20;
  unsigned int prec = 18;
  int seed = 7;

  if (argc < 2) print_usage(bench_name);

  for (int i=1 ; i < argc ;  i++) {
    bool random_flag = true;
    if (strcmp(argv[i],"insert")==0 && do_insert)
      do_remove = false;
    else if (strcmp(argv[i],"remove")==0 && do_remove)
      do_insert = false;
    else if (strcmp(argv[i],"--nospatial")==0)
      do_spatial_sort = false;
    else if (strcmp(argv[i],"--noread")==0)
      do_read = false;
    else if (strcmp(argv[i],"--nowrite")==0)
      do_write = false;
    else if (strcmp(argv[i],"-i")==0 && random_flag) {
      dfile = argv[++i];
      random_data = false;
      random_flag = false;
    }
    else if (strcmp(argv[i],"-c")==0 && random_flag) {
      no_cube = atoi(argv[++i]);
      no_sphe = atoi(argv[++i]);
      random_flag = false;
    }
    else if (strcmp(argv[i],"-r")==0) {
      no_rem = atoi(argv[++i]);
    }
    else if (strcmp(argv[i],"-s")==0)
      seed = atoi(argv[++i]);
    else if (strcmp(argv[i],"-f")==0)
      ifile = argv[++i];
    else if (strcmp(argv[i],"-w")==0)
      ofile = argv[++i];
    else
      print_usage(bench_name);
  }

  if (do_insert) {
    Timer t_data;
    if (random_data) {
      std::cout << "creating random points: " ;
      t_data.start();
      create_random_points<Point>( std::back_inserter(pts),
          no_cube, no_sphe, seed );
      t_data.stop();
    } else {
      std::cout << "reading test data: " ;
      std::ifstream data_in(dfile.c_str());
      assert(data_in.good());
      t_data.start();
      read_point_set_from_skel<Triang>( data_in, std::back_inserter(pts) );
      t_data.stop();
    }
    std::cout << t_data.time() << " sec" <<std::endl;
  }

  if (do_read) {
    Timer t_read;
    std::cout<<"reading Triangulation: " ;
    std::ifstream in(ifile.c_str());
    in.precision(prec);
    assert(in.good());
    t_read.start();
    in >> T;
    t_read.stop();
    std::cout << t_read.time() << " sec" <<std::endl;
  }

  if (do_insert) {
    Timer t_insert;
    std::cout<<"benchmarking insert: " ;
    if (do_spatial_sort) {
      t_insert.start();
      T.insert( pts.begin(), pts.end(), false );
      t_insert.stop();
    } else {
      t_insert.start();
      for (typename std::vector<Point>::iterator pit=pts.begin();
           pit!=pts.end(); ++pit)
        T.insert(*pit);
      t_insert.stop();
    }
    std::cout << t_insert.time() << " sec" <<std::endl;
  }

  if (do_remove) {
    Timer t_remove;
    std::cout<<"benchmarking remove: ";
    t_remove.start();
    test_remove( T, no_rem ) ;
    t_remove.stop();
    std::cout << t_remove.time() << " sec" <<std::endl;
  }

  if (do_write) {
    Timer t_write;
    std::cout<<"writing Triangulation: ";
    std::ofstream out(ofile.c_str());
    out.precision(prec);
    t_write.start();
    out << T;
    t_write.stop();
    std::cout << t_write.time() << " sec" <<std::endl;
  }

  std::cout<<"finished"<<std::endl;

  return 0;
}

#endif // TYPES_H
