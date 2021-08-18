#define CGAL_GENERIC_EXTRACT 1

#include <cstdio>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <CGAL/Timer.h>
#include <CGAL/IO/io.h>
#include <boost/lexical_cast.hpp>



int main(int argc, char* argv[])
{
  CGAL::Timer t;
  int repeat = boost::lexical_cast<int>(argv[1]);
  int choice = boost::lexical_cast<int>(argv[2]);

   double sum=0;
  for(int j=0; j<repeat;++j){
    std::ifstream in(argv[3]);
    std::string off;
    int n,f,z;
    double d;
    in >> off >> n >> f >> z;
    n*=3;

  if(choice == 0){
    std::cerr << "operator"<< std::endl;
    t.start();
    for(int i=0; i<n; i++){
      in >> d;
      sum+= d;
    }
    t.stop();
  }

  if(choice == 1){
    std::cerr << "lexical_cast"<< std::endl;
    t.start();
    std::string sd;
    for(int i=0; i<n; i++){
      in >> sd;
      d = boost::lexical_cast<double>(sd);
      sum+= d;
    }
    t.stop();
  }

  if(choice == 2){
    std::cerr << "strtod"<< std::endl;
    t.start();
    std::string sd;
    for(int i=0; i<n; i++){
      in >> sd;
      d = strtod(sd.c_str(),NULL);
      sum+= d;
    }
    t.stop();
  }

  if(choice == 3){
    std::cerr << "sscanf"<< std::endl;
    t.start();
    std::string sd;
    for(int i=0; i<n; i++){
      in >> sd;
      sscanf(sd.c_str(),"%lf", &d);
      sum+= d;
    }
    t.stop();
  }

  if(choice == 4){
    std::cerr << "iformat"<< std::endl;
    t.start();
    for(int i=0; i<n; i++){
      in >> CGAL::IO::iformat(d);
      sum+= d;
    }
    t.stop();
  }
  }
  std::cerr.precision(17);
  std::cerr << "sum = " << sum << std::endl;
  std::cerr << t.time() << "sec."<< std::endl;

  return 0;
}
