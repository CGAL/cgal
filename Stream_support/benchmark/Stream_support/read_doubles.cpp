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
  int choice = 0;
  if(argc>1){
    choice = boost::lexical_cast<int>(argv[1]);
  }
  std::string off;
  int n,f,z;
  double d, sum=0;
  std::cin >> off >> n >> f >> z;
  n*=3;

  if(choice == 0){
    std::cerr << "operator"<< std::endl;
    t.start();
    for(int i=0; i<n; i++){
      std::cin >> d;
      sum+= d;
    }
    t.stop();
  }

  if(choice == 1){
    std::cerr << "lexical_cast"<< std::endl;
    t.start();
    std::string sd;
    for(int i=0; i<n; i++){
      std::cin >> sd;
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
      std::cin >> sd;
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
      std::cin >> sd;
      sscanf(sd.c_str(),"%lf", &d);
      sum+= d;
    }
    t.stop();
  } 

  if(choice == 4){
    std::cerr << "iformat"<< std::endl;
    t.start();
    for(int i=0; i<n; i++){
      std::cin >> CGAL::iformat(d);
      sum+= d;
    }
    t.stop();
  }
  std::cerr.precision(17);
  std::cerr << "sum = " << sum << std::endl;
  std::cerr << t.time() << "sec."<< std::endl;
    
  return 0;
}
