#include <fstream>
#include <CGAL/IO/io.h>
#include <CGAL/Timer.h>

int main()
{
  std::cout.precision(17);

  std::ofstream out("data.xyz");
  out.precision(17);

  std::cout << "Write 3 million doubles" << std::endl;
  for(int i=0; i < 1000000; i++){
    out << i*0.000001 << " " << i *0.000001 << " " << i << std::endl;
  }
  out.close();

  double x,y,z;
  {
    CGAL::Timer t;
    t.start();
    double sum=0;
    
    std::ifstream in("data.xyz");
    while(in >> x >> y >> z){
      sum += (x + y + z);
    }
    
    std::cout << "operator>>(istream&,double&) Sum = " <<  sum << std::endl << t.time() << " sec\n";
  }    

 {
    CGAL::Timer t;
    t.start();
    double sum=0;
    
    std::ifstream in("data.xyz");
    while(CGAL::extract(in,x), CGAL::extract(in,y), CGAL::extract(in,z)){
      sum += (x + y + z);
    }
    
    std::cout << "extract(istream&, double&)   Sum = " <<  sum << std::endl << t.time() << " sec\n";
  }    
 return 0;
}
