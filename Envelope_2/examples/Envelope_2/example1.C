#include <CGAL/Envelope_2.h>

#include <iostream>
#include <CGAL/enum.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Arr_segment_traits_2.h>

#include <vector>


//typedef CGAL::Gmpq                                      NT;
typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits_2;

typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::Curve_2                               Curve_2;
typedef Traits_2::X_monotone_curve_2                    X_monotone_curve_2;

typedef CGAL::Envelope_2<Traits_2>                      Env_2;

typedef std::vector<Curve_2>                            CurveVector;



int main(int argc,char **argv)
{
  if (argc<2)
  {
    std::cerr<< "Missing file name\n";
    return -1;
  }
  std::ifstream in_file(argv[1]);
  if (!in_file.is_open()) {
    std::cerr << "Cannot open file: " << argv[1] << "!" << std::endl;
    return -1;
  }

  CurveVector curves;
  int num_of_curves;

  in_file >> num_of_curves;
  curves.resize(num_of_curves);
  for(int i = 0 ; i < num_of_curves ; i++)
  {
    NT x1 , y1 , x2 , y2;
#ifdef MACHINE_NATIVE
			read_native(in_file, x1);	
			read_native(in_file, y1);
			read_native(in_file, x2);
			read_native(in_file, y2);
#else
    in_file >> x1 >> y1 >> x2 >> y2;
#endif
    
    curves[i] = Curve_2(Point_2(x1,y1), Point_2(x2,y2));
  }
  

  Env_2 env;
  
  env.insert (curves.begin(), curves.end(), Env_2::LOWER);
  Env_2::Vertex_iterator vi = env.begin();
  Env_2::Vertex_iterator vi_end = env.end();
  

  while (vi != vi_end)
  {
    std::cout<<"( " <<*vi<<" )" << " => ";
    if(vi.right() == NULL)
      std::cout<<"NULL => ";
    else
      std::cout<<(*(vi.right()))<<" => ";
    ++vi;
  }
  std::cout<<"\n";


  return 0;
}

