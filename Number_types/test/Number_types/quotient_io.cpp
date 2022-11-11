#include <cassert>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <CGAL/Quotient.h>

int
main()
{
  std::ifstream S;
  std::ofstream T;
  std::ifstream U;
  S.open("quotient_in");
  assert( S );

  std::vector<CGAL::Quotient<double> >   V1;
  std::vector<CGAL::Quotient<double> >   V2;
  std::copy( std::istream_iterator<CGAL::Quotient<double> >(S),
             std::istream_iterator<CGAL::Quotient<double> >(),
             std::back_inserter(V1) );
  S.close();

  T.open("quotient_out");
  assert( T );
  std::copy( V1.begin(),
             V1.end(),
             std::ostream_iterator<CGAL::Quotient<double> >(T,"\n") );
  T.close();

  U.open("quotient_out");
  assert( T );
  std::copy( std::istream_iterator<CGAL::Quotient<double> >(U),
             std::istream_iterator<CGAL::Quotient<double> >(),
             std::back_inserter(V2) );
  U.close();


  if ( V1 != V2 )
  {
    std::cout << "V1" << std::endl;
    std::copy( V1.begin(), V1.end(), std::ostream_iterator<CGAL::Quotient<double> >(std::cout,"\n") );
    std::cout << "V2" << std::endl;
    std::copy( V2.begin(), V2.end(), std::ostream_iterator<CGAL::Quotient<double> >(std::cout,"\n") );
  }

  assert( V1 == V2 );

  return 0;
}
