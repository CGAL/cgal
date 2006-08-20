#include <CGAL/basic.h>

#include <CGAL/Nef_2/Polynomial.h>

#include <CGAL/test_macros.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer Integer;
#else
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz Integer;
#else
typedef int Integer;
#endif
#endif

using namespace CGAL;




#define PRT(t1,t2) std::cout<<"testing instances "<<#t1<<" "<<#t2<<std::endl

int main()
{
  //CGAL_NEF_SETDTHREAD(3); CGAL::set_pretty_mode ( std::cerr );
  CGAL_TEST_START;
  
  { PRT(Integer,Integer);
    typedef Integer NT; typedef Nef::Polynomial<Integer> RP;
      RP::NT seq[4] = { 0, 1, 2, 0 };
      RP p1, p2(NT(1)), p3(NT(1),NT(1)), p4(5,2), p5(-2,5), p6(4,1), 
         p7(3,0), p8(std::make_pair(&seq[0], &seq[4]));
      RP p10(-1,0,1), p11(-1,1), p12(1,1);
      NT r1(2), r2(-2);
      CGAL_TEST(p1.degree()==-1);
      CGAL_TEST(p2.degree()==0);
      CGAL_TEST(p4.degree()==1);
      CGAL_TEST(p7.degree()==0);
      CGAL_TEST(p8.degree()==2);
      CGAL_TEST((-(-p4)) == p4);
      CGAL_TEST((-(-p7)) == p7);
      CGAL_TEST((p4+p5) == RP(3,7));
      CGAL_TEST((p4-p5) == RP(7,-3));
      RP::NT prod[3] = { -10, 21, 10 };
      CGAL_TEST((p4*p5) == RP(std::make_pair(&prod[0],&prod[3])));
      CGAL_TEST((p2*p3) == p3);
      r1+p3;
      p3+r1;
      CGAL_TEST((r1+p3) == RP(3,1));
      CGAL_TEST((r1-p3) == RP(1,-1));
      CGAL_TEST((r1*p3) == RP(2,2));
      CGAL_TEST((p3+r1) == RP(3,1));
      CGAL_TEST((p3-r1) == RP(-1,1));
      CGAL_TEST((p3*r1) == RP(2,2));
      CGAL_TEST(p2 != p3);
      CGAL_TEST(p2 < p3);
      CGAL_TEST(p2 <= p3);
      CGAL_TEST(p5 > p4);
      CGAL_TEST(p5 >= p4);

      CGAL_TEST(r1 != p2);
      CGAL_TEST(r2 < p2);
      CGAL_TEST(r2 <= p2);
      CGAL_TEST(r1 > p2);
      CGAL_TEST(r1 >= p2);
      CGAL_TEST(p2 != r1);
      CGAL_TEST(p2 > r2);
      CGAL_TEST(p2 >= r2);
      CGAL_TEST(p2 < r1);
      CGAL_TEST(p2 <= r1);
      
      CGAL_TEST(CGAL_NTS sign(p5)==+1);
      CGAL_TEST(CGAL_NTS sign(-p5)==-1);
      CGAL_TEST(CGAL_NTS sign(p2)==+1);
      CGAL_TEST(CGAL_NTS sign(-p2)==-1);
      p3 += p2;
      p3 -= p2;
      p3 *= p5;
      p3 += r1;
      p3 -= r1;
      p3 *= r2;
      
      RP::NT D;
      RP q1(17),q2(5),q3,q4; 
      RP::pseudo_div(q1,q2,q3,q4,D); 
      CGAL_TEST(D*q1==q2*q3+q4);
      RP::pseudo_div(-q1,q2,q3,q4,D); 
      CGAL_TEST(D*-q1==q2*q3+q4);
      RP::pseudo_div(q1,-q2,q3,q4,D); 
      CGAL_TEST(D*q1==-q2*q3+q4);
      RP::pseudo_div(-q1,-q2,q3,q4,D); 
      CGAL_TEST(D*-q1==-q2*q3+q4);
      RP qq1(5),qq2(17),qq3,qq4;
      RP::pseudo_div(qq1,qq2,qq3,qq4,D); 
      CGAL_TEST(D*qq1==qq2*qq3+qq4);
      RP::pseudo_div(-qq1,qq2,qq3,qq4,D); 
      CGAL_TEST(D*-qq1==qq2*qq3+qq4);
      RP::pseudo_div(qq1,-qq2,qq3,qq4,D); 
      CGAL_TEST(D*qq1==-qq2*qq3+qq4);
      RP::pseudo_div(-qq1,-qq2,qq3,qq4,D); 
      CGAL_TEST(D*-qq1==-qq2*qq3+qq4);
      CGAL_TEST(p10/p11 == p12);

      q3 = RP::gcd(q1,q2);
      CGAL_TEST(q3 == 1);
      //CGAL_IO_TEST(p4,p1,CGAL::IO::BINARY);
      CGAL_IO_TEST(p4,p1,CGAL::IO::ASCII);
      CGAL_IO_TEST(p4,p1,CGAL::IO::PRETTY);
      CGAL::to_double(p6);
      CGAL::is_finite(p6);
      CGAL::is_valid(p6);
  }
  { PRT(int,Integer);
    typedef int NT; typedef Nef::Polynomial<Integer> RP;
      RP::NT seq[4] = { 0, 1, 2, 0 };
      RP p1, p2(NT(1)), p3(NT(1),NT(1)), p4(5,2), p5(-2,5), p6(4,1), 
         p7(3,0), p8(std::make_pair(&seq[0],&seq[4]));
      RP p10(-1,0,1), p11(-1,1), p12(1,1);
      NT r1(2), r2(-2);

      CGAL_TEST(p1.degree()==-1);
      CGAL_TEST(p2.degree()==0);
      CGAL_TEST(p4.degree()==1);
      CGAL_TEST(p7.degree()==0);
      CGAL_TEST(p8.degree()==2);
      CGAL_TEST((-(-p4)) == p4);
      CGAL_TEST((-(-p7)) == p7);
      CGAL_TEST((p4+p5) == RP(3,7));
      CGAL_TEST((p4-p5) == RP(7,-3));
      RP::NT prod[3] = { -10, 21, 10 };
      CGAL_TEST((p4*p5) == RP(std::make_pair(&prod[0],&prod[3])));
      CGAL_TEST((p2*p3) == p3);
      r1+p3;
      p3+r1;
      CGAL_TEST((r1+p3) == RP(3,1));
      CGAL_TEST((r1-p3) == RP(1,-1));
      CGAL_TEST((r1*p3) == RP(2,2));
      CGAL_TEST((p3+r1) == RP(3,1));
      CGAL_TEST((p3-r1) == RP(-1,1));
      CGAL_TEST((p3*r1) == RP(2,2));
      CGAL_TEST(p2 != p3);
      CGAL_TEST(p2 < p3);
      CGAL_TEST(p2 <= p3);
      CGAL_TEST(p5 > p4);
      CGAL_TEST(p5 >= p4);

      CGAL_TEST(r1 != p2);
      CGAL_TEST(r2 < p2);
      CGAL_TEST(r2 <= p2);
      CGAL_TEST(r1 > p2);
      CGAL_TEST(r1 >= p2);
      CGAL_TEST(p2 != r1);
      CGAL_TEST(p2 > r2);
      CGAL_TEST(p2 >= r2);
      CGAL_TEST(p2 < r1);
      CGAL_TEST(p2 <= r1);
      
      
      CGAL_TEST(CGAL_NTS sign(p5)==+1);
      CGAL_TEST(CGAL_NTS sign(-p5)==-1);
      CGAL_TEST(CGAL_NTS sign(p2)==+1);
      CGAL_TEST(CGAL_NTS sign(-p2)==-1);
      p3 += p2;
      p3 -= p2;
      p3 *= p5;
      p3 += r1;
      p3 -= r1;
      p3 *= r2;
      
      RP::NT D;
      RP q1(17),q2(5),q3,q4; 
      RP::pseudo_div(q1,q2,q3,q4,D); 
      CGAL_TEST(D*q1==q2*q3+q4);
      RP::pseudo_div(-q1,q2,q3,q4,D); 
      CGAL_TEST(D*-q1==q2*q3+q4);
      RP::pseudo_div(q1,-q2,q3,q4,D); 
      CGAL_TEST(D*q1==-q2*q3+q4);
      RP::pseudo_div(-q1,-q2,q3,q4,D); 
      CGAL_TEST(D*-q1==-q2*q3+q4);
      RP qq1(5),qq2(17),qq3,qq4;
      RP::pseudo_div(qq1,qq2,qq3,qq4,D); 
      CGAL_TEST(D*qq1==qq2*qq3+qq4);
      RP::pseudo_div(-qq1,qq2,qq3,qq4,D); 
      CGAL_TEST(D*-qq1==qq2*qq3+qq4);
      RP::pseudo_div(qq1,-qq2,qq3,qq4,D); 
      CGAL_TEST(D*qq1==-qq2*qq3+qq4);
      RP::pseudo_div(-qq1,-qq2,qq3,qq4,D); 
      CGAL_TEST(D*-qq1==-qq2*qq3+qq4);
      CGAL_TEST(p10/p11 == p12);

      q3 = RP::gcd(q1,q2);
      CGAL_TEST(q3 == 1);
      //CGAL_IO_TEST(p4,p1,CGAL::IO::BINARY);
      CGAL_IO_TEST(p4,p1,CGAL::IO::ASCII);
      CGAL_IO_TEST(p4,p1,CGAL::IO::PRETTY);
      CGAL::to_double(p6);
      CGAL::is_finite(p6);
      CGAL::is_valid(p6);

  }
  
  { PRT(double,Integer);
    typedef double NT; typedef Nef::Polynomial<Integer> RP;
      RP::NT seq[4] = { 0, 1, 2, 0 };
      RP p1, p2(NT(1)), p3(NT(1),NT(1)), p4(5,2), p5(-2,5), p6(4,1), 
         p7(3,0), p8(std::make_pair(&seq[0],&seq[4]));
      RP p10(-1,0,1), p11(-1,1), p12(1,1);

      NT r1(2), r2(-2);

      CGAL_TEST(p1.degree()==-1);
      CGAL_TEST(p2.degree()==0);
      CGAL_TEST(p4.degree()==1);
      CGAL_TEST(p7.degree()==0);
      CGAL_TEST(p8.degree()==2);
      CGAL_TEST((-(-p4)) == p4);
      CGAL_TEST((-(-p7)) == p7);
      CGAL_TEST((p4+p5) == RP(3,7));
      CGAL_TEST((p4-p5) == RP(7,-3));
      RP::NT prod[3] = { -10, 21, 10 };
      CGAL_TEST((p4*p5) == RP(std::make_pair(&prod[0],&prod[3])));
      CGAL_TEST((p2*p3) == p3);
      r1+p3;
      p3+r1;
      CGAL_TEST((r1+p3) == RP(3,1));
      CGAL_TEST((r1-p3) == RP(1,-1));
      CGAL_TEST((r1*p3) == RP(2,2));
      CGAL_TEST((p3+r1) == RP(3,1));
      CGAL_TEST((p3-r1) == RP(-1,1));
      CGAL_TEST((p3*r1) == RP(2,2));
      CGAL_TEST(p2 != p3);
      CGAL_TEST(p2 < p3);
      CGAL_TEST(p2 <= p3);
      CGAL_TEST(p5 > p4);
      CGAL_TEST(p5 >= p4);

      CGAL_TEST(r1 != p2);
      CGAL_TEST(r2 < p2);
      CGAL_TEST(r2 <= p2);
      CGAL_TEST(r1 > p2);
      CGAL_TEST(r1 >= p2);
      CGAL_TEST(p2 != r1);
      CGAL_TEST(p2 > r2);
      CGAL_TEST(p2 >= r2);
      CGAL_TEST(p2 < r1);
      CGAL_TEST(p2 <= r1);
      
      
      CGAL_TEST(CGAL_NTS sign(p5)==+1);
      CGAL_TEST(CGAL_NTS sign(-p5)==-1);
      CGAL_TEST(CGAL_NTS sign(p2)==+1);
      CGAL_TEST(CGAL_NTS sign(-p2)==-1);
      p3 += p2;
      p3 -= p2;
      p3 *= p5;
      p3 += r1;
      p3 -= r1;
      p3 *= r2;
      
      RP::NT D;
      RP q1(17),q2(5),q3,q4; 
      RP::pseudo_div(q1,q2,q3,q4,D); 
      CGAL_TEST(D*q1==q2*q3+q4);
      RP::pseudo_div(-q1,q2,q3,q4,D); 
      CGAL_TEST(D*-q1==q2*q3+q4);
      RP::pseudo_div(q1,-q2,q3,q4,D); 
      CGAL_TEST(D*q1==-q2*q3+q4);
      RP::pseudo_div(-q1,-q2,q3,q4,D); 
      CGAL_TEST(D*-q1==-q2*q3+q4);
      RP qq1(5),qq2(17),qq3,qq4;
      RP::pseudo_div(qq1,qq2,qq3,qq4,D); 
      CGAL_TEST(D*qq1==qq2*qq3+qq4);
      RP::pseudo_div(-qq1,qq2,qq3,qq4,D); 
      CGAL_TEST(D*-qq1==qq2*qq3+qq4);
      RP::pseudo_div(qq1,-qq2,qq3,qq4,D); 
      CGAL_TEST(D*qq1==-qq2*qq3+qq4);
      RP::pseudo_div(-qq1,-qq2,qq3,qq4,D); 
      CGAL_TEST(D*-qq1==-qq2*qq3+qq4);
      CGAL_TEST(p10/p11 == p12);

      q3 = RP::gcd(q1,q2);
      CGAL_TEST(q3 == 1);
      //CGAL_IO_TEST(p4,p1,CGAL::IO::BINARY);
      CGAL_IO_TEST(p4,p1,CGAL::IO::ASCII);
      CGAL_IO_TEST(p4,p1,CGAL::IO::PRETTY);
      CGAL::to_double(p6);
      CGAL::is_finite(p6);
      CGAL::is_valid(p6);

  }
  
  { PRT(int,int);
    typedef int NT; typedef Nef::Polynomial<int> RP;
      RP::NT seq[4] = { 0, 1, 2, 0 };
      RP p1, p2(NT(1)), p3(NT(1),NT(1)), p4(5,2), p5(-2,5), p6(4,1), 
         p7(3,0), p8(std::make_pair(&seq[0],&seq[4]));
      RP p10(-1,0,1), p11(-1,1), p12(1,1);
      NT r1(2), r2(-2);

      CGAL_TEST(p1.degree()==-1);
      CGAL_TEST(p2.degree()==0);
      CGAL_TEST(p4.degree()==1);
      CGAL_TEST(p7.degree()==0);
      CGAL_TEST(p8.degree()==2);
      CGAL_TEST((-(-p4)) == p4);
      CGAL_TEST((-(-p7)) == p7);
      CGAL_TEST((p4+p5) == RP(3,7));
      CGAL_TEST((p4-p5) == RP(7,-3));
      RP::NT prod[3] = { -10, 21, 10 };
      CGAL_TEST((p4*p5) == RP(std::make_pair(&prod[0],&prod[3])));
      CGAL_TEST((p2*p3) == p3);
      r1+p3;
      p3+r1;
      CGAL_TEST((r1+p3) == RP(3,1));
      CGAL_TEST((r1-p3) == RP(1,-1));
      CGAL_TEST((r1*p3) == RP(2,2));
      CGAL_TEST((p3+r1) == RP(3,1));
      CGAL_TEST((p3-r1) == RP(-1,1));
      CGAL_TEST((p3*r1) == RP(2,2));
      CGAL_TEST(p2 != p3);
      CGAL_TEST(p2 < p3);
      CGAL_TEST(p2 <= p3);
      CGAL_TEST(p5 > p4);
      CGAL_TEST(p5 >= p4);

      CGAL_TEST(r1 != p2);
      CGAL_TEST(r2 < p2);
      CGAL_TEST(r2 <= p2);
      CGAL_TEST(r1 > p2);
      CGAL_TEST(r1 >= p2);
      CGAL_TEST(p2 != r1);
      CGAL_TEST(p2 > r2);
      CGAL_TEST(p2 >= r2);
      CGAL_TEST(p2 < r1);
      CGAL_TEST(p2 <= r1);
      
      
      CGAL_TEST(CGAL_NTS sign(p5)==+1);
      CGAL_TEST(CGAL_NTS sign(-p5)==-1);
      CGAL_TEST(CGAL_NTS sign(p2)==+1);
      CGAL_TEST(CGAL_NTS sign(-p2)==-1);
      p3 += p2;
      p3 -= p2;
      p3 *= p5;
      p3 += r1;
      p3 -= r1;
      p3 *= r2;
      
      RP::NT D;
      RP q1(17),q2(5),q3,q4; 
      RP::pseudo_div(q1,q2,q3,q4,D); 
      CGAL_TEST(D*q1==q2*q3+q4);
      RP::pseudo_div(-q1,q2,q3,q4,D); 
      CGAL_TEST(D*-q1==q2*q3+q4);
      RP::pseudo_div(q1,-q2,q3,q4,D); 
      CGAL_TEST(D*q1==-q2*q3+q4);
      RP::pseudo_div(-q1,-q2,q3,q4,D); 
      CGAL_TEST(D*-q1==-q2*q3+q4);
      RP qq1(5),qq2(17),qq3,qq4;
      RP::pseudo_div(qq1,qq2,qq3,qq4,D); 
      CGAL_TEST(D*qq1==qq2*qq3+qq4);
      RP::pseudo_div(-qq1,qq2,qq3,qq4,D); 
      CGAL_TEST(D*-qq1==qq2*qq3+qq4);
      RP::pseudo_div(qq1,-qq2,qq3,qq4,D); 
      CGAL_TEST(D*qq1==-qq2*qq3+qq4);
      RP::pseudo_div(-qq1,-qq2,qq3,qq4,D); 
      CGAL_TEST(D*-qq1==-qq2*qq3+qq4);
      CGAL_TEST(p10/p11 == p12);

      q3 = RP::gcd(q1,q2);
      CGAL_TEST(q3 == 1);
      //CGAL_IO_TEST(p4,p1,CGAL::IO::BINARY);
      CGAL_IO_TEST(p4,p1,CGAL::IO::ASCII);
      CGAL_IO_TEST(p4,p1,CGAL::IO::PRETTY);
      CGAL::to_double(p6);
      CGAL::is_finite(p6);
      CGAL::is_valid(p6);

  }
  
  { PRT(double,int);
    typedef double NT; typedef Nef::Polynomial<int> RP;
      RP::NT seq[4] = { 0, 1, 2, 0 };
      RP p1, p2(NT(1)), p3(NT(1),NT(1)), p4(5,2), p5(-2,5), p6(4,1), 
         p7(3,0), p8(std::make_pair(&seq[0],&seq[4]));
      RP p10(-1,0,1), p11(-1,1), p12(1,1);
      NT r1(2), r2(-2);

      CGAL_TEST(p1.degree()==-1);
      CGAL_TEST(p2.degree()==0);
      CGAL_TEST(p4.degree()==1);
      CGAL_TEST(p7.degree()==0);
      CGAL_TEST(p8.degree()==2);
      CGAL_TEST((-(-p4)) == p4);
      CGAL_TEST((-(-p7)) == p7);
      CGAL_TEST((p4+p5) == RP(3,7));
      CGAL_TEST((p4-p5) == RP(7,-3));
      RP::NT prod[3] = { -10, 21, 10 };
      CGAL_TEST((p4*p5) == RP(std::make_pair(&prod[0],&prod[3])));
      CGAL_TEST((p2*p3) == p3);
      r1+p3;
      p3+r1;
      CGAL_TEST((r1+p3) == RP(3,1));
      CGAL_TEST((r1-p3) == RP(1,-1));
      CGAL_TEST((r1*p3) == RP(2,2));
      CGAL_TEST((p3+r1) == RP(3,1));
      CGAL_TEST((p3-r1) == RP(-1,1));
      CGAL_TEST((p3*r1) == RP(2,2));
      CGAL_TEST(p2 != p3);
      CGAL_TEST(p2 < p3);
      CGAL_TEST(p2 <= p3);
      CGAL_TEST(p5 > p4);
      CGAL_TEST(p5 >= p4);

      CGAL_TEST(r1 != p2);
      CGAL_TEST(r2 < p2);
      CGAL_TEST(r2 <= p2);
      CGAL_TEST(r1 > p2);
      CGAL_TEST(r1 >= p2);
      CGAL_TEST(p2 != r1);
      CGAL_TEST(p2 > r2);
      CGAL_TEST(p2 >= r2);
      CGAL_TEST(p2 < r1);
      CGAL_TEST(p2 <= r1);
      
      
      CGAL_TEST(CGAL_NTS sign(p5)==+1);
      CGAL_TEST(CGAL_NTS sign(-p5)==-1);
      CGAL_TEST(CGAL_NTS sign(p2)==+1);
      CGAL_TEST(CGAL_NTS sign(-p2)==-1);
      p3 += p2;
      p3 -= p2;
      p3 *= p5;
      p3 += r1;
      p3 -= r1;
      p3 *= r2;
      
      RP::NT D;
      RP q1(17),q2(5),q3,q4; 
      RP::pseudo_div(q1,q2,q3,q4,D); 
      CGAL_TEST(D*q1==q2*q3+q4);
      RP::pseudo_div(-q1,q2,q3,q4,D); 
      CGAL_TEST(D*-q1==q2*q3+q4);
      RP::pseudo_div(q1,-q2,q3,q4,D); 
      CGAL_TEST(D*q1==-q2*q3+q4);
      RP::pseudo_div(-q1,-q2,q3,q4,D); 
      CGAL_TEST(D*-q1==-q2*q3+q4);
      RP qq1(5),qq2(17),qq3,qq4;
      RP::pseudo_div(qq1,qq2,qq3,qq4,D); 
      CGAL_TEST(D*qq1==qq2*qq3+qq4);
      RP::pseudo_div(-qq1,qq2,qq3,qq4,D); 
      CGAL_TEST(D*-qq1==qq2*qq3+qq4);
      RP::pseudo_div(qq1,-qq2,qq3,qq4,D); 
      CGAL_TEST(D*qq1==-qq2*qq3+qq4);
      RP::pseudo_div(-qq1,-qq2,qq3,qq4,D); 
      CGAL_TEST(D*-qq1==-qq2*qq3+qq4);
      CGAL_TEST(p10/p11 == p12);

      q3 = RP::gcd(q1,q2);
      CGAL_TEST(q3 == 1);
      //CGAL_IO_TEST(p4,p1,CGAL::IO::BINARY);
      CGAL_IO_TEST(p4,p1,CGAL::IO::ASCII);
      CGAL_IO_TEST(p4,p1,CGAL::IO::PRETTY);
      CGAL::to_double(p6);
      CGAL::is_finite(p6);
      CGAL::is_valid(p6);

  }


  CGAL_TEST_END;
}


