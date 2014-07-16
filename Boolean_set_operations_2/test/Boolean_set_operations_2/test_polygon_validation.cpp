/*test file for polygon validation. Intended for testing the global functions defined at Gps_polygon_validation.h*/

#include <CGAL/basic.h>
#include <CGAL/assertions_behaviour.h>

#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Gps_segment_traits_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Boolean_set_operations_2/Gps_polygon_validation.h> 
#include <iterator>
#include <string>
#include <sstream>
#include <iostream>

// leda_rational, or Gmpq, or Quotient<MP_float>
typedef CGAL::Exact_rational                       Number_type;
typedef CGAL::Cartesian<Number_type>               Kernel;
typedef CGAL::Gps_segment_traits_2<Kernel>         Traits_2;
typedef Traits_2::Polygon_2                        Polygon_2;
typedef Traits_2::Polygon_with_holes_2             Polygon_with_holes_2;

/*test files:

1. val_test1.dat - invalid polygon with holes. The hole is relatively simple instead of strictly simple.
2. val_test2.dat - invalid polygon with holes. The hole overlaps the outer boundary (the intersection results in a polygon).
3. val_test3.dat - invalid polygon with holes. Two holes intersect (the intersection results in a polygon). 
4. val_test4.dat - invalid polygon with holes. Two holes intersect (one contains the other).
5. val_test5.dat -  invalid polygon with holes. Two holes share an edge. (non regularized intersection 
results in an edge).
6. val_test6.dat - invalid polygon with holes. A hole and the outer boundary share an edge. (non regularized intersection 
results in an edge).
7. val_test7.dat - invalid polygon with holes. The outer boundary is not relatively simple because a "crossover" occurs
    at an intersection
8. val_test8.dat - valid polygon with holes. Outer boundary is relatively simple.
9. val_test9.dat - valid polygon with holes. Outer Boundary and holes are pairwise disjoint except on vertices
*/


/*test an input file. isValid indicates the input polygon is valid. ErrorMsg is displayed if the validation result does't equal isValid */
bool testValidationForFile(const char * infilename, std::ofstream & outfile ,
                           bool isValid)
{
   std::ifstream input_file (infilename);

  if (! input_file.is_open()) {
    std::cerr << "Failed to open the " << infilename <<std::endl;
    return (false);
  }
  // Read a polygon with holes from a file.
  Polygon_2               outerP;
  unsigned int            num_holes;

  input_file >> outerP;
  input_file >> num_holes;

  std::vector<Polygon_2>  holes (num_holes);
  unsigned int            k;

  for (k = 0; k < num_holes; k++)
    input_file >> holes[k];
  Polygon_with_holes_2    P (outerP, holes.begin(), holes.end());
  Traits_2 tr;  
  bool testValid = CGAL::is_valid_polygon_with_holes(P,tr);
  bool res = true;
  if (testValid != isValid) {
    res=false;
    outfile<< "Error validating " << infilename <<std::endl;
    //outfile << "P = " ;
    //print_polygon_with_holes (P);
    outfile<<std::endl;  
  }  	  
  input_file.close();
  return res;
}

void
special_warnings(const char *,
                 const char* expr,
                 const char* file,
                 int         line,
                 const char* msg )
{
  std::cerr << "  // CGAL: check violation! THIS MESSAGE IS PROBABLY WANTED." << std::endl
            << "  // Expression : " << expr << std::endl
            << "  // File       : " << file << std::endl
            << "  // Line       : " << line << std::endl
            << "  // Explanation: " << msg << std::endl
            << "  // Refer to the bug-reporting instructions at http://www.cgal.org/bug_report.html"
            << std::endl;
}

int main()
{
  std::cerr << "Modify the w-a-r-n-i-n-g-s handler...\n";
  CGAL::set_warning_handler(special_warnings);
  std::string testfilePrefix = "data/validation/val_test";
  std::string testfileSuffix = ".dat";   
  const char* outputfilename = "data/validation/validation_test_output.txt"; 
  std::ofstream output_file (outputfilename);
  if (! output_file.is_open()) {
    std::cerr << "Failed to open the " << outputfilename <<std::endl;
    return (1);
  }
  
  int result = 0;  
  for (int i=1;i<10;i++) {
    std::stringstream strs;
    std::string si;
    strs << i;
    strs >> si;     
    std::string filename = testfilePrefix + si + testfileSuffix;
    const char *cfilename = filename.c_str();
    bool isValidPgn=false;
    if (i>7)
      isValidPgn=true;    
    bool res =  testValidationForFile(cfilename, output_file, isValidPgn);
    if (!res) {
        std::cout << "test " << i << " failed" << std::endl;
        result=1;
    }  
    else {
      std::cout <<"test " << i << " succeeded" << std::endl;      
    }
  }
  if (result == 0)
    std::cout <<"ALL TESTS SUCCEEDED!" << std::endl; 
  else
  {
    std::cout <<"SOME TESTS FAILED" << std::endl; 
    return 1;
  }
  return (0);
}

