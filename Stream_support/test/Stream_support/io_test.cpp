#include <CGAL/Bbox_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/OFF.h>
#include <CGAL/IO/io.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <array>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

int return_code = 0;

void state(const char* label, std::istream& is, bool expected)
{
  std::cout << label << " good=" << is.good() << " fail=" << is.fail()
            << " eof=" << is.eof() << "\n";
  if(!is.fail() != expected) {
    std::cerr << "ERROR: expected good=" << expected << " but got good=" << is.good() << "\n";
    return_code = 1;
  }
}

int main()
{
#ifdef CGAL_DISABLE_FROM_CHARS
  return 0;
#endif
  {
    std::istringstream is("+1.5 2 3");
    K::Point_3 p;
    is >> p;
    std::cout << "[leading plus] p=" << p << "  ";
    state("", is, true);
  }
  {
    std::istringstream is("1 2 3");
    K::Point_3 p;
    is >> p;
    std::cout << "[eof after last number] p=" << p << "  ";
    state("", is, true);
  }
  {
    std::istringstream is("1-2 3");
    K::Point_2 p;
    is >> p;
    std::cout << "[1-2 3 compact] p=" << p << "  ";
    state("", is, false);
  }
  {
    std::istringstream is("1 2 3 end");
    K::Point_3 p; std::string tok; double d;
    is >> p;
    is >> CGAL::IO::iformat(d);  // fails
    is.clear(); is >> tok;
    std::cout << "[failed number read, then token] token='" << tok << "'  ";
    state("", is, true);
    if(tok != "end") {
      std::cerr << "ERROR: [failed number read, then token] expected 'end', got '" << tok << "'\n";
      return_code = 1;
    }
  }
  {
    std::string off =
      "OFF\n3 1 0\n"
      "+0.0 +0.0 +0.0\n"
      "1.0 0.0 0.0\n"
      "0.0 1.0 0.0\n"
      "3 0 1 2\n";
    std::istringstream is(off);
    std::vector<K::Point_3> pts;
    std::vector<std::array<std::size_t,3> > faces;
    bool ok = CGAL::IO::read_OFF(is, pts, faces);
    std::cout << "[OFF with '+' coords] ok=" << ok << " npts=" << pts.size()
              << " nfaces=" << faces.size();
    if(!pts.empty()) std::cout << " p0=" << pts[0];
    std::cout << "\n";
  }
  return return_code;
}
