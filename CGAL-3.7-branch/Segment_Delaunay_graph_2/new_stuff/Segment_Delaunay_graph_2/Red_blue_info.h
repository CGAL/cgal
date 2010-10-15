#ifndef REB_BLUE_INFO_H
#define REB_BLUE_INFO_H 1

#include <CGAL/basic.h>
#include <iostream>
#include <CGAL/Random.h>

enum Red_blue {
  RED = 1,
  BLUE = 2,
  PURPLE = 3
};


std::ostream&
operator<<(std::ostream& os, const Red_blue& rb)
{
  if ( rb == RED ) { os << "Red"; }
  else if ( rb == BLUE ) { os << "Blue"; }
  else if ( rb == PURPLE ) { os << "Purple"; }
  else {
    std::cerr << "INFO: " << static_cast<int>(rb) << std::endl;
    CGAL_error();
  }
  return os;
}


struct Red_blue_convert_info
{
  typedef Red_blue   Info;
  typedef Info       result_type;

  inline
  Info operator()(const Info& info0, bool) const {
    return info0;
  }

  inline
  Info operator()(const Info& info0, const Info& info1, bool) const {
    return info0;
  }
};


struct Red_blue_merge_info
{
  typedef Red_blue   Info;
  typedef Info       result_type;

  inline
  Info operator()(const Info& info0, const Info& info1) const {
    if ( info0 == info1 ) { return info0; }
    return PURPLE;
  }
};


class Random_red_blue
{
public:
  Random_red_blue(int seed = 0) : r_(seed) {}

  Red_blue operator*() {
    double d = r_.get_double(0, 1);
    return (d < 0.5) ? RED : BLUE;
  }
private:
  CGAL::Random r_;
};

#endif // REB_BLUE_INFO_H
