#ifndef PROGRESS_INDICATOR_H
#define PROGRESS_INDICATOR_H

#include <iostream>
#include <CGAL/assertions.h>

class Progress_indicator
{
 protected:
  long total, current;
 public:
  Progress_indicator( long n) : total(n), current(0) {
    CGAL_warning( total > 0);
  }
  void operator++(int) { 
    CGAL_warning( current != total);
    ++current; 
  }
  void operator++() {
    operator++(0);
  }
  float percentage() { 
    CGAL_assertion( total > 0);
    return 100.0*current/total; 
  }
};

class Progress_indicator_ostream : public Progress_indicator
{
  typedef Progress_indicator Base;
 protected:
  std::ostream& os;
  char separator;
 public:
  Progress_indicator_ostream( std::ostream& o, long n, char *msg, char s='\n') 
    : Base(n), os(o), separator(s) {
    os<<msg<<std::endl;
    os.precision(2);
    os.setf( std::ios::fixed | std::ios::right);    
  }
  void operator++(int) {
    Base::operator++();
    os.width(6);
    os<<percentage()<<'%'<<separator;
    os<<separator;
    os.flush();
  }
  void operator++() {
    operator++(0);
  }
};

class Progress_indicator_cout 
: public Progress_indicator_ostream
{
  typedef Progress_indicator_ostream Base;
 public:
  Progress_indicator_cout( long n, char *msg, char s='\r') 
    : Base( std::cout, n, msg, s) {}
};

#endif // PROGRESS_INDICATOR_H
