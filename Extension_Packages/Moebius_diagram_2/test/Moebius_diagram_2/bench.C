#include <iostream>
#include <fstream>

#include "timer.h"

#include "mobius_cgal_types.h"

MD md;


void print_delay (std::ostream &os, double delay, struct tms *_tms)
{
  os << timer::convert (_tms->tms_utime) << " user, "
     << timer::convert (_tms->tms_stime) << " sys, "
     << delay << " total";
}

void load (std::ifstream &fin)
{
  std::istream_iterator<WPoint> start (fin);
  std::istream_iterator<WPoint> stop;

  struct tms _tms;
  double _delay;
  timer _timer, _total;

  _total.restart ();
  _timer.restart ();
  md.init (start, stop);
  _delay = _timer.elapsed (&_tms);
  std::cout << "       init: ";
  print_delay (std::cout, _delay, &_tms);
  std::cout << "\n";
  
  _timer.restart ();
  md.build ();
  _delay = _timer.elapsed (&_tms);
  std::cout << "      build: ";
  print_delay (std::cout, _delay, &_tms);
  std::cout << "\n";

  _delay = _total.elapsed (&_tms);
  std::cout << "  pre-total: ";
  print_delay (std::cout, _delay, &_tms);
  std::cout << "\n";


  _timer.restart ();
  md.construct ();
  _delay = _timer.elapsed (&_tms);
  std::cout << "  construct: ";
  print_delay (std::cout, _delay, &_tms);
  std::cout << "\n";

  _delay = _total.elapsed (&_tms);
  std::cout << "      total: ";
  print_delay (std::cout, _delay, &_tms);
  std::cout << "\n";

#if 0
  printf ("  stats: %5d/%5d edges\n"
	  "         %5d/%5d facets\n"
	  "         %5d/%5d cells\n",
	  md.number_of_edges (), md.rt().number_of_edges (),
	  md.number_of_facets (), md.rt().number_of_facets (),
	  md.number_of_cells (), md.rt().number_of_cells ());
#endif
}

int main(int argc, char **argv)
{
  
  for (int i = 1; i < argc; ++i) {
    md.clear ();
    std::cout << "loading " << argv[i] << "...\n";
    std::ifstream fin (argv[i]);
    load (fin);
  }
  return 0;
}
