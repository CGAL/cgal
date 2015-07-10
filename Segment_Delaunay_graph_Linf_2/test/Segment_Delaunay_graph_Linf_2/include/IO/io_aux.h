#ifndef CGAL_SDG_IO_AUX_H
#define CGAL_SDG_IO_AUX_H


#include <iostream>
#include <cassert>
#include <string>

bool assert_no_warning(bool b)
{
  assert(b);
  return b;
}

//---------------------------------------------------------------

std::ostream& start_testing(std::ostream& os,const std::string& msg) {
  os << "testing " << msg << "... " << std::flush;
  return os;
}

void start_testing(const std::string& msg) {
  start_testing(std::cout, msg);
}

//---------------------------------------------------------------

std::ostream& end_testing(std::ostream& os, const std::string& ) {
  os << "done." << std::endl;
  return os;
}

void end_testing(const std::string& msg) {
  end_testing(std::cout, msg);
}

//---------------------------------------------------------------

std::ostream& print_separator(std::ostream& os)
{
  os << std::endl;
  char separator[] = "----------------------";
  os << separator << separator << separator << std::endl;
  os << std::endl;
  return os;
}

void print_separator()
{
  print_separator(std::cout);
}

//---------------------------------------------------------------

template<class SDG>
std::ostream& print_all_input_sites(std::ostream& os, const SDG& sdg)
{
  typename SDG::Input_sites_iterator it = sdg.input_sites_begin();
  os << "Input sites: " << std::endl;
  for (; it != sdg.input_sites_end(); ++it) {
    os << *it << std::endl;
  }
  return os;
}

template<class SDG>
void print_all_input_sites(const SDG& sdg)
{
  print_all_input_sites(std::cout, sdg);
}

//---------------------------------------------------------------

template<class SDG>
std::ostream& print_all_output_sites(std::ostream& os, const SDG& sdg)
{
  typename SDG::Output_sites_iterator it = sdg.output_sites_begin();
  os << "Output sites: " << std::endl;
  for (; it != sdg.output_sites_end(); ++it) {
    os << *it << std::endl;
  }
  return os;
}

template<class SDG>
void print_all_output_sites(const SDG& sdg)
{
  print_all_output_sites(std::cout, sdg);
}

#endif // CGAL_SDG_IO_AUX_H
