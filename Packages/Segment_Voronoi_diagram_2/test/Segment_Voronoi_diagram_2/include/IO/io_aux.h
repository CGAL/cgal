#ifndef CGAL_SVD_IO_AUX_H
#define CGAL_SVD_IO_AUX_H


#include <iostream>
#include <string>


std::ostream& start_testing(std::ostream& os,const std::string& msg) {
  os << "testing " << msg << "... " << std::flush;
  return os;
}

void start_testing(const std::string& msg) {
  start_testing(std::cout, msg);
}

//---------------------------------------------------------------

std::ostream& end_testing(std::ostream& os, const std::string& msg) {
  os << "done." << std::endl;
  return os;
}

void end_testing(const std::string& msg) {
  end_testing(std::cout, msg);
}

//---------------------------------------------------------------

std::ostream& print_separator(std::ostream& os)
{
  char separator[] = "----------------------";
  os << separator << separator << separator << std::endl;
  return os;
}

void print_separator()
{
  print_separator(std::cout);
}

//---------------------------------------------------------------

template<class SVD>
std::ostream& print_all_input_sites(std::ostream& os, const SVD& svd)
{
  typename SVD::Input_sites_iterator it = svd.input_sites_begin();
  os << "Input sites: " << std::endl;
  for (; it != svd.input_sites_end(); ++it) {
    os << *it << std::endl;
  }
  return os;
}

template<class SVD>
void print_all_input_sites(const SVD& svd)
{
  print_all_input_sites(std::cout, svd);
}

//---------------------------------------------------------------

template<class SVD>
std::ostream& print_all_output_sites(std::ostream& os, const SVD& svd)
{
  typename SVD::Output_sites_iterator it = svd.output_sites_begin();
  os << "Output sites: " << std::endl;
  for (; it != svd.output_sites_end(); ++it) {
    os << *it << std::endl;
  }
  return os;
}

template<class SVD>
void print_all_output_sites(const SVD& svd)
{
  print_all_output_sites(std::cout, svd);
}

#endif // CGAL_SVD_IO_AUX_H
