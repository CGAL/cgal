#ifndef CGAL_IO_TEST_H
#define CGAL_IO_TEST_H

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "IO_base_test.h"

template <typename T_Traits>
class IO_test : public IO_base_test<T_Traits> {
public:
  typedef T_Traits                                      Traits;
  typedef IO_base_test<Traits>                          Base;

  typedef typename Base::Point_2                        Point_2;
  typedef typename Base::X_monotone_curve_2             X_monotone_curve_2;
  typedef typename Base::Curve_2                        Curve_2;

  typedef typename std::vector<Point_2>                 Points_vector;
  typedef typename std::vector<X_monotone_curve_2>      Xcurves_vector;
  typedef typename std::vector<Curve_2>                 Curves_vector;

  /*! Constructor */
  IO_test();

  /*! Destructor */
  virtual ~IO_test();
  
  bool read_xcurves(const char* filename, Xcurves_vector& xcurves);
  bool read_points(const char* filename, Points_vector& points);
  bool read_curves(const char* filename, Curves_vector& curves);

  /*! Set the file names */
  void set_filenames(const char* points_filename,
                     const char* xcurves_filename,
                     const char* curves_filename);
  
  /*! Parse the command line */
  virtual bool parse(int argc, char* argv[]);

  /*! Initialize the data structure */
  virtual bool init();

  /*! Clear the data structures */
  virtual void clear();

protected:
  /*! Skip comments */
  std::istream& skip_comments(std::istream& is, std::string& line);  

  /*! Remove blanks */
  std::string remove_blanks(char* str);

  /*! Print an error message */
  void print_error(const std::string& msg)
  { std::cerr << "Error: " << msg.c_str() << std::endl; }

  /*! Print the end-of-line */
  void print_eol()
  {
    std::cout << std::endl;
    m_eol_printed = true;
  }

  /*! Print information */
  void print_info(std::string& info,
                  bool start_line = true, bool end_line = true)
  {
    if (start_line && !m_eol_printed) print_eol();
    std::cout << info.c_str();
    if (end_line) print_eol();
  }

  /*! Print final results */
  void print_result(bool result)
  {
    std::string result_str((result) ? "Passed" : "Failed");
    print_info(result_str, false);
  }

  /*! Print expected answer and real answer */
  void print_answer(const std::string& exp, const std::string& real,
                    const std::string& str)
  {
    print_info(std::string("Expected ").append(str).append(": ").append(exp));
    print_info(std::string("Obtained ").append(str).append(": ").append(real));
  }

  /*! Indicates whether the end-of-line has been printed */
  bool m_eol_printed;

  /*! The input data file of points*/
  std::string m_filename_points;

  /*! The input data file of curves*/
  std::string m_filename_curves;

  /*! The input data file of xcurves*/
  std::string m_filename_xcurves;

  /*! The container of input points */
  Points_vector m_points;

  /*! The container of input curves */
  Curves_vector m_curves;

  /*! The container of x-monotone curves */
  Xcurves_vector m_xcurves;
};

/*!
 * Constructor. 
 * Accepts test data file name.
 */
template <typename T_Traits>
IO_test<T_Traits>::IO_test() : m_eol_printed(true) {}

/*!
 * Destructor. 
 */
template <typename T_Traits>
IO_test<T_Traits>::~IO_test() { clear(); }

/*! Set the file names */
template <typename T_Traits>
void IO_test<T_Traits>::set_filenames(const char* points_filename,
                                      const char* xcurves_filename,
                                      const char* curves_filename)
{
  m_filename_points.assign(points_filename);
  m_filename_xcurves.assign(xcurves_filename);
  m_filename_curves.assign(curves_filename);
}

template <typename T_Traits>
bool IO_test<T_Traits>::parse(int argc, char* argv[])
{
  if (argc < 4) {
    print_info(std::string("Usage: ").append(argv[0]).
               append(" points_file xcurves_file curves_file"));
    return false;
  }

  m_filename_points.assign(argv[1]);
  m_filename_xcurves.assign(argv[2]);
  m_filename_curves.assign(argv[3]);
  return true;
}

/*! Initialize the data structure */
template <typename T_Traits>
bool IO_test<T_Traits>::init()
{
  if (!read_points(m_filename_points.c_str(), m_points)) return false;
  if (!read_xcurves(m_filename_xcurves.c_str(), m_xcurves)) return false;
  if (!read_curves(m_filename_curves.c_str(), m_curves)) return false;
  return true;
}

/*! Clear the data structures */
template <typename T_Traits>
void IO_test<T_Traits>::clear()
{
  m_filename_points.clear();
  m_filename_xcurves.clear();
  m_filename_curves.clear();

  m_points.clear();
  m_curves.clear();
  m_xcurves.clear();
}

/*!
 * Skip comments. Comments start with the '#' character and extend to the
 * end of the line
 */
template <typename T_Traits>
std::istream& IO_test<T_Traits>::skip_comments(std::istream& is,
                                               std::string& line)
{
  while (std::getline(is, line))
    if ((line[0] != '#') && !line.empty()) break;
  return is;
}

/*!
 */
template <typename T_Traits>
std::string IO_test<T_Traits>::remove_blanks(char* str)
{
  std::string result = "";
  bool flag = false;
  //only alphanumeric characters and underscores are allowed
  for (; *str != '\0'; ++str) {
    if ((*str >= '0' && *str <= '9') || //digits
        (*str >= 'A' && *str <= 'Z') || //upper case letters
        (*str >= 'a' && *str <= 'z') || //lower case letters
         *str == '_') //underscores
    {
      if (!flag)
        flag = true;
      result += *str;
    }
    if (*str == ' ' && flag) break;
  }
  return result;
}

/*! */
template <typename T_Traits>
bool IO_test<T_Traits>::read_points(const char* filename, Points_vector& points)
{
  typedef T_Traits                              Traits;

  // read points from file into associative container
  std::ifstream p_stream(filename);
  if (!p_stream.is_open()) {
    std::cerr << "Cannot open file " << filename << "!" << std::endl;
    return false;
  }

  std::string line;
  while (skip_comments(p_stream, line)) {
    std::istringstream line_stream(line);
    typename Traits::Point_2 p;
    this->read_point(line_stream, p);
    points.push_back(p);
    line_stream.clear();
  }
  return true;
}

/*! */
template <typename T_Traits>
bool
IO_test<T_Traits>::read_xcurves(const char* filename, Xcurves_vector& xcurves)
{
  typedef T_Traits                              Traits;

  // read x-monotone curves from file into associative container
  std::ifstream xcv_stream(filename);
  if (!xcv_stream.is_open()) {
    std::cerr << "Cannot open file " << filename << "!" << std::endl;
    return false;
  }

  std::string line;
  while (skip_comments(xcv_stream, line)) {
    std::istringstream line_stream(line);
    typename Traits::X_monotone_curve_2 xcv;
    this->read_xcurve(line_stream, xcv);
    xcurves.push_back(xcv);
    line_stream.clear();
  }
  return true;
}

/*! */
template <typename T_Traits>
bool
IO_test<T_Traits>::read_curves(const char* filename, Curves_vector& curves)
{
  typedef T_Traits                              Traits;

  // Read curves from file into associative container
  std::ifstream cv_stream(filename);
  if (!cv_stream.is_open()) {
    std::cerr << "Cannot open file " << filename << "!" << std::endl;
    return false;
  }

  std::string line;
  while (skip_comments(cv_stream, line)) {
    std::istringstream line_stream(line);
    typename Traits::Curve_2 cv;
    this->read_curve(line_stream, cv);
    curves.push_back(cv);
    line_stream.clear();
  }
  return true;
}

#endif
