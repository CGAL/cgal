#ifndef CGAL_IO_TEST_H
#define CGAL_IO_TEST_H

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "IO_base_test.h"

template <typename Geom_traits_T>
class IO_test : public IO_base_test<Geom_traits_T> {
public:
  typedef Geom_traits_T                                 Geom_traits;
  typedef IO_base_test<Geom_traits>                     Base;

  typedef typename Base::Point_2                        Point_2;
  typedef typename Base::X_monotone_curve_2             X_monotone_curve_2;
  typedef typename Base::Curve_2                        Curve_2;

  typedef typename std::vector<Point_2>                 Points_vector;
  typedef typename std::vector<X_monotone_curve_2>      Xcurves_vector;
  typedef typename std::vector<Curve_2>                 Curves_vector;

#if TEST_GEOM_TRAITS == POLYCURVE_CONIC_GEOM_TRAITS ||          \
  TEST_GEOM_TRAITS == POLYCURVE_CIRCULAR_ARC_GEOM_TRAITS ||     \
  TEST_GEOM_TRAITS == POLYCURVE_BEZIER_GEOM_TRAITS ||           \
  TEST_GEOM_TRAITS == POLYLINE_GEOM_TRAITS ||                   \
  TEST_GEOM_TRAITS == NON_CACHING_POLYLINE_GEOM_TRAITS

  typedef typename Geom_traits_T::X_monotone_subcurve_2 X_monotone_subcurve_2;
  typedef typename Geom_traits_T::Subcurve_2            Subcurve_2;
  typedef typename std::vector<X_monotone_subcurve_2>   X_segment_vector;
  typedef typename std::vector<Subcurve_2>              Subcurve_vector;

  //vector containers for segments and xsegments
  X_segment_vector m_xsegments;
  Subcurve_vector m_segments;

#endif

  /*! Constructor */
  IO_test(const Geom_traits& traits);

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
  std::string remove_blanks(const char* str);

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
template <typename Geom_traits_T>
IO_test<Geom_traits_T>::IO_test(const Geom_traits_T& traits) :
  Base(traits),
  m_eol_printed(true)
{}

/*!
 * Destructor.
 */
template <typename Geom_traits_T>
IO_test<Geom_traits_T>::~IO_test() { clear(); }

/*! Set the file names */
template <typename Geom_traits_T>
void IO_test<Geom_traits_T>::set_filenames(const char* points_filename,
                                           const char* xcurves_filename,
                                           const char* curves_filename)
{
  m_filename_points.assign(points_filename);
  m_filename_xcurves.assign(xcurves_filename);
  m_filename_curves.assign(curves_filename);
}

template <typename Geom_traits_T>
bool IO_test<Geom_traits_T>::parse(int argc, char* argv[])
{
  /*
  The arguments are
  argv 0 is ./test_traits (string)
  argv 1 is data/polycurves_conics/compare_y_at_x.pt
  argv 2 is data/polycurves_conics/compare_y_at_x.xcv
  argv 3 is data/polycurves_conics/compare_y_at_x.cv
  argv 4 is data/polycurves_conics/compare_y_at_x
  argv 5 is polycurve_conic_traits (string)
  */
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
template <typename Geom_traits_T>
bool IO_test<Geom_traits_T>::init()
{
  if (!read_points(m_filename_points.c_str(), m_points)) return false;
  if (!read_xcurves(m_filename_xcurves.c_str(), m_xcurves)) return false;
  if (!read_curves(m_filename_curves.c_str(), m_curves)) return false;

  // for (int i = 0; i < m_points.size(); ++i)
  // {
  //   std::cout<< m_points[i] <<  " " ;
  // }
  return true;
}

/*! Clear the data structures */
template <typename Geom_traits_T>
void IO_test<Geom_traits_T>::clear()
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
template <typename Geom_traits_T>
std::istream& IO_test<Geom_traits_T>::skip_comments(std::istream& is,
                                                    std::string& line)
{
  while (std::getline(is, line))
    if (!line.empty() && (line[0] != '#')) break;
  return is;
}

/*!
 */
template <typename Geom_traits_T>
std::string IO_test<Geom_traits_T>::remove_blanks(const char* str)
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
      if (!flag) flag = true;
      result += *str;
    }
    if (*str == ' ' && flag) break;
  }
  return result;
}

/*! */
template <typename Geom_traits_T>
bool IO_test<Geom_traits_T>::read_points(const char* filename,
                                         Points_vector& points)
{
  typedef Geom_traits_T                              Geom_traits;

  // read points from file into associative container
  std::ifstream p_stream(filename);
  if (!p_stream.is_open()) {
    std::cerr << "Cannot open file " << filename << "!" << std::endl;
    return false;
  }

  std::string line;
  while (skip_comments(p_stream, line)) {
    std::istringstream line_stream(line);
    typename Geom_traits::Point_2 p;
    this->read_point(line_stream, p);
    points.push_back(p);
    line_stream.clear();
  }
  p_stream.close();
  return true;
}

/*! */
template <typename Geom_traits_T>
bool IO_test<Geom_traits_T>::read_xcurves(const char* filename,
                                          Xcurves_vector& xcurves)
{
  typedef Geom_traits_T                              Geom_traits;

  // read x-monotone curves from file into associative container
  std::ifstream xcv_stream(filename);

  if (!xcv_stream.is_open()) {
    std::cerr << "Cannot open file " << filename << "!" << std::endl;
    return false;
  }

  std::string line;

  while (skip_comments(xcv_stream, line)) {
#if TEST_GEOM_TRAITS == POLYCURVE_CONIC_GEOM_TRAITS ||            \
  TEST_GEOM_TRAITS == POLYCURVE_CIRCULAR_ARC_GEOM_TRAITS ||       \
  TEST_GEOM_TRAITS == POLYCURVE_BEZIER_GEOM_TRAITS ||             \
  TEST_GEOM_TRAITS == POLYLINE_GEOM_TRAITS ||                     \
  TEST_GEOM_TRAITS == NON_CACHING_POLYLINE_GEOM_TRAITS

    if (line[0] == 's') { //segment (see segment in 'Arr_polyline_traits.h')
      std::istringstream line_stream(line);
      typename Geom_traits::X_monotone_subcurve_2 xseg;
      this->read_xsegment(line_stream, xseg);
      m_xsegments.push_back(xseg);
      line_stream.clear();
    }
    else {
      std::istringstream line_stream(line);
      typename Geom_traits::X_monotone_curve_2 xcv;
      this->read_xcurve(line_stream, xcv);
      xcurves.push_back(xcv);
      line_stream.clear();
    }

#else

    std::istringstream line_stream(line);
    typename Geom_traits::X_monotone_curve_2 xcv;
    this->read_xcurve(line_stream, xcv);
    xcurves.push_back(xcv);
    line_stream.clear();

#endif
  }

  xcv_stream.close();
  return true;
}

/*! */
template <typename Geom_traits_T>
bool
IO_test<Geom_traits_T>::read_curves(const char* filename,
                                    Curves_vector& curves)
{
  typedef Geom_traits_T                         Geom_traits;

  // Read curves from file into associative container
  std::ifstream cv_stream(filename);

  if (!cv_stream.is_open()) {
    std::cerr << "Cannot open file " << filename << "!" << std::endl;
    return false;
  }

  std::string line;

  while (skip_comments(cv_stream, line)) {
#if TEST_GEOM_TRAITS == POLYCURVE_CONIC_GEOM_TRAITS ||             \
  TEST_GEOM_TRAITS == POLYCURVE_CIRCULAR_ARC_GEOM_TRAITS ||        \
  TEST_GEOM_TRAITS == POLYCURVE_BEZIER_GEOM_TRAITS ||              \
  TEST_GEOM_TRAITS == POLYLINE_GEOM_TRAITS ||                      \
  TEST_GEOM_TRAITS == NON_CACHING_POLYLINE_GEOM_TRAITS

    if (line[0] == 's') { //segment (see segment in 'Arr_polyline_traits.h')
      std::istringstream line_stream(line);
      typename Geom_traits::Subcurve_2 seg;
      this->read_segment(line_stream, seg);
      m_segments.push_back(seg);
      line_stream.clear();
    }
    else {
      std::istringstream line_stream(line);
      typename Geom_traits::Curve_2 cv;
      this->read_curve(line_stream, cv);
      curves.push_back(cv);
      line_stream.clear();
    }

#else
    std::istringstream line_stream(line);
    typename Geom_traits::Curve_2 cv;
    this->read_curve(line_stream, cv);
    curves.push_back(cv);
    line_stream.clear();
#endif
  }

  cv_stream.close();
  return true;
}

#endif
