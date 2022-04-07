#ifndef CGAL_POINT_LOCATION_DYNAMIC_TEST_H
#define CGAL_POINT_LOCATION_DYNAMIC_TEST_H


#include "Point_location_test.h"

/*! Point location test */
template <typename T_Geom_traits, typename T_Topol_traits>
class Point_location_dynamic_test :
  public Point_location_test<T_Geom_traits, T_Topol_traits>
{
private:
  typedef T_Geom_traits                                  Geom_traits;
  typedef T_Topol_traits                                 Topol_traits;
  typedef Point_location_test<Geom_traits, Topol_traits> Base;

public:
  typedef typename Base::Point_2                         Point_2;
  typedef typename Base::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Base::Curve_2                         Curve_2;

  typedef typename Base::Arrangement                     Arrangement;

  /*! Constructor */
  Point_location_dynamic_test(const Geom_traits& geom_traits);

  /*! Destructor */
  virtual ~Point_location_dynamic_test()
  {
    this->deallocate_arrangement();
    clear();
    this->deallocate_pl_strategies();
  }

  /*! Clear the data structures */
  virtual void clear();

  void set_filenames(const char* points_filename, const char* xcurves_filename,
                     const char* curves_filename, const char* commands_filename,
                     const char* queries_filename)
  {
    Base::set_filenames(points_filename, xcurves_filename,
                        curves_filename, queries_filename);
    m_filename_commands.assign(commands_filename);
  }

  bool construct_arrangement();

private:
  bool read_perform_opts(std::istream& is);

  bool remove(const X_monotone_curve_2& xcv);

  /*! The input data file of commands*/
  std::string m_filename_commands;
};

/*!
 * Constructor.
 */
template <typename T_Geom_traits, typename T_Topol_traits>
Point_location_dynamic_test<T_Geom_traits, T_Topol_traits>::
Point_location_dynamic_test(const Geom_traits& geom_traits) :
  Base(geom_traits) {}

/*! Clear the data structures */
template <typename T_Geom_traits, typename T_Topol_traits>
void Point_location_dynamic_test<T_Geom_traits, T_Topol_traits>::clear()
{
  Base::clear();
  m_filename_commands.clear();
}

template <typename T_Geom_traits, typename T_Topol_traits>
bool Point_location_dynamic_test<T_Geom_traits, T_Topol_traits>::
construct_arrangement()
{
  std::ifstream in_com(this->m_filename_commands.c_str());
  if (!in_com.is_open()) {
    this->print_error(std::string("cannot open file ").
                      append(this->m_filename_commands));
    return false;
  }

  if (!read_perform_opts(in_com)) {
    in_com.close();
    return false;
  }
  in_com.close();

  // Print the size of the arrangement.
  std::cout << "V = " << this->m_arr->number_of_vertices()
            << ",  E = " << this->m_arr->number_of_edges()
            << ",  F = " << this->m_arr->number_of_faces() << std::endl;

  return true;
}

template <typename T_Geom_traits, typename T_Topol_traits>
bool Point_location_dynamic_test<T_Geom_traits, T_Topol_traits>::
read_perform_opts(std::istream& is)
{
  bool rc = true;

  CGAL::Timer timer;
  timer.reset();
  timer.start();

  std::string sline;
  while (this->skip_comments(is, sline)) {
    std::istringstream line(sline);
    char cmd;
    line >> cmd;

    if (cmd == 'a') {
      // Insert all into the arrangement
      CGAL::insert(*(this->m_arr),
                   this->m_xcurves.begin(), this->m_xcurves.end());
      // insert(*(this->m_arr), m_points.begin(), m_points.end());
      CGAL::insert(*(this->m_arr), this->m_curves.begin(),
                   this->m_curves.end());
      continue;
    }

    size_t id;
    line >> id;
    if (id >= this->m_xcurves.size()) {
      std::cerr << "Index of x-monotone curve " << id << " is out of range ("
                << this->m_xcurves.size() << ") in " << "m_filename_commands"
                << "!" << std::endl;
      rc = false;
      continue;
    }
    if (cmd == 'i') CGAL::insert(*(this->m_arr), this->m_xcurves[id]);

    if (cmd == 'd') {
      if (!remove(this->m_xcurves[id])) rc = false;
    }
  }
  timer.stop(); ///END
  std::cout << "Arrangement aggregate construction took "
            << timer.time() << std::endl;

  return rc;
}

template <typename T_Geom_traits, typename T_Topol_traits>
bool Point_location_dynamic_test<T_Geom_traits, T_Topol_traits>::
remove(const X_monotone_curve_2& xcv)
{
  typedef T_Geom_traits                                 Geom_traits;
  bool rc = false;          // be pasimistic, assume nothing is removed.

  const Geom_traits* traits = this->m_arr->geometry_traits();
  typename Geom_traits::Equal_2 equal = traits->equal_2_object();

  typename Arrangement::Edge_iterator eit = this->m_arr->edges_begin();
  for (; eit != this->m_arr->edges_end(); ++eit) {
    const X_monotone_curve_2& xcv_arr = eit->curve();
    if (equal(xcv, xcv_arr)) {
      this->m_arr->remove_edge(eit);
      rc = true;           // found a curve to remove.
      break;
    }
  }
  return rc;
}

#endif
