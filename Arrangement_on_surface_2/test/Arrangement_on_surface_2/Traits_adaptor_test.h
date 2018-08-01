#ifndef CGAL_TRAITS_ADAPTOR_TEST_H
#define CGAL_TRAITS_ADAPTOR_TEST_H

#include <string>
#include <map>

#include "Traits_base_test.h"

template <typename GeomTraits, typename BaseGeomTraits>
class Traits_adaptor_test :
  public Traits_base_test<BaseGeomTraits>
{
public:
  typedef GeomTraits                                          Geom_traits;
  typedef BaseGeomTraits                                      Base_geom_traits;
  typedef Traits_base_test<Base_geom_traits>                  Base;

private:
  typedef Traits_adaptor_test<Geom_traits, Base_geom_traits>  Self;

  /*! A map between (strings) commands and (member functions) operations */
  typedef bool (Traits_adaptor_test::* Wrapper)(std::istringstream&);
  typedef std::map<std::string, Wrapper> Wrapper_map;
  typedef typename Wrapper_map::iterator Wrapper_iter;
  Wrapper_map m_wrappers;

  virtual bool exec(std::istringstream& str_stream,
                    const std::string& str_command,
                    bool& result)
  {
    Wrapper_iter wi = m_wrappers.find(str_command);
    str_stream.clear();
    if (wi == m_wrappers.end()) return true;
    Wrapper wrapper = (*wi).second;
    result = (this->*wrapper)(str_stream);
    return false;
  }

  //@{
  bool ta_compare_xy(std::istringstream&);
  bool ta_compare_xy_imp(std::istringstream&);

  bool ta_compare_y_at_x_left_wrapper(std::istringstream&);
  bool ta_compare_y_at_x_left_wrapper_imp(std::istringstream&,
                                          CGAL::Tag_false);
  bool ta_compare_y_at_x_left_wrapper_imp(std::istringstream&,
                                          CGAL::Tag_true);

  bool ta_is_in_x_range_wrapper(std::istringstream&);
  bool ta_compare_y_position_wrapper(std::istringstream&);
  bool ta_is_between_cw_wrapper(std::istringstream&);
  bool ta_compare_cw_around_point_wrapper(std::istringstream&);

  bool ta_are_mergeable_wrapper(std::istringstream&);
  bool ta_are_mergeable_wrapper_imp(std::istringstream&, CGAL::Tag_false);
  bool ta_are_mergeable_wrapper_imp(std::istringstream&, CGAL::Tag_true);

  bool ta_merge_wrapper(std::istringstream&);
  bool ta_merge_wrapper_imp(std::istringstream&, CGAL::Tag_false);
  bool ta_merge_wrapper_imp(std::istringstream&, CGAL::Tag_true);

  //@}

protected:
  /*! An instance of the traits */
  const Geom_traits& m_geom_traits;

public:
  /*! Construct.
   */
  Traits_adaptor_test(const Geom_traits& geom_traits);

  /*! Destruct.
   */
  virtual ~Traits_adaptor_test();
};

/*! Construct.
 * Accepts test data file name.
 */
template <typename GeomTraits, typename BaseGeomTraits>
Traits_adaptor_test<GeomTraits, BaseGeomTraits>::
Traits_adaptor_test(const GeomTraits& geom_traits) :
  Base(geom_traits),
  m_geom_traits(geom_traits)
{
  m_wrappers[std::string("compare_xy")] = &Self::ta_compare_xy;
  m_wrappers[std::string("compare_y_at_x_left")] =
    &Self::ta_compare_y_at_x_left_wrapper;
  m_wrappers[std::string("is_in_x_range")] = &Self::ta_is_in_x_range_wrapper;
  m_wrappers[std::string("compare_y_position")] =
    &Self::ta_compare_y_position_wrapper;
  m_wrappers[std::string("is_between_cw")] = &Self::ta_is_between_cw_wrapper;
  m_wrappers[std::string("compare_cw_around_point")] =
    &Self::ta_compare_cw_around_point_wrapper;
  m_wrappers[std::string("are_mergeable")] = &Self::ta_are_mergeable_wrapper;
  m_wrappers[std::string("merge")] = &Self::ta_merge_wrapper;
}

/*! Destruct.
 * Declares as virtual.
 */
template <typename GeomTraits, typename BaseGeomTraits>
Traits_adaptor_test<GeomTraits, BaseGeomTraits>::~Traits_adaptor_test() {}

template <typename GeomTraits, typename BaseGeomTraits>
bool Traits_adaptor_test<GeomTraits, BaseGeomTraits>::
ta_compare_xy(std::istringstream& str_stream)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  unsigned int exp_answer = this->get_expected_enum(str_stream);
  std::cout << "Test: compare_xy( " << this->m_xcurves[id1]
            << "," << this->m_xcurves[id2] << " ) ? " << exp_answer << " ";

  unsigned int real_answer =
    m_geom_traits.compare_xy_2_object()(this->m_xcurves[id1] ,
                                        this->m_xcurves[id2]);
  return this->compare(exp_answer, real_answer);
}

template <typename GeomTraits, typename BaseGeomTraits>
bool Traits_adaptor_test<GeomTraits, BaseGeomTraits>::
ta_compare_y_at_x_left_wrapper(std::istringstream& str_stream)
{
  typedef typename Geom_traits::Has_left_category          Has_left_category;
  return ta_compare_y_at_x_left_wrapper_imp(str_stream, Has_left_category());
}

template <typename GeomTraits, typename BaseGeomTraits>
bool Traits_adaptor_test<GeomTraits, BaseGeomTraits>::
ta_compare_y_at_x_left_wrapper_imp(std::istringstream&, CGAL::Tag_false)
{
  CGAL_error();
  return false;
}

template <typename GeomTraits, typename BaseGeomTraits>
bool Traits_adaptor_test<GeomTraits, BaseGeomTraits>::
ta_compare_y_at_x_left_wrapper_imp(std::istringstream& str_stream,
                                   CGAL::Tag_true)
{
  unsigned int id1, id2, id3;
  str_stream >> id1 >> id2 >> id3;
  unsigned int exp_answer = this->get_expected_enum(str_stream);
  std::cout << "Test: compare_y_at_x_left( " << this->m_xcurves[id1] << ","
            << this->m_xcurves[id2] << ", " << this->m_points[id3]
            << " ) ? " << exp_answer << " ";

  unsigned int real_answer =
    this->m_geom_traits.compare_y_at_x_left_2_object()(this->m_xcurves[id1],
                                                       this->m_xcurves[id2],
                                                       this->m_points[id3]);
  return this->compare(exp_answer, real_answer);
}

template <typename GeomTraits, typename BaseGeomTraits>
bool Traits_adaptor_test<GeomTraits, BaseGeomTraits>::
ta_is_in_x_range_wrapper(std::istringstream& str_stream)
{
  unsigned int id1, id2;
  char c;
  str_stream >> c >> id1 >> id2;
  bool exp_answer = this->get_expected_boolean(str_stream);
  std::cout << "Test: is_in_x_range( " << this->m_xcurves[id1] << ",";
  if (c == 'p') std::cout << this->m_points[id2];
  else if (c == 'x') std::cout << this->m_xcurves[id2];
  else CGAL_error();
  std::cout << " ) ? " << " ";

  bool real_answer = (c == 'p') ?
    m_geom_traits.is_in_x_range_2_object()(this->m_xcurves[id1],
                                           this->m_points[id2]) :
    m_geom_traits.is_in_x_range_2_object()(this->m_xcurves[id1],
                                           this->m_xcurves[id2]);
  return this->compare(exp_answer, real_answer);
}

template <typename GeomTraits, typename BaseGeomTraits>
bool Traits_adaptor_test<GeomTraits, BaseGeomTraits>::
ta_compare_y_position_wrapper(std::istringstream& str_stream)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  unsigned int exp_answer = this->get_expected_enum(str_stream);
  std::cout << "Test: compare_y_position( " << this->m_xcurves[id1]
            << "," << this->m_xcurves[id2] << " ) ? " << exp_answer << " ";

  unsigned int real_answer =
    m_geom_traits.compare_y_position_2_object()(this->m_xcurves[id1] ,
                                                this->m_xcurves[id2]);
  return this->compare(exp_answer, real_answer);
}

template <typename GeomTraits, typename BaseGeomTraits>
bool Traits_adaptor_test<GeomTraits, BaseGeomTraits>::
ta_is_between_cw_wrapper(std::istringstream& str_stream)
{
  unsigned int xcv , b , xcv1 , b1 , xcv2 , b2 , p;
  //note that b_ref1 b_ref2 are outputs so they can be tested also
  bool b_ref1,b_ref2;
  str_stream >> xcv >> b >> xcv1 >> b1 >> xcv2 >> b2 >> p;
  bool exp_answer = this->get_expected_boolean(str_stream);
  std::cout << "Test: is_between_cw( " << this->m_xcurves[xcv] << " , "
            << (b == 0 ? "false" : "true") << " , " << this->m_xcurves[xcv1]
            << " , "
            << (b1 == 0 ? "false" : "true") << " , " << this->m_xcurves[xcv2]
            << " , "
            << (b2 == 0 ? "false" : "true") << " , " << this->m_points[p]
            << " ) ? ";
  std::cout << exp_answer << " ";

  bool real_answer =
    m_geom_traits.is_between_cw_2_object()(this->m_xcurves[xcv], (b != 0),
                                           this->m_xcurves[xcv1], (b1 != 0),
                                           this->m_xcurves[xcv2], (b2 != 0),
                                           this->m_points[p], b_ref1, b_ref2);
  return this->compare(exp_answer, real_answer);
}

template <typename GeomTraits, typename BaseGeomTraits>
bool Traits_adaptor_test<GeomTraits, BaseGeomTraits>::
ta_compare_cw_around_point_wrapper(std::istringstream& str_stream)
{
  unsigned int xcv1 , b1 , xcv2 , b2 , p , b3;
  str_stream >> xcv1 >> b1 >> xcv2 >> b2 >> p >> b3;
  unsigned int exp_answer = this->get_expected_enum(str_stream);
  std::cout << "Test: compare_cw_around_point( " << this->m_xcurves[xcv1]
            << " , "
            << (b1 == 0 ? "false" : "true") << " , " << this->m_xcurves[xcv2]
            << " , "
            << (b2 == 0 ? "false" : "true") << " , " << this->m_points[p]
            << " , "
            << (b3 == 0 ? "false" : "true") << " ) ? ";
  std::cout << exp_answer << " ";

  unsigned int real_answer =
    m_geom_traits.compare_cw_around_point_2_object()(this->m_xcurves[xcv1],
                                                     (b1 != 0),
                                                     this->m_xcurves[xcv2],
                                                     (b2 != 0),
                                                     this->m_points[p],
                                                     (b3 != 0));
  return this->compare(exp_answer, real_answer);
}

template <typename GeomTraits, typename BaseGeomTraits>
bool Traits_adaptor_test<GeomTraits, BaseGeomTraits>::
ta_are_mergeable_wrapper(std::istringstream& str_stream)
{
  typedef typename GeomTraits::Has_merge_category      Has_merge_category;
  return ta_are_mergeable_wrapper_imp(str_stream, Has_merge_category());
}

template <typename GeomTraits, typename BaseGeomTraits>
bool Traits_adaptor_test<GeomTraits, BaseGeomTraits>::
ta_are_mergeable_wrapper_imp(std::istringstream&, CGAL::Tag_false)
{
  CGAL_error();
  return false;
}

template <typename GeomTraits, typename BaseGeomTraits>
bool Traits_adaptor_test<GeomTraits, BaseGeomTraits>::
ta_are_mergeable_wrapper_imp(std::istringstream& str_stream, CGAL::Tag_true)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  bool exp_answer = this->get_expected_boolean(str_stream);
  std::cout << "Test: are_mergeable( " << this->m_xcurves[id1] << ", "
            << this->m_xcurves[id2] << " ) ? " << exp_answer << " ";

  bool real_answer =
    this->m_geom_traits.are_mergeable_2_object()(this->m_xcurves[id1],
                                                 this->m_xcurves[id2]);
  return this->compare(exp_answer, real_answer);
}

template <typename GeomTraits, typename BaseGeomTraits>
bool Traits_adaptor_test<GeomTraits, BaseGeomTraits>::ta_merge_wrapper
(std::istringstream& str_stream)
{
  typedef typename GeomTraits::Has_merge_category      Has_merge_category;
  return ta_merge_wrapper_imp(str_stream, Has_merge_category());
}

template <typename GeomTraits, typename BaseGeomTraits>
bool Traits_adaptor_test<GeomTraits, BaseGeomTraits>::
ta_merge_wrapper_imp(std::istringstream&, CGAL::Tag_false)
{
  CGAL_error();
  return false;
}

template <typename GeomTraits, typename BaseGeomTraits>
bool Traits_adaptor_test<GeomTraits, BaseGeomTraits>::
ta_merge_wrapper_imp(std::istringstream& str_stream, CGAL::Tag_true)
{
  typedef typename Geom_traits::X_monotone_curve_2  X_monotone_curve_2;
  typedef typename Geom_traits::Equal_2             Equal_2;
  CGAL_USE_TYPE(Equal_2);

  unsigned int id1, id2, id;
  str_stream >> id1 >> id2 >> id;
  X_monotone_curve_2 cv;
  std::cout << "Test: merge( " << this->m_xcurves[id1] << ", "
            << this->m_xcurves[id2] << " ) ? " << this->m_xcurves[id] << " ";

  this->m_geom_traits.merge_2_object()(this->m_xcurves[id1],
                                       this->m_xcurves[id2],
                                       cv);
  return this->compare_curves(this->m_xcurves[id], cv);
}

#endif
