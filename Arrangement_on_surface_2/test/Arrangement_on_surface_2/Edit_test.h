#ifndef CGAL_EDIT_TEST_H
#define CGAL_EDIT_TEST_H

#include <boost/tokenizer.hpp>

#include <CGAL/basic.h>
#include <CGAL/Timer.h>
#include <CGAL/Arrangement_2.h>

#include "IO_test.h"

/*! Edit test */
template <typename T_Traits>
class Edit_test : public IO_test<T_Traits> {
private:
  typedef T_Traits                                      Traits;
  typedef IO_test<Traits>                               Base;

public:
  typedef typename Base::Point_2                        Point_2;
  typedef typename Base::X_monotone_curve_2             X_monotone_curve_2;
  typedef typename Base::Curve_2                        Curve_2;

  typedef typename Base::Points_vector                  Points_vector;
  typedef typename Base::Xcurves_vector                 Xcurves_vector;
  typedef typename Base::Curves_vector                  Curves_vector;
  
  typedef CGAL::Arrangement_2<Traits>                   Arrangement_2;

  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Edge_const_iterator   Edge_const_iterator;
  typedef typename Arrangement_2::Vertex_const_iterator Vertex_const_iterator;

  typedef typename Points_vector::iterator              Point_iterator;
  typedef std::vector<CGAL::Object>                     Objects_vector;
  typedef Objects_vector::iterator                      Object_iterator;
  
protected:
  bool read_perform_opts(std::istream& is);

public:
  /*! Constructor */
  Edit_test();

  /*! Destructor */
  virtual ~Edit_test()
  {
    deallocate_arrangement();
    clear();
  }

  void set_filenames(const char* points_filename, const char* xcurves_filename,
                     const char* curves_filename, const char* cmds_filename);
  
  /*! Initialize the data structures */
  virtual bool init();
  
  /*! Perform the test */
  virtual bool perform();

  /*! Clear the data structures */
  virtual void clear();

  bool allocate_arrangement();

  void deallocate_arrangement();
  
  bool construct_arrangement();

  void clear_arrangement();
  
  /*! The arrangement */
  Arrangement_2* m_arr;  

private:
  /*! The input data file of the query points*/
  std::string m_filename_cmds;

  /*! The query points */
  Points_vector m_cmds;
};

/*!
 * Constructor. 
 * Accepts test data file name.
 */
template <typename T_Traits>
Edit_test<T_Traits>::Edit_test() :
  m_arr(NULL)
{}

/*! Set the file names */
template <typename T_Traits>
void Edit_test<T_Traits>::set_filenames(const char* points_filename,
                                        const char* xcurves_filename,
                                        const char* curves_filename,
                                        const char* cmds_filename)
{
  Base::set_filenames(points_filename, xcurves_filename, curves_filename);
  m_filename_cmds.assign(cmds_filename);
}

/*! Initialize the data structures */
template <typename T_Traits>
bool Edit_test<T_Traits>::init()
{
  if (!Base::init()) return false;
  return true;
}

/*! Clear the data structures */
template<class T_Traits>
void Edit_test<T_Traits>::clear()
{
  Base::clear();
  m_cmds.clear();
  m_filename_cmds.clear();
}

template <typename T_Traits>
bool Edit_test<T_Traits>::allocate_arrangement()
{
  if (!(m_arr = new Arrangement_2())) return false;
  return true;
}

template <typename T_Traits>
void Edit_test<T_Traits>::deallocate_arrangement()
{
  if (m_arr) {
    delete m_arr;
    m_arr = NULL;
  }
}

template <typename T_Traits>
void Edit_test<T_Traits>::clear_arrangement()
{
  if (m_arr) m_arr->clear();
}

// Perform the test
template <typename T_Traits>
bool Edit_test<T_Traits>::perform()
{
  CGAL::Timer timer;
  timer.reset(); timer.start();

  std::cout << std::endl;
  
  std::ifstream in_cmds(this->m_filename_cmds.c_str());
  if (!in_cmds.is_open()) {
    this->print_error(std::string("cannot open file ").append(this->m_filename_cmds));
    return false;
  }
  
  if (!read_perform_opts(in_cmds)) {
    in_cmds.close();
    return false;
  }
  in_cmds.close(); 

  timer.stop(); ///END
  std::cout << "Arrangement aggregate construction took " 
            << timer.time() << std::endl;  
  
  // Print the size of the arrangement.
  std::cout << "V = " << this->m_arr->number_of_vertices()
            << ",  E = " << this->m_arr->number_of_edges() 
            << ",  F = " << this->m_arr->number_of_faces() << std::endl;

  return true;
}

template <typename T_Traits>
bool Edit_test<T_Traits>::read_perform_opts(std::istream& is)
{
  typename Arrangement_2::Vertex_const_handle    vh_ref, vh_curr;
  typename Arrangement_2::Halfedge_const_handle  hh_ref, hh_curr;
  typename Arrangement_2::Face_const_handle      fh_ref, fh_curr;

  bool rc = true;
  
  std::string sline;
  while (this->skip_comments(is, sline)) {
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep(", \t\n\r");
    tokenizer tokens(sline, sep);
    for (tokenizer::iterator it = tokens.begin(); it != tokens.end(); ++it)
      std::cout << *it << std::endl;
    
    // std::istringstream line(sline);
    // char cmd;
    // line >> cmd;
  
    // if (cmd == 'a') {
    //   // Insert all into the arrangement
    //   CGAL::insert(*(this->m_arr), this->m_xcurves.begin(), this->m_xcurves.end());
    //   // insert(*(this->m_arr), m_points.begin(), m_points.end());
    //   CGAL::insert(*(this->m_arr), this->m_curves.begin(), this->m_curves.end());
    //   continue;
    // }
    
    // unsigned int id;
    // line >> id;
    // if (id >= this->m_xcurves.size()) {
    //   std::cerr << "Index of x-monotone curve " << id << " is out of range ("
    //             << this->m_xcurves.size() << ") in " << "m_filename_commands"
    //             << "!" << std::endl;
    //   rc = false;
    //   continue;
    // }
    // if (cmd == 'i') CGAL::insert(*(this->m_arr), this->m_xcurves[id]);

    // if (cmd == 'd') {
    //   if (!remove(this->m_xcurves[id])) rc = false;
    // }
  }

  return rc;
}

#endif
