#include <iostream>
#include <string>

// Base class for all arrangement traits test classes
// Class is virtual. Define build_curve_list in your derivative.
template <class arr_traits>
class Arr_traits_test {

public:
  typedef arr_traits              Traits;
  typedef typename Traits::Curve  Curve;

protected:
  Traits tr;

  // a data structure to store test info connected with each curve
  struct Curve_with_info {
    Curve         cv;
    bool          is_x_monotone;
    unsigned      x_monotone_num;
    std::string   description;
    
    // Construction
    Curve_with_info() {};
    
    Curve_with_info(Curve cv_in, bool is_x, int parts, std::string desc) :
      cv(cv_in), is_x_monotone(is_x), x_monotone_num(parts), 
      description(desc)
    {};
  };
 
  // Builds list of curves, depends on template parameter arr_traits
  virtual void build_curve_list(std::list<Curve_with_info>& curve_list) = 0;

  // Checks is_x_monotone and make_x_monotone
  bool check_monotony_functions(std::list<Curve_with_info> curves_list, 
				bool flip)
  {
    bool test_success = true,
      is_x_monotone;
    Curve_with_info cv_w_info;
    std::list<Curve> x_monotone_parts;

    typename std::list<Curve_with_info>::iterator cit = curves_list.begin();
    for (;
	 cit != curves_list.end();
	 cit++) {
      if (flip) 
	cit->cv = tr.curve_flip(cit->cv);
      std::cout << cit->description.c_str();
      if (flip)
	std::cout << " flipped";
      std::cout << std::endl;
    
      is_x_monotone = tr.is_x_monotone(cit->cv);
      if (is_x_monotone != cit->is_x_monotone) {
	std::cout << "  Should ";
	std::cout << ((cit->is_x_monotone) ? "be " : "not be ");
	std::cout << "x-monotone but traits tells the opposite." << std::endl;

	test_success = false;
      }

      if ( ! is_x_monotone ) {
	x_monotone_parts.clear();
	tr.make_x_monotone(cit->cv, x_monotone_parts);
      
	if (x_monotone_parts.size() != cit->x_monotone_num) {
	  std::cout << "  cut into wrong number of parts." << std::endl;
	  test_success = false;
	}
      }
    } // for

    return test_success;
  }	       

public:

  // main test routine
  bool start() 
  {
    bool result = true;

    std::list<Curve_with_info> curves_list;
    
    build_curve_list(curves_list);

    std::cout << "Checking is_x_monotone and make_x_monotone" << std::endl;
    std::cout << "------------------------------------------" << std::endl;
    if ( ! check_monotony_functions(curves_list, false) ) // don't flip
      result = false;
    if ( !  check_monotony_functions(curves_list, true) )  // flip curves
      result = false;
     
    // if ( ! check_bla_function(this, that) )
    //  result = false;

    return result; 
  }

}; // Arr_traits_test
