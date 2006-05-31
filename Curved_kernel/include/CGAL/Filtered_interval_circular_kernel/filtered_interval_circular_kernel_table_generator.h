#ifndef FILTERED_INTERVAL_CIRCULAR_KERNEL_TABLE_GENERATOR_H_
#define FILTERED_INTERVAL_CIRCULAR_KERNEL_TABLE_GENERATOR_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <CGAL/basic.h>
#include <CGAL/Filtered_interval_circular_kernel/filtered_interval_circular_kernel_defines.h>

CGAL_BEGIN_NAMESPACE

class Filtered_interval_circular_kernel_table_generator {
public:

  typedef enum {
    BUILD_X_AND_Y = 0,
  	BUILD_X
  } KIND_OF_OUTPUT;

  typedef std::vector< double > Axe;

private:
  	
  double eps; 
  double max_r2_div_d2;
  std::string file_name;
  mutable Axe paraboloid_axe_xs; 
  mutable Axe paraboloid_axe_ys; 
  KIND_OF_OUTPUT ko;
  
  inline double f(double x) const {
	return std::sqrt(1 - x*x);
  }

  void build_table_x(double eps, double xmin, double xmax) const {
	if((f(xmin) - f(xmax)) <= eps) {
	  paraboloid_axe_xs.push_back(xmin);
	} else {
	  double xmid = (xmin + xmax)/2.0;
      build_table_x(eps, xmin, xmid);
      build_table_x(eps, xmid, xmax);        
	}
  }

  inline double g(double y) const {
	return std::sqrt(y);
  }

  void build_table_y(double eps, double ymin, double ymax) const {
	if((g(ymax) - g(ymin)) <= eps) {
	  paraboloid_axe_ys.push_back(ymin);
	} else {
	  double ymid = (ymin + ymax)/2.0;
      build_table_y(eps, ymin, ymid);
      build_table_y(eps, ymid, ymax);        
	}
 }
	
public:

  static inline const double default_epsilon() {
  	return 0.05;
  }
  
  static inline const double default_max_r_div_by_d() {
  	return 100.0;
  }
  
  static inline const std::string default_file_name() {
  	return std::string("table.ftck");
  }
  
  static inline const KIND_OF_OUTPUT default_kind_of_output() {
  	return BUILD_X_AND_Y;
  }

  Filtered_interval_circular_kernel_table_generator(
  	double epsilon = default_epsilon(), 
  	double max_r_div_by_d = default_max_r_div_by_d(),
  	std::string file = default_file_name(),
  	KIND_OF_OUTPUT kind_of_output = default_kind_of_output())
  	: eps(epsilon), max_r2_div_d2(max_r_div_by_d * max_r_div_by_d),
  	  file_name(file)  {
  	ko = kind_of_output;
    paraboloid_axe_xs.clear();
    build_table_x(eps/2.0,0.0,1.0);
    paraboloid_axe_xs.push_back(1.0);
  	if(ko == BUILD_X_AND_Y) { 
  	  paraboloid_axe_ys.clear();
  	  build_table_y(eps/2.0,0.0,
  	    max_r2_div_d2 + FCKTG_EPSILON);
  	  paraboloid_axe_ys.push_back(max_r2_div_d2);
  	}
  }
  
  ~Filtered_interval_circular_kernel_table_generator() {
  	paraboloid_axe_xs.clear();
  	paraboloid_axe_ys.clear();
  }

  friend inline std::ostream& operator<<(std::ostream &os, 
    const Filtered_interval_circular_kernel_table_generator &fcktp);
  
  void save() {
  	std::ofstream file;
    file.open (file_name.c_str(), std::ios::out|std::ios::binary);
    CGAL_assertion(file.is_open() && !file.bad());
    file << *this;
    file.close();
  }
  
  void save(std::string file_name_) {
  	std::ofstream file;
    file.open (file_name.c_str(), std::ios::out|std::ios::binary);
    CGAL_assertion(file.is_open() && !file.bad());
    file << *this;
    file.close();
  }
  	
}; // end of class: Filtered_interval_circular_kernel_table_generator

inline std::ostream& operator<<(std::ostream &os, 
  const Filtered_interval_circular_kernel_table_generator &fcktp) {
  	
  double r = fcktp.max_r2_div_d2;
  Filtered_interval_circular_kernel_table_generator::KIND_OF_OUTPUT ko = 
    fcktp.ko;  	

  os.write(reinterpret_cast<char *>(&ko),sizeof(ko));
  os.write(reinterpret_cast<char *>(&r),sizeof(r));
  
  unsigned nx = fcktp.paraboloid_axe_xs.size();
  os.write(reinterpret_cast<char *>(&nx),sizeof(nx));
  for(Filtered_interval_circular_kernel_table_generator::Axe::const_iterator 
    it = fcktp.paraboloid_axe_xs.begin(); 
    it != fcktp.paraboloid_axe_xs.end(); it++) {
  	double x = *it;
    os.write(reinterpret_cast<char *>(&x),sizeof(x));
  }
  
  if(ko == Filtered_interval_circular_kernel_table_generator::BUILD_X_AND_Y) {
    unsigned ny = fcktp.paraboloid_axe_ys.size();
    os.write(reinterpret_cast<char *>(&ny),sizeof(ny));
    for(Filtered_interval_circular_kernel_table_generator::Axe::const_iterator 
      it = fcktp.paraboloid_axe_ys.begin(); 
      it != fcktp.paraboloid_axe_ys.end(); it++) {
  	  double y = *it;
      os.write(reinterpret_cast<char *>(&y),sizeof(y));
    }  	
  }

  return os;
}

CGAL_END_NAMESPACE

#endif /*FILTERED_INTERVAL_CIRCULAR_KERNEL_TABLE_GENERATOR_H_*/
