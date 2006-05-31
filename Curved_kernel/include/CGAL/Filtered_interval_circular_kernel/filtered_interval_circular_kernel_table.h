#ifndef FILTERED_INTERVAL_CIRCULAR_KERNEL_TABLE_H_
#define FILTERED_INTERVAL_CIRCULAR_KERNEL_TABLE_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <CGAL/basic.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Filtered_interval_circular_kernel/filtered_interval_circular_kernel_defines.h>
#include <CGAL/Filtered_interval_circular_kernel/filtered_interval_circular_kernel_table_generator.h>

CGAL_BEGIN_NAMESPACE

class Filtered_interval_circular_kernel_table {
private:

  typedef Interval_nt<false> Interval;
  typedef CGAL::Interval_nt<false>::Protector IntervalProtector;
 
public:

  typedef Filtered_interval_circular_kernel_table_generator::Axe Axe;
  typedef Filtered_interval_circular_kernel_table_generator::KIND_OF_OUTPUT
    KIND_OF_OUTPUT;
  typedef std::vector< Interval > TableX;
  typedef std::vector< TableX > TableXY;
  
  typedef enum {
  	WITH_TABLE,
  	WITH_X_VECTOR
  } KIND_OF_SEARCH;
  
  typedef enum {
     INTERSECTS=0,
     DONT_KNOW,
     DONT_INTERSECT
   } SEARCH_RESULT;
  
private:

  KIND_OF_OUTPUT ko;
  KIND_OF_SEARCH ks;
  std::string file_name;
  double max_r2_div_d2;
  mutable Axe xs_paraboloid;
  mutable Axe ys_paraboloid;
  mutable TableX table_x;
  mutable TableXY table_xy;
  unsigned capacity;

  friend std::ifstream& operator>>(std::ifstream &is,
    Filtered_interval_circular_kernel_table &table);
  
#ifdef DEBUG_FILTERED_INTERVAL_CIRCULAR_KERNEL_TABLE
  
    friend std::ofstream& operator<<(std::ofstream &os,
    Filtered_interval_circular_kernel_table &table);
  
#endif

  void open() {
  	std::ifstream file;
    file.open (file_name.c_str(),  
      std::ios::in|std::ios::binary);
    CGAL_assertion(!file.fail() && !file.bad());
    file >> *this;
    file.close();      
  }
  
  void make_table() const {
  	table_xy.clear();
  	for(unsigned j=0; j<ys_paraboloid.size()-1;j++) {
  	  TableX tx; 
  	  tx.clear();
  	  for(unsigned i=0; i<xs_paraboloid.size()-1;i++) {
        double ys0 = ys_paraboloid[j];
    	double ys1 = ys_paraboloid[j+1];
    	double xs0 = xs_paraboloid[i];
    	double xs1 = xs_paraboloid[i+1];
    	double z2 = std::sqrt(ys0 - ys0*xs1*xs1);
    	double z3 = std::sqrt(ys1 - ys1*xs0*xs0);      
    	tx.push_back(Interval(z2,z3));
      }
      table_xy.push_back(tx);
    }
  }
  
  void make_table_x() const {
  	table_x.clear();
  	for(unsigned i=0; i<xs_paraboloid.size()-1;i++) {
      double z2 = sqrt(1.0 - xs_paraboloid[i]*xs_paraboloid[i]);
      double z1 = sqrt(1.0 - xs_paraboloid[i+1]*xs_paraboloid[i+1]);
      table_x.push_back(Interval(z1,z2));
  	}
  }
  
  // I think it is a litte bit more efficient than the STL lower_bound()
  // but not sure, maybe needs tests
  inline void binary_search_in_x(double v, int &i) const {
    int min = 0, max = xs_paraboloid.size(); 
  	i = ((min+max)>>1);
  	do {
  		if(xs_paraboloid[i] <= v) min = i;
  		else max = i;
  		i = ((min+max)>>1);
  	} while(i != min); 
  }

  // I think it is a litte bit more efficient than the STL lower_bound()
  // but not sure, maybe needs tests
  inline void binary_search_in_y(double v, int &i) const {
    int min = 0, max = ys_paraboloid.size();
  	i = ((min+max)>>1);
  	do {
  		if(ys_paraboloid[i] <= v) min = i;
  		else max = i;
  		i = ((min+max)>>1);
  	} while(i != min); 
  }
  
  inline void search_table_x(Interval a_r, Interval &res) const {
  	int imin, imax;
  	binary_search_in_x(a_r.inf(), imin);
  	binary_search_in_x(a_r.sup(), imax); // very often imax = imin + 1
  	res = Interval(table_x[imin].inf(),table_x[imax].sup());
  }
  
  inline void search_table_xy(Interval r2_d2, 
    Interval a_d, Interval &res) const {
  	int imin, imax, jmin, jmax;
  	binary_search_in_x(r2_d2.inf(), imin);
  	binary_search_in_x(a_d.inf(), jmin);
  	binary_search_in_x(r2_d2.sup(), imax);
  	binary_search_in_x(a_d.sup(), jmax);
  	if(imin == imax && jmin == jmax) {
  	  res = table_xy[imin][jmin];
  	} else {
  	  res = Interval(table_xy[imin][jmax].inf(), 
  	    table_xy[imax][jmin].sup());	
  	} 
  }  

public:
  
  static inline SEARCH_RESULT test_intersection (
    Interval r2_d2, Interval a2_d2) {
    if(a2_d2.inf() > r2_d2.sup()) return DONT_INTERSECT;
    if(a2_d2.sup() >= r2_d2.inf()) return DONT_KNOW;
    return INTERSECTS;
  }
  
  static inline SEARCH_RESULT get_h_div_d_with_excellent_percision(
    Interval r2_d2, Interval a_d, Interval &res) {    	
    SEARCH_RESULT t;
    IntervalProtector ip;
    Interval a2_d2 = square(a_d);
    if((t = test_intersection(r2_d2, a2_d2)) == DONT_INTERSECT) return t;
    res = sqrt(r2_d2 - a2_d2);
    return t;
  }
  
  inline SEARCH_RESULT get_h_div_d_from_table(
    Interval r2_d2, Interval a_d, Interval &res) const {
    SEARCH_RESULT t;
    IntervalProtector ip;
    Interval a2_d2 = square(a_d);
    if((t = test_intersection(r2_d2, a2_d2)) == DONT_INTERSECT) return t;
    if(ks == WITH_X_VECTOR) {
      Interval r_d = sqrt(r2_d2);
      Interval a_r = a_d / r_d; 
      search_table_x(a_r, res);
      res = res * r_d;
    } else {
      if(r2_d2.sup() >= max_r2_div_d2) {
      	res = sqrt(r2_d2 - a2_d2);
      }	else {
        search_table_xy(r2_d2, a_d, res);	
      }
    }
    return t;
  }

  static inline const std::string default_file_name() {
  	return std::string("table.ftck");
  }
  
  static inline const unsigned default_maximum_capacity() {
  	return 1<<24;
  }

  Filtered_interval_circular_kernel_table(
    std::string file = default_file_name(),
    int max_capacity = default_maximum_capacity())
    : file_name(file), capacity(max_capacity) {
    open();
  }
  
  ~Filtered_interval_circular_kernel_table() {
  	xs_paraboloid.clear();
  	table_x.clear();
    ys_paraboloid.clear();
    for(TableXY::iterator it = table_xy.begin(); 
      it != table_xy.end(); it++) {
    	(*it).clear();
    } table_xy.clear();
  }

}; // end of class: Filtered_interval_circular_kernel_table

std::ifstream& operator>>(std::ifstream &is,
  Filtered_interval_circular_kernel_table &table) {
      	
  	is.read(reinterpret_cast<char *>(&table.ko), sizeof(table.ko));
  	is.read(reinterpret_cast<char *>(&table.max_r2_div_d2),
  	  sizeof(table.max_r2_div_d2));
  	CGAL_assertion(!is.fail());

    unsigned n;

    table.xs_paraboloid.clear();
    is.read(reinterpret_cast<char *>(&n),sizeof(n));  
    CGAL_assertion(!is.fail());
    CGAL_assertion((2*n-1) < table.capacity);
    for(unsigned i=0; i<n; i++) {
      double x;
      is.read(reinterpret_cast<char *>(&x), sizeof(x));
  	  CGAL_assertion(!is.fail());
  	  table.xs_paraboloid.push_back(x);
    } table.ks = Filtered_interval_circular_kernel_table::WITH_X_VECTOR;
    
    if(table.ko ==  Filtered_interval_circular_kernel_table_generator::
      BUILD_X_AND_Y) {	
      table.ys_paraboloid.clear();
      is.read(reinterpret_cast<char *>(&n),sizeof(n));  
      CGAL_assertion(!is.fail());
      if((((table.xs_paraboloid.size() - 1) * (n-1)) +
        2*table.xs_paraboloid.size() + n - 1) < table.capacity) {
        table.ks = Filtered_interval_circular_kernel_table::WITH_TABLE;
      }
      for(unsigned i=0; i<n; i++) {
    	double y;
    	is.read(reinterpret_cast<char *>(&y), sizeof(y));
  		CGAL_assertion(!is.fail());
  		if(table.ks == Filtered_interval_circular_kernel_table::WITH_TABLE) 
  		  table.ys_paraboloid.push_back(y);
      }  	
    }
    
    if(table.ks == Filtered_interval_circular_kernel_table::WITH_TABLE) 
      table.make_table();
    else table.make_table_x();
    
  	// Must reach the eof otherwise we have a parser error
  	char c;
  	is.read(reinterpret_cast<char *>(&c), sizeof(c));
  	CGAL_assertion(is.eof());
    
    return is;
}

#ifdef DEBUG_FILTERED_CURVED_KERNEL_TABLE_PARABOLOID

std::ofstream& operator<<(std::ofstream &os,
    Filtered_interval_circular_kernel_table &table) {
  os.write(reinterpret_cast<char *>(&table.ko),sizeof(table.ko));	
  os.write(reinterpret_cast<char *>(&table.max_r2_div_d2),
    sizeof(table.max_r2_div_d2));
  
  unsigned nx = table.xs_paraboloid.size();
  os.write(reinterpret_cast<char *>(&nx),sizeof(nx));
  for(Filtered_interval_circular_kernel_table_generator::Axe::const_iterator 
    it = table.xs_paraboloid.begin(); 
    it != table.xs_paraboloid.end(); it++) {
  	double x = *it;
    os.write(reinterpret_cast<char *>(&x),sizeof(x));
  }
  
  if(table.ko == 
    Filtered_interval_circular_kernel_table_generator::BUILD_X_AND_Y) {
    unsigned ny = table.ys_paraboloid.size();
    os.write(reinterpret_cast<char *>(&ny),sizeof(ny));
    for(Filtered_interval_circular_kernel_table_generator::Axe::const_iterator 
      it = table.ys_paraboloid.begin(); 
      it != table.ys_paraboloid.end(); it++) {
  	  double y = *it;
      os.write(reinterpret_cast<char *>(&y),sizeof(y));
    }  	
  }

  return os;
}
  
#endif

CGAL_END_NAMESPACE

#endif /*FILTERED_INTERVAL_CIRCULAR_KERNEL_TABLE_H_*/
