#include <iostream>
#include <list>
#include <vector>

#include <CGAL/Nested_iterator.h>

int main()
{
  std::vector< std::vector<int> > v;

  std::vector<int> empty;
  std::vector<int> v1;
  std::vector<int> v2;
  std::vector<int> v3;
  std::vector<int> v4;


  int n = 5;
  for (int i = 0; i < n; i++) {
    v1.push_back(i);
    v2.push_back(i + n * 10);
    v3.push_back(i + n * 100);
    v4.push_back(i + n * 1000);
  }

  v.push_back(empty);
  v.push_back(empty);
  v.push_back(v1);
  v.push_back(v2);
  v.push_back(empty);
  v.push_back(v3);
  v.push_back(empty);
  v.push_back(v4);
  v.push_back(empty);
  v.push_back(empty);
  v.push_back(empty);

  {
    for (unsigned int i = 0; i < v.size(); i++) {
      std::cout << "List " << (i+1) << ": ";
      for (unsigned int j = 0; j < v[i].size(); j++) {
	std::cout << " " << v[i][j];
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  typedef std::vector< std::vector<int> >::iterator b_iterator;

  typedef CGAL::Nested_iterator<b_iterator>  n_iterator;

  //------------------------------------------------------------------
  {
    std::cout << "================================================"
	      << std::endl << std::endl;
    n_iterator ni_begin(v.begin(), v.end(), v.begin());
    n_iterator ni_end(v.begin(), v.end(), v.end());


    std::cout << "Nested iterator:" << std::endl;
    n_iterator ni_it;
    for (ni_it = ni_begin; ni_it != ni_end; ++ni_it) {
      std::cout << " " << (*ni_it);
    }
    std::cout << std::endl << std::endl;


    std::cout << "Nested iterator in reverse order:" << std::endl;
    for (ni_it = --ni_end; ni_it != ni_begin; --ni_it) {
      std::cout << " " << (*ni_it);
    }
    std::cout << " " << (*ni_it);
    std::cout << std::endl << std::endl;
  }
  //------------------------------------------------------------------
  {
    std::cout << "================================================"
	      << std::endl << std::endl;
    v.clear();
    v.push_back(empty);
    v.push_back(empty);
    v.push_back(empty);
    v.push_back(empty);


    {
      for (unsigned int i = 0; i < v.size(); i++) {
	std::cout << "List " << (i+1) << ": ";
	for (unsigned int j = 0; j < v[i].size(); j++) {
	  std::cout << " " << v[i][j];
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;
    }

    n_iterator ni_begin(v.begin(), v.end(), v.begin());
    n_iterator ni_end(v.begin(), v.end(), v.end());
    
    std::cout << "Nested iterator:" << std::endl;
    n_iterator ni_it;
    for (ni_it = ni_begin; ni_it != ni_end; ++ni_it) {
      std::cout << " " << (*ni_it);
    }
    std::cout << std::endl << std::endl;


    std::cout << "Nested iterator in reverse order:" << std::endl;
    if ( ni_begin != ni_end ) {
      for (ni_it = --ni_end; ni_it != ni_begin; --ni_it) {
	std::cout << " " << (*ni_it);
      }
      std::cout << " " << (*ni_it);
    }
    std::cout << std::endl << std::endl;
  }

  
  //------------------------------------------------------------------
  {
    std::cout << "================================================"
	      << std::endl << std::endl;
    v.clear();
    v1.clear();

    v1.push_back(1);

    v.push_back(empty);
    v.push_back(v1);
    v.push_back(empty);
    v.push_back(empty);


    {
      for (unsigned int i = 0; i < v.size(); i++) {
	std::cout << "List " << (i+1) << ": ";
	for (unsigned int j = 0; j < v[i].size(); j++) {
	  std::cout << " " << v[i][j];
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;
    }

    n_iterator ni_begin(v.begin(), v.end(), v.begin());
    n_iterator ni_end(v.begin(), v.end(), v.end());
    
    std::cout << "Nested iterator:" << std::endl;
    n_iterator ni_it;
    for (ni_it = ni_begin; ni_it != ni_end; ++ni_it) {
      std::cout << " " << (*ni_it);
    }
    std::cout << std::endl << std::endl;


    std::cout << "Nested iterator in reverse order:" << std::endl;
    if ( ni_begin != ni_end ) {
      for (ni_it = --ni_end; ni_it != ni_begin; --ni_it) {
	std::cout << " " << (*ni_it);
      }
      std::cout << " " << (*ni_it);
    }
    std::cout << std::endl << std::endl;
  }
 
  //------------------------------------------------------------------
  {
    std::cout << "================================================"
	      << std::endl << std::endl;
    v.clear();
    v1.clear();

    v1.push_back(1);
    v1.push_back(2);

    v.push_back(v1);
    v.push_back(empty);
    v.push_back(empty);
    v.push_back(empty);
    v.push_back(v1);


    {
      for (unsigned int i = 0; i < v.size(); i++) {
	std::cout << "List " << (i+1) << ": ";
	for (unsigned int j = 0; j < v[i].size(); j++) {
	  std::cout << " " << v[i][j];
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;
    }

    n_iterator ni_begin(v.begin(), v.end(), v.begin());
    n_iterator ni_end(v.begin(), v.end(), v.end());
    
    std::cout << "Nested iterator:" << std::endl;
    for (n_iterator ni_it = ni_begin; ni_it != ni_end; ++ni_it) {
      std::cout << " " << (*ni_it);
    }
    std::cout << std::endl << std::endl;


    std::cout << "Nested iterator in reverse order:" << std::endl;
    if ( ni_begin != ni_end ) {
      n_iterator ni_it;
      for (ni_it = --ni_end; ni_it != ni_begin; --ni_it) {
	std::cout << " " << (*ni_it);
      }
      std::cout << " " << (*ni_it);
    }
    std::cout << std::endl << std::endl;
  }
}
