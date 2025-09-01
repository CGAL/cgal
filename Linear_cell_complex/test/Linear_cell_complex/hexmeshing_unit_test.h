#include <CGAL/Hexmeshing_for_linear_cell_complex_sequential.h>
#include <iostream>

class TestFramework {
private:
  int total_tests = 0;
  int passed_tests = 0;
  int failed_tests = 0;

public:
  void test(const std::string test_name, std::function<bool()> test_func) {
    total_tests++;
    std::cout << "Running test: " << test_name << " ... ";

    try {
      bool result = test_func();
      if (result) {
        std::cout << "PASSED" << std::endl;
        passed_tests++;
      } else {
        std::cout << "FAILED" << std::endl;
        failed_tests++;
      }
    } catch (const std::exception& e) {
      std::cout << "FAILED (Exception: " << e.what() << ")" << std::endl;
      failed_tests++;
    } catch (...) {
      std::cout << "FAILED (Unknown exception)" << std::endl;
      failed_tests++;
    }
  }

  template<typename T>
  void test_exact_value(const std::string test_name, const T expected_value, const std::function<T()> test_func) {
    test(test_name, [expected_value, test_func]() {
      T actual_value = test_func();
      assert(actual_value == expected_value);
      return true;
    });
  }

  template<typename T>
  void test_approximate_value(const std::string test_name, const T expected_value, const std::function<T()> test_func, const std::function<bool(T, T)> is_approximately_equal) {
    test(test_name, [expected_value, test_func, is_approximately_equal]() {
      T actual_value = test_func();
      return is_approximately_equal(actual_value, expected_value);
    });
  }

  void print_summary() {
    std::cout << "\n=== Test Summary ===" << std::endl;
    std::cout << "Total tests: " << total_tests << std::endl;
    std::cout << "Passed: " << passed_tests << std::endl;
    std::cout << "Failed: " << failed_tests << std::endl;
    std::cout << "Success rate: " << (total_tests > 0 ? (passed_tests * 100.0 / total_tests) : 0) << "%" << std::endl;
  }
};


using namespace CGAL::internal::Hexmeshing;

bool is_approximately_equal_double(const double a, const double b, const double tolerance = 1e-10) {
  return std::abs(a - b) < tolerance;
}

bool is_approximately_equal_Point(const Point a, const Point b, const double tolerance = 1e-10) {
  return std::abs(a.x() - b.x()) < tolerance &&
         std::abs(a.y() - b.y()) < tolerance &&
         std::abs(a.z() - b.z()) < tolerance;
}

bool is_approximately_equal_Vector(const Vector a, const Vector b, const double tolerance = 1e-10) {
  return std::abs(a.x() - b.x()) < tolerance &&
         std::abs(a.y() - b.y()) < tolerance &&
         std::abs(a.z() - b.z()) < tolerance;
}

bool is_approximately_equal_vector_Point(const std::vector<Point> a, const std::vector<Point> b, const double tolerance = 1e-10) {
  if (a.size() != b.size()) {
    return false;
  }
  for (size_t i = 0; i < a.size(); i++) {
    // std::cout << a[i] << ' ' << b[i] << std::endl;
    if (!is_approximately_equal_Point(a[i], b[i], tolerance)) {
      return false;
    }
  }
  return true;
}

bool is_approximately_equal_Points(std::vector<Point> a, std::vector<Point> b, const double tolerance = 1e-10) {
  auto cmp = [tolerance](Point a, Point b) {
    if(is_approximately_equal_double(a.x(), b.x(), tolerance)) {
      if(is_approximately_equal_double(a.y(), b.y(), tolerance)) {
        return a.z() > b.z();
      }
      return a.y() > b.y();
    }
    return a.x() > b.x();
  };
  std::sort(a.begin(), a.end(), cmp);
  std::sort(b.begin(), b.end(), cmp);

  return is_approximately_equal_vector_Point(a, b, tolerance);
}

bool is_approximately_equal_Segment(const Segment a, const Segment b, const double tolerance = 1e-10) {
  Point a1 = a.source(), a2 = a.target();
  Point b1 = b.source(), b2 = b.target();
  return (is_approximately_equal_Point(a1, b1, tolerance) and is_approximately_equal_Point(a2, b2, tolerance))
          or (is_approximately_equal_Point(a1, b2, tolerance) and is_approximately_equal_Point(a2, b1, tolerance));
}

LCC create_grid_mesh(int nx, int ny, int nz) {
  Grid grid(Point(0, 0, 0), {1., 1., 1.}, {nx, ny, nz});

  LCC lcc;
  Grid::generate_grid(lcc, grid);

  return lcc;
}

template<int FACE=1>
LCC create_grid_mesh_with_volume_fraction() {
  LCC lcc = create_grid_mesh(11, 11, 11);

  set_centroids(lcc);
  auto volumes = lcc.one_dart_per_cell<3>();
  for(auto it = volumes.begin(); it != volumes.end(); it++) {
    auto &vol_attr = lcc.attribute<3>(it)->info();
    Point c = vol_attr.centroid;
    double f = 1.;
    std::array<double, 3> xyz = {c.x(), c.y(), c.z()};
    if constexpr(FACE == 6) {
      for(auto x:xyz) {
        if(2 < x and x < 9) f *= 1.;
        else if(x < 1 or 10 < x) f *= 0.;
        else f *= 0.7;
      }
    }
    else if constexpr(FACE == 1) {
      const double maku = 5.5;
      double s = 0;
      for(auto x:xyz) {
        s += int(x);
      }
      if(maku < s) f = 0.;
      else if(maku < s+1) {
        double t = maku-s;
        f = t*t*t/6.;
      }
      else if(maku < s+2) {
        double t = maku-s, t1 = t-1.;
        f = t*t*t/6. - t1*t1*t1/2.;
      }
      else if(maku < s+3) {
        double t = s+3-maku;
        f = 1. - t*t*t/6.;
      }
    }
    vol_attr.fraction = f;
  }
  return lcc;
}

// set x+y+z < c_xyz and x+y+z > d_xyz by fraction
void set_plane(LCC &lcc, double c_xyz, double d_xyz, size_type inner_mark) {
  set_centroids(lcc);

  auto volumes = lcc.one_dart_per_cell<3>();

  for(auto volume = volumes.begin(); volume != volumes.end(); volume++) {
    Point p = lcc.attribute<3>(volume)->info().centroid;
    double &frac = lcc.attribute<3>(volume)->info().fraction = 0.;
    double x = p.x(), y = p.y(), z = p.z();
    double sum_xyz = x+y+z;
    if(sum_xyz < c_xyz-1.5) {
      frac = 1.;
    }
    else if(sum_xyz < c_xyz-0.5) {
      double l = int(c_xyz)+1. - c_xyz;
      frac += 1. - l*l*l/6.;
    }
    else if(sum_xyz < c_xyz+0.5) {
      double l = c_xyz - int(c_xyz);
      frac += (l+1)*(l+1)*(l+1)/6. - l*l*l/2.;
    }
    else if(sum_xyz < c_xyz+1.5) {
      double l = c_xyz - int(c_xyz);
      frac += l*l*l/6.;
    }
    if (sum_xyz > d_xyz+1.5) {
      frac = 1.;
    }
    else if(sum_xyz > d_xyz+0.5) {
      double l = d_xyz - int(d_xyz);
      frac += 1. - l*l*l/6.;
    }
    else if(sum_xyz > d_xyz-0.5) {
      double l = int(d_xyz)+1 - d_xyz;
      frac += (l+1)*(l+1)*(l+1)/6. - l*l*l/2.;
    }
    else if(sum_xyz > d_xyz-1.5) {
      double l = int(d_xyz)+1 - d_xyz;
      frac += l*l*l/6.;
    }

    if(frac > 0.5) {
      lcc.mark_cell<3>(volume, inner_mark);
    }
  }
}