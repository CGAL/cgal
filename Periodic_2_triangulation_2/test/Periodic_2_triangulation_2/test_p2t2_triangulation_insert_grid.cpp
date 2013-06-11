// Author(s)     : Nico Kruithof  <Nico@nghk.nl>

#include "./types.h"
#include <map>

const int N = 10;

void test_insertion_xy(int x_order, int y_order)
{
  std::cout << "xy: " << x_order << " " << y_order << std::endl;
  Triangulation t;

  Triangulation::Face_handle fh;

  // Insert the first point
  for (int x = -N; x < N; x++)
    {
      for (int y = -N; y < N; y++)
        {
          Point p(0.5 + x_order * (0.4999 / N)*x, 0.5 + y_order * (0.4999 / N)*y);
          t.insert(p);
        }
    }
  CGAL_assertion(t.is_valid());
}

void test_insertion_yx(int x_order, int y_order)
{
  std::cout << "yx: " << x_order << " " << y_order << std::endl;
  Triangulation t;

  // Insert the first point
  for (int y = -N; y < N; y++)
    {
      for (int x = -N; x < N; x++)
        {
          Point p(0.5 + x_order * (0.4999 / N)*x, 0.5 + y_order * (0.4999 / N)*y);
          t.insert(p);
        }
    }
  CGAL_assertion(t.is_valid());
}

int main()
{
  test_insertion_xy( 1,  1);
  test_insertion_xy( 1, -1);
  test_insertion_xy(-1,  1);
  test_insertion_xy(-1, -1);

  test_insertion_yx( 1,  1);
  test_insertion_yx( 1, -1);
  test_insertion_yx(-1,  1);
  test_insertion_yx(-1, -1);

  return 0;
}
