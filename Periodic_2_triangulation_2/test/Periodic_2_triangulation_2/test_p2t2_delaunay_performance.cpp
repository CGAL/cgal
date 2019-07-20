// Author(s)     : Nico Kruithof  <Nico@nghk.nl>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Timer.h>
#include <CGAL/point_generators_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<K> Gt;
typedef Gt::Point_2                                         Point;
typedef Gt::Vector_2                                        Vector;

typedef CGAL::Periodic_2_Delaunay_triangulation_2<Gt>       P2DT2;
typedef CGAL::Delaunay_triangulation_2<Gt>                  DT2;

template <class Dt>
class DT2_inserter
{
  Dt t;
public:
  template <class Iterator>
  void insert(Iterator begin, Iterator end)
  {
    t.insert(begin, end);
  }
};


template <class PT, bool large>
class P2DT2_inserter
{
  PT t;
public:
  template <class Iterator>
  void insert(Iterator begin, Iterator end)
  {
    t.insert(begin, end, large);
  }
};

template <class Inserter>
void test_performance(const std::string &name, int maximum = 1e5)
{
  // Create point sets
  typedef CGAL::Creator_uniform_2<double, Point>  Creator;
  CGAL::Random rnd(7);
  CGAL::Random_points_in_square_2<Point, Creator> in_square(0.5, rnd);


  for (int n = 0; n <= maximum; n += 5000)
    {
      CGAL::Timer timer;

      std::vector<Point> pts;
      for (int i = 0 ; i < n ; i++)
        {
          pts.push_back(*in_square++ + Vector(0.5, 0.5));
        }

      Inserter inserter;

      timer.start();
      inserter.insert(pts.begin(), pts.end());
      timer.stop();
      std::cout << name << "; " << pts.size() << "; " << timer.time() << std::endl;
    }
}

int main(int argc, char *argv[])
{
  int maximum = 100000;
  if (argc > 1) maximum = atoi(argv[1]);
  test_performance<DT2_inserter<DT2> >("Euclidean Delaunay", maximum);
  test_performance<P2DT2_inserter<P2DT2, false> >("Periodic Delaunay", maximum);
  test_performance<P2DT2_inserter<P2DT2, true>  >("Periodic Delaunay, large point set", maximum);

  return 0;
}

// For generating the plot:

/*
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html>
 <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
    <title>P2T2 performance</title>
    <script language="javascript" type="text/javascript" src="js/jquery.js"></script>
    <script language="javascript" type="text/javascript" src="js/jquery.flot.js"></script>
    <script language="javascript" type="text/javascript" src="js/jquery.csv-0.71.js"></script>
 </head>
    <body>
    <div id="placeholder" style="width:600px;height:300px;"></div>

    <script type="text/javascript">
      function createPlot(csv_data) {
        var data = $.csv.toArrays(csv_data, {separator:";"});
        var processed_data = {};
        data.map(function (elem) {
          if (!(elem[0] in processed_data)) processed_data[elem[0]] = []
          return processed_data[elem[0]].push([ elem[1], elem[2] ]);
        });

        var processed_data2 = {};
        for (var key in processed_data) {
          var data0 = processed_data["Euclidean Delaunay"];
          var data1 = processed_data[key];
          var data = []
          for (var i=0; i<.25*data1.length; ++i) {
            data.push([data0[i][0], data1[i][1]/data0[i][1]]);
          }
          processed_data2[key] = data;
        }

        var plot_data = [];
        for (var key in processed_data) {
          plot_data.push({
            label : key,
            data  : processed_data[key],
            lines : { lineWidth: 1 },
            shadowSize: 0
          });
        }

        $.plot($("#placeholder"), plot_data, { legend: { position: "nw" } } );
      }

      $.get('performance.txt', createPlot);

    </script>
    </body>
</html>
*/
