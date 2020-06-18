#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Timer.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>

#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Segment_2 Segment_2;
typedef Kernel::Point_2 Point_2;
typedef CGAL::Random_points_in_square_2< Point_2 > Point_generator;
typedef CGAL::Polygon_2<Kernel>  Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel> Polygon_with_holes_2;
typedef CGAL::Polygon_set_2<Kernel> Polygon_set_2;

static int numLines;
static int numArrays;

void get_rotated_line_array( double x0,
                             double y0,
                             double angle,
                             double scale,
                             std::vector< Segment_2 > &results )
{
    double thickness = .1 / numLines;
    double cc = cos( angle );
    double ss = sin( angle );
    double xl,xr,yb,yt;
    xl = -scale * thickness;
    xr = scale * thickness;
    yb = -scale;
    yt = scale;
    Point_2 p1(( x0 + xl * cc - yb * ss ), ( y0 + xl * ss + yb * cc ));
    Point_2 p2(( x0 + xr * cc - yb * ss ), ( y0 + xr * ss + yb * cc ));
    Point_2 p3(( x0 + xr * cc - yt * ss ), ( y0 + xr * ss + yt * cc ));
    Point_2 p4(( x0 + xl * cc - yt * ss ), ( y0 + xl * ss + yt * cc ));
    results.push_back( Segment_2( p1, p2 ));
    results.push_back( Segment_2( p2, p3 ));
    results.push_back( Segment_2( p3, p4 ));
    results.push_back( Segment_2( p4, p1 ));
    for( int i = 1; i <= numLines; i++ )
    {
        double delta = (double)i / (double)numLines;
        xl = scale * ( delta - thickness );
        xr = scale * ( delta + thickness );
        Point_2 p1(( x0 + xl * cc - yb * ss ), ( y0 + xl * ss + yb * cc ));
        Point_2 p2(( x0 + xr * cc - yb * ss ), ( y0 + xr * ss + yb * cc ));
        Point_2 p3(( x0 + xr * cc - yt * ss ), ( y0 + xr * ss + yt * cc ));
        Point_2 p4(( x0 + xl * cc - yt * ss ), ( y0 + xl * ss + yt * cc ));
        results.push_back( Segment_2( p1, p2 ));
        results.push_back( Segment_2( p2, p3 ));
        results.push_back( Segment_2( p3, p4 ));
        results.push_back( Segment_2( p4, p1 ));
        xl = -scale * ( delta + thickness );
        xr = -scale * ( delta - thickness );
        Point_2 p5(( x0 + xl * cc - yb * ss ), ( y0 + xl * ss + yb * cc ));
        Point_2 p6(( x0 + xr * cc - yb * ss ), ( y0 + xr * ss + yb * cc ));
        Point_2 p7(( x0 + xr * cc - yt * ss ), ( y0 + xr * ss + yt * cc ));
        Point_2 p8(( x0 + xl * cc - yt * ss ), ( y0 + xl * ss + yt * cc ));
        results.push_back( Segment_2( p5, p6 ));
        results.push_back( Segment_2( p6, p7 ));
        results.push_back( Segment_2( p7, p8 ));
        results.push_back( Segment_2( p8, p5 ));
    }
};


void build_segments(std::vector< Segment_2 >& all_segments)
{
    double x0 = 358 * (double)rand() / RAND_MAX - 179;
    double y0 = 178 * (double)rand() / RAND_MAX - 89;

    for( int i = 0; i < numArrays; ++i )
    {
        double angle = CGAL_PI * (double)rand() / RAND_MAX;
        double scale = 1 + (double)rand() / RAND_MAX;

        get_rotated_line_array( x0, y0, angle, scale, all_segments );
    }
}

void print_polygon(std::ostream& out, const Polygon_2& polygon)
{
  out << polygon.size()+1 << " ";
  for (Polygon_2::Vertex_const_iterator vit=polygon.vertices_begin();
                                        vit!=polygon.vertices_end();++vit)
    out << *vit <<" 0 ";
  out << *polygon.vertices_begin() << " 0\n";
}

void print_polygons(std::ostream& out, const Polygon_set_2& polygon_set)
{
  std::vector<Polygon_with_holes_2> polygons_wh(polygon_set.number_of_polygons_with_holes());
  polygon_set.polygons_with_holes(&polygons_wh[0]);
  for(Polygon_with_holes_2& polygon_wh : polygons_wh)
  {
    print_polygon(out, polygon_wh.outer_boundary());
    for(Polygon_with_holes_2::Hole_const_iterator it=polygon_wh.holes_begin();
                                                  it!=polygon_wh.holes_end();++it)
    {
      print_polygon(out, *it);
    }
  }
}

int main( int  argc , char ** argv )
{
  double scale_factor = argc>1 ? atof(argv[1]):1;
  double sqrt_scale_factor = std::sqrt(scale_factor);

  //testing grid rectangle intersection
  numLines  = 100 * sqrt_scale_factor;
  numArrays = 2;

  srand( 0 );

  std::cout.precision( std::numeric_limits< double >::digits10 + 1 );

{
  std::vector<Segment_2> all_segments;
  build_segments(all_segments);

  std::ofstream out("polygons_grid.cgal");
  for(const Segment_2& s : all_segments)
  {
    out   << "2 " << s.source() << " 0"
          << "  " << s.target() << " 0\n";
  }

  CGAL::Timer time;

  std::vector<Polygon_2> polygons(all_segments.size()/4);
  for (std::size_t i=0; i< polygons.size(); ++i)
  {
    polygons[i].push_back( all_segments[ 4*i   ].source() );
    polygons[i].push_back( all_segments[ 4*i+1 ].source() );
    polygons[i].push_back( all_segments[ 4*i+2 ].source() );
    polygons[i].push_back( all_segments[ 4*i+3 ].source() );
  }
  out.close();

  time.reset(); time.start();
  Polygon_set_2 polygon_set;
  polygon_set.join(polygons.begin(), polygons.end());
  time.stop();
  std::cout << "Time for " << polygons.size() << " grid crossing polygons " << time.time() << "\n";

  out.open("polygons_grid_out.cgal");
  print_polygons(out, polygon_set);
}

{
  // testing random polygons
  int nb_poly=1000 * scale_factor;
  int ngon=55;

  std::ofstream out("polygons_random.cgal");
  std::vector<Polygon_2> polygons(nb_poly);
  for (int i=0;i<nb_poly;++i)
  {
    Polygon_2& polygon=polygons[i];
    CGAL::random_polygon_2(ngon, std::back_inserter(polygon),
                           Point_generator(0.5));
    print_polygon(out, polygon);
  }
  out.close();

  CGAL::Timer time; time.start();
  Polygon_set_2 polygon_set;
  polygon_set.join(polygons.begin(), polygons.end());
  time.stop();

  std::cout << "Time for " << polygons.size() << " random polygons " << time.time() << "s" << std::endl;
  out.open("polygons_random_out.cgal");
  print_polygons(out, polygon_set);
}

{
  // testing nested polygons
  int nb_poly=100000 * scale_factor;
  std::ofstream out("polygons_nested.cgal");
  std::vector<Polygon_2> polygons(nb_poly);
  const double epsilon = 1./(nb_poly+1);
  for (int i=0;i<nb_poly;++i)
  {
    Polygon_2& polygon=polygons[i];

    polygon.push_back( Point_2(0+epsilon*i,0+epsilon*i) );
    polygon.push_back( Point_2(0+epsilon*i,1-epsilon*i) );
    polygon.push_back( Point_2(1-epsilon*i,1-epsilon*i) );
    polygon.push_back( Point_2(1-epsilon*i,0+epsilon*i) );
    print_polygon(out, polygon);
  }
  out.close();

  CGAL::Timer time; time.start();
  Polygon_set_2 polygon_set;
  polygon_set.join(polygons.begin(), polygons.end());
  time.stop();

  std::cout << "Time for " << polygons.size() << " nested polygons " << time.time() << "s" << std::endl;

  out.open("polygons_nested_out.cgal");
  print_polygons(out, polygon_set);
}

{
  // testing disjoint polygons
  int grid_size=200 * scale_factor; // grid_size^2 polygons
  std::ofstream out("polygons_disjoint.cgal");
  std::vector<Polygon_2> polygons;
  polygons.reserve(grid_size*grid_size);

  double epsilon = 1./(2* grid_size);

  for (int i=0;i<grid_size;++i)
    for (int j=0;j<grid_size; ++j)
    {
      polygons.push_back( Polygon_2() );
      Polygon_2& polygon=polygons.back();
      polygon.push_back( Point_2(0+epsilon*2*i    ,0+epsilon*2*j) );
      polygon.push_back( Point_2(0+epsilon*(2*i+1),0+epsilon*2*j) );
      polygon.push_back( Point_2(0+epsilon*(2*i+1),0+epsilon*(2*j+1)) );
      polygon.push_back( Point_2(0+epsilon*2*i    ,0+epsilon*(2*j+1)) );
      print_polygon(out, polygon);
    }
  out.close();

  CGAL::Timer time; time.start();
  Polygon_set_2 polygon_set;
  polygon_set.join(polygons.begin(), polygons.end());
  time.stop();

  std::cout << "Time for " << polygons.size() << " disjoint polygons " << time.time() << "s" << std::endl;

  out.open("polygons_disjoint_out.cgal");
  print_polygons(out, polygon_set);
}
    return 0;
}
