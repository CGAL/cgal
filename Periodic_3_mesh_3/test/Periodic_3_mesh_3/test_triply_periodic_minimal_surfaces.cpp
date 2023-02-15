#include <CGAL/Periodic_3_mesh_3/config.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/make_periodic_3_mesh_3.h>
#include <CGAL/Periodic_3_mesh_3/IO/File_medit.h>
#include <CGAL/Periodic_3_mesh_triangulation_3.h>
#include <CGAL/Periodic_3_function_wrapper.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Implicit_to_labeling_function_wrapper.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/functional.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::FT                                               FT;
typedef K::Point_3                                          Point;
typedef K::Segment_3                                        Segment;
typedef K::Iso_cuboid_3                                     Iso_cuboid;

// Domain
typedef FT (Function)(const Point&);
typedef CGAL::Periodic_3_function_wrapper<Function, K>      Periodic_function;
typedef CGAL::Implicit_multi_domain_to_labeling_function_wrapper<Periodic_function>
                                                            Labeling_function;
typedef CGAL::Labeled_mesh_domain_3<K>                      Periodic_mesh_domain;

// Triangulation
typedef CGAL::Periodic_3_mesh_triangulation_3<Periodic_mesh_domain>::type  Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr>                        C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr>                           Periodic_mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

const FT PI = CGAL_PI;

FT schwarz_p(const Point& p)
{
  assert(p.x() >= 1 && p.y() >= 1 && p.z() >= 1 && p.x() < 3 && p.y() < 3 && p.z() < 3);
  const FT x2 = std::cos( p.x() * 2*PI ),
           y2 = std::cos( p.y() * 2*PI ),
           z2 = std::cos( p.z() * 2*PI );
  return x2 + y2 + z2;
}

FT schwarz_p_transl (const Point& p)
{
  assert(p.x() >= 1 && p.y() >= 1 && p.z() >= 1 && p.x() < 3 && p.y() < 3 && p.z() < 3);
  const FT x2 = std::cos(p.x() * 2 * PI + PI / 2.0),
           y2 = std::cos(p.y() * 2 * PI + PI / 2.0),
           z2 = std::cos(p.z() * 2 * PI + PI / 2.0);
  return x2 + y2 + z2;
}

FT gyroid (const Point& p)
{
  assert(p.x() >= 1 && p.y() >= 1 && p.z() >= 1 && p.x() < 3 && p.y() < 3 && p.z() < 3);
  const FT cx = std::cos(p.x() * 2 * PI),
           cy = std::cos(p.y() * 2 * PI),
           cz = std::cos(p.z() * 2 * PI);
  const FT sx = std::sin(p.x() * 2 * PI),
           sy = std::sin(p.y() * 2 * PI),
           sz = std::sin(p.z() * 2 * PI);
  return cx * sy + cy * sz + cz * sx;
}

FT diamond (const Point& p)
{
  assert(p.x() >= 1 && p.y() >= 1 && p.z() >= 1 && p.x() < 3 && p.y() < 3 && p.z() < 3);
  const FT cx = std::cos(p.x() * 2 * PI),
           cy = std::cos(p.y() * 2 * PI),
           cz = std::cos(p.z() * 2 * PI);
  const FT sx = std::sin(p.x() * 2 * PI),
           sy = std::sin(p.y() * 2 * PI),
           sz = std::sin(p.z() * 2 * PI);
  return sx * sy * sz + sx * cy * cz + cx * sy * cz + cx * cy * sz;
}

FT double_p (const Point& p)
{
  assert(p.x() >= 1 && p.y() >= 1 && p.z() >= 1 && p.x() < 3 && p.y() < 3 && p.z() < 3);
  const FT cx = std::cos(p.x() * 2 * PI),
           cy = std::cos(p.y() * 2 * PI),
           cz = std::cos(p.z() * 2 * PI);
  const FT c2x = std::cos(2 * p.x() * 2 * PI),
           c2y = std::cos(2 * p.y() * 2 * PI),
           c2z = std::cos(2 * p.z() * 2 * PI);
  return 0.5 * (cx * cy  + cy * cz + cz * cx ) + 0.2*(c2x + c2y + c2z);
}

FT G_prime (const Point& p)
{
  assert(p.x() >= 1 && p.y() >= 1 && p.z() >= 1 && p.x() < 3 && p.y() < 3 && p.z() < 3);
  const FT cx = std::cos(p.x() * 2 * PI),
           cy = std::cos(p.y() * 2 * PI),
           cz = std::cos(p.z() * 2 * PI);
  const FT c2x = std::cos(2 * p.x() * 2 * PI),
           c2y = std::cos(2 * p.y() * 2 * PI),
           c2z = std::cos(2 * p.z() * 2 * PI);
  const FT sx = std::sin(p.x() * 2 * PI),
           sy = std::sin(p.y() * 2 * PI),
           sz = std::sin(p.z() * 2 * PI);
  const FT s2x = std::sin(2 * p.x() * 2 * PI),
           s2y = std::sin(2 * p.y() * 2 * PI),
           s2z = std::sin(2 * p.z() * 2 * PI);
  return 5 * (s2x * sz * cy  + s2y * sx * cz  + s2z * sy * cx)
           + 1*(c2x * c2y  + c2y * c2z  + c2z * c2x);
}

FT lidinoid (const Point& p)
{
  assert(p.x() >= 1 && p.y() >= 1 && p.z() >= 1 && p.x() < 3 && p.y() < 3 && p.z() < 3);
  const FT cx = std::cos(p.x() * 2 * PI),
           cy = std::cos(p.y() * 2 * PI),
           cz = std::cos(p.z() * 2 * PI);
  const FT c2x = std::cos(2 * p.x() * 2 * PI),
           c2y = std::cos(2 * p.y() * 2 * PI),
           c2z = std::cos(2 * p.z() * 2 * PI);
  const FT sx = std::sin(p.x() * 2 * PI),
           sy = std::sin(p.y() * 2 * PI),
           sz = std::sin(p.z() * 2 * PI);
  const FT s2x = std::sin(2 * p.x() * 2 * PI),
           s2y = std::sin(2 * p.y() * 2 * PI),
           s2z = std::sin(2 * p.z() * 2 * PI);
  return 1 * (s2x * sz * cy  + s2y * sx * cz  + s2z * sy * cx)
           - 1 * (c2x * c2y  + c2y * c2z  + c2z * c2x) + 0.3;
}

FT D_prime (const Point& p)
{
  assert(p.x() >= 1 && p.y() >= 1 && p.z() >= 1 && p.x() < 3 && p.y() < 3 && p.z() < 3);
  const FT cx = std::cos(p.x() * 2 * PI),
           cy = std::cos(p.y() * 2 * PI),
           cz = std::cos(p.z() * 2 * PI);
  const FT c2x = std::cos(2 * p.x() * 2 * PI),
           c2y = std::cos(2 * p.y() * 2 * PI),
           c2z = std::cos(2 * p.z() * 2 * PI);
  const FT sx = std::sin(p.x() * 2 * PI),
           sy = std::sin(p.y() * 2 * PI),
           sz = std::sin(p.z() * 2 * PI);
  return 1 * ( sx * sy * sz) + 1 * ( cx * cy * cz)
           - 1 * ( c2x * c2y  + c2y * c2z  + c2z * c2x) - 0.4;
}

FT split_p (const Point& p)
{
  assert(p.x() >= 1 && p.y() >= 1 && p.z() >= 1 && p.x() < 3 && p.y() < 3 && p.z() < 3);
  const FT cx = std::cos(p.x() * 2 * PI),
           cy = std::cos(p.y() * 2 * PI),
           cz = std::cos(p.z() * 2 * PI);
  const FT c2x = std::cos(2 * p.x() * 2 * PI),
           c2y = std::cos(2 * p.y() * 2 * PI),
           c2z = std::cos(2 * p.z() * 2 * PI);
  const FT sx = std::sin(p.x() * 2 * PI),
           sy = std::sin(p.y() * 2 * PI),
           sz = std::sin(p.z() * 2 * PI);
  const FT s2x = std::sin(2 * p.x() * 2 * PI),
           s2y = std::sin(2 * p.y() * 2 * PI),
           s2z = std::sin(2 * p.z() * 2 * PI);
  return  1.1 * (s2x * sz * cy  + s2y * sx * cz  + s2z * sy * cx)
              - 0.2 * (c2x * c2y  + c2y * c2z  + c2z * c2x)
              - 0.4 * (cx + cy + cz);
}

struct Segments_function
  : public CGAL::cpp98::unary_function<Point, FT>
{
  typedef std::vector<Segment> Segments;

  Segments_function()
    : segments(), nb_evals(0)
  {
    const FT vmin = 1, vmax = 3;
    const FT mid = 0.5 * (vmin + vmax);
    const Point pmid(mid, mid, mid);

    segments.push_back(Segment(Point(vmin, mid, vmin), pmid));
    segments.push_back(Segment(Point(vmax, mid, vmin), pmid));
    segments.push_back(Segment(Point(vmin, mid, vmax), pmid));
    segments.push_back(Segment(Point(vmax, mid, vmax), pmid));
    segments.push_back(Segment(Point(mid, vmin, vmin), pmid));
    segments.push_back(Segment(Point(mid, vmax, vmin), pmid));
    segments.push_back(Segment(Point(mid, vmin, vmax), pmid));
    segments.push_back(Segment(Point(mid, vmax, vmax), pmid));
  }

  FT operator()(const Point& p)
  {
    ++nb_evals;

    FT min_distance = 1000000;
    for (Segments::const_iterator si = segments.begin(); si != segments.end(); ++si)
      min_distance = (std::min)(CGAL::squared_distance(p, *si), min_distance);

    return min_distance - 0.01; // Change the squared beam radius here
  }

  Segments segments;
  int nb_evals;
};

FT segments(const Point& p) {
  assert(p.x() >= 1 && p.y() >= 1 && p.z() >= 1 && p.x() < 3 && p.y() < 3 && p.z() < 3);
  static Segments_function instance;
  return instance(p);
}

std::vector<std::string> make_vps_in ()
{
  return std::vector<std::string>(1, "-");
}

std::vector<std::string> make_vps_in_out ()
{
  std::vector<std::string> vps;
  vps.push_back("-");
  vps.push_back("+");
  return vps;
}

C3t3 make_mesh(const Labeling_function& labeling_function, const Iso_cuboid& canonical_cube)
{
  Periodic_mesh_domain domain(labeling_function, canonical_cube);

  Periodic_mesh_criteria criteria(facet_angle = 30.,
                                  facet_size = 0.03 * 2 /*domain's edge length*/,
                                  facet_distance = 0.03 * 2 /*domain's edge length*/,
                                  cell_radius_edge_ratio = 2.,
                                  cell_size = 0.05);

  return CGAL::make_periodic_3_mesh_3<C3t3>(domain, criteria);
}

typedef std::vector<std::string> Position_vector;

int main(int, char**)
{
  Iso_cuboid canonical_cube(1, 1, 1, 3, 3, 3);

  std::map<std::string, Periodic_function> functions;
#ifdef CGAL_NDEBUG
  // Only test those when not in debug (otherwise it takes too long)
  functions["D_prime"] = Periodic_function(D_prime, canonical_cube);
  functions["G_prime"] = Periodic_function(G_prime, canonical_cube);
  functions["diamond"] = Periodic_function(diamond, canonical_cube);
  functions["double_p"] = Periodic_function(double_p, canonical_cube);
  functions["gyroid"] = Periodic_function(gyroid, canonical_cube);
  functions["lidinoid"] = Periodic_function(lidinoid, canonical_cube);
#endif
  functions["schwarz_p"] = Periodic_function(schwarz_p, canonical_cube);
  functions["schwarz_p_transl"] = Periodic_function(schwarz_p_transl, canonical_cube);
  functions["segments"] = Periodic_function(segments, canonical_cube);
  functions["split_p"] = Periodic_function(split_p, canonical_cube);

  std::map<std::string, Position_vector> v_vps;
  v_vps["in"] = make_vps_in();
  v_vps["in-out"] = make_vps_in_out();

  std::vector<unsigned> v_ncopy;
  v_ncopy.push_back(1);
//  v_ncopy.push_back(2); // if you wish to print more copies
//  v_ncopy.push_back(4);
//  v_ncopy.push_back(8);

  for (std::map<std::string, Periodic_function>::iterator iter = functions.begin();
       iter != functions.end(); ++iter)
  {
    for (std::map<std::string, Position_vector>::iterator it = v_vps.begin();
         it != v_vps.end(); ++it)
    {
      std::ostringstream oss;
      oss << iter->first << "__" << it->first;
      std::string mesh_id = oss.str();

      std::cout << "Meshing " << mesh_id << "..." << std::flush;

      std::vector<Periodic_function> funcs;
      funcs.push_back(iter->second);

      Labeling_function labeling_function(funcs, it->second);
      C3t3 c3t3 = make_mesh(labeling_function, canonical_cube);

      assert(c3t3.is_valid());
      assert(c3t3.triangulation().is_valid());

      std::cout << " Done" << std::flush;

      for (std::vector<unsigned>::iterator i = v_ncopy.begin(); i != v_ncopy.end(); ++i)
      {
        std::ostringstream oss_2;
        oss_2 << iter->first << "__" << it->first << "__" << *i  << ".mesh";
        std::string output_filename = oss_2.str();
        std::ofstream medit_file( output_filename.data() );
        CGAL::IO::output_periodic_mesh_to_medit(medit_file, c3t3, *i);
        medit_file.close();

        std::cout << ", " << *i << "-copy Saved" << std::flush;
      }
      std::cout << std::endl;
    }
  }

  std::cout << "EXIT SUCCESS" << std::endl;
  return 0;
}
