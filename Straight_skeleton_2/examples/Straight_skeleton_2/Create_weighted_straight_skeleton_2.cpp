// @todo convert taper angle to edge weight
// @todo convert height to offset
// @todo exterior skeleton
// @todo filtering optimization with weights
// @todo extracting tops with holes

#if 1

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

bool lAppToLog = false ;

void Straight_skeleton_external_trace ( std::string m )
{
  std::ofstream out("/home/mrouxell/sls_log.txt", ( lAppToLog ? std::ios::app | std::ios::ate : std::ios::trunc | std::ios::ate ) );
  out << std::setprecision(19) << m << std::endl << std::flush ;
  lAppToLog = true ;
}
void Straight_skeleton_traits_external_trace ( std::string m )
{
  std::ofstream out("/home/mrouxell/sls_log.txt", ( lAppToLog ? std::ios::app | std::ios::ate : std::ios::trunc | std::ios::ate ) ) ;
  out << std::setprecision(19) << m << std::endl << std::flush ;
  lAppToLog = true ;
}

void error_handler ( char const* what, char const* expr, char const* file, int line, char const* msg )
{
  std::cerr << "CGAL error: " << what << " violation!" << std::endl
       << "Expr: " << expr << std::endl
       << "File: " << file << std::endl
       << "Line: " << line << std::endl;
  if ( msg != nullptr)
      std::cerr << "Explanation:" << msg << std::endl;

  std::exit(1);
}

#define CGAL_STRAIGHT_SKELETON_ENABLE_TRACE 4
#define CGAL_STRAIGHT_SKELETON_TRAITS_ENABLE_TRACE
#define CGAL_POLYGON_OFFSET_ENABLE_TRACE 4
#define CGAL_STSKEL_BUILDER_TRACE 4

#endif // if 0|1

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <CGAL/create_weighted_straight_skeleton_2.h>
#include <CGAL/create_weighted_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/create_offset_polygons_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/repair.h> // for remove_isolated_vertices()
#include <CGAL/Random.h>
#include <CGAL/Real_timer.h>

#include <CGAL/draw_straight_skeleton_2.h>
#include <CGAL/draw_polygon_2.h>
#include "print.h" // @fixme should be included first

#include <boost/container/flat_map.hpp>
#include <boost/shared_ptr.hpp>

#include <deque>
#include <iostream>
#include <unordered_map>

// @todo For robustness, use SLS predicates and comparisons:
//       - Compare_offset_against_event_time
//       - Construct_offset_point_2
// @todo: limit the SS construction using the offset value


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::FT                                               FT;
typedef K::Point_3                                          Point_3;

typedef K::Point_2                                          Point_2;
typedef CGAL::Polygon_2<K>                                  Polygon_2;
typedef CGAL::Polygon_with_holes_2<K>                       Polygon_with_holes_2;

typedef CGAL::Straight_skeleton_2<K>                        Ss;
typedef boost::shared_ptr<Ss>                               SsPtr;

typedef Ss::Base                                            HDS;

// @fixme for partial skeletons, the geometry is currently wrong because we interpolate
// between the time at the polygon vertex and a vertex at infinity (usually something like 10^308)
// so it creates pretty much vertical segments
const bool only_use_full_skeletons = true;
const FT def_offset = std::numeric_limits<double>::max();

bool read_input_polygon(const char* filename,
                        Polygon_with_holes_2& p)
{
  std::ifstream in(filename);
  if(!in)
  {
    std::cerr << "Error: Count not read " << filename << std::endl;
    return false;
  }

  bool is_number_of_CC_in_input = false;
  if(CGAL::IO::internal::get_file_extension(filename) == "poly")
  {
    is_number_of_CC_in_input = true;
  }

  std::vector<Polygon_2> polys;

  auto read_polygon = [&in, &polys](int i) -> void
  {
    std::vector<Point_2> poly;

    int v_count = 0;
    in >> v_count;
    for(int j=0; j<v_count && in; ++j)
    {
      double x = 0., y = 0.;
      in >> x >> y;
      poly.push_back(Point_2(x, y));
    }

    if(poly.size() >= 3)
    {
      bool is_simple = CGAL::is_simple_2(poly.begin(), poly.end(), K());
      if(!is_simple)
        std::cerr << "Input polygon not simple (hopefully it is strictly simple...)" << std::endl;

      CGAL::Orientation expected = (i == 0 ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE);

      const double area = CGAL::to_double(CGAL::polygon_area_2(poly.begin(), poly.end(), K()));
      CGAL::Orientation orientation = area > 0 ? CGAL::COUNTERCLOCKWISE : area < 0 ? CGAL::CLOCKWISE : CGAL::COLLINEAR;

      if(orientation == expected)
        polys.emplace_back(poly.begin(), poly.end());
      else
        polys.emplace_back(poly.rbegin(), poly.rend());
    }
  };

  if(is_number_of_CC_in_input)
  {
    int ccb_count = 0;
    in >> ccb_count;
    for(int i=0; i<ccb_count && in; ++i)
      read_polygon(i);
  }
  else
  {
    int i = 0;
    while(in)
      read_polygon(i++);
  }

  if(polys.empty())
  {
    std::cerr << "Error: empty input?" << std::endl;
    return false;
  }

  std::cout <<"Polygon with border of size: " << polys[0].size() << std::endl;
  if(polys.size() > 1)
    std::cout << polys.size() - 1 << " hole(s)" << std::endl;

  p = Polygon_with_holes_2(polys[0]);
  for(std::size_t i=0; i<polys.size()-1; ++i)
    p.add_hole(polys[i+1]);

  std::cout << p.outer_boundary().size() << " outer vertices" << std::endl;

  return true;
}

bool read_input_weights(const char* filename,
                        std::vector<std::vector<FT> >& weights)
{
  std::ifstream in(filename);
  if(!in)
  {
    std::cerr << "Error: Count not read " << filename << std::endl;
    return false;
  }

  std::vector<FT> border_weights;

  std::string line;
  while(getline(in, line))
  {
    if(line.empty())
    {
      weights.push_back(border_weights);
      border_weights.clear();
    }

    std::istringstream iss(line);
    double w;
    if(iss >> w)
      border_weights.push_back(w);
  }

  // in case the last line is not empty
  if(!border_weights.empty())
    weights.push_back(border_weights);

  return true;
}

// @todo for speed, one could abuse the fact that input vertices are first in the HDS
// @fixme bottom with holes not supported yet
template <typename Graph>
void add_bottom_face(const Polygon_with_holes_2& p,
                     Graph& g)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor       vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::face_descriptor         face_descriptor;

  std::unordered_map<Point_2, vertex_descriptor> p2v;
  for(vertex_descriptor v : vertices(g))
    p2v[v->point()] = v;

  std::deque<vertex_descriptor> vs;
  for(auto pit = CGAL::CGAL_SS_i::vertices_begin(p.outer_boundary());
            pit != CGAL::CGAL_SS_i::vertices_end(p.outer_boundary()); ++pit)
  {
    vs.push_front(p2v.at(*pit)); // flip the face so that it is facing outside
  }

  face_descriptor f = CGAL::Euler::add_face(vs, g);
  if(f == boost::graph_traits<HDS>::null_face())
    std::cerr << "Error: Failed to add bottom face?" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

struct HDS_to_3D_VPM
{
  typedef HDS_to_3D_VPM                                     Self;

  typedef boost::graph_traits<HDS>::vertex_descriptor       vertex_descriptor;
  typedef vertex_descriptor                                 key_type;
  typedef Point_3                                           value_type;
  typedef Point_3                                           reference;
  typedef boost::readable_property_map_tag                  category;

  friend value_type get(const Self&, const key_type v)
  {
    const auto& p = v->point();
    return {p.x(), p.y(), v->time()};
  }
};

void draw_3D_mesh(const Polygon_with_holes_2& p,
                  const Ss& ss,
                  const bool add_bot_face)
{
  // Convert to a 3D mesh, with time giving altitude
  HDS& hds = const_cast<HDS&>(static_cast<const HDS&>(ss));

  if(add_bot_face)
    add_bottom_face(p, hds);

  std::cout << "write" << std::endl;
  CGAL::IO::write_OFF("sm_3D.off", hds,
                      CGAL::parameters::vertex_point_map(HDS_to_3D_VPM())
                                       .stream_precision(17));
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename Mesh>
void add_horizontal_faces(Mesh& sm)
{
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor      vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor    halfedge_descriptor;
  typedef typename boost::graph_traits<Mesh>::face_descriptor        face_descriptor;

  std::vector<halfedge_descriptor> cycle_reps;
  CGAL::Polygon_mesh_processing::extract_boundary_cycles(sm, std::back_inserter(cycle_reps));

  for(halfedge_descriptor bh : cycle_reps)
  {
    std::vector<vertex_descriptor> vs;
    for(halfedge_descriptor h : halfedges_around_face(bh, sm))
      vs.push_back(target(h, sm));

    face_descriptor f = CGAL::Euler::add_face(vs, sm);
    if(f == boost::graph_traits<Mesh>::null_face())
      std::cerr << "Error while adding bot/top faces" << std::endl;
  }
}

void draw_3D_mesh(const Polygon_with_holes_2& p,
                  const Ss& ss,
                  const FT offset_time,
                  const bool add_bot_face)
{
  typedef CGAL::Surface_mesh<Point_3>                       Mesh;

  typedef boost::graph_traits<Mesh>::vertex_descriptor      SM_vertex_descriptor;
  typedef boost::graph_traits<Mesh>::face_descriptor        SM_face_descriptor;

  typedef boost::graph_traits<HDS>::vertex_descriptor       HDS_vertex_descriptor;
  typedef boost::graph_traits<HDS>::halfedge_descriptor     HDS_halfedge_descriptor;
  typedef boost::graph_traits<HDS>::face_descriptor         HDS_face_descriptor;

  const bool with_offset = (offset_time != std::numeric_limits<double>::max());

  const HDS& hds = static_cast<const HDS&>(ss);

  Mesh sm;
  auto sm_vpm = get(CGAL::vertex_point, sm);

  std::unordered_map<HDS_vertex_descriptor, SM_vertex_descriptor> v2v;
  std::unordered_map<HDS_halfedge_descriptor, SM_vertex_descriptor> offset_points;

  for(const HDS_vertex_descriptor hds_v : vertices(hds))
  {
    auto& hds_p = hds_v->point();

    SM_vertex_descriptor sm_v = add_vertex(sm);
    put(sm_vpm, sm_v, Point_3(hds_p.x(), hds_p.y(), hds_v->time()));

    v2v[hds_v] = sm_v;
  }

  for(const HDS_face_descriptor hds_f : faces(hds))
  {
//    std::cout << "\ncheck face" << std::endl;

    // Loop around, find if this face is intersected by the offset
    std::vector<SM_vertex_descriptor> sm_vs;

    HDS_halfedge_descriptor hds_h = halfedge(hds_f, hds), done = hds_h;
    do
    {
      HDS_vertex_descriptor hds_sv = source(hds_h, hds);
      HDS_vertex_descriptor hds_tv = target(hds_h, hds);

//      std::cout << "check halfedge\n"
//                << "\t" << hds_sv->point() << " t: " << hds_sv->time() << "\n"
//                << "\t" << hds_tv->point() << " t: " << hds_tv->time() << std::endl;

      if(!with_offset)
      {
        sm_vs.push_back(v2v.at(hds_tv)); // @speed 'at' -> []
      }
      else
      {
        const CGAL::Comparison_result sc = CGAL::compare(hds_sv->time(), offset_time);
        const CGAL::Comparison_result tc = CGAL::compare(hds_tv->time(), offset_time);

//        std::cout << "offset_time = " << offset_time << std::endl;
//        std::cout << "sc/tc " << sc << " " << tc << std::endl;

        // if the offset is crossing at the source, it will be added when seen as a target
        // from the previous halfedge

        if(sc != tc && sc != CGAL::EQUAL && tc != CGAL::EQUAL)
        {
          CGAL_assertion(sc != CGAL::EQUAL && tc != CGAL::EQUAL);
//          std::cout << "sc != tc" << std::endl;

          HDS_halfedge_descriptor hds_off_h = hds_h;
          if(hds_sv->time() > hds_tv->time()) // ensure a canonical representation
            hds_off_h = opposite(hds_off_h, hds);

          SM_vertex_descriptor null_sm_v = boost::graph_traits<Mesh>::null_vertex();
          auto insert_res = offset_points.emplace(hds_off_h, null_sm_v);
          if(insert_res.second) // never computed that offset point before
          {
            HDS_vertex_descriptor hds_off_sv = source(hds_off_h, hds);
            HDS_vertex_descriptor hds_off_tv = target(hds_off_h, hds);
            CGAL_assertion(hds_off_tv->time() > hds_off_sv->time());

            const FT lambda = (offset_time - hds_off_sv->time()) / (hds_off_tv->time() - hds_off_sv->time());
            const Point_2 off_p = hds_off_sv->point() + lambda * (hds_off_tv->point() - hds_off_sv->point());
//            std::cout << "Add offset point at " << off_p << std::endl;

            SM_vertex_descriptor sm_off_v = add_vertex(sm);
            insert_res.first->second = sm_off_v;
            put(sm_vpm, sm_off_v, Point_3(off_p.x(), off_p.y(), offset_time));
          }

          sm_vs.push_back(insert_res.first->second);
        }

        if(tc != CGAL::LARGER)
          sm_vs.push_back(v2v.at(hds_tv)); // @speed 'at' -> []
      }

      hds_h = next(hds_h, hds);
    }
    while(hds_h != done);

// --
    // std::cout << "add face of size " << sm_vs.size() << std::endl;
    // for(SM_vertex_descriptor sm_v : sm_vs)
    //   std::cout << sm_v << " ";
    // std::cout << std::endl;
// --

    if(sm_vs.size() < 3)
    {
      std::cerr << "Warning: sm_vs has size 1 or 2: offset crossing face at a single point?" << std::endl;
    }
    else
    {
      SM_face_descriptor sm_f = CGAL::Euler::add_face(sm_vs, sm);
      if(sm_f == boost::graph_traits<Mesh>::null_face())
        std::cerr << "Error while adding face to mesh" << std::endl;
    }
  }

  add_horizontal_faces(sm);

  // some vertices might not be used if the offset is smaller than the max offset value
  CGAL::Polygon_mesh_processing::remove_isolated_vertices(sm);

  CGAL::IO::write_OFF("sm_3D.off", sm, CGAL::parameters::stream_precision(17));
}

int main(int argc, char** argv)
{
  const int argc_check = argc - 1;

  char* poly_filename = nullptr;
  char* weights_filename = nullptr;
  FT offset = def_offset;
  bool add_bottom_face = false;
  bool inward = true;
  std::size_t seed = std::time(nullptr);

  for(int i = 1; i < argc; ++i)
  {
    if(!strcmp("-h", argv[i]) || !strcmp("--help", argv[i]) || !strcmp("-?", argv[i]))
    {
      std::cout << "Usage: " << argv[0] << "[options].\n"
        "Options:\n"
        "   -i <input_filename>: input polygon filename.\n"
        "   -w <weights_filename>: weights. Format: one weight per line, a space to separate borders.\n"
        "   -t <value>: max time. Must be strictly positive.\n"
        "   -in: grow inward.\n"
        "   -out: grow outward.\n"
                << std::endl;

      return EXIT_FAILURE;
    } else if(!strcmp("-i", argv[i]) && i < argc_check) {
      poly_filename = argv[++i];
    } else if(!strcmp("-w", argv[i])) {
      weights_filename = argv[++i];
    } else if(!strcmp("-t", argv[i]) && i < argc_check) {
      offset = std::stod(argv[++i]);
    } else if(!strcmp("-out", argv[i]) && i < argc_check) {
      inward = false;
    } else if(!strcmp("-in", argv[i]) && i < argc_check) {
      inward = true;
    } else if(!strcmp("-s", argv[i]) && i < argc_check) {
      seed = std::stoi(argv[++i]);
    }
  }

  CGAL::Real_timer timer;
  timer.start();

  Polygon_with_holes_2 p;
  if(poly_filename == nullptr)
  {
    Polygon_2 poly;
    poly.push_back(Point_2(0,0));
    poly.push_back(Point_2(1,0));
    poly.push_back(Point_2(1,1));
    poly.push_back(Point_2(0,1));

    // poly.push_back(Point_2(0,0));
    // poly.push_back(Point_2(10,0));
    // poly.push_back(Point_2(10,5));
    // poly.push_back(Point_2(5,5));
    // poly.push_back(Point_2(5,1));
    // poly.push_back(Point_2(0,1));

    assert(poly.is_counterclockwise_oriented());
    p = Polygon_with_holes_2(poly);
  }
  else
  {
    if(!read_input_polygon(poly_filename, p))
      return EXIT_FAILURE;
  }

  // read weights
  std::vector<std::vector<FT> > weights;

  if(weights_filename == nullptr)
  {
    // random weights
    CGAL::Random rnd(seed);
    std::cout << "Seed is " << rnd.get_seed() << std::endl;

    std::vector<FT> border_weights;
    for(std::size_t i=0; i<p.outer_boundary().size(); ++i)
    {
      border_weights.push_back(rnd.get_int(1, 10));
      std::cout << border_weights.back() << std::endl;
    }

    std::cout << border_weights.size() << " outer contour weights" << std::endl;
    weights.push_back(border_weights);
    border_weights.clear();

    for(auto hit=p.holes_begin(); hit!=p.holes_end(); ++hit)
    {
      for(std::size_t i=0; i<hit->size(); ++i)
      {
        border_weights.push_back(rnd.get_int(1, 10));
        std::cout << border_weights.back() << std::endl;
      }

      std::cout << border_weights.size() << " hole weights" << std::endl;
      weights.push_back(border_weights);
      border_weights.clear();
    }
  }
  else
  {
    read_input_weights(weights_filename, weights);
  }

  timer.stop();
  std::cout << "Read took " << timer.time() << " s." << std::endl;

  // CGAL::draw(p.outer_boundary());

  timer.reset();
  timer.start();

  // You can pass the polygon via an iterator pair
  SsPtr ss;

  if(inward)
  {
    if(only_use_full_skeletons || offset == def_offset)
    {
      ss = CGAL::create_interior_weighted_straight_skeleton_2(p, weights, K());
    }
    else
    {
#if 1
      std::cerr << "Warning: partial interior SS not yet supported with weights" << std::endl;
#else
      std::vector<Polygon_2> no_holes;
      ss = CGAL::CGAL_SS_i::create_partial_interior_straight_skeleton_2(
              2 * offset,
              CGAL::CGAL_SS_i::vertices_begin(p),
              CGAL::CGAL_SS_i::vertices_end(p),
              no_holes.begin(),
              no_holes.end(),
              K());
#endif
    }
  }
  else // outward
  {
#if 1
    std::cerr << "Warning: exterior SS not yet supported with weights" << std::endl;
#else
    if(offset == def_offset)
    {
      std::cerr << "You cannot use the default offset value with outward offset" << std::endl;
    }
    else if(only_use_full_skeletons)
    {
      ss = CGAL::create_exterior_straight_skeleton_2(2 * offset, p, K());
    }
    else
    {
      std::vector<Polygon_2> no_holes;
      ss = CGAL::CGAL_SS_i::create_partial_exterior_straight_skeleton_2(
              2 * offset,
              CGAL::CGAL_SS_i::vertices_begin(p),
              CGAL::CGAL_SS_i::vertices_end(p),
              K());
    }
#endif
  }

  if(!ss)
  {
    std::cerr << "Error: encounter an issue during SS computation" << std::endl;
    return EXIT_FAILURE;
  }

  timer.stop();
  std::cout << "Straight skeleton computation took " << timer.time() << " s." << std::endl;

  FT max_time = 0;
  for(auto v : vertices(static_cast<const HDS&>(*ss)))
    if(max_time < v->time())
      max_time = v->time();
  std::cout << "Max time in skeleton is " << max_time << std::endl;

  // print_straight_skeleton(*ss);
  CGAL::draw(*ss);

  timer.reset();
  timer.start();

  draw_3D_mesh(p, *ss, offset, add_bottom_face);

  timer.stop();
  std::cout << "Conversion to 3D took " << timer.time() << " s." << std::endl;

  return EXIT_SUCCESS;
}