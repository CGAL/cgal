#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>

#include <CGAL/Bbox_2.h>
#include <CGAL/number_utils.h>

#include <fstream>
#include <iostream>
#include <list>
#include <set>
#include <vector>

#define SDG_DRAW_DEBUG // debug log
#define SDG_DRAW_DUMP_FILES // print input / ouput
// #define SINGLE_INPUT_FILE // if not defined, each segment of the input has its own file

#ifdef SDG_DRAW_DUMP_FILES_PP // also print parabolas
#define SDG_DRAW_DUMP_FILES
#endif

typedef CGAL::Simple_cartesian<double>                                            CK;
typedef CK::Point_2                                                               Point_2;
typedef CK::Segment_2                                                             Segment_2;
typedef CGAL::Field_with_sqrt_tag                                                 CM;
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt               EK;
typedef CGAL::Field_with_sqrt_tag                                                 EM;
typedef CGAL::Simple_cartesian<CGAL::Interval_nt<false> >                         FK;
typedef CGAL::Field_with_sqrt_tag                                                 FM;
typedef CGAL::Segment_Delaunay_graph_filtered_traits_2<CK, CM, EK, EM, FK, FM>    Gt;
typedef CGAL::Segment_Delaunay_graph_2<Gt>                                        SDG;

///////////////////////////////// CODE ABOUT EXACT DUALS ///////////////////////////////////////////

template < typename ExactSite, typename ExactKernel,
           typename SDGSite,   typename InputKernel >
ExactSite convert_site_to_exact(const SDGSite &site,
                                const InputKernel & /*k*/,
                                const ExactKernel & /*ek*/)
{
  using To_exact = CGAL::Cartesian_converter<InputKernel, ExactKernel>;
  To_exact to_exact;

  // Note: in theory, a site can be constructed from more than just one or two points
  // (e.g. 4 points for the segment defined by the intersection of two segments). Thus, it
  // would be better to convert the input points at the very beginning and just maintain
  // a type of map between the base and exact sites.
  ExactSite es;
  if(site.is_point())
    es = ExactSite::construct_site_2(to_exact(site.point()));
  else
    es = ExactSite::construct_site_2(to_exact(site.segment().source()), to_exact(site.segment().target()));

  return es;
}

// Dual (Voronoi site) of an SDG face
template < typename FiniteFacesIterator, typename InputKernel, typename ExactKernel >
typename ExactKernel::Point_2 exact_primal(const FiniteFacesIterator sdg_f,
                                           const InputKernel& k,
                                           const ExactKernel& ek)
{
  using Exact_SDG_traits = CGAL::Segment_Delaunay_graph_traits_2<ExactKernel>;
  using Exact_site_2 = typename Exact_SDG_traits::Site_2;

  static Exact_SDG_traits e_sdg_gt;
  const Exact_site_2 es0 = convert_site_to_exact<Exact_site_2>(sdg_f->vertex(0)->site(), k, ek);
  const Exact_site_2 es1 = convert_site_to_exact<Exact_site_2>(sdg_f->vertex(1)->site(), k, ek);
  const Exact_site_2 es2 = convert_site_to_exact<Exact_site_2>(sdg_f->vertex(2)->site(), k, ek);

  return e_sdg_gt.construct_svd_vertex_2_object()(es0, es1, es2);
}

// Dual (Voronoi edge) of an SDG edge
// this function is identical 'SDG::primal()', but with a conversion to exact sites
template < typename Edge, typename InputKernel, typename ExactKernel >
CGAL::Object exact_primal(const Edge& e,
                          const SDG& sdg,
                          const InputKernel& k,
                          const ExactKernel& ek)
{
  using Exact_SDG_traits = CGAL::Segment_Delaunay_graph_traits_2<ExactKernel>;
  using Exact_site_2 = typename Exact_SDG_traits::Site_2;

  using DT = CGAL::Field_with_sqrt_tag;
  using Construct_sdg_bisector_2 = CGAL::SegmentDelaunayGraph_2::Construct_sdg_bisector_2<Exact_SDG_traits, DT>;
  using Construct_sdg_bisector_ray_2 = CGAL::SegmentDelaunayGraph_2::Construct_sdg_bisector_ray_2<Exact_SDG_traits, DT>;
  using Construct_sdg_bisector_segment_2 = CGAL::SegmentDelaunayGraph_2::Construct_sdg_bisector_segment_2<Exact_SDG_traits, DT>;

  CGAL_precondition(!sdg.is_infinite(e));

  if(sdg.dimension() == 1)
  {
    Exact_site_2 p = convert_site_to_exact<Exact_site_2>((e.first)->vertex(sdg.cw(e.second))->site(), k, ek);
    Exact_site_2 q = convert_site_to_exact<Exact_site_2>((e.first)->vertex(sdg.ccw(e.second))->site(), k, ek);

    return make_object(Construct_sdg_bisector_2()(p, q));
  }

  // dimension == 2
  // neither of the two adjacent faces is infinite
  if((!sdg.is_infinite(e.first)) && (!sdg.is_infinite(e.first->neighbor(e.second))))
  {
    Exact_site_2 p = convert_site_to_exact<Exact_site_2>((e.first)->vertex(sdg.ccw(e.second))->site(), k, ek);
    Exact_site_2 q = convert_site_to_exact<Exact_site_2>((e.first)->vertex(sdg.cw(e.second))->site(), k, ek);
    Exact_site_2 r = convert_site_to_exact<Exact_site_2>((e.first)->vertex(e.second)->site(), k, ek);
    Exact_site_2 s = convert_site_to_exact<Exact_site_2>(sdg.tds().mirror_vertex(e.first, e.second)->site(), k, ek);

    return Construct_sdg_bisector_segment_2()(p, q, r, s);
  }

  // both of the adjacent faces are infinite
  if(sdg.is_infinite(e.first) && sdg.is_infinite(e.first->neighbor(e.second)))
  {
    Exact_site_2 p = convert_site_to_exact<Exact_site_2>((e.first)->vertex(sdg.cw(e.second))->site(), k, ek);
    Exact_site_2 q = convert_site_to_exact<Exact_site_2>((e.first)->vertex(sdg.ccw(e.second))->site(), k, ek);

    return make_object(Construct_sdg_bisector_2()(p, q));
  }

  // only one of the adjacent faces is infinite
  CGAL_assertion(sdg.is_infinite(e.first) || sdg.is_infinite(e.first->neighbor(e.second)));
  CGAL_assertion(!(sdg.is_infinite(e.first) && sdg.is_infinite(e.first->neighbor(e.second))));
  CGAL_assertion(sdg.is_infinite(e.first->vertex(e.second)) || sdg.is_infinite(sdg.tds().mirror_vertex(e.first, e.second)));

  Edge ee = e;
  if(sdg.is_infinite(e.first->vertex(e.second)))
  {
    ee = Edge(e.first->neighbor(e.second),
              e.first->neighbor(e.second)->index(sdg.tds().mirror_vertex(e.first, e.second)));
  }

  Exact_site_2 p = convert_site_to_exact<Exact_site_2>(ee.first->vertex(sdg.ccw(ee.second))->site(), k, ek);
  Exact_site_2 q = convert_site_to_exact<Exact_site_2>(ee.first->vertex(sdg.cw(ee.second))->site(), k, ek);
  Exact_site_2 r = convert_site_to_exact<Exact_site_2>(ee.first->vertex(ee.second)->site(), k, ek);

  return make_object(Construct_sdg_bisector_ray_2()(p, q, r));
}

///////////////////////////////// CODE TO DRAW A CROPPED DIAGRAM ///////////////////////////////////

// Split a Voronoi edge that is a parabola (one site is a point, one site is a segment) into small segments
template <typename OutputKernel, typename SegmentContainer>
void segment_parabola(const CGAL::Parabola_segment_2<OutputKernel>& p,
                      const CGAL::Bbox_2& scaled_bbox,
                      SegmentContainer& segment_list)
{
  using FT = typename OutputKernel::FT;
  using Point_2 = typename OutputKernel::Point_2;

  const Point_2& o = p.origin();
  const Point_2& c = p.center();

  // @todo could be cached
  const Point_2 mm(scaled_bbox.xmin(), scaled_bbox.ymin());
  const Point_2 mM(scaled_bbox.xmin(), scaled_bbox.ymax());
  const Point_2 Mm(scaled_bbox.xmax(), scaled_bbox.ymin());
  const Point_2 MM(scaled_bbox.xmax(), scaled_bbox.ymax());

  FT s = CGAL::squared_distance(mm, c);
  s = (std::max)(s, CGAL::squared_distance(mM, c));
  s = (std::max)(s, CGAL::squared_distance(Mm, c));
  s = (std::max)(s, CGAL::squared_distance(MM, c));

  s = CGAL::sqrt(s) - CGAL::sqrt(CGAL::squared_distance(c, o));

  const double lx = scaled_bbox.xmax() - scaled_bbox.xmin();
  const double ly = scaled_bbox.ymax() - scaled_bbox.ymin();

  const double max_length = 0.001 * (std::min)(lx, ly); // max length in the discretization of a parabola

#ifdef SDG_DRAW_DEBUG
  std::cout << "origin & center: " << o << " [[]] " << c << std::endl;
  std::cout << "max length " << max_length << std::endl;
#endif

  std::vector<Point_2> points;
  p.generate_points(points, max_length, -s, s);

#ifdef SDG_DRAW_DEBUG
  std::cout << "Discretize parabola into " << points.size() << " points" << std::endl;
  for(const auto& pt : points)
    std::cout << CGAL::to_double(pt.x()) << " " << CGAL::to_double(pt.y()) << " 0" << std::endl;
  std::cout << "end of parabola" << std::endl;
#endif

  if(points.size() < 2)
    return;

  for(std::size_t i=0, ps=points.size()-1; i<ps; ++i)
    segment_list.emplace_back(points[i], points[i+1]);

#ifdef SDG_DRAW_DUMP_FILES_PP
  std::cout << "filled segment list. Dumping parabola..." << std::endl;
  static int parabola_id = 0;
  std::stringstream oss;
  oss << "parabola_" << parabola_id++ << ".cgal" << std::ends;
  std::ofstream out(oss.str().c_str());
  out.precision(17);

  for(std::size_t i=1, ps=(points.size()-1); i<ps; ++i)
  {
    out << "2 " << CGAL::to_double(points[i].x()) << " ";
    out << CGAL::to_double(points[i].y()) << " 0 ";
    out << CGAL::to_double(points[i+1].x()) << " ";
    out << CGAL::to_double(points[i+1].y()) << " 0\n";
  }
#endif
}

template< typename OutputKernel, typename LineContainer, typename RayContainer, typename SegmentContainer>
void fill_Voronoi_structure(const SDG& sdg,
                            const CGAL::Bbox_2& scaled_bbox,
                            LineContainer& line_list,
                            RayContainer& ray_list,
                            SegmentContainer& segment_list)
{
  using Line_2 = typename OutputKernel::Line_2;
  using Ray_2 = typename OutputKernel::Ray_2;
  using Segment_2 = typename OutputKernel::Segment_2;
  using SDG_traits = CGAL::Segment_Delaunay_graph_traits_2<OutputKernel>;

  CK k;
  OutputKernel ek;

  Line_2 l;
  Ray_2 r;
  Segment_2 s;
  CGAL::Parabola_segment_2<SDG_traits> p;

  int nl = 0, ns = 0, nr = 0, np = 0;

  typename SDG::Finite_edges_iterator eit = sdg.finite_edges_begin(),
                                      eend = sdg.finite_edges_end();
  for (; eit != eend; ++eit)
  {
    CGAL::Object o = exact_primal(*eit, sdg, k, ek);

    if(CGAL::assign(l, o)) { line_list.push_back(l); ++nl; }
    if(CGAL::assign(s, o)) { segment_list.push_back(s); ++ns; }
    if(CGAL::assign(r, o)) { ray_list.push_back(r); ++nr; }
    if(CGAL::assign(p, o)) { segment_parabola(p, scaled_bbox, segment_list); ++np; }
  }

#ifdef SDG_DRAW_DEBUG
  std::cout << nl << " lines" << std::endl;
  std::cout << ns << " segments" << std::endl;
  std::cout << nr << " rays" << std::endl;
  std::cout << np << " parabolas" << std::endl;
#endif
}

template <typename OutputKernel, typename T, typename OutputIterator>
struct Box_clipper
{
  bool operator()(const T& obj,
                  const typename OutputKernel::Iso_rectangle_2& bbox,
                  OutputIterator oit) const
  {
    CGAL::Object obj_cgal = CGAL::intersection(obj, bbox);

    typename OutputKernel::Segment_2 s;
    bool ret = CGAL::assign(s, obj_cgal);
    if(ret)
      oit++ = s;

    return ret;
  }
};

template <typename OutputKernel, typename OutputIterator>
struct Box_clipper<OutputKernel, typename OutputKernel::Segment_2, OutputIterator>
{
  bool operator()(const typename OutputKernel::Segment_2& is,
                  const typename OutputKernel::Iso_rectangle_2& bbox,
                  OutputIterator oit) const
  {
    if(bbox.has_on_unbounded_side(is.source()) && bbox.has_on_unbounded_side(is.target()))
      return false;

    if(!bbox.has_on_unbounded_side(is.source()) && !bbox.has_on_unbounded_side(is.target()))
    {
      oit++ = is;
      return true;
    }

    CGAL::Object obj_cgal = CGAL::intersection(is, bbox);

    typename OutputKernel::Segment_2 s;
    bool ret = CGAL::assign(s, obj_cgal);
    if(ret)
      oit++ = s;

    return ret;
  }
};

template <typename OutputKernel,
          typename LineContainer, typename RayContainer, typename SegmentContainer,
          typename OutputIterator>
void draw_dual(const LineContainer& lines,
               const RayContainer& rays,
               const SegmentContainer& segments,
               const typename OutputKernel::Iso_rectangle_2& bbox,
               OutputIterator oit)
{
  using Line_2 = typename OutputKernel::Line_2;
  using Ray_2 = typename OutputKernel::Ray_2;
  using Segment_2 = typename OutputKernel::Segment_2;

  typedef Box_clipper<OutputKernel, Line_2, OutputIterator> Line_clipper;
  typedef Box_clipper<OutputKernel, Ray_2, OutputIterator> Ray_clipper;
  typedef Box_clipper<OutputKernel, Segment_2, OutputIterator> Segment_clipper;

  Line_clipper lc;
  Ray_clipper rc;
  Segment_clipper sc;

  for(const Line_2& line : lines)
    lc(line, bbox, oit);
  for(const Ray_2& ray : rays)
    rc(ray, bbox, oit);
  for(const Segment_2& segment : segments)
    sc(segment, bbox, oit);
}

template< typename OutputKernel, typename OutputIterator >
OutputIterator draw_SDG(const SDG& sdg,
                        const typename OutputKernel::Iso_rectangle_2& bbox,
                        const CGAL::Bbox_2& scaled_bbox,
                        OutputIterator oit)
{
  using Line_2 = typename OutputKernel::Line_2;
  using Ray_2 = typename OutputKernel::Ray_2;
  using Segment_2 = typename OutputKernel::Segment_2;

  std::list<Line_2> line_list;
  std::list<Ray_2> ray_list;
  std::list<Segment_2> segment_list;

  fill_Voronoi_structure<OutputKernel>(sdg, scaled_bbox, line_list, ray_list, segment_list);
  draw_dual<OutputKernel>(line_list, ray_list, segment_list, bbox, oit);

  return oit;
}

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cout << std::fixed;

  if(argc < 2)
    std::cout << "Usage: " << argv[0] << " input.cin" << std::endl;

#ifdef SDG_DRAW_DUMP_FILES
  // dump the input
  std::ifstream input((argc > 1) ? argv[1] : "./data/sites.cin");
  assert(input);

  std::string str;
  int i = 0;

# ifdef SINGLE_INPUT_FILE
  std::ofstream in_out("input.cgal");
  in_out.precision(17);
# endif

  // read line by line, and distinguish between point and segment
  std::string line;
  while(std::getline(input, line))
  {
# ifndef SINGLE_INPUT_FILE
    std::stringstream oss;
    oss << "input_" << i++ << ".cgal" << std::ends;
    std::ofstream in_out(oss.str().c_str());
    in_out.precision(17);
# endif

    std::stringstream ss(line);
    ss >> str;
    if(str == "s")
    {
      double x0, y0, x1, y1;
      ss >> x0 >> y0 >> x1 >> y1;

      in_out << "2 " << x0 << " " << y0 << " 0 ";
      in_out << x1 << " " << y1 << " 0\n";
    }
    else if(str == "p")
    {
      double x, y;
      ss >> x >> y;

      in_out << "2 " << x << " " << y << " 0\n";
      in_out << x << " " << y << " 0\n"; // abusive 0-length polyline
    }
    else
    {
      std::cerr << "Error: Unknown input: " << str << std::endl;
      return EXIT_FAILURE;
    }

# ifndef SINGLE_INPUT_FILE
    in_out.close();
# endif
  }
# ifdef SINGLE_INPUT_FILE
    in_out.close();
# endif
#endif

  std::ifstream ifs((argc > 1) ? argv[1] : "./data/sites.cin");
  assert(ifs);

  // polygon points
  std::set<Gt::Point_2> all_points;

  // segments of the polygon as a pair of point indices
  std::vector<Gt::Point_2> points;
  std::vector<Gt::Segment_2> segments;
  SDG::Site_2 site;

  // read segment input, format:
  // s [x0 y0 x1 y1]
  // p [x y]

  while (ifs >> site)
  {
    if(site.is_segment())
    {
      all_points.insert(site.source_of_supporting_site());
      all_points.insert(site.target_of_supporting_site());
      segments.push_back(site.segment());
    }
    else
    {
      all_points.insert(site.point());
      points.push_back(site.point());
    }
  }

  if(points.empty() && segments.empty())
  {
    std::cerr << "Nothing in input..." << std::endl;
    return EXIT_SUCCESS;
  }

  // insert the sites all at once using spatial sorting to speed the insertion
  SDG sdg;
  sdg.insert_points(points.begin(), points.end());
  sdg.insert_segments(segments.begin(), segments.end());
  assert(sdg.is_valid(true, 1)); // validate the segment Delaunay graph

  typedef EK OK; // output kernel
  std::list<OK::Segment_2> svd_edges;

  // Get the bbox of the input points, and grow it a bit
  const CGAL::Bbox_2 bbox = bbox_2(all_points.begin(), all_points.end());
  const double xmin = bbox.xmin(), xmax = bbox.xmax();
  const double ymin = bbox.ymin(), ymax = bbox.ymax();
  const double xmid = 0.5 * (xmin + xmax), ymid = 0.5 * (ymin + ymax);
  const double scaling_factor = 3.; // '0.5' gives the identity
  const double lx = scaling_factor * (xmax - xmin),
               ly = scaling_factor * (ymax - ymin);
  const CGAL::Bbox_2 scaled_bbox(xmid - lx, ymid - ly, xmid + lx, ymid + ly);
  const OK::Iso_rectangle_2 bounding_iso_rec(scaled_bbox);

#ifdef SDG_DRAW_DEBUG
  std::cout << "bbox: " << bbox.xmin() << " " << bbox.ymin() << std::endl;
  std::cout << "bbox: " << bbox.xmax() << " " << bbox.ymax() << std::endl;
  std::cout << "lx/y: " << lx << " " << ly << std::endl;
  std::cout << "Scaled bbox: " << scaled_bbox.xmin() << " " << scaled_bbox.ymin() << std::endl;
  std::cout << "Scaled bbox: " << scaled_bbox.xmax() << " " << scaled_bbox.ymax() << std::endl;
#endif

  draw_SDG<OK>(sdg, bounding_iso_rec, scaled_bbox, std::back_inserter(svd_edges));

#ifdef SDG_DRAW_DEBUG
  std::cout << "Edges of the diagram (exact form):" << std::endl;
  for(const OK::Segment_2& edge : svd_edges)
    std::cout << edge << "\n";

  std::cout << "Now in the double, approximate format..." << std::endl;
  for(const OK::Segment_2& edge : svd_edges)
  {
    std::cout << CGAL::to_double(edge.source().x()) << " ";
    std::cout << CGAL::to_double(edge.source().y()) << " ";
    std::cout << CGAL::to_double(edge.target().x()) << " ";
    std::cout << CGAL::to_double(edge.target().y()) << "\n";
  }
#endif

#ifdef SDG_DRAW_DUMP_FILES
  // This file can be visualized with the CGAL 3D Polyhedron Demo
  std::ofstream out("dual.cgal");
  out.precision(17);
  for(const OK::Segment_2& edge : svd_edges)
  {
    out << "2 " << CGAL::to_double(edge.source().x()) << " ";
    out << CGAL::to_double(edge.source().y()) << " 0 ";
    out << CGAL::to_double(edge.target().x()) << " ";
    out << CGAL::to_double(edge.target().y()) << " 0\n";
  }
#endif

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
