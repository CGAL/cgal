#ifndef CGAL_ORIENT_POLYHEDRON_3
#define CGAL_ORIENT_POLYHEDRON_3
#include <CGAL/IO/generic_print_polyhedron.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Modifier_base.h>

#include <set>
#include <stack>
#include <algorithm>
#include <iostream>

#include <boost/array.hpp>

namespace CGAL {

namespace internal {
// going to be removed

//typedef internal::Polygon_soup_orienter<Polyhedron> Polygon_soup_orienter;
//typedef internal::Polyhedron_to_polygon_soup_writer<Polygon_soup_orienter> Soup_writer;

//Polygon_soup_orienter soup;
//Soup_writer writer(&soup);
//generic_print_polyhedron(std::cerr, polyhedron, writer);

template<class Polygon_soup_orienter>
struct Polyhedron_to_polygon_soup_writer {

  Polygon_soup_orienter* soup;
  typename Polygon_soup_orienter::Polygon_3 polygon;

  Polyhedron_to_polygon_soup_writer(Polygon_soup_orienter* soup) : soup(soup), polygon() {
  }

  void write_header( std::ostream&,std::size_t /* vertices */,std::size_t /* halfedges */,std::size_t /* facets */,bool /* normals */ = false ) 
  { }

  void write_footer() { }

  void write_vertex( const double& x, const double& y, const double& z) {
    soup->points.push_back(typename Polygon_soup_orienter::Point_3(x, y, z));
  }

  void write_normal( const double& /* x */, const double& /* y */, const double& /* z */) { }

  void write_facet_header() { }

  void write_facet_begin( std::size_t no) {
    polygon.clear();
    polygon.reserve(no);
  }

  void write_facet_vertex_index( std::size_t index) {
    polygon.push_back(index);
  }

  void write_facet_end() {
    soup->polygons.push_back(polygon);
    polygon.clear();
  }
}; // end struct Polyhedron_to_soup_writer
}// namespace internal

/** 
 * Class providing the functionality of orienting a polygon soup.
 * @tparam Polyhedron a %CGAL polyhedron
 */ 
template<class Polyhedron>
class Polygon_soup_orienter :
  public Modifier_base<typename Polyhedron::HalfedgeDS> 
{
  template<class PS>
  friend class Polyhedron_to_polygon_soup_writer;
public:
  typedef typename Polyhedron::Point_3 Point_3;

private:
  typedef Polygon_soup_orienter<Polyhedron> Self;
  typedef std::vector<std::size_t> Polygon_3;
  typedef std::vector<Point_3> Points;
  typedef std::map<std::pair<std::size_t, std::size_t>, std::set<std::size_t> > Edges_map;
  typedef boost::array<std::size_t, 2> Edge;
  typedef std::vector<Polygon_3> Polygons;
  typedef std::set<Edge> Edges;
  typedef Polygons::size_type size_type;

  Points points;
  Polygons polygons;
  Edges_map edges;
  Edges non_manifold_edges;

public:
  friend std::istream& operator>>(std::istream& stream, Self& soup)
  {
#if CGAL_VERSION_NR >= 1030700091
    typedef std::size_t indices_t;
#else
    typedef boost::int32_t indices_t;
#endif

    CGAL::File_scanner_OFF scanner(stream);
    soup.clear();
    soup.points.resize(scanner.size_of_vertices());
    soup.polygons.resize(scanner.size_of_facets());
    for (indices_t i = 0; i < scanner.size_of_vertices(); ++i) {
      double x, y, z, w;
      scanner.scan_vertex( x, y, z, w);
      soup.points[i] = Self::Point_3(x, y, z, w);
      scanner.skip_to_next_vertex( i);
    }
    if(!stream)
      return stream;

    for (indices_t i = 0; i < scanner.size_of_facets(); ++i) {
      indices_t no;

      scanner.scan_facet( no, i);
      soup.polygons[i].resize(no);
      for(indices_t j = 0; j < no; ++j) {
        indices_t id;
        scanner.scan_facet_vertex_index(id, i);
        if(id < scanner.size_of_vertices())
        {
          soup.polygons[i][j] = id;
        }
        else
          return stream;
      }
    }
    soup.fill_edges();
    return stream;
  }
  friend std::ostream& operator>>(std::ostream& stream, const Self& soup)
  {
    //CGAL::File_writer_OFF writer(stream);
    //writer.write_header(stream,
    //  soup.points.size(),
    //  0,
    //  soup.polygons.size());
    //for(Self::size_type i = 0, end = soup.points.size();
    //  i < end; ++i)
    //{
    //  const Point_3& p = soup->points[i];
    //  writer.write_vertex( p.x(), p.y(), p.z() );
    //}
    //writer.write_facet_header();
    //for(size_type i = 0, end = soup->polygons.size();
    //  i < end; ++i)
    //{
    //  const Polygon_soup_orienter::Polygon_3& polygon = soup->polygons[i]; 
    //  const size_type size = polygon.size();
    //  writer.write_facet_begin(size);
    //  for(size_type j = 0; j < size; ++j) {
    //    writer.write_facet_vertex_index(polygon[j]);
    //  }
    //  writer.write_facet_end();
    //}
    //writer.write_footer();

    //return out;
    return stream;
  }

  /** 
   * Default constructor which creates an empty polygon soup.
   */ 
  Polygon_soup_orienter() { }

  /** 
   * Constructor for creating polygon soup from points and polygons.
   * @param points points inside polygon soup
   * @param polygons each internal vector defines a polygon with indices of @a points
   */ 
  Polygon_soup_orienter(const std::vector<Point_3>& points,
    const std::vector<std::vector<std::size_t> >& polygons)
    : points(points), polygons(polygons)
  {
    fill_edges();
  }

private:
  void fill_edges() {
    // Fill edges
    edges.clear();
    for(size_type i = 0; i < polygons.size(); ++i)
    {
      const size_type size = polygons[i].size();
      for(size_type j = 0; j < size; ++j) {
        const std::size_t& i0 = polygons[i][j];
        const std::size_t& i1 = polygons[i][ j+1 < size ? j+1: 0];
        edges[std::make_pair(i0, i1)].insert(i);
//         qDebug() << tr("edges[std::make_pair(%1, %2)].insert(%3). Size=%4")
//           .arg(i0).arg(i1).arg(i).arg(edges[std::make_pair(i0, i1)].size());
      }
    }

    // Fill non-manifold edges
    non_manifold_edges.clear();
    for(size_type i = 0; i < polygons.size(); ++i)
    {
      const size_type size = polygons[i].size();
      for(size_type j = 0; j < size; ++j) {
        const std::size_t& i0 = polygons[i][j];
        const std::size_t& i1 = polygons[i][ j+1 < size ? j+1: 0];
        if( (i0 < i1) && 
            (edges[std::make_pair(i0, i1)].size() +
             edges[std::make_pair(i1, i0)].size() > 2) )
        {
          Edge edge;
          edge[0] = i0;
          edge[1] = i1;
          if(i0 > i1) std::swap(edge[0], edge[1]);
          non_manifold_edges.insert(edge);
        }
      }
    }
  }

  void clear() {
    points.clear();
    polygons.clear();
    edges.clear();
    non_manifold_edges.clear();
  }

  void inverse_orientation(const std::size_t index) {
    std::reverse(polygons[index].begin(), polygons[index].end());
  }

  void operator()(typename Polyhedron::HalfedgeDS& out_hds)
  {
    typedef typename Polyhedron::HalfedgeDS Output_HDS;

    Polyhedron_incremental_builder_3<Output_HDS> builder(out_hds);

    builder.begin_surface(points.size(),
      polygons.size(),
      edges.size() * 2);
    for(size_type i = 0, end = points.size();
      i < end; ++i)
    {
      builder.add_vertex(points[i]);
    }
    for(size_type i = 0, end = polygons.size();
      i < end; ++i)
    {
      const Polygon_3& polygon = polygons[i]; 
      const size_type size = polygon.size();
      builder.begin_facet();
      for(size_type j = 0; j < size; ++j) {
        builder.add_vertex_to_facet(polygon[j]);
      }
      builder.end_facet();
    }
    builder.end_surface();
  }

public:
  /** 
   * Orients polygon soup and fills @a polyhedron if orientation is successful.
   * @param[out] polyhedron polyhedron to be filled
   * @return true if orientation is successful
   */ 
  bool orient(Polyhedron& polyhedron)
  {
    std::vector<bool> oriented;
    std::stack<std::size_t> stack;
    using std::make_pair;

    // no polygon is oriented
    oriented.resize(polygons.size());

    size_type polygon_index = 0;
    bool success = true;

    while (polygon_index != polygons.size()) 
    {
      while ( polygon_index != polygons.size() && oriented[polygon_index] ) {
        ++polygon_index;
      }
      if(polygon_index == polygons.size()) break;

      //     qDebug() << tr("Seed %1...\n").arg(polygon_index);
      oriented[polygon_index] = true;
      stack.push(polygon_index);
      while(! stack.empty() )
      {
        const size_type to_be_oriented_index = stack.top();
        //       qDebug() << tr("polygon #%1").arg(to_be_oriented_index);
        stack.pop();
        const size_type size = polygons[to_be_oriented_index].size();
        for(size_type ih = 0 ; ih < size ; ++ih) {
          size_type ihp1 = ih+1;
          if(ihp1>=size) ihp1 = 0;
          const std::size_t& i1 = polygons[to_be_oriented_index][ih];
          const std::size_t& i2 = polygons[to_be_oriented_index][ihp1];

          Edge edge;
          edge[0] = i1;
          edge[1] = i2;
          if(i1 > i2) std::swap(edge[0], edge[1]);

          if(non_manifold_edges.count(edge) > 0) {
            continue;
          }

          //         qDebug() << tr("edge %3-%4 (%1,%2)").arg(i1).arg(i2).arg(ih).arg(ihp1);
          // edge (i1,i2)
          Edges_map::iterator it_same_orient = edges.find(make_pair(i1, i2));
          // edges (i2,i1)
          Edges_map::iterator it_other_orient = edges.find(make_pair(i2, i1));

          CGAL_assertion(it_same_orient != edges.end());
          if(it_same_orient->second.size() > 1) {
            if((it_other_orient != edges.end() && it_other_orient->second.size() > 0) ||
              it_same_orient->second.size() > 2) {
                // three polygons at the edge
                //             qDebug() << "three polygons at the edge";
                success = false; // non-orientable
            }
            {
              // one neighbor polyhedron, opposite orientation
              size_type index = *(it_same_orient->second.begin());
              if(index == to_be_oriented_index)
                index = *(++it_same_orient->second.begin());
              if(oriented[index]) {
                //               qDebug() << tr("neighbor polygon #%1 is already oriented, but in opposite orientation").arg(index);
                success = false; // non-orientable
                continue; // next edge
              }

              // reverse the orientation
              const size_type size = polygons[index].size();
              for(size_type j = 0; j < size; ++j) {
                const std::size_t& i0 = polygons[index][j];
                const std::size_t& i1 = polygons[index][ j+1 < size ? j+1: 0];
                CGAL_assertion_code(const bool r = )
                  edges[std::make_pair(i0, i1)].erase(index);
                CGAL_assertion(r);
              }
              inverse_orientation(index);
              for(size_type j = 0; j < size; ++j) {
                const std::size_t& i0 = polygons[index][j];
                const std::size_t& i1 = polygons[index][ j+1 < size ? j+1: 0];
                edges[std::make_pair(i0, i1)].insert(index);
              }
              //             qDebug() << tr("inverse the orientation of polygon #%1\n").arg(index);
              oriented[index] = true;
              stack.push(index);
            }
          }
          else if(it_other_orient != edges.end() && it_other_orient->second.size() == 1) {
            // one polygon, same orientation
            const size_type index = *(it_other_orient->second.begin());
            if(oriented[index])
              continue;
            oriented[index] = true;
            //           qDebug() << tr("keep the orientation of polygon #%1\n").arg(index);
            stack.push(index);
          }
          else {
            // qDebug() << "else" << it_same_orient->second.size() << 
            //   (it_other_orient == edges.end() ? 0 : it_other_orient->second.size());
            if(it_same_orient->second.size() != 1 || 
              (it_other_orient != edges.end() && it_other_orient->second.size() > 0)) 
            {
              // qDebug() << tr("non orientable");
              success = false; // non-orientable
            }
          }
        } // end for on all edges of one 
      } // end while loop on the polygons of the connected component
    } // end while loop on all non-oriented polygons remaining 

    if(success) {
      polyhedron.delegate(*this);
    }
    return success;
  }
};

}// namespace CGAL

#endif // CGAL_ORIENT_POLYHEDRON_3
