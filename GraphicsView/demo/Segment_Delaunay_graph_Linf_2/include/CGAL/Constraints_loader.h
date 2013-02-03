#ifndef CONSTRAINTS_LOADER
#define CONSTRAINTS_LOADER

#include <CGAL/config.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Timer.h>
#include <vector>
#include <utility>
#include <iostream>

#include <boost/progress.hpp>

namespace CGAL {

template <class CDT>
class Constraints_loader {
  typedef typename CDT::Point_2 Point;

  typedef std::vector<Point> Points_container;

  typedef typename Points_container::size_type Index_type;
  typedef std::pair<Index_type,Index_type> Constraint;

  typedef std::vector<Constraint> Constraints_container;

  Points_container points;
  Constraints_container constraints;
  CDT& cdt;

  Bbox_2 points_bbox;

  bool insert_all_points_at_once_;

  template <typename Kernel, typename Iterator>
  struct Sort_traits_2 {
    Kernel k;
    Sort_traits_2 (const Kernel &kernel = Kernel())
      : k (kernel)
    {}

    typedef Iterator Point_2;

    struct Less_x_2 {
      Kernel k;
      Less_x_2 (const Kernel &kernel = Kernel())
        : k (kernel)
      {}
      bool operator() (const Point_2 &p, const Point_2 &q) const
      {
        return k.less_x_2_object() (*p, *q);
      }
    };

    Less_x_2
    less_x_2_object() const
    {
      return Less_x_2(k);
    }

    struct Less_y_2 {
      Kernel k;
      Less_y_2 (const Kernel &kernel = Kernel())
        : k (kernel)
      {}
      bool operator() (const Point_2 &p, const Point_2 &q) const
      {
        return k.less_y_2_object() (*p, *q);
      }
    };

    Less_y_2
    less_y_2_object() const
    {
      return Less_y_2(k);
    }
  };

  void update_bbox() {
    if(points.empty())
      return;
    points_bbox = points[0].bbox();
    for(Index_type i = 1; i < points.size(); ++i) {
      points_bbox = points_bbox + points[i].bbox();
    }
  }

  void insert_constaints_using_spatial_sort() const {
    typedef typename Points_container::const_iterator Points_iterator;
    typedef std::vector<Points_iterator> Indices;
    typedef std::vector<typename CDT::Vertex_handle> Vertices;
    
    Sort_traits_2<typename CDT::Geom_traits, Points_iterator> sort_traits;

    Indices indices;
    indices.reserve(points.size());
    for(Points_iterator it = points.begin(); it != points.end(); ++it) {
      indices.push_back(it);
    }
    std::random_shuffle(indices.begin(), indices.end());
    CGAL::spatial_sort(indices.begin(), indices.end(),
                       sort_traits);

    std::cerr << "Inserting points...";
    CGAL::Timer timer;
    timer.start();
    Vertices vertices;
    vertices.resize(points.size());
    typename CDT::Vertex_handle hint;
    for(typename Indices::const_iterator 
          pt_it_it = indices.begin(), end = indices.end();
        pt_it_it != end; ++pt_it_it) {
      typename CDT::Vertex_handle vh = cdt.insert(**pt_it_it, hint);
      hint = vh;
      vertices[*pt_it_it - points.begin()] = vh;
    }
    timer.stop();
    std::cerr << " done (" << timer.time() << "s)\n";

    std::cerr << "Inserting constraints...\n";
    boost::progress_display show_progress(constraints.size(), 
                                          std::cerr,
                                          "");
    timer.reset();
    timer.start();
    for(typename Constraints_container::const_iterator 
          cit = constraints.begin(), end = constraints.end();
        cit != end; ++cit) {
      ++show_progress;
      const typename CDT::Vertex_handle& v1 = vertices[cit->first];
      const typename CDT::Vertex_handle& v2 = vertices[cit->second];
      if(v1 != v2)
         {
           cdt.insert(v1, v2);
      }
    }
    timer.stop();
    std::cerr << " done (" << timer.time() << "s)\n";
  }

public:
  bool insert_all_points_at_once() const {
    return insert_all_points_at_once_;
  }

  void set_insert_all_points_at_once(const bool b) {
    insert_all_points_at_once_ = b;
  }

  Constraints_loader(CDT& _cdt, bool _insert_all_points_at_once = true)
    : cdt(_cdt), insert_all_points_at_once_(_insert_all_points_at_once)
  {

  }

  void clear() {
    points.clear();
    constraints.clear();
  }

  bool load_edg_file(std::istream& ifs) {
    std::cerr << "Loading edg file... ";
    CGAL::Timer timer;
    timer.start();
    bool not_first = false;
    int n;
    ifs >> n;
    if(!ifs.good())
      return false;
    Point p, q, qold;
    int point_counter = 0;
    while(ifs >> p >> q) {
      if(not_first && p == q) {
        continue;
      }
      if(p == qold) {
        points.push_back(q);
        constraints.push_back(std::make_pair(point_counter-1, point_counter));
        ++point_counter;
      }
      else {
        points.push_back(p);
        points.push_back(q);
        constraints.push_back(std::make_pair(point_counter, point_counter+1));
        point_counter += 2;
      }
      qold = q;
      not_first = true;
    }
    std::cerr << "done (" << timer.time() << "s)" << std::endl;
    update_bbox();
    insert_constaints_using_spatial_sort();
    return true;
  }

  bool load_plg_file(std::istream& ifs) {
    std::cerr << "Loading plg file... ";
    CGAL::Timer timer;
    timer.start();
    int n;
    int points_counter = 0;
    while(ifs >> n){
      Point first, p;
      ifs >> first;
      p = first;
      points.push_back(p);
      ++points_counter;
      while(--n){
        Point q;
        ifs >> q;
        if(p != q){
          points.push_back(q);
          constraints.push_back(std::make_pair(points_counter-1, points_counter));
          ++points_counter;
        }
        p = q;
      }
    }
    std::cerr << "done (" << timer.time() << "s)" << std::endl;
    update_bbox();
    insert_constaints_using_spatial_sort();
    return true;
  }

}; // end class Constraints_loader

} // namespace CGAL

#endif // CONSTRAINTS_LOADER
