#ifndef CGAL_BINOP_INTERSECTION_TESTS_H
#define CGAL_BINOP_INTERSECTION_TESTS_H

#include <bbox/segment_tree.h>
#include <bbox/box_traits.h>
#include <vector>
#include <iostream>
#include <CGAL/Timer.h>

CGAL_BEGIN_NAMESPACE

template< class SNC_decorator, class Selection >
class binop_intersection_tests {
protected:
  template< bool benchmark >
  struct auto_benchmark {
    CGAL::Timer timer;
    auto_benchmark() { timer.start(); }
    ~auto_benchmark() {
      timer.stop();
      if( benchmark ) {
        std::cout << std::endl
                  << "number of intersection tests: "
                  <<  intersection_counter << std::endl
                  << "time spent in intersection tests: "
                  << timer.time() << std::endl;
        intersection_counter = 0;
      }
    }
  };

  // why doesnt this work?
  //template<> struct auto_benchmark_impl<false> {};

  typedef typename SNC_decorator::SNC_structure          SNC_structure;
  typedef typename SNC_decorator::SNC_intersection       SNC_intersection;

  static unsigned int intersection_counter;

public:
  virtual void operator()( SNC_decorator d,
                           const Selection& BOP,
                           SNC_structure& sncp,
                           SNC_structure& snc1i,
                           SNC_structure& result ) = 0;
};


template< class SNC_decorator, class Selection >
unsigned int
binop_intersection_tests< SNC_decorator, Selection >::intersection_counter = 0;

template< class SNC_decorator, class Selection, bool benchmark = false >
struct binop_intersection_tests_allpairs :
    public binop_intersection_tests< SNC_decorator, Selection >
{
  typedef binop_intersection_tests< SNC_decorator,
                                    Selection>           parent;
  typedef typename parent::auto_benchmark<benchmark>     auto_benchmark;

  typedef typename SNC_decorator::SNC_structure          SNC_structure;
  typedef typename SNC_decorator::SNC_intersection       SNC_intersection;

  typedef typename SNC_decorator::Vertex_iterator        Vertex_iterator;
  typedef typename SNC_decorator::Vertex_handle          Vertex_handle;
  typedef typename SNC_decorator::Vertex_const_handle    Vertex_const_handle;
  typedef typename SNC_decorator::Halfedge_iterator      Halfedge_iterator;
  typedef typename SNC_decorator::Halfedge_handle        Halfedge_handle;
  typedef typename SNC_decorator::Halfedge_const_handle
                                  Halfedge_const_handle;
  typedef typename SNC_decorator::Halffacet_iterator     Halffacet_iterator;
  typedef typename SNC_decorator::Halffacet_handle       Halffacet_handle;
  typedef typename SNC_decorator::Halffacet_const_handle
                                  Halffacet_const_handle;
  typedef typename SNC_decorator::Halffacet_cycle_iterator
                                  Halffacet_cycle_iterator;
  typedef typename SNC_decorator::Infi_box               Infi_box;
  typedef typename SNC_decorator::Point_3                Point_3;
  typedef typename SNC_decorator::SHalfedge_handle       SHalfedge_handle;
  typedef typename SNC_decorator::SHalfedge_iterator     SHalfedge_iterator;
  typedef typename SNC_decorator::SHalfedge_around_facet_circulator
                                  SHalfedge_around_facet_circulator;

  virtual void operator()( SNC_decorator d,
                           const Selection& BOP,
                           SNC_structure& sncp,
                           SNC_structure& snc1i,
                           SNC_structure& result )
  {
    auto_benchmark bench;
    SNC_intersection is( sncp );
    TRACEN("start edge0 face1");
    Halfedge_iterator e0, e1;
    Halffacet_iterator f0, f1;
    Vertex_iterator v0, v1;

    Unique_hash_map<Halfedge_handle, bool> Ignore_halfedge(false);
    CGAL_nef3_forall_halfedges( e0, sncp) {
      TRACEN(PH(e0));
      if(!Ignore_halfedge[e0]) {
        CGAL_nef3_forall_facets( f1, snc1i) { 
          if(Infi_box::degree(d.plane(f1).d())>0) continue;
          Point_3 ip;
          ++intersection_counter;
          if( is.does_intersect_internally( d.segment(e0), f1, ip )) {
            TRACEN(" edge0 face1 intersection...");
            ip = normalized(ip);
            v0 = d.qualify_with_respect( ip, sncp, result);
            v1 = d.qualify_with_respect( ip, snc1i, result);
            d.binop_local_views( v0, v1, BOP, result);
            result.delete_vertex(v0);
            result.delete_vertex(v1);
            Ignore_halfedge[d.twin(e0)]=true;
          }
        }
      }
    }

    TRACEN("start edge1 face0");
    CGAL_nef3_forall_edges( e1, snc1i) {
        CGAL_nef3_forall_facets( f0,  sncp ) {
        if(Infi_box::degree(d.plane(f0).d())>0) continue;
        Point_3 ip;
        ++intersection_counter;
        if( is.does_intersect_internally( d.segment(e1), f0, ip )) {
          TRACEN(" edge1 face0 intersection...");
          ip = normalized(ip);
          Halffacet_cycle_iterator it; 
          CGAL_nef3_forall_facet_cycles_of(it,f0){ 
            TRACEN("facet cycle");
            SHalfedge_handle es;
            if ( assign(es,it)) {
              SHalfedge_around_facet_circulator start(es), end(es);
              CGAL_For_all(start,end) {
                TRACEN("vertex " << PH(d.source(d.previous(start))));
              }
            }
          }
          
          v1 = d.qualify_with_respect( ip, snc1i, result);
          v0 = d.qualify_with_respect( ip, sncp, result);

          d.binop_local_views( v0, v1, BOP, result);
          result.delete_vertex(v0);
          result.delete_vertex(v1);
        }
      }
    }

    //        SETDTHREAD(19*37);

    TRACEN("start edge0 edge1");
    TRACEN("=> edge edge intersection");
    CGAL_nef3_forall_edges( e0, sncp) {
      CGAL_nef3_forall_edges( e1, snc1i) {
        ++intersection_counter;
        Point_3 ip;
        if( is.does_intersect_internally( d.segment(e0),
                                          d.segment(e1), ip ) )
        {
          TRACEN(" edge0 edge1 intersection..." << ip);
          ip = normalized(ip);
          Vertex_handle v0, v1;
          v0 = d.qualify_with_respect( ip, sncp, result);
          v1 = d.qualify_with_respect( ip, snc1i, result);

          d.binop_local_views( v0, v1, BOP, result);
          result.delete_vertex(v0);
          result.delete_vertex(v1);
        }
      }
    }
  }
};

template< class SNC_decorator, class Selection, bool benchmark = false >
struct binop_intersection_tests_segment_tree :
  public binop_intersection_tests< SNC_decorator, Selection >
{
  typedef binop_intersection_tests< SNC_decorator,
                                    Selection >          parent;
  typedef typename parent::auto_benchmark<benchmark>     auto_benchmark;

  typedef typename SNC_decorator::SNC_structure          SNC_structure;
  typedef typename SNC_decorator::SNC_intersection       SNC_intersection;

  typedef typename SNC_decorator::Vertex_iterator        Vertex_iterator;
  typedef typename SNC_decorator::Vertex_handle          Vertex_handle;
  typedef typename SNC_decorator::Vertex_const_handle    Vertex_const_handle;
  typedef typename SNC_decorator::Halfedge_iterator      Halfedge_iterator;
  typedef typename SNC_decorator::Halfedge_handle        Halfedge_handle;
  typedef typename SNC_decorator::Halfedge_const_handle
                                  Halfedge_const_handle;
  typedef typename SNC_decorator::Halffacet_iterator
                                  Halffacet_iterator;
  typedef typename SNC_decorator::Halffacet_handle       Halffacet_handle;
  typedef typename SNC_decorator::Halffacet_const_handle
                                  Halffacet_const_handle;
  typedef typename SNC_decorator::Halffacet_cycle_iterator
                                  Halffacet_cycle_iterator;
  typedef typename SNC_decorator::Infi_box               Infi_box;
  typedef typename SNC_decorator::Point_3                Point_3;
  typedef typename SNC_decorator::SHalfedge_handle       SHalfedge_handle;
  typedef typename SNC_decorator::SHalfedge_iterator     SHalfedge_iterator;
  typedef typename SNC_decorator::SHalfedge_around_facet_circulator
                                  SHalfedge_around_facet_circulator;

  template< class T >
  class Nef_3_Bbox
  {
    T __min[3], __max[3];
    unsigned int __num;
    enum Type { FACET, EDGE };
    Type type;

    Halffacet_handle halffacet;
    Halfedge_handle  halfedge;


    void init ( bool complete = false ) {
        if( complete ) {
            // initialize to full box
            for( unsigned int dim = 0; dim < 3; ++dim ) {
                // this is tailored for double
                // for integers, -max() should be min()
                // TODO generalize
                __min[ dim ] = -std::numeric_limits< NumberType >::max();
                __max[ dim ] =  std::numeric_limits< NumberType >::max();
            }
        }
        else
            // initialize to empty box
            for( unsigned int dim = 0; dim < 3; ++dim ) {
                // this is tailored for double
                // for integers, -max() should be min()
                __min[ dim ] =  std::numeric_limits< NumberType >::max();
                __max[ dim ] = -std::numeric_limits< NumberType >::max();
            }
        __num = getCounter();
        //dump();
    }

    static unsigned int getCounter( bool reset = false ) {
        static unsigned int counter = 0;
        if( reset )
            counter = 0;
        else
            ++counter;
        return counter;
    }

    void extend( const Point_3& p ) {
        //dump( p );
        double x = CGAL::to_double( p.x() );
        double y = CGAL::to_double( p.y() );
        double z = CGAL::to_double( p.z() );
        if( __min[0] > x )
            __min[0] = x;
        if( __min[1] > y )
            __min[1] = y;
        if( __min[2] > z )
            __min[2] = z;

        if( __max[0] < x )
            __max[0] = x;
        if( __max[1] < y )
            __max[1] = y;
        if( __max[2] < z )
            __max[2] = z;
    }

    void finish() {
      //const double small = std::numeric_limits< double > :: min();
      const double small = 0.01;
      for( unsigned int dim = 0; dim < 3; ++dim ) {
          __max[dim] += small;
          CGAL_nef3_assertion( __min[dim] < __max[dim] );
      }
      //for( unsigned int dim = 0; dim < 3; ++dim ) {
      //    __min[dim] = 1;
      //    __max[dim] = 2;
      //}
    }

    void dump() {
      for( unsigned int dim = 0; dim < 3; ++dim )
          std::cout << " min " << __min[dim]
                    << " max " << __max[dim] << std::endl;
      std::cout << std::endl;
    }

    void dump( const Point_3& p ) {
      double x = CGAL::to_double( p.x() );
      double y = CGAL::to_double( p.y() );
      double z = CGAL::to_double( p.z() );
      std::cout << "p: " << x << " " << y << " " << z << std::endl;
    }
  public:
    typedef T NumberType;

    Nef_3_Bbox() { init(); }
    Nef_3_Bbox( Halffacet_handle f ) : type( FACET ), halffacet( f )
    {
      if( SNC_decorator::is_infbox_facet( f ) ) {
        //std::cout << "inf facet!" << std::endl;
        init( true );
      } else {
        init();
        Halffacet_cycle_iterator cycle_it = f->facet_cycles_begin();
        SHalfedge_iterator edge_it;
        if( assign( edge_it, cycle_it ) ) {
          SHalfedge_around_facet_circulator
            start( edge_it ), end( edge_it );
          CGAL_For_all( start, end ) {
            const Point_3& p =
              SNC_decorator::point(
                SNC_decorator::source(
                  SNC_decorator::previous( start ) ) );
            extend( p );
          }
        } else
          CGAL_nef3_assertion_msg(0, "is facet first cycle a SHalfloop?");
        finish();
      }
    }

    Nef_3_Bbox( Halfedge_handle e ) : type( EDGE ), halfedge( e )
    {
      if( SNC_decorator::is_infbox_vertex( SNC_decorator::source( e ) ) ||
          SNC_decorator::is_infbox_vertex( SNC_decorator::target( e ) ) )
      {
          init( true );
      } else {
          init();
          extend( SNC_decorator::point( SNC_decorator::source( e ) ) );
          extend( SNC_decorator::point( SNC_decorator::target( e ) ) );
          finish();
      }
    }

    bool has_infinite_component() {
      for( unsigned int dim = 0; dim < 3; ++dim )
          if( __min[dim] ==  std::numeric_limits< NumberType >::max() ||
              __max[dim] == -std::numeric_limits< NumberType >::max()   )
              return true;
      return false;
    }

    Halffacet_handle get_halffacet() {
      assert( type == FACET );
      return halffacet;
    }

    const Halfedge_handle get_halfedge() {
      assert( type == EDGE );
      return halfedge;
    }

    T min( unsigned int dim ) const { return __min[dim]; }
    T max( unsigned int dim ) const { return __max[dim]; }
    unsigned int num()       const { return __num;     }
  }; // end Nef_3_Bbox

  template< class _Box >
  struct Nef_3_Bbox_Adapter {
    typedef _Box Box;
    typedef typename _Box::NumberType NumberType;

    static NumberType get_lo( const Box& b, unsigned int dim )
    { return b.min( dim ); }

    static NumberType get_hi( const Box& b, unsigned int dim )
    { return b.max( dim ); }

    static unsigned int get_num( const Box& b )
    { return b.num();     }
  };


  typedef Nef_3_Bbox< double > Nef_3_Box;
  typedef Nef_3_Bbox_Adapter< Nef_3_Box > Nef_3_Box_Adapter;
  typedef Default_Box_Traits< Nef_3_Box_Adapter > Default_Box_Traits;
  typedef std::vector< Nef_3_Box > Box_Container;

  struct BOP_Callback_Environment {
      Halfedge_iterator e0, e1;
      Halffacet_iterator f0, f1;

      SNC_decorator             &d;
      const Selection  &BOP;
      Vertex_iterator  v0, v1;
      SNC_intersection &is;
      // traditional names. ask Peter Hachenbacher
      SNC_structure    &_sncp, &snc1i, &result;
      BOP_Callback_Environment( SNC_decorator &d,
                                const Selection &BOP,
                                SNC_intersection &is,
                                SNC_structure &sncp,
                                SNC_structure &snc1i,
                                SNC_structure &result)
          : d( d ), BOP( BOP ),
            is( is ), _sncp( sncp ), snc1i( snc1i ), result( result )
      {}
  };

  struct BOP_Edge0_Face1_Callback : public BOP_Callback_Environment {
    struct Pair_hash_function {
      typedef std::size_t result_type;

      template <class H>
      std::size_t
      operator() (const H& h) const {
        return
          std::size_t(&*(h.first)) / sizeof
          (typename std::iterator_traits<typename H::first_type>::value_type)
              +
          std::size_t(&*(h.second)) / sizeof
          (typename std::iterator_traits<typename H::second_type>::value_type);
      }
    };

    typedef std::pair<Halfedge_handle,
                      Halffacet_handle> Hash_key_type;


    Unique_hash_map< Hash_key_type,
                     bool,
                     Pair_hash_function > Ignore;

    BOP_Edge0_Face1_Callback( SNC_decorator& d,
                              const Selection &BOP,
                              SNC_intersection &is,
                              SNC_structure &sncp,
                              SNC_structure &snc1i,
                              SNC_structure &result)
    : BOP_Callback_Environment( d, BOP, is, sncp, snc1i, result ),
      Ignore( false )
    {}

    void operator()( Nef_3_Box& box0, Nef_3_Box& box1 ) {
      e0 = box0.get_halfedge();
      f1 = box1.get_halffacet();
      if(Ignore[ std::make_pair( e0, f1 ) ])
          return;
      if( Infi_box::degree( SNC_decorator::plane( f1 ).d() ) > 0 )
          return;
      Point_3 ip;
      ++intersection_counter;
      if( is.does_intersect_internally( SNC_decorator::segment(e0), f1, ip )) {
        //std::cout << " edge0 face1 intersection..." << std::endl;
        TRACEN(" edge0 face1 intersection...");
        ip = normalized(ip);
        v0 = d.qualify_with_respect( ip, _sncp, result );
        v1 = d.qualify_with_respect( ip, snc1i, result );
        d.binop_local_views( v0, v1, BOP, result );
        result.delete_vertex( v0 );
        result.delete_vertex( v1 );
        Ignore[ std::make_pair( SNC_decorator::twin( e0 ), f1 ) ] = true;
      }
    }
  };
  
  struct BOP_Edge1_Face0_Callback : public BOP_Callback_Environment {
      BOP_Edge1_Face0_Callback( SNC_decorator &d,
                                const Selection &BOP,
                                SNC_intersection &is,
                                SNC_structure &sncp,
                                SNC_structure &snc1i,
                                SNC_structure &result)
          : BOP_Callback_Environment( d, BOP, is, sncp, snc1i, result )
      {}
  
      void operator()( Nef_3_Box& box0, Nef_3_Box& box1 ) {
          f0 = box0.get_halffacet();
          e1 = box1.get_halfedge();
          if( Infi_box::degree( SNC_decorator::plane( f0 ).d() ) > 0 )
              return;
          Point_3 ip;
          ++intersection_counter;
          if( is.does_intersect_internally( SNC_decorator::segment( e1 ),
                                            f0, ip ) )
          {
              //std::cout << " edge1 face0 intersection..." << std::endl;
              TRACEN(" edge1 face0 intersection...");
              ip = normalized(ip);
              v1 = d.qualify_with_respect( ip, snc1i, result);
              v0 = d.qualify_with_respect( ip, _sncp,  result);
  
              d.binop_local_views( v0, v1, BOP, result);
              result.delete_vertex(v0);
              result.delete_vertex(v1);
          }
      }
  };
  
  struct BOP_Edge0_Edge1_Callback : public BOP_Callback_Environment {
    BOP_Edge0_Edge1_Callback( SNC_decorator& d,
                              const Selection& BOP,
                              SNC_intersection &is,
                              SNC_structure &sncp,
                              SNC_structure &snc1i,
                              SNC_structure &result)
        : BOP_Callback_Environment( d, BOP, is, sncp, snc1i, result )
    {}

    void operator()( Nef_3_Box& box0, Nef_3_Box& box1 ) {
      ++intersection_counter;
      Point_3 ip;
      e0 = box0.get_halfedge();
      e1 = box1.get_halfedge();
      if( is.does_intersect_internally( SNC_decorator::segment( e0 ),
                                        SNC_decorator::segment( e1 ), ip ))
      {
          //std::cout << " edge0 edge1 intersection..." << std::endl;
          TRACEN(" edge0 edge1 intersection..." << ip);
          ip = normalized(ip);
          Vertex_handle v0, v1;
          v0 = d.qualify_with_respect( ip, _sncp,  result );
          v1 = d.qualify_with_respect( ip, snc1i, result );

          d.binop_local_views( v0, v1, BOP, result );
          result.delete_vertex( v0 );
          result.delete_vertex( v1 );
      }
    }
  };

  static void dump_bboxes( Box_Container& container, char* filename ) {
    std::ofstream out( filename );
    out << container.size() << " 3" << std::endl;
    typedef typename Box_Container::iterator IT;
    for( IT it = container.begin(); it != container.end(); ++it ) {
        Nef_3_Box box = *it;
        for( unsigned int dim = 0; dim < 3; ++dim )
            out << "[" << box.min( dim ) << "," << box.max( dim ) << ") ";
        out << std::endl;
    }
  }

  virtual void operator()( SNC_decorator d,
                          const Selection& BOP,
                          SNC_structure& sncp,
                          SNC_structure& snc1i,
                          SNC_structure& result )
  {
    auto_benchmark bench;
    Halfedge_iterator e0, e1;
    Halffacet_iterator f0, f1;
    SNC_intersection is( sncp );

    Unique_hash_map<Halfedge_handle, bool> Ignore_halfedge(false);

    Box_Container
        halfedge_container, halffacet_container,
        edge_container, edge_container_2;
    BOP_Edge0_Face1_Callback
        callback_edge0_face1( d, BOP, is, sncp, snc1i, result );
    CGAL_nef3_forall_halfedges( e0, sncp )
        halfedge_container.push_back( Nef_3_Box( e0 ) );
    CGAL_nef3_forall_facets( f1, snc1i)
        halffacet_container.push_back( Nef_3_Box( f1 ) );


    //Default_Box_Traits::cutoff =
    //   (halfedge_container.size() + halffacet_container.size()) / 20;

    Default_Box_Traits::cutoff = 1;
    segment_tree( halfedge_container.begin(),  halfedge_container.end(),
                  halffacet_container.begin(), halffacet_container.end(),
                  callback_edge0_face1, Default_Box_Traits(), 2 );


    TRACEN("start edge1 face0");
    halffacet_container.clear();
    edge_container.clear();

    BOP_Edge1_Face0_Callback
        callback_edge1_face0( d, BOP, is, sncp, snc1i, result  );
    CGAL_nef3_forall_edges( e1, snc1i)
        edge_container.push_back( Nef_3_Box( e1 ) );
    CGAL_nef3_forall_facets( f0, sncp )
        halffacet_container.push_back( Nef_3_Box( f0 ) );

    segment_tree( halffacet_container.begin(), halffacet_container.end(),
                  edge_container.begin(),  edge_container.end(),
                  callback_edge1_face0, Default_Box_Traits(), 2 );
    //        SETDTHREAD(19*37);

    TRACEN("start edge0 edge1");
    TRACEN("=> edge edge intersection");
    edge_container.clear();
    edge_container_2.clear();

    BOP_Edge0_Edge1_Callback
        callback_edge0_edge1( d, BOP, is, sncp, snc1i, result );
    CGAL_nef3_forall_edges( e0, sncp)
        edge_container.push_back( Nef_3_Box( e0 ) );
    CGAL_nef3_forall_edges( e1, snc1i)
        edge_container_2.push_back( Nef_3_Box( e1 ) );

    segment_tree( edge_container.begin(),  edge_container.end(),
                  edge_container_2.begin(), edge_container_2.end(),
                  callback_edge0_edge1, Default_Box_Traits(), 2 );
  }
};

CGAL_END_NAMESPACE


#endif
