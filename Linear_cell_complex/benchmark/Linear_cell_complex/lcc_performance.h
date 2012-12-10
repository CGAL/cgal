//== INCLUDES =================================================================


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <iostream>
#include <fstream>
#include "performance.h"

#define HAS_NORMALS 1

//== CLASS DEFINITION =========================================================

typedef CGAL::Simple_cartesian<double>  MyKernelLCC;
typedef MyKernelLCC::Point_3            Point_3;
typedef MyKernelLCC::Vector_3           Vector_3;

class MyItemsLCC
{
public:

  template <class LCC>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<2, LCC> Dart; // Dim 2 == Polyhedron

#if HAS_NORMALS
    typedef CGAL::Cell_attribute_with_point<LCC, Vector_3, CGAL::Tag_true> Vertex;
    typedef CGAL::Cell_attribute<LCC, Vector_3, CGAL::Tag_true> Face;
#else
    typedef CGAL::Cell_attribute_with_point<LCC, void, CGAL::Tag_true> Vertex;
    typedef CGAL::Cell_attribute<LCC, void, CGAL::Tag_true> Face;
#endif
    typedef CGAL::cpp0x::tuple<Vertex,void,Face> Attributes;
  };
};

typedef CGAL::Linear_cell_complex_traits<3, MyKernelLCC> MyTraitsLCC;
typedef CGAL::Linear_cell_complex<2, 3, MyTraitsLCC, MyItemsLCC> LCC;
                                  
//== CLASS DEFINITION =========================================================

class LCC_performance : public Performance_test
{

private:

  LCC lcc;

private:

  virtual bool read_mesh(const char* _filename)
  {
    std::ifstream ifs(_filename);
    CGAL::load_off(lcc, ifs);
    // lcc.display_characteristics(std::cout);

    for ( LCC::Dart_range::iterator dit=lcc.darts().begin(),
            dend=lcc.darts().end(); dit!=dend; ++dit )
    {
      if ( dit->attribute<2>()==NULL )
        lcc.set_attribute<2>(dit, lcc.create_attribute<2>());        
    }
    
    return true;
  }


  virtual bool write_mesh(const char* _filename)
  {
    std::ofstream ofs(_filename);
    //    CGAL::write_off(lcc, ofs);
    return true;
  }


  virtual int circulator_test()
  {
    LCC::Vertex_attribute_range::iterator vit, vend = lcc.vertex_attributes().end();
    LCC::Attribute_range<2>::type::iterator fit, fend = lcc.attributes<2>().end();

    int counter = 0;
      
    for (vit = lcc.vertex_attributes().begin(); vit != vend; ++vit)
    {
      for( LCC::Dart_of_cell_range<0>::iterator
             vhit  = lcc.darts_of_cell<0>(vit->dart()).begin(),
             vhend = lcc.darts_of_cell<0>(vit->dart()).end();
           vhit!=vhend; ++vhit )          
      {
        ++counter;
      }
    }
      
    for (fit = lcc.attributes<2>().begin(); fit != fend; ++fit)
    {
      for( LCC::Dart_of_cell_range<2>::iterator
             fhit  = lcc.darts_of_cell<2>(fit->dart()).begin(),
             fhend = lcc.darts_of_cell<2>(fit->dart()).end();
           fhit!=fhend; ++fhit )          
      {
        --counter;
      }
    }
      
    return counter;
  }


  virtual void barycenter_test()
  {
    LCC::Vertex_attribute_range::iterator vit, vend = lcc.vertex_attributes().end();
    Vector_3 v(CGAL::NULL_VECTOR);
    for (vit = lcc.vertex_attributes().begin(); vit != vend; ++vit)
    {
      v = v + (vit->point() - CGAL::ORIGIN);
    }
    v = v / lcc.number_of_vertex_attributes();
    for (vit = lcc.vertex_attributes().begin(); vit != vend; ++vit)
    {
      vit->point() = vit->point() - v;
    }
  }


  virtual void normal_test()
  {
#if HAS_NORMALS
    LCC::Vertex_attribute_range::iterator vit, vend = lcc.vertex_attributes().end();
    LCC::Attribute_range<2>::type::iterator fit, fend = lcc.attributes<2>().end();

    for (fit = lcc.attributes<2>().begin(); fit != fend; ++fit)
    {
      LCC::Dart_handle dh = fit->dart();
      Point_3 p0 = LCC::point(dh);
      dh = dh->beta(1);
      Point_3 p1 = LCC::point(dh);
      dh = dh->beta(1);
      Point_3 p2 = LCC::point(dh);
      Vector_3 n = cross_product(p0-p1, p2-p1);
      n = n / sqrt(n.squared_length());
      fit->info() = n; // info is here the normal
    }
      
    for (vit = lcc.vertex_attributes().begin(); vit != vend; ++vit)
    {
      Vector_3 n(0,0,0);
      for( LCC::Dart_of_cell_range<0>::iterator
             vhit  = lcc.darts_of_cell<0>(vit->dart()).begin(),
             vhend = lcc.darts_of_cell<0>(vit->dart()).end();
           vhit!=vhend; ++vhit )          
      {
        n = n + vhit->attribute<2>()->info();
      }
      n = n / sqrt(n.squared_length());
      vit->info() = n; // info is here the normal
    }
#endif
  }

  virtual void smoothing_test()
  {
    LCC::Vertex_attribute_range::iterator vit, vend = lcc.vertex_attributes().end();

    for (vit = lcc.vertex_attributes().begin(); vit != vend; ++vit)
    {
      bool vertex_is_border = false;
      Vector_3 v(0,0,0);
      float c(0);
      for( LCC::Dart_of_cell_range<0>::iterator
             vhit  = lcc.darts_of_cell<0>(vit->dart()).begin(),
             vhend = lcc.darts_of_cell<0>(vit->dart()).end();
           vhit!=vhend; ++vhit )          
      {
        if (vhit->is_free(2))
        {
          vertex_is_border = true;
          break;
        }
        v = v + (LCC::point(vhit->other_extremity()) - CGAL::ORIGIN);
        ++c;
      }
      if (!vertex_is_border)
        vit->point() = CGAL::ORIGIN + (v / c);
    }
     
  }

  void flip_edge(LCC::Dart_handle d)
  {
    LCC::Dart_handle d1 = d->beta(1);
    LCC::Dart_handle d2 = d->beta(2)->beta(0);

    CGAL_assertion ( !d1->is_free(1) && !d2->is_free(0) );

    LCC::Dart_handle d3 = d1->beta(1);
    LCC::Dart_handle d4 = d2->beta(0);

    // We isolated the edge
    lcc.basic_link_beta_1(d->beta(0), d->beta(2)->beta(1));
    lcc.basic_link_beta_0(d->beta(1), d->beta(2)->beta(0));

    // Then we push the two extremities.
    lcc.basic_link_beta_0(d3, d);
    lcc.basic_link_beta_0(d2, d->beta(2));
    lcc.link_beta_1(d4, d);
    lcc.link_beta_1(d1, d->beta(2));
  }
  
  virtual void subdivision_test()
  {
    // for vector-storage we *must* reserve memory first!
    int nv = lcc.number_of_vertex_attributes();
    /*int nh = lcc.number_of_darts();
      int nf = lcc.number_of_attributes<2>();*/
    //  lcc.reserve(nv+nf, nh+6*nf, 3*nf); /: TODO implement in CMap ?

    // iterators
    LCC::Vertex_attribute_range::iterator vit, vend = lcc.vertex_attributes().end();
    LCC::Attribute_range<2>::type::iterator fit, fend = lcc.attributes<2>().end();
    LCC::Dart_range::iterator dit, dend =lcc.darts().end();
         
    // compute new positions of old vertices
    int i;
    std::vector<Point_3> new_pos(nv);
    for (vit = lcc.vertex_attributes().begin(), i=0; vit != vend; ++vit, ++i)
    {
      bool      is_border = false;
      Vector_3  v(0,0,0);
      float     n = 0;

      for( LCC::Dart_of_cell_range<0>::iterator
             vhit  = lcc.darts_of_cell<0>(vit->dart()).begin(),
             vhend = lcc.darts_of_cell<0>(vit->dart()).end();
           vhit!=vhend; ++vhit )     
      {
        if (vhit->is_free(2))
        {
          is_border = true;
          break;
        }

        v = v + (LCC::point(vhit->other_extremity()) - CGAL::ORIGIN);
        ++n;
      }

      if ( !is_border )
      {
        float alpha = (4.0 - 2.0*cos(2.0*M_PI/n)) / 9.0;
        v = (1.0f-alpha)*(vit->point() - CGAL::ORIGIN) + alpha/n*v;
        new_pos[i] = CGAL::ORIGIN + v;
      }
      else
      {
        new_pos[i] = vit->point();
      }
    }
        
    // adjust end iterators
    --vend; --fend; --dend;

    // split faces
    fit = lcc.attributes<2>().begin();
    do
    {
      Vector_3  v(CGAL::NULL_VECTOR);
      float     c(0);

      for( LCC::Dart_of_cell_range<2>::iterator
             vhit  = lcc.darts_of_cell<2>(fit->dart()).begin(),
             vhend = lcc.darts_of_cell<2>(fit->dart()).end();
           vhit!=vhend; ++vhit )          
      {
        v = v + (LCC::point(vhit) - CGAL::ORIGIN);
        ++c;
      }
      v = v / c;

      { // Here we insert a point in cell_2: equivalent to
        // lcc.insert_point_in_cell<2>(fit->dart(), CGAL::ORIGIN + v);
        // but optimized, and we create faces.
        typename LCC::Dart_handle firstdart = fit->dart();
        typename LCC::Dart_handle currentdart = firstdart;
        typename LCC::Dart_handle nextdart = currentdart->beta(1);

        typename LCC::Attribute_handle<2>::type newf;

        typename LCC::Vertex_attribute_handle
          newv = lcc.create_vertex_attribute(CGAL::ORIGIN + v);
        typename LCC::Dart_handle newdart = lcc.create_dart(newv);
        typename LCC::Dart_handle newdart2 = lcc.create_dart(nextdart->attribute<0>());

        lcc.basic_link_beta(newdart2, newdart, 2);
        lcc.basic_link_beta(newdart, nextdart, 1);
        lcc.basic_link_beta(currentdart, newdart2, 1);

        for ( currentdart=nextdart, nextdart=nextdart->beta(1);
              currentdart!=firstdart;
              currentdart=nextdart, nextdart=nextdart->beta(1) )
        {
          newf = lcc.create_attribute<2>();
          
          newdart = lcc.create_dart(newv);
          newdart2 = lcc.create_dart(nextdart->attribute<0>());
          
          lcc.set_attribute_of_dart<2>(newdart2, newf);
          lcc.set_attribute_of_dart<2>(currentdart, newf);
          lcc.set_attribute_of_dart<2>(currentdart->beta(0), newf);
          
          lcc.basic_link_beta(newdart2, newdart, 2);
          lcc.basic_link_beta(newdart, nextdart, 1);
          lcc.basic_link_beta(currentdart, newdart2, 1);
          lcc.basic_link_beta(newdart2, currentdart->beta(0), 1);
        }

        lcc.basic_link_beta(firstdart->beta(1), firstdart->beta(0), 1);
        lcc.set_attribute_of_dart<2>(firstdart->beta(0),
                                     firstdart->attribute<2>());
        lcc.set_attribute_of_dart<2>(firstdart->beta(1),
                                     firstdart->attribute<2>());
      }

    }
    while (fit++ != fend);

    // adjust end iterators
    ++vend; ++fend; ++dend;

    // set new positions of old vertices
    for (vit = lcc.vertex_attributes().begin(), i=0; vit != vend; ++vit, ++i)
      vit->point() = new_pos[i];

    // flip old edges
    for (dit = lcc.darts().begin(); dit != dend; ++dit)
    {
      if (!dit->is_free(2) && dit<dit->beta(2) )
        flip_edge(dit);
    }

    // write_mesh("reslcc.off");
  }

  virtual void collapse_test()
  {

    // reserve memory
    /*        int nv = P.size_of_vertices();
              int nh = P.size_of_halfedges();
              int nf = P.size_of_facets();
              P.reserve(nv+nf, nh+6*nf, 3*nf);

              // iterators
              Polyhedron::Vertex_iterator vit, vend = P.vertices_end();
              Polyhedron::Face_iterator   fit, fend = P.facets_end();


              // adjust end iterators
              --vend; --fend;


              // split faces (a for loop does not work for list-kernel!)
              Point_3  p(0,0,0);
              fit = P.facets_begin();
              do
              {
              Polyhedron::Halfedge_handle h = P.create_center_vertex(fit->halfedge());
              h->vertex()->point() = p;
              }
              while (++fit != fend);


              // adjust end iterators
              ++vend; ++fend;


              // collapse new edges
              vit=vend; vend=P.vertices_end();
              for (; vit!=vend; ++vit)
              halfedge_collapse(vit->halfedge()->opposite());*/
  }


  void halfedge_collapse(LCC::Dart_handle pq)
  {
    // this code is copied from the CGAL surface simplification package
    /*
      Polyhedron::Halfedge_handle qp = pq->opposite();
      Polyhedron::Halfedge_handle pt = pq->prev()->opposite();
      Polyhedron::Halfedge_handle qb = qp->prev()->opposite();

      bool lTopFaceExists         = !pq->is_border() ;
      bool lBottomFaceExists      = !qp->is_border() ;
      bool lTopLeftFaceExists     = lTopFaceExists    && !pt->is_border() ;
      bool lBottomRightFaceExists = lBottomFaceExists && !qb->is_border() ;

      Polyhedron::Vertex_handle q = pq->vertex();
      Polyhedron::Vertex_handle p = pq->opposite()->vertex();

      bool lP_Erased = false, lQ_Erased = false ;

      if ( lTopFaceExists )
      {
      if ( lTopLeftFaceExists )
      {
      P.join_facet (pt);
      }
      else
      {
      P.erase_facet(pt->opposite());

      if ( !lBottomFaceExists )
      {
      lP_Erased = true ;
      }
      }
      }

      if ( lBottomFaceExists )
      {
      if ( lBottomRightFaceExists )
      {
      P.join_facet (qb);
      }
      else
      {
      P.erase_facet(qb->opposite());

      if ( !lTopFaceExists )
      {
      lQ_Erased = true ;
      }
      }
      }

      if ( !lP_Erased && !lQ_Erased )
      {
      P.join_vertex(pq);
      lP_Erased = true ;
      }
    */
  }

};


//=============================================================================
