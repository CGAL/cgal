//== INCLUDES =================================================================

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <iostream>
#include <fstream>
#include "performance_2.h"

// Comment to use cmap with handle
#define CMAP_WITH_INDEX 1

//== CLASS DEFINITION =========================================================

typedef CGAL::Simple_cartesian<double>  MyKernelLCC;
typedef MyKernelLCC::Point_3            LCCPoint_3;
typedef MyKernelLCC::Vector_3           LCCVector_3;

typedef CGAL::Linear_cell_complex_traits<3, MyKernelLCC> MyTraitsLCC;

#ifndef CMAP_WITH_INDEX

class MyItemsLCC
{
public:

  template <class LCC>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<2, LCC> Dart; // Dim 2 == Polyhedron

    typedef CGAL::Cell_attribute_with_point<LCC, LCCVector_3, CGAL::Tag_true> Vertex;
    typedef CGAL::Cell_attribute<LCC, LCCVector_3, CGAL::Tag_true> Face;
    typedef std::tuple<Vertex,void,Face> Attributes;
  };
};

typedef CGAL::Linear_cell_complex<2, 3, MyTraitsLCC, MyItemsLCC> LCC;

#else

class MyItemsLCCWithIndex
{
public:

  template <class LCC>
  struct Dart_wrapper
  {
    typedef CGAL::Dart_for_index<2, LCC> Dart; // Dim 2 == Polyhedron

    typedef CGAL::Cell_attribute_for_index_with_point<LCC, LCCVector_3, CGAL::Tag_true> Vertex;
    typedef CGAL::Cell_attribute_for_index<LCC, LCCVector_3, CGAL::Tag_true> Face;
    typedef std::tuple<Vertex,void,Face> Attributes;
  };
};
typedef CGAL::Linear_cell_complex_for_index<2, 3, MyTraitsLCC, MyItemsLCCWithIndex> LCC;

#endif

//== CLASS DEFINITION =========================================================

class LCC_performance_2 : public Performance_test_2
{

private:

  LCC lcc;

private:
  void display_info()
  {
    lcc.display_characteristics(std::cout);
    std::cout << "\t" << std::endl;
  }

  virtual bool read_mesh(const char* _filename)
  {
    std::ifstream ifs(_filename);
    CGAL::load_off(lcc, ifs);

    for ( LCC::Dart_range::iterator dit=lcc.darts().begin(),
            dend=lcc.darts().end(); dit!=dend; ++dit )
    {
      if ( lcc.attribute<2>(dit)==NULL )
        lcc.set_attribute<2>(dit, lcc.create_attribute<2>());
    }
    assert( lcc.is_valid());
    return true;
  }


  virtual bool write_mesh(const char* _filename)
  {
    std::ofstream ofs(_filename);
    CGAL::write_off(lcc, ofs);
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


  virtual void barycenter_test(bool draw)
  {
    LCC::Vertex_attribute_range::iterator vit, vend = lcc.vertex_attributes().end();
    LCCVector_3 v(CGAL::NULL_VECTOR);
    for (vit = lcc.vertex_attributes().begin(); vit != vend; ++vit)
    {
      v = v + (vit->point() - CGAL::ORIGIN);
    }
    v = v / lcc.number_of_vertex_attributes();
    for (vit = lcc.vertex_attributes().begin(); vit != vend; ++vit)
    {
      vit->point() = vit->point() - v;
    }

    if ( draw ) std::cout<<"Barycenter: "<<v<<std::endl;
  }


  virtual void normal_test()
  {
    LCC::Vertex_attribute_range::iterator vit, vend = lcc.vertex_attributes().end();
    LCC::Attribute_range<2>::type::iterator fit, fend = lcc.attributes<2>().end();

    for (fit = lcc.attributes<2>().begin(); fit != fend; ++fit)
    {
      LCC::Dart_handle dh = fit->dart();
      LCCPoint_3& p0 = lcc.point(dh);
      dh = lcc.beta<1>(dh);
      LCCPoint_3& p1 = lcc.point(dh);
      dh = lcc.beta<1>(dh);
      LCCPoint_3& p2 = lcc.point(dh);
      LCCVector_3 n = cross_product(p0-p1, p2-p1);
      n = n / sqrt(n.squared_length());
      fit->info() = n;
    }

    for (vit = lcc.vertex_attributes().begin(); vit != vend; ++vit)
    {
      LCCVector_3 n(0,0,0);
      for( LCC::Dart_of_cell_range<0>::iterator
             vhit  = lcc.darts_of_cell<0>(vit->dart()).begin(),
             vhend = lcc.darts_of_cell<0>(vit->dart()).end();
           vhit!=vhend; ++vhit )
      {
        n = n + lcc.info<2>(vhit);
      }
      n = n / sqrt(n.squared_length());
      vit->info() = n;
    }
  }

  virtual void smoothing_test()
  {
    LCC::Vertex_attribute_range::iterator vit, vend = lcc.vertex_attributes().end();

    for (vit = lcc.vertex_attributes().begin(); vit != vend; ++vit)
    {
      bool vertex_is_border = false;
      LCCVector_3 v(0,0,0);
      float c(0);
      for( LCC::Dart_of_cell_range<0>::iterator
             vhit  = lcc.darts_of_cell<0>(vit->dart()).begin(),
             vhend = lcc.darts_of_cell<0>(vit->dart()).end();
           vhit!=vhend; ++vhit )
      {
        if (lcc.is_free<2>(vhit))
        {
          vertex_is_border = true;
          break;
        }

        v = v + (lcc.point(lcc.other_extremity(vhit)) - CGAL::ORIGIN);
        ++c;
      }
      if (!vertex_is_border)
        vit->point() = CGAL::ORIGIN + (v / c);
    }
    assert( lcc.is_valid());
  }

  void flip_edge(LCC::Dart_handle d)
  {
    LCC::Dart_handle d1 = lcc.beta<1>(d);
    LCC::Dart_handle d2 = lcc.beta<2, 1>(d);

    CGAL_assertion ( !lcc.is_free<1>(d1) && !lcc.is_free<1>(d2) );

    LCC::Dart_handle d3 = lcc.beta<1>(d1);
    LCC::Dart_handle d4 = lcc.beta<1>(d2);

    // We isolated the edge
    lcc.link_beta_1(lcc.beta<0>(d), d2);
    lcc.link_beta_1(lcc.beta<2,0>(d), d1);

    // Then we push the two extremities.
    lcc.basic_link_beta_0(d3, d);
    lcc.basic_link_beta_1(d1, lcc.beta<2>(d));
    lcc.basic_link_beta_1(d2, d);
    lcc.basic_link_beta_0(d4, lcc.beta<2>(d));

    // And we update the vertex attribute
    lcc.set_vertex_attribute_of_dart(d, lcc.vertex_attribute(d4));
    lcc.set_vertex_attribute_of_dart(lcc.beta<2>(d), lcc.vertex_attribute(d3));
  }

  virtual void subdivision_test()
  {
    int nv = lcc.number_of_vertex_attributes();

    // iterators
    LCC::Vertex_attribute_range::iterator vit, vend = lcc.vertex_attributes().end();
    LCC::Attribute_range<2>::type::iterator fit, fend = lcc.attributes<2>().end();
    LCC::Dart_range::iterator dit, dend =lcc.darts().end();

    // compute new positions of old vertices
    int i;
    std::vector<LCCPoint_3> new_pos(nv);
    for (vit = lcc.vertex_attributes().begin(), i=0; vit != vend; ++vit, ++i)
    {
      bool      vertex_is_border = false;
      LCCVector_3  v(0,0,0);
      float     n = 0;

      for( LCC::Dart_of_cell_range<0>::iterator
             vhit  = lcc.darts_of_cell<0>(vit->dart()).begin(),
             vhend = lcc.darts_of_cell<0>(vit->dart()).end();
           vhit!=vhend; ++vhit )
      {
        if (lcc.is_free<2>(vhit))
        {
          vertex_is_border = true;
          break;
        }

        v = v + (lcc.point(lcc.other_extremity(vhit)) - CGAL::ORIGIN);
        ++n;
      }
      if (!vertex_is_border)
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
      lcc.insert_barycenter_in_cell<2>(fit->dart());
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
      if (!lcc.is_free<2>(dit) && dit<lcc.beta<2>(dit) )
        flip_edge(dit);
    }

    assert( lcc.is_valid());
  }

  virtual void collapse_test()
  {
    LCC::Attribute_range<0>::type::iterator vit, vend = lcc.attributes<0>().end();
    LCC::Attribute_range<2>::type::iterator fit, fend = lcc.attributes<2>().end();

    // adjust end iterators
    --vend; --fend;

    fit = lcc.attributes<2>().begin();
    do
    {
      lcc.insert_point_in_cell<2>(fit->dart(), CGAL::ORIGIN);
    }
    while (fit++ != fend);

    // adjust end iterators
    ++vend;

    // collapse new vertices
    vit=vend; vend=lcc.attributes<0>().end();
    for (; vit!=vend; )
    {
      LCC::Dart_handle cur = vit->dart();
      ++vit;
      collapse_edge(cur);
    }
    assert( lcc.is_valid());
  }

  void contract_face(LCC::Dart_handle dh)
  {
    CGAL_assertion( dh!=lcc.null_dart_handle );
    LCC::Dart_handle d1=lcc.beta<2>(dh);
    LCC::Dart_handle d2=lcc.beta<1,2>(dh);
    CGAL_assertion(d1!=lcc.null_dart_handle &&
                   d2!=lcc.null_dart_handle);

    lcc.basic_link_beta<2>(d1, d2);
    lcc.set_dart_of_attribute<0>(lcc.vertex_attribute(d1), d1);
    lcc.set_dart_of_attribute<0>(lcc.vertex_attribute(d2), d2);

    lcc.erase_dart(lcc.beta<1>(dh));
    lcc.erase_dart(dh);
  }

  void collapse_edge(LCC::Dart_handle dh)
  {
    CGAL_assertion( dh!=lcc.null_dart_handle );
    CGAL_assertion(!lcc.is_free<2>(dh));

    LCC::Dart_handle h1=lcc.beta<0>(dh);
    LCC::Dart_handle o1=lcc.beta<2,1>(dh);

    lcc.set_attribute<0>(dh, lcc.attribute<0>(lcc.beta<2>(dh)));

    lcc.basic_link_beta_1(lcc.beta<2,0>(dh),lcc.beta<2,1>(dh));
    lcc.basic_link_beta_1(lcc.beta<0>(dh),lcc.beta<1>(dh));

    lcc.erase_dart(lcc.beta<2>(dh));
    lcc.erase_dart(dh);

    if (lcc.beta<1,1>(h1)==h1)
    {
      contract_face(h1);
    }

    if (lcc.beta<1,1>(o1)==o1)
    {
      contract_face(o1);
    }
  }
};
//=============================================================================
