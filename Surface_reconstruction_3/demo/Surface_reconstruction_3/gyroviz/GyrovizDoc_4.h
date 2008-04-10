// GyrovizDoc_4.h : interface of the CGyrovizDoc_4 class
#pragma once

// STL
#include <list>
#include <iostream>

// CGAL
#include "GyrovizKernel.h"
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include "Gyroviz_constrained_triangle_soup.h"

// CImg
#include <CImg.h>
using namespace cimg_library;

// This demo 
//#include "Gyroviz_dt2.h"
#include "Gyroviz_cdt2.h"

// Gyroviz's Delaunay triangulation 2-3
//typedef CGAL::Triangulation_vertex_base_with_info_2<Gyroviz_info_for_dt2,K> Vb;
//typedef CGAL::Triangulation_face_base_2<K> Fb;
//typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
//typedef Gyroviz_dt2<K,Tds>         Dt2;
//typedef Dt2::Face_handle Face_handle;
//typedef Dt2::Finite_faces_iterator Finite_faces_iterator;
//typedef Dt2::Finite_edges_iterator Finite_edges_iterator;
//typedef Dt2::Finite_vertices_iterator Finite_vertices_iterator;
typedef CGAL::Triangulation_vertex_base_with_info_2<Gyroviz_info_for_cdt2,K> CVb;
typedef CGAL::Constrained_triangulation_face_base_2<K> CFb;
typedef CGAL::Triangulation_data_structure_2<CVb,CFb>  CTds;
typedef CGAL::Exact_predicates_tag                     Itag;
typedef Gyroviz_cdt2<K,CTds, Itag>					   CDt2;
typedef Gyroviz_triangle_with_cam<CDt2>                Triangle_with_cam;
typedef Gyroviz_constrained_triangle_soup<CDt2>        Gyroviz_const_triangle_soup;


class CGyrovizDoc_4 : public CDocument
{
protected: // create from serialization only
  CGyrovizDoc_4();
  DECLARE_DYNCREATE(CGyrovizDoc_4)

  // Attributes
public:

  // Triangulation
  //CDt2 m_gyroviz_dt; // The Gyroviz equation is solved on the vertices of m_gyroviz_dt
  std::list<CDt2> m_list_cdt2;
  Gyroviz_const_triangle_soup m_vector_triangle3_wc;

// Public methods
public:

  //// Get triangulation.
  //CDt2& get_cdt2()
  //{
  //  return m_gyroviz_dt;
  //}
  //const CDt2& get_cdt2() const
  //{
  //  return m_gyroviz_dt;
  //}

	// Get m_vector_triangle...
	Gyroviz_const_triangle_soup& get_cts()
	{
		return m_vector_triangle3_wc;
	}
	const Gyroviz_const_triangle_soup& get_cts() const
	{
		return m_vector_triangle3_wc;
	}


  // Private methods
private:

  bool read_sequence(const std::string& first_image, 
  					 const std::string& first_feature_file);

  // misc status stuff
  void update_status();
  void status_message(char* fmt,...);
  double duration(const double time_init);

  // Overrides
public:
  virtual BOOL OnNewDocument();
  virtual void Serialize(CArchive& ar);

  // Implementation
public:
  virtual ~CGyrovizDoc_4();
#ifdef _DEBUG
  virtual void AssertValid() const;
  virtual void Dump(CDumpContext& dc) const;
#endif

protected:

  // Generated message map functions
protected:
  DECLARE_MESSAGE_MAP()
public:
  virtual BOOL OnOpenDocument(LPCTSTR lpszPathName);
  virtual BOOL OnSaveDocument(LPCTSTR lpszPathName);


  //unsigned char* cimg_image_multiplexer_char(const CImg <unsigned char>& image)
  //{
  //  int ix,iy,i=0;

  //  unsigned char* result = new unsigned char[image.dimx()*image.dimy()*3];

  //  for(iy=image.dimy() - 1; iy >= 0; --iy) 
  //  {
  //    for(ix=0; ix < image.dimx(); ++ix)
  //    {
  //      result[i++] = *(image.ptr() + image.dimx()*iy + ix);
  //      result[i++] = *(image.ptr() + image.dimx()*iy + ix + image.dimx()*image.dimy());
  //      result[i++] = *(image.ptr() + image.dimx()*iy + ix + image.dimx()*image.dimy()*2); 
  //    }
  //  }
  //  return result;
  //}
  unsigned char* cimg_image_multiplexer_char(const CImg <unsigned char>& image)
  {
    int ix,iy,i=0;

    unsigned char* result = new unsigned char[image.dimx()*image.dimy()*3];

    for(iy=0; iy < image.dimy(); ++iy) 
    {
      for(ix=0; ix < image.dimx(); ++ix)
      {
        result[i++] = *(image.ptr() + image.dimx()*iy + ix);
        result[i++] = *(image.ptr() + image.dimx()*iy + ix + image.dimx()*image.dimy());
        result[i++] = *(image.ptr() + image.dimx()*iy + ix + image.dimx()*image.dimy()*2); 
      }
    }
    return result;
  }

};


