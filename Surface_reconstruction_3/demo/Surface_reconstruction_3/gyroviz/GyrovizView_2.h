// GyrovizView_2.h : interface of the CGyrovizView_2 class


#pragma once
class CGyrovizDoc_2;

#include "GyrovizKernel.h"
class CGyrovizDoc_2;


class CGyrovizView_2 : public CView
{
protected: // create from serialization only
  CGyrovizView_2();
  DECLARE_DYNCREATE(CGyrovizView_2)

  // Attributes
public:
  CGyrovizDoc_2* GetDocument() const;

  // Operations
public:

  // OpenGL
  HGLRC m_hGLContext;
  int m_GLPixelIndex;
  BOOL SetWindowPixelFormat(HDC hDC);
  BOOL CreateViewGLContext(HDC hDC);

  // mouse
  bool m_LeftButtonDown;
  bool m_RightButtonDown;
  CPoint m_LeftDownPos;
  CPoint m_RightDownPos;
  Point_2 m_mouse_point;
  bool m_Moving;
  bool m_has_moved;
  bool m_add_constraint;

  Point_2 m_query;
  double m_xmin,m_xmax;
  double m_ymin,m_ymax;
  bool m_view_zoom_in;


  Point_2 convert(const CPoint& point);
  Point_2 m_corner1,m_corner2;
  void resize(UINT nType, int cx, int cy);

  // options
  bool m_view_all_2D_vertices;
  bool m_view_border_2D_vertices;
  bool m_view_constraints;
  bool m_view_2D_delaunay_triangulation;
  bool m_view_delaunay_edges;
  bool m_view_constrained_edges;

  bool m_view_all_cdt_edges;

  bool m_view_original;
  bool m_view_grayscale;
  bool m_view_segmentation;
  bool m_view_filtering;

  // Overrides
public:
  virtual void OnDraw(CDC* pDC);  // overridden to draw this view
  virtual BOOL PreCreateWindow(CREATESTRUCT& cs);

  // Implementation
public:
  virtual ~CGyrovizView_2();
#ifdef _DEBUG
  virtual void AssertValid() const;
  virtual void Dump(CDumpContext& dc) const;
#endif

protected:

  // Generated message map functions
protected:
  DECLARE_MESSAGE_MAP()
public:
  afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
  afx_msg void OnSize(UINT nType, int cx, int cy);
  afx_msg void OnDestroy();
  afx_msg void OnPaint();
  afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
  afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
  afx_msg void OnMouseMove(UINT nFlags, CPoint point);
  afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
  afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
  afx_msg void OnViewDelaunayedges();
  afx_msg void OnViewGenerators();
  afx_msg void OnUpdateViewGenerators(CCmdUI *pCmdUI);
  afx_msg BOOL OnEraseBkgnd(CDC* pDC);
protected:
  virtual void OnActivateView(BOOL bActivate, CView* pActivateView, CView* pDeactiveView);
public:
  afx_msg void OnEditCopy();
  afx_msg void OnViewZoom();
  afx_msg void OnUpdateViewZoom(CCmdUI *pCmdUI);
  afx_msg void OnViewZoom1();
  afx_msg void OnPslgDraw();
  afx_msg void OnViewZoomin();
  afx_msg void OnViewZoomout();
  afx_msg void OnViewSeeds();
  afx_msg void OnViewComponents32960();
  afx_msg void OnViewVoronoiedges32965();
  afx_msg void OnConstraintsAdd();
  afx_msg void OnSeedsInsert();
  afx_msg void OnUpdateSeedsInsert(CCmdUI *pCmdUI);
  virtual void OnInitialUpdate();
  afx_msg void OnInsertionSpray();
  afx_msg void OnUpdateInsertionSpray(CCmdUI *pCmdUI);
  afx_msg void OnPanLeft();
  afx_msg void OnPanRight();
  afx_msg void OnPanTop();
  afx_msg void OnPanBottom();
  afx_msg void OnViewBlindtriangles();
  afx_msg void OnUpdateViewBlindtriangles(CCmdUI *pCmdUI);
  afx_msg void OnUpdateViewDelaunayedges(CCmdUI *pCmdUI);
  afx_msg void OnViewConstrainededges();
  afx_msg void OnUpdateViewConstrainededges(CCmdUI *pCmdUI);
  afx_msg void OnUpdateViewCentroids(CCmdUI *pCmdUI);
  afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);
  afx_msg void OnViewAllcdtedges();
  afx_msg void OnUpdateViewAllcdtedges(CCmdUI *pCmdUI);
  afx_msg void OnBvdViewConstructedBVD();
  afx_msg void OnUpdateBvdViewConstructedBVD(CCmdUI *pCmdUI);
  afx_msg void OnRenderGrayscaled();
  afx_msg void OnUpdateRenderGrayscaled(CCmdUI *pCmdUI);
  afx_msg void OnRenderFiltered();
  afx_msg void OnUpdateRenderFiltered(CCmdUI *pCmdUI);
  afx_msg void OnRenderSegmented();
  afx_msg void OnUpdateRenderSegmented(CCmdUI *pCmdUI);
  afx_msg void OnRenderOriginalimage();
  afx_msg void OnUpdateRenderOriginalimage(CCmdUI *pCmdUI);
  afx_msg void OnRenderPoints();
  afx_msg void OnUpdateRenderPoints(CCmdUI *pCmdUI);
  afx_msg void OnRenderTriangles();
  afx_msg void OnUpdateRenderTriangles(CCmdUI *pCmdUI);
  afx_msg void OnRenderBordervertices();
  afx_msg void OnUpdateRenderBordervertices(CCmdUI *pCmdUI);
  afx_msg void OnRenderConstraints();
  afx_msg void OnUpdateRenderConstraints(CCmdUI *pCmdUI);
};

#ifndef _DEBUG  // debug version in GyrovizView_2.cpp
inline CGyrovizDoc_2* CGyrovizView_2::GetDocument() const
{ return reinterpret_cast<CGyrovizDoc_2*>(m_pDocument); }
#endif

