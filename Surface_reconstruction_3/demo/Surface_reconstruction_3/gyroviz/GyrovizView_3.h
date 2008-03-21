// GyrovizView_3.h : interface of the CGyrovizView_3 class
//

#pragma once

#include "arcball.h"
class CGyrovizDoc_3;


class CGyrovizView_3 : public CView
{
protected: // create from serialization only
	CGyrovizView_3();
	DECLARE_DYNCREATE(CGyrovizView_3)

// Attributes
public:
	CGyrovizDoc_3* GetDocument() const;

// Operations
public:

	// mouse
	bool m_LeftButtonDown;
	bool m_MiddleButtonDown;
	bool m_RightButtonDown;
	CPoint m_prev_mouse_pos;
	bool m_moving;
	
	// keyboard
	UINT m_KeyboardFlags;

	// arcball
	Arcball *m_pArcball;
	
	bool first_paint;

	// OpenGL
	HGLRC m_hGLContext;
	int m_GLPixelIndex;

	// Drawing options
	bool m_view_delaunay_triangles;
	bool m_view_delaunay_edges;	
	bool m_view_delaunay_vertices;
	bool m_view_coff_rays;
	bool m_view_coff_inside;
	
	bool m_view_arcball;

	// OpenGL specific
	BOOL SetWindowPixelFormat(HDC hDC);
	BOOL CreateViewGLContext(HDC hDC);

// Overrides
	public:
	virtual void OnDraw(CDC* pDC);  // overridden to draw this view
  virtual BOOL PreCreateWindow(CREATESTRUCT& cs);

// Implementation
public:
	virtual ~CGyrovizView_3();
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
	afx_msg void OnPaint();
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnMButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnMButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnRenderDelaunayedges();
	afx_msg void OnUpdateRenderPoints(CCmdUI *pCmdUI);
	afx_msg void OnRenderPoints();
	afx_msg void OnUpdateRenderDelaunayedges(CCmdUI *pCmdUI);
	afx_msg void OnArcballReset();
  afx_msg void OnRenderArcball();
  afx_msg void OnUpdateRenderArcball(CCmdUI *pCmdUI);
  afx_msg void OnRenderTriangles();
  afx_msg void OnUpdateRenderTriangles(CCmdUI *pCmdUI);
  afx_msg void OnRendercoffRays();
  afx_msg void OnUpdateRendercoffRays(CCmdUI *pCmdUI);
  afx_msg void OnRendercoffInside();
  afx_msg void OnUpdateRendercoffInside(CCmdUI *pCmdUI);
};

#ifndef _DEBUG  // debug version in GyrovizView_3.cpp
inline CGyrovizDoc_3* CGyrovizView_3::GetDocument() const
   { return reinterpret_cast<CGyrovizDoc_3*>(m_pDocument); }
#endif

