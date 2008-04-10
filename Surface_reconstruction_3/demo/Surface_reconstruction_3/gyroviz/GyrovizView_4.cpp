// GyrovizView_4.cpp : implementation of the CGyrovizView_4 class
//

#include "stdafx.h"
#include "GyrovizView_4.h"
#include "Gyroviz.h"
#include "GyrovizDoc_4.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CGyrovizView_4

IMPLEMENT_DYNCREATE(CGyrovizView_4, CView)

BEGIN_MESSAGE_MAP(CGyrovizView_4, CView)
	ON_WM_CREATE()
	ON_WM_PAINT()
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_MBUTTONDOWN()
	ON_WM_MBUTTONUP()
	ON_WM_RBUTTONDOWN()
	ON_WM_RBUTTONUP()
	ON_WM_MOUSEMOVE()
	ON_WM_MOUSEWHEEL()
	ON_WM_SIZE()
	ON_COMMAND(ID_RENDER_DELAUNAYEDGES, OnRenderDelaunayedges)
	ON_UPDATE_COMMAND_UI(ID_RENDER_POINTS, OnUpdateRenderPoints)
	ON_COMMAND(ID_RENDER_POINTS, OnRenderPoints)
	ON_UPDATE_COMMAND_UI(ID_RENDER_DELAUNAYEDGES, OnUpdateRenderDelaunayedges)
	ON_COMMAND(ID_ARCBALL_RESET, OnArcballReset)
	ON_COMMAND(ID_RENDER_ARCBALL, OnRenderArcball)
	ON_UPDATE_COMMAND_UI(ID_RENDER_ARCBALL, OnUpdateRenderArcball)
	ON_COMMAND(ID_RENDER_TRIANGLES, &CGyrovizView_4::OnRenderTriangles)
	ON_UPDATE_COMMAND_UI(ID_RENDER_TRIANGLES, &CGyrovizView_4::OnUpdateRenderTriangles)
	ON_COMMAND(ID_RENDERCOFF_RAYS, &CGyrovizView_4::OnRendercoffRays)
	ON_UPDATE_COMMAND_UI(ID_RENDERCOFF_RAYS, &CGyrovizView_4::OnUpdateRendercoffRays)
	ON_COMMAND(ID_RENDERCOFF_INSIDE, &CGyrovizView_4::OnRendercoffInside)
	ON_UPDATE_COMMAND_UI(ID_RENDERCOFF_INSIDE, &CGyrovizView_4::OnUpdateRendercoffInside)
	ON_COMMAND(ID_RENDER_BORDERVERTICES32845, &CGyrovizView_4::OnRenderBordervertices32845)
	ON_UPDATE_COMMAND_UI(ID_RENDER_BORDERVERTICES32845, &CGyrovizView_4::OnUpdateRenderBordervertices32845)
END_MESSAGE_MAP()

// CGyrovizView_4 construction/destruction

CGyrovizView_4::CGyrovizView_4()
{
	// OpenGL
	m_hGLContext = NULL;
	m_GLPixelIndex = 0;

	// mouse
	m_LeftButtonDown   = false;
	m_MiddleButtonDown = false;
	m_RightButtonDown  = false;
	m_moving           = false;

	// arcballl
	m_pArcball  = NULL;
	first_paint = true;

	// options
	m_view_delaunay_vertices  = true;
	m_view_border_delaunay_vertices = false;
	m_view_delaunay_edges     = false;
	m_view_delaunay_triangles = false;
	m_view_arcball            = true;

	m_view_coff_rays          = false;
	m_view_coff_inside        = false;

}

CGyrovizView_4::~CGyrovizView_4()
{
	delete m_pArcball;
}

BOOL CGyrovizView_4::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CView::PreCreateWindow(cs);
}

// CGyrovizView_4 drawing

void CGyrovizView_4::OnDraw(CDC* /*pDC*/)
{
	CGyrovizDoc_4* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	// TODO: add draw code for native data here
}


// CGyrovizView_4 diagnostics

#ifdef _DEBUG
void CGyrovizView_4::AssertValid() const
{
	CView::AssertValid();
}

void CGyrovizView_4::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CGyrovizDoc_4* CGyrovizView_4::GetDocument() const // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CGyrovizDoc_4)));
	return (CGyrovizDoc_4*)m_pDocument;
}
#endif 


//********************************************
// SetWindowPixelFormat
//********************************************
BOOL CGyrovizView_4::SetWindowPixelFormat(HDC hDC)
{
	PIXELFORMATDESCRIPTOR pixelDesc;

	pixelDesc.nSize = sizeof(PIXELFORMATDESCRIPTOR);
	pixelDesc.nVersion = 1;

	pixelDesc.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL |
		PFD_DOUBLEBUFFER | PFD_STEREO_DONTCARE;

	pixelDesc.iPixelType = PFD_TYPE_RGBA;
	pixelDesc.cColorBits = 32;
	pixelDesc.cRedBits = 8;
	pixelDesc.cRedShift = 16;
	pixelDesc.cGreenBits = 8;
	pixelDesc.cGreenShift = 8;
	pixelDesc.cBlueBits = 8;
	pixelDesc.cBlueShift = 0;
	pixelDesc.cAlphaBits = 0;
	pixelDesc.cAlphaShift = 0;
	pixelDesc.cAccumBits = 64;
	pixelDesc.cAccumRedBits = 16;
	pixelDesc.cAccumGreenBits = 16;
	pixelDesc.cAccumBlueBits = 16;
	pixelDesc.cAccumAlphaBits = 0;
	pixelDesc.cDepthBits = 32;
	pixelDesc.cStencilBits = 8;
	pixelDesc.cAuxBuffers = 0;
	pixelDesc.iLayerType = PFD_MAIN_PLANE;
	pixelDesc.bReserved = 0;
	pixelDesc.dwLayerMask = 0;
	pixelDesc.dwVisibleMask = 0;
	pixelDesc.dwDamageMask = 0;

	m_GLPixelIndex = ChoosePixelFormat(hDC,&pixelDesc);
	if(m_GLPixelIndex == 0) // Choose default
	{
		m_GLPixelIndex = 1;
		if(DescribePixelFormat(hDC,m_GLPixelIndex,
			sizeof(PIXELFORMATDESCRIPTOR),&pixelDesc)==0)
			return FALSE;
	}

	if(!SetPixelFormat(hDC,m_GLPixelIndex,&pixelDesc))
		return FALSE;

	return TRUE;
}

//********************************************
// CreateViewGLContext
// Create an OpenGL rendering context
//********************************************
BOOL CGyrovizView_4::CreateViewGLContext(HDC hDC)
{
	m_hGLContext = wglCreateContext(hDC);

	if(m_hGLContext==NULL)
		return FALSE;

	if(wglMakeCurrent(hDC,m_hGLContext)==FALSE)
		return FALSE;

	return TRUE;	
}


// CGyrovizView_4 message handlers

int CGyrovizView_4::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if(CView::OnCreate(lpCreateStruct) == -1)
		return -1;

	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);

	if(SetWindowPixelFormat(hDC)==FALSE)
		return -1;

	if(CreateViewGLContext(hDC)==FALSE)
		return -1;

	// Default mode
	glEnable(GL_DEPTH_TEST);
	glClearColor(1.0f,1.0f,1.0f,1.0f);

	// activate lighting
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);
	float lpos[4] = { -.2f, .2f, .9797958971f, 0.0f };
	glLightfv(GL_LIGHT0,GL_POSITION,lpos);
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0f);

	::glEnable(GL_NORMALIZE);

	::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

	m_pArcball = new Arcball(100,100,1.0f,0.0f,0.0f,0.0f);

	return 0;
}


void CGyrovizView_4::OnPaint()
{
	// device context for painting
	CPaintDC dc(this); 

	// model is stored in Document
	CGyrovizDoc_4 *pDoc = GetDocument();

	// useful in multidoc templates
	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);
	wglMakeCurrent(dc.m_ps.hdc,m_hGLContext);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//if(first_paint)
	//{
	//	// allocate new arcball
	//	CRect rect;
	//	GetClientRect(&rect);
	//	int width = rect.Width();
	//	int height = rect.Height();

	//	Sphere arcball_sphere; 
	//	arcball_sphere = pDoc->get_cdt2().region_of_interest();

	//	float size = (float)sqrt(arcball_sphere.squared_radius());
	//	float cx = (float)arcball_sphere.center().x();
	//	float cy = (float)arcball_sphere.center().y();
	//	float cz = (float)arcball_sphere.center().z();
	//	delete m_pArcball;
	//	m_pArcball = new Arcball(width,height,size,cx,cy,cz);
	//	first_paint = false;
	//}

	// shading option
	glShadeModel(GL_FLAT);

	// test
	glMatrixMode(GL_MODELVIEW);

	// drawing
	glPushMatrix();

	::glDisable(GL_LIGHTING);

	// Draw arcball's "scene" sphere
	if ( m_view_arcball && (m_KeyboardFlags & MK_SHIFT) )
	{
		m_pArcball->setCamera_no_offset_no_scale();
		m_pArcball->draw_sphere_scene();
	}

	// setup viewpoint from current arcball
	m_pArcball->setCamera();

	// Draw arcball's "object" sphere
	if ( m_view_arcball && !(m_KeyboardFlags & MK_SHIFT) )
		m_pArcball->draw_sphere();

	//// draw Delaunay vertices
	//if(m_view_delaunay_vertices)
	//{
	//	pDoc->get_cdt2().gl_draw_soup_vertices(0,0,0,2);
	//}

	////draw border vertices
	//if(m_view_border_delaunay_vertices)
	//{
	//	pDoc->get_cdt2().gl_draw_on_border_3D_vertices(255,0,255,4);
	//}

	//// draw Delaunay edges
	//if(m_view_delaunay_edges)
	//{
	//	pDoc->get_cdt2().gl_draw_soup_constrained_edges(0,0,255,1.0f);
	//}

	// draw Delaunay triangles
	if(m_view_delaunay_triangles)
	{
		pDoc->get_cts().gl_draw_soup_triangles(0,255,0);
	}

	glPopMatrix();

	// Double buffer
	SwapBuffers(dc.m_ps.hdc);
	glFlush();

	// release hdc
	::ReleaseDC(hWnd,hDC);
}

void CGyrovizView_4::OnLButtonDown(UINT nFlags, CPoint point)
{
	SetCapture();
	m_LeftButtonDown = TRUE;
	m_prev_mouse_pos = point;
	if (m_KeyboardFlags != nFlags) // update keyboard status and redraw arcball sphere
	{
		m_KeyboardFlags = nFlags;
		InvalidateRect(NULL,FALSE);
	}
	CView::OnLButtonDown(nFlags, point);
}

void CGyrovizView_4::OnLButtonUp(UINT nFlags, CPoint point)
{
	ReleaseCapture();
	m_LeftButtonDown = FALSE;
	m_moving = false;
	if (m_KeyboardFlags != nFlags) // update keyboard status and redraw arcball sphere
	{
		m_KeyboardFlags = nFlags;
		InvalidateRect(NULL,FALSE);
	}
	CView::OnLButtonUp(nFlags, point);
}

void CGyrovizView_4::OnMButtonDown(UINT nFlags, CPoint point)
{
	SetCapture();
	m_MiddleButtonDown = TRUE;
	m_prev_mouse_pos = point;
	if (m_KeyboardFlags != nFlags) // update keyboard status and redraw arcball sphere
	{
		m_KeyboardFlags = nFlags;
		InvalidateRect(NULL,FALSE);
	}
	CView::OnMButtonDown(nFlags, point);
}

void CGyrovizView_4::OnMButtonUp(UINT nFlags, CPoint point)
{
	ReleaseCapture();
	m_MiddleButtonDown = FALSE;
	m_moving = false;
	if (m_KeyboardFlags != nFlags) // update keyboard status and redraw arcball sphere
	{
		m_KeyboardFlags = nFlags;
		InvalidateRect(NULL,FALSE);
	}
	CView::OnMButtonUp(nFlags, point);
}

void CGyrovizView_4::OnRButtonDown(UINT nFlags, CPoint point)
{
	SetCapture();
	m_RightButtonDown = TRUE;
	m_prev_mouse_pos = point;
	if (m_KeyboardFlags != nFlags) // update keyboard status and redraw arcball sphere
	{
		m_KeyboardFlags = nFlags;
		InvalidateRect(NULL,FALSE);
	}
	CView::OnRButtonDown(nFlags, point);
}

void CGyrovizView_4::OnRButtonUp(UINT nFlags, CPoint point)
{
	ReleaseCapture();
	m_RightButtonDown = FALSE;
	m_moving = false;
	if (m_KeyboardFlags != nFlags) // update keyboard status and redraw arcball sphere
	{
		m_KeyboardFlags = nFlags;
		InvalidateRect(NULL,FALSE);
	}
	CView::OnRButtonUp(nFlags, point);
}

// Mouse move callback
void CGyrovizView_4::OnMouseMove(UINT nFlags, CPoint point)
{
	if (m_KeyboardFlags != nFlags) // update keyboard status and redraw arcball sphere
	{
		m_KeyboardFlags = nFlags;
		InvalidateRect(NULL,FALSE);
	}

	// Drag with middle button (or both left and right buttons): zoom
	if( m_MiddleButtonDown || (m_LeftButtonDown && m_RightButtonDown) )
	{
		m_moving = true;
		int dy = point.y - m_prev_mouse_pos.y;
		float factor = (float)(1.0 - dy * 0.01);
		m_pArcball->zoom(factor);
		m_prev_mouse_pos = point;
		InvalidateRect(NULL,FALSE);
	}
	// Drag with left button: rotate
	else if(m_LeftButtonDown)
	{
		m_moving = true;
		if (m_KeyboardFlags & MK_SHIFT)
			m_pArcball->rotate_scene(m_prev_mouse_pos.x,m_prev_mouse_pos.y,point.x,point.y);
		else
			m_pArcball->rotate_obj(m_prev_mouse_pos.x,m_prev_mouse_pos.y,point.x,point.y);
		m_prev_mouse_pos = point;
		InvalidateRect(NULL,FALSE);
	}
	// Drag withn right button: translate
	else if(m_RightButtonDown)
	{
		m_moving = true;
		m_pArcball->translate(m_prev_mouse_pos.x,m_prev_mouse_pos.y,point.x,point.y);
		m_prev_mouse_pos = point;
		InvalidateRect(NULL,FALSE);
	}
	else
		m_moving = false;

	CView::OnMouseMove(nFlags, point);
}

// Mouse wheel callback: zoom
BOOL CGyrovizView_4::OnMouseWheel(UINT nFlags, short zDelta, CPoint pt)
{
	m_KeyboardFlags = nFlags;
	float nb = (float)zDelta / 120.0f;
	float factor = (float)(1.0 + 0.3 * nb);
	m_pArcball->zoom(factor);
	InvalidateRect(NULL,FALSE);
	return CView::OnMouseWheel(nFlags, zDelta, pt);
}


void CGyrovizView_4::OnSize(UINT nType, int cx, int cy)
{
	CView::OnSize(nType, cx, cy);

	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);
	wglMakeCurrent(hDC,m_hGLContext);

	glMatrixMode(GL_MODELVIEW);
	m_pArcball->resize(cx,cy);
	glDrawBuffer(GL_BACK);
	::ReleaseDC(hWnd,hDC);
}

void CGyrovizView_4::OnRenderTriangles()
{
	m_view_delaunay_triangles = !m_view_delaunay_triangles;
	InvalidateRect(NULL,FALSE);
}
void CGyrovizView_4::OnUpdateRenderTriangles(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_view_delaunay_triangles);
}

void CGyrovizView_4::OnRenderDelaunayedges()
{
	m_view_delaunay_edges = !m_view_delaunay_edges;
	InvalidateRect(NULL,FALSE);
}
void CGyrovizView_4::OnUpdateRenderDelaunayedges(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_view_delaunay_edges);
}

void CGyrovizView_4::OnRenderPoints()
{
	m_view_delaunay_vertices = !m_view_delaunay_vertices;
	InvalidateRect(NULL,FALSE);
}
void CGyrovizView_4::OnUpdateRenderPoints(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_view_delaunay_vertices);
}

void CGyrovizView_4::OnArcballReset()
{
	m_pArcball->reset();
	InvalidateRect(NULL,FALSE);
}

void CGyrovizView_4::OnRenderArcball()
{
	m_view_arcball = !m_view_arcball; 
	InvalidateRect(NULL,FALSE);
}
void CGyrovizView_4::OnUpdateRenderArcball(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_view_arcball);
}


void CGyrovizView_4::OnRendercoffRays()
{
	m_view_coff_rays = !m_view_coff_rays; 
	InvalidateRect(NULL,FALSE);
}

void CGyrovizView_4::OnUpdateRendercoffRays(CCmdUI *pCmdUI)
{      
	//if (!pDoc->is_dt2())
	pCmdUI->SetCheck(m_view_coff_rays);

	//else 
	//pCmdUI->Enable(!m_view_coff_rays);
}



void CGyrovizView_4::OnRendercoffInside()
{
	m_view_coff_inside = !m_view_coff_inside; 
	InvalidateRect(NULL,FALSE);
}

void CGyrovizView_4::OnUpdateRendercoffInside(CCmdUI *pCmdUI)
{
	//  if (!pDoc->is_dt2())
	pCmdUI->SetCheck(m_view_coff_inside);
	//else 
	//pCmdUI->Enable(!m_view_coff_inside);
}

void CGyrovizView_4::OnRenderBordervertices32845()
{
	m_view_border_delaunay_vertices = !m_view_border_delaunay_vertices;
	InvalidateRect(NULL,FALSE);
}

void CGyrovizView_4::OnUpdateRenderBordervertices32845(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_view_border_delaunay_vertices);
}
