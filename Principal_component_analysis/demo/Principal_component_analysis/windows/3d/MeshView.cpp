// MeshView.cpp : implementation of the CMeshView class
//

#include "stdafx.h"
#include "Mesh.h"
#include "MeshDoc.h"
#include "MeshView.h"
#include ".\meshview.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CMeshView

IMPLEMENT_DYNCREATE(CMeshView, CView)

BEGIN_MESSAGE_MAP(CMeshView, CView)
	// Standard printing commands
	ON_COMMAND(ID_FILE_PRINT, CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, CView::OnFilePrintPreview)
	ON_WM_CREATE()
	ON_WM_PAINT()
	ON_WM_SIZE()
	ON_WM_ACTIVATE()
	ON_WM_ERASEBKGND()
	ON_WM_DESTROY()
	ON_WM_MOUSEMOVE()
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_RBUTTONDOWN()
	ON_WM_RBUTTONUP()
END_MESSAGE_MAP()

// CMeshView construction/destruction

CMeshView::CMeshView()
{
	// OpenGL
	m_hGLContext = NULL;
	m_GLPixelIndex = 0;

	// mouse
	m_LeftButtonDown = false;
	m_RightButtonDown = false;
	m_Moving = false;
}

CMeshView::~CMeshView()
{
}

BOOL CMeshView::PreCreateWindow(CREATESTRUCT& cs)
{
	return CView::PreCreateWindow(cs);
}

// CMeshView drawing
// OnPaint message is preferred
void CMeshView::OnDraw(CDC* /*pDC*/)
{
}


// CMeshView printing
BOOL CMeshView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// default preparation
	return DoPreparePrinting(pInfo);
}

void CMeshView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
}
void CMeshView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
}

// CMeshView diagnostics

#ifdef _DEBUG
void CMeshView::AssertValid() const
{
	CView::AssertValid();
}

void CMeshView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CMeshDoc* CMeshView::GetDocument() const // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CMeshDoc)));
	return (CMeshDoc*)m_pDocument;
}
#endif //_DEBUG

//////////////////////////////////////////////
// OPENGL
//////////////////////////////////////////////

//********************************************
// SetWindowPixelFormat
//********************************************
BOOL CMeshView::SetWindowPixelFormat(HDC hDC)
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
BOOL CMeshView::CreateViewGLContext(HDC hDC)
{
	m_hGLContext = wglCreateContext(hDC);

	if(m_hGLContext==NULL)
		return FALSE;

	if(wglMakeCurrent(hDC,m_hGLContext)==FALSE)
		return FALSE;

	return TRUE;	
}

//********************************************
// Create an OpenGL framework
//********************************************
int CMeshView::OnCreate(LPCREATESTRUCT lpCreateStruct)
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

	// lighting
	glDisable(GL_LIGHTING);

	// init camera
	InitCamera();

	return 0;
}

// init camera parameters
void CMeshView::InitCamera()
{
	// perspective
	CRect rect;
	GetClientRect(&rect);

	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);
	wglMakeCurrent(hDC,m_hGLContext);

	// set viewport and camera
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	m_Viewport.SetOrigin(0,0);
	m_Viewport.SetSize(rect.Width(),rect.Height());
	m_Camera.SetHeightAngle(45.);
	m_Camera.SetPosition(0.,0.,5.);
	m_Camera.SetOrientation(0.,1.,0.,0);
	m_Camera.SetNearDistance(.1);
	m_Camera.SetFarDistance(1000.);
	m_Viewport.glDraw();
	m_Camera.glDraw(m_Viewport);
	glDrawBuffer(GL_BACK);

	::ReleaseDC(hWnd,hDC);
}

void CMeshView::OnPaint()
{
	// device context for painting
	CPaintDC dc(this); 

	// model is stored in Document
	CMeshDoc *pDoc = (CMeshDoc *)GetDocument();

	// useful in multidoc templates
	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);
	wglMakeCurrent(dc.m_ps.hdc,m_hGLContext);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// shading option
	glShadeModel(GL_FLAT);

	// culling option
	glDisable(GL_CULL_FACE);

	// polygon mode (point, line or fill)
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

	// drawing
	glPushMatrix();

	  // setup viewpoint from current arcball
		m_Arcball.glDraw();

  	pDoc->m_mesh.gl_draw_sites(0,0,100,2.0f);
    pDoc->gl_draw_fitting_primitives();

	glPopMatrix();

	// Double buffer
	SwapBuffers(dc.m_ps.hdc);
	glFlush();

	// release
	::ReleaseDC(hWnd,hDC);
}

// resize client
void CMeshView::OnSize(UINT nType, int cx, int cy)
{
	CView::OnSize(nType, cx, cy);
	
	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);
	wglMakeCurrent(hDC,m_hGLContext);

	// set OpenGL viewport and perspective
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	m_Viewport.SetSize(cx,cy);
	m_Viewport.glDraw();
	m_Camera.glDraw(m_Viewport);
	glDrawBuffer(GL_BACK);
	::ReleaseDC(hWnd,hDC);
}

void CMeshView::OnActivateView(BOOL bActivate, 
															 CView* pActivateView, 
															 CView* pDeactiveView) 
{
	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);
	wglMakeCurrent(hDC,m_hGLContext);
	::ReleaseDC(hWnd,hDC);
	CView::OnActivateView(bActivate, 
		                    pActivateView, 
												pDeactiveView);
}

// avoid flickering
BOOL CMeshView::OnEraseBkgnd(CDC*)
{
	return TRUE;
}

// destroy rendering context
void CMeshView::OnDestroy()
{
	if(wglGetCurrentContext() != NULL)
		wglMakeCurrent(NULL,NULL);

	if(m_hGLContext != NULL)
	{
		wglDeleteContext(m_hGLContext);
		m_hGLContext = NULL;
	}
}

// handle mouse button
void CMeshView::HandleMouseButton(int x, 
																	int y)
{
  CVector3d vec = m_Arcball.Intersect(x,m_Viewport.yRes()-y,
		                                  m_Camera,m_Viewport);
  m_Arcball.EndDrag(vec);
	m_Arcball.SetMode(m_LeftButtonDown+2*m_RightButtonDown);
  vec = m_Arcball.Intersect(x,m_Viewport.yRes()-y,m_Camera,m_Viewport);
	m_Arcball.BeginDrag(vec);
}

// mouse move
void CMeshView::OnMouseMove(UINT nFlags, CPoint point)
{
	if(m_LeftButtonDown || m_RightButtonDown)
	{
		CVector3d vec = m_Arcball.Intersect(point.x,m_Viewport.yRes()-point.y,
										                    m_Camera,m_Viewport);
		m_Arcball.Motion(vec);
		m_Moving = true;
		InvalidateRect(NULL,FALSE);
	}
	else
		m_Moving = false;
}

void CMeshView::OnLButtonDown(UINT nFlags, CPoint point)
{
	SetCapture();
	m_LeftButtonDown = TRUE;
	HandleMouseButton(point.x,point.y);
	CView::OnLButtonDown(nFlags, point);
}

void CMeshView::OnLButtonUp(UINT nFlags, CPoint point)
{
	ReleaseCapture();
	m_LeftButtonDown = FALSE;
	m_Moving = false;
	HandleMouseButton(point.x,point.y);
	CView::OnLButtonUp(nFlags, point);
	InvalidateRect(NULL,FALSE);
}

void CMeshView::OnRButtonDown(UINT nFlags, CPoint point)
{
	SetCapture();
	m_RightButtonDown = TRUE;
	HandleMouseButton(point.x,point.y);
	CView::OnRButtonDown(nFlags, point);
}

void CMeshView::OnRButtonUp(UINT nFlags, CPoint point)
{
	ReleaseCapture();
	m_RightButtonDown = FALSE;
	m_Moving = false;
	HandleMouseButton(point.x,point.y);
	CView::OnRButtonUp(nFlags, point);
	InvalidateRect(NULL,FALSE);
}


void CMeshView::OnInitialUpdate()
{
	CView::OnInitialUpdate();
}

