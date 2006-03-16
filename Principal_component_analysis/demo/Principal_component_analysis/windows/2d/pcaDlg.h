// pcaDlg.h : header file
//

#pragma once


#include <CGAL/Basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_2.h>



// kernel
typedef CGAL::Cartesian<double> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Line_2 Line_2;


#include <vector>

// CpcaDlg dialog
class CpcaDlg : public CDialog
{
// Construction
public:
	CpcaDlg(CWnd* pParent = NULL);	// standard constructor

  // data set
	std::vector<Point_2> m_points;
  Point_2 m_centroid;
  Line_2 m_fitting_line;
  FT m_quality;

	// OpenGL
	HGLRC m_hGLContext;
	int m_GLPixelIndex;
	BOOL SetWindowPixelFormat(HDC hDC);
	BOOL CreateViewGLContext(HDC hDC);
	void render();

	Point_2 convert(const CPoint& point);
  bool inside(const CPoint& point);

  void draw_disc(const Point_2& c,const FT radius,unsigned char r,
                 unsigned char g,unsigned char b);

	// mouse
	CPoint m_RightDownPos;
	CPoint m_LeftDownPos;
	BOOL m_RightButtonDown;
	BOOL m_LeftButtonDown;

// Dialog Data
	enum { IDD = IDD_PCA_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support


// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnDestroy();
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnSetClear();
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
  afx_msg void OnFitLine();
  afx_msg void OnRandomHorizontalline();
  afx_msg void OnRandomVerticalline();
  afx_msg void OnFitDebug();
  afx_msg void OnDebugManytestsfornumerouspointsonaline();
};
