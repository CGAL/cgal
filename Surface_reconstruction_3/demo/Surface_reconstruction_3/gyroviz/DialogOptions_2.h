#pragma once

#include "resource.h"


// Boîte de dialogue CDialogOptions_2

class CDialogOptions_2 : public CDialog
{
	DECLARE_DYNAMIC(CDialogOptions_2)

public:
	CDialogOptions_2(CWnd* pParent = NULL);   // constructeur standard
	virtual ~CDialogOptions_2();

// Données de boîte de dialogue
	enum { IDD = IDD_DIALOG_OPTIONS_3 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // Prise en charge DDX/DDV

	DECLARE_MESSAGE_MAP()

public:
	double m_point_size; // OpenGL point size
};
