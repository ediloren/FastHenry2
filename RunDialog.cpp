// Copyright (c) 2013 FastFieldSolvers S.r.l.
// http://www.fastfieldsolvers.com
// All Rights reserved
//
// Usage is subject to the License that you should have received
// as a text file together with these source files. The License text is also available
// in the relevant source code distribution from http://www.fastfieldsolvers.com
// or contacting FastFieldSolvers S.r.l.


// RunDialog.cpp : implementation file
//

#include "stdafx.h"
#include "FastHenry2.h"
#include "RunDialog.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CRunDialog dialog


CRunDialog::CRunDialog(CWnd* pParent /*=NULL*/)
	: CDialog(CRunDialog::IDD, pParent)
{
	//{{AFX_DATA_INIT(CRunDialog)
	m_bAutoRefinement = TRUE;
	m_sComputePorts = _T("");
	m_bDebugInfo = FALSE;
	m_iDumpFileType = 0;
	m_iDumpMatrices = 0;
	m_bExitAfterROM = FALSE;
	m_sFilenameSuffix = _T("");
	m_iGndPlane = 0;
	m_iInitialRefinement = 0;
	m_iMatrixSolution = 0;
	m_iMtxVctProduct = 0;
	m_iMultipoleOrder = 2;
	m_sPartitioningLevels = _T("auto");
	m_iPrecondMethod = 4;
	m_bRegurgitate = FALSE;
	m_iROMOrder = 0;
	m_iShellRadius = 0;
	m_iSolveIterations = 200;
	m_fATol = 0.01f;
	m_fRTol = 0.001f;
	m_iVisualizationMode = 0;
	m_sInputFileName = _T("");
	m_bPrintResults = TRUE;
	//}}AFX_DATA_INIT
}


void CRunDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CRunDialog)
	DDX_Control(pDX, IDC_INPUT_FILE_NAME, m_ctrlInputFileName);
	DDX_Check(pDX, IDC_AUTOMATIC_REFINE, m_bAutoRefinement);
	DDX_Text(pDX, IDC_COMPUTE_PORTS, m_sComputePorts);
	DDX_Check(pDX, IDC_DEBUG_INFO, m_bDebugInfo);
	DDX_CBIndex(pDX, IDC_DUMP_FILE_TYPE, m_iDumpFileType);
	DDX_CBIndex(pDX, IDC_DUMP_MATRICES, m_iDumpMatrices);
	DDX_Check(pDX, IDC_EXIT_AFTER_ROM, m_bExitAfterROM);
	DDX_Text(pDX, IDC_FILENAME_SUFFIX, m_sFilenameSuffix);
	DDX_CBIndex(pDX, IDC_GROUND_PLANE, m_iGndPlane);
	DDX_Text(pDX, IDC_INITIAL_REFINEMENT, m_iInitialRefinement);
	DDX_CBIndex(pDX, IDC_MATRIX_SOLUTION_METHOD, m_iMatrixSolution);
	DDX_CBIndex(pDX, IDC_MATRIX_VECTOR_PRODUCT, m_iMtxVctProduct);
	DDX_Text(pDX, IDC_MULTIPOLE_ORDER, m_iMultipoleOrder);
	DDX_CBString(pDX, IDC_PARTITIONING_LEVELS, m_sPartitioningLevels);
	DDX_CBIndex(pDX, IDC_PRECONDITIONER, m_iPrecondMethod);
	DDX_Check(pDX, IDC_REGURGITATE, m_bRegurgitate);
	DDX_Text(pDX, IDC_ROM_ORDER, m_iROMOrder);
	DDX_Text(pDX, IDC_SHELLS_RADIUS, m_iShellRadius);
	DDX_Text(pDX, IDC_SOLVE_ITERATIONS, m_iSolveIterations);
	DDX_Text(pDX, IDC_TOLERANCE_ATOL, m_fATol);
	DDX_Text(pDX, IDC_TOLERANCE_RTOL, m_fRTol);
	DDX_CBIndex(pDX, IDC_VISUALIZATION_MODE, m_iVisualizationMode);
	DDX_Text(pDX, IDC_INPUT_FILE_NAME, m_sInputFileName);
	DDX_Check(pDX, IDC_PRINT_RESULTS_ON_SCREEN, m_bPrintResults);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CRunDialog, CDialog)
	//{{AFX_MSG_MAP(CRunDialog)
	ON_BN_CLICKED(IDC_RESET, OnReset)
	ON_BN_CLICKED(IDC_BROWSE_INPUT_FILE, OnBrowseInputFile)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CRunDialog message handlers

BOOL CRunDialog::OnInitDialog() 
{
	int len;

	CDialog::OnInitDialog();
	
	// position caret to end of line
	UpdateData(TRUE);
	len = m_ctrlInputFileName.LineLength();
	m_ctrlInputFileName.SetSel(len, len);
	UpdateData(FALSE);

	// set controls to default values
//	OnReset();
	
	return TRUE;  // return TRUE unless you set the focus to a control
	              // EXCEPTION: OCX Property Pages should return FALSE
}

// reset dialog box controls to defaults, except input file name
void CRunDialog::OnReset() 
{
	m_bAutoRefinement = TRUE;
	m_sComputePorts = _T("");
	m_bDebugInfo = FALSE;
	m_iDumpFileType = 0;
	m_iDumpMatrices = 0;
	m_bExitAfterROM = FALSE;
	m_sFilenameSuffix = _T("");
	m_iGndPlane = 0;
	m_iInitialRefinement = 0;
	m_iMatrixSolution = 0;
	m_iMtxVctProduct = 0;
	m_iMultipoleOrder = 2;
	m_sPartitioningLevels = _T("auto");
	m_iPrecondMethod = 4;
	m_bRegurgitate = FALSE;
	m_iROMOrder = 0;
	m_iShellRadius = 0;
	m_iSolveIterations = 200;
	m_fATol = 0.01f;
	m_fRTol = 0.001f;
	m_iVisualizationMode = 0;
	m_bPrintResults = TRUE;

	UpdateData(FALSE);
}

void CRunDialog::OnBrowseInputFile() 
{
	int len;
	DWORD firstTime, reglen;
	LONG ret;
	bool isFirstTime;
	HKEY rkey;

	CFileDialog	dlgFile(TRUE, NULL, NULL, OFN_HIDEREADONLY,
		"FastHenry files (*.inp)|*.inp|All files (*.*)|*.*|");

	isFirstTime = false;

	// Load current (CURRENT_USER) FasterCap variables from the registry
	ret = RegOpenKeyEx(HKEY_CURRENT_USER, _T("Software\\FastFieldSolvers\\FastHenry2\\Settings"), 0, KEY_READ | KEY_SET_VALUE, &rkey);
	if ( ret == ERROR_SUCCESS) {
		// get the value to understand if this is the first time FasterCap is launched;
		// if the key does not exist (for the current user CU), then this is the first time
		reglen = sizeof(firstTime);
		if ( RegQueryValueEx( rkey, _T("FirstLaunch"), NULL, NULL, (unsigned char *)&firstTime, &reglen) != ERROR_SUCCESS ) {
			isFirstTime = true;
			// and record in the registry that the user already launched FasterCap once
			firstTime = 1;
			if ( RegSetValueEx( rkey, _T("FirstLaunch"), NULL, REG_DWORD,  (unsigned char *)&firstTime, sizeof(firstTime)) != ERROR_SUCCESS ) {
				AfxMessageBox("Cannot write to the Registry\nPlease check user privileges", MB_ICONSTOP);
			}
		}
	}
	RegCloseKey(rkey);

	if(isFirstTime == true ) { 
		dlgFile.m_ofn.lpstrInitialDir = m_strSamplePath;
	}

	if( dlgFile.DoModal() != IDCANCEL ) {
		// must copy into class vars the values of control,
		// not to loose their status when UpdateData(FALSE)
		UpdateData(TRUE);
		m_sInputFileName = dlgFile.GetPathName();
		UpdateData(FALSE);

		// position caret to end of line
		UpdateData(TRUE);
		len = m_ctrlInputFileName.LineLength();
		m_ctrlInputFileName.SetSel(len, len);
		UpdateData(FALSE);
	}	
}

char CRunDialog::generateCmdLine(int *argc, char ***argv)
{
	char argStr[256], *pointers[256], *strPnt;
	int len, skip, i;
	
	// initalize first argument with program name
	pointers[0] = new char[10];
	strcpy(pointers[0], "fasthenry");
	*argc = 1;
	
	// then generate all arguments
	
	// filename
	//
	// check that file name is not void or composed only of spaces
	len = sscanf((LPCTSTR)m_sInputFileName, "%s", argStr);
	if(sscanf((LPCTSTR)m_sInputFileName, "%s", argStr) != EOF) {
		pointers[*argc] = new char[ m_sInputFileName.GetLength() + 1];
		strcpy(pointers[*argc], (LPCTSTR)m_sInputFileName);
		(*argc)++;
	}
	else {
		MessageBox( "No file name given!", "Warning", MB_ICONEXCLAMATION);
		delete pointers[0];
		return 0;
	}


	if( m_iMatrixSolution == 1) {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-s");
		(*argc)++;
		pointers[*argc] = new char[10];
		strcpy(pointers[*argc], "ludecomp");
		(*argc)++;
	}

	if( m_iMtxVctProduct == 1) {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-m");
		(*argc)++;
		pointers[*argc] = new char[10];
		strcpy(pointers[*argc], "direct");
		(*argc)++;
	}

	if( m_iPrecondMethod != 4) {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-p");
		(*argc)++;

		pointers[*argc] = new char[10];
		switch(m_iPrecondMethod) {
			case 0:
				strcpy(pointers[*argc], "on");
				break;
			case 1:
				strcpy(pointers[*argc], "off");
				break;
			case 2:
				strcpy(pointers[*argc], "loc");
				break;
			case 3:
				strcpy(pointers[*argc], "posdef");
				break;
			case 5:
				strcpy(pointers[*argc], "seg");
				break;
			case 6:
				strcpy(pointers[*argc], "diag");
				break;
			case 7:
				strcpy(pointers[*argc], "shells");
				break;
		}
		(*argc)++;
	}

	if( m_iMultipoleOrder != 2) {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-o");
		(*argc)++;
		pointers[*argc] = new char[10];
		sprintf(pointers[*argc], "%i", m_iMultipoleOrder);
		(*argc)++;
	}

	if( m_sPartitioningLevels != "auto") {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-l");
		(*argc)++;
		pointers[*argc] = new char[10];
		strcpy(pointers[*argc], (LPCTSTR)m_sPartitioningLevels);
		(*argc)++;
	}

	if( m_iVisualizationMode != 0) {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-f");
		(*argc)++;

		pointers[*argc] = new char[10];
		switch(m_iVisualizationMode) {
			case 1:
				strcpy(pointers[*argc], "simple");
				break;
			case 2:
				strcpy(pointers[*argc], "refined");
				break;
			case 3:
				strcpy(pointers[*argc], "both");
				break;
			case 4:
				strcpy(pointers[*argc], "hierarchy");
				break;
		}
		(*argc)++;
	}

	if( m_iGndPlane != 0) {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-g");
		(*argc)++;

		pointers[*argc] = new char[6];
		switch(m_iGndPlane) {
			case 1:
				strcpy(pointers[*argc], "on");
				break;
			case 2:
				strcpy(pointers[*argc], "thin");
				break;
			case 3:
				strcpy(pointers[*argc], "thick");
				break;
		}
		(*argc)++;
	}

	if( m_bAutoRefinement != TRUE) {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-a");
		(*argc)++;

		pointers[*argc] = new char[4];
		strcpy(pointers[*argc], "off");
		(*argc)++;
	}

	if( m_iInitialRefinement != 0) {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-i");
		(*argc)++;

		pointers[*argc] = new char[10];
		sprintf(pointers[*argc], "%d", m_iInitialRefinement);
		(*argc)++;
	}

	if( m_iDumpMatrices != 0) {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-d");
		(*argc)++;

		pointers[*argc] = new char[10];
		switch(m_iDumpMatrices) {
			case 1:
				strcpy(pointers[*argc], "on");
				break;
			case 2:
				strcpy(pointers[*argc], "mrl");
				break;
			case 3:
				strcpy(pointers[*argc], "mzmt");
				break;
			case 4:
				strcpy(pointers[*argc], "grids");
				break;
			case 5:
				strcpy(pointers[*argc], "meshes");
				break;
			case 6:
				strcpy(pointers[*argc], "pre");
				break;
			case 7:
				strcpy(pointers[*argc], "a");
				break;
			case 8:
				strcpy(pointers[*argc], "m");
				break;
			case 9:
				strcpy(pointers[*argc], "rl");
				break;
			case 10:
				strcpy(pointers[*argc], "ls");
				break;
		}
		(*argc)++;
	}

	if( m_iDumpFileType != 0) {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-k");
		(*argc)++;

		pointers[*argc] = new char[7];
		switch(m_iDumpFileType) {
			case 1:
				strcpy(pointers[*argc], "matlab");
				break;
			case 2:
				strcpy(pointers[*argc], "text");
				break;
			case 3:
				strcpy(pointers[*argc], "both");
				break;
		}
		(*argc)++;
	}

	if( m_fRTol != 0.001f) {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-t");
		(*argc)++;

		pointers[*argc] = new char[20];
		sprintf(pointers[*argc], "%f", m_fRTol);
		(*argc)++;
	}

	if( m_fATol != 0.01f) {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-b");
		(*argc)++;

		pointers[*argc] = new char[20];
		sprintf(pointers[*argc], "%f", m_fATol);
		(*argc)++;
	}

	if( m_iSolveIterations != 200) {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-c");
		(*argc)++;

		pointers[*argc] = new char[10];
		sprintf(pointers[*argc], "%d", m_iSolveIterations);
		(*argc)++;
	}

	if( m_bDebugInfo != FALSE) {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-D");
		(*argc)++;

		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "on");
		(*argc)++;
	}

	if( m_sComputePorts != _T("")  ) {

		len = 0;
		strPnt = (char *) ((LPCTSTR)m_sComputePorts);
		while (sscanf(strPnt, "%s%n", argStr, &skip) != EOF) {
			pointers[*argc] = new char[3];
			strcpy(pointers[*argc], "-x");
			(*argc)++;

			pointers[*argc] = new char[strlen(argStr)+1];
			strcpy(pointers[*argc], argStr);
			(*argc)++;

			strPnt += skip;
		}
	}

	if( m_sFilenameSuffix != _T("") ) {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-S");
		(*argc)++;

		pointers[*argc] = new char[40];
		strcpy(pointers[*argc], (LPCTSTR)m_sFilenameSuffix);
		(*argc)++;
	}

	if( m_iROMOrder != 0) {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-r");
		(*argc)++;

		pointers[*argc] = new char[10];
		sprintf(pointers[*argc], "%d", m_iROMOrder);
		(*argc)++;
	}

	if( m_bExitAfterROM != FALSE) {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-M");
		(*argc)++;
	}

	if( m_iShellRadius != 0) {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-R");
		(*argc)++;

		pointers[*argc] = new char[10];
		sprintf(pointers[*argc], "%d", m_iShellRadius);
		(*argc)++;
	}

	if( m_bRegurgitate != FALSE) {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-v");
		(*argc)++;
	}

	if( m_bPrintResults == TRUE) {
		pointers[*argc] = new char[3];
		strcpy(pointers[*argc], "-O");
		(*argc)++;
	}

	// in the end, copy all the options to the program arguments
	*argv = new char*[*argc];
		
	for (i=0; i<*argc; i++)
		(*argv)[i] = pointers[i];

	return 1;
}
