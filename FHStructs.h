// Copyright (c) 2013 FastFieldSolvers S.r.l.
// http://www.fastfieldsolvers.com
// All Rights reserved
//
// Usage is subject to the License that you should have received
// as a text file together with these source files. The License text is also available
// in the relevant source code distribution from http://www.fastfieldsolvers.com
// or contacting FastFieldSolvers S.r.l.


// C Header file with definitions common to both FastHenry and Windows
// I/O window


#ifndef FHSTRUCTS_H
#define FHSTRUCTS_H

#define FHV_BLACK	RGB(0,0,0)
#define FHV_RED		RGB(255,0,0)

// defines used with FHOnClosing function defined in FastHenryWindow.h/.cpp
#define FH_NORMAL_END 0
// flags a normal termination but with no matrices computation / output
// (e.g. only ROM model is generated)
#define FH_NO_MTX_NORMAL_END 1
#define FH_USER_BREAK 2
#define FH_GENERIC_ERROR 1000

struct impMatrix {
    long   m_lRowNum;
    long   m_lColNum;
	char   **m_sRowNames;
	char   **m_sColNames;
	long   m_lFreqNum;
	double *m_daFrequencies;
	double ***m_daMatricesReal;
	double ***m_daMatricesImag;
};

#endif
