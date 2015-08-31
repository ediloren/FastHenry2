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
