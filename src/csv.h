/*******************
 * Vanilla real valued CSV IO functions
 *
 * Work to be done asap:
 *  Header support
 *  Other types (not only doubles)
 * 
 * $Date: 2015-04-20 17:09:57 +0200 (lun 20 de abr de 2015) $
 * $Revision: 14 $
 * $Author: lrodriguez $
 * $HeadUrl:$
 * $Id: stats.h 14 2015-04-20 15:09:57Z lrodriguez $
 **/

#ifndef _CSV_H__
#define _CSV_H__

#include <stdio.h>

/**
 * Writes a real valued CSV stream
 *
 * @param [in] vals real values matrix
 * @param [in] n Number of rows
 * @param [in] m Number of cols
 * @param [in] f Field delimiter character
 * @param [in] rd Row delimiter character
 * @param [in] FILE Output stream
 *
 * @return number of characters written
 */
int writeCSV(double *vals, int n, int m, char f, char rd, FILE* s);

/**
 * Reads a real valued CSV stream
 *
 * @param [in] FILE input stream
 * @param [in] delim field delimiters
 * @param [out] m Number of cols
 * @param [out] n Number of fields
 *
 * @return real valued row-leading matrix
 */
double *readCSV(FILE* stream, char* delim, int *m, int *n);
#endif
