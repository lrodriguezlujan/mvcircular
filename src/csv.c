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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "csv.h"

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
int writeCSV(double *vals, int n, int m, char f, char rd, FILE *s){

    int i,j;

    /* Print as tab separated CSV */
    for(i=0;i<n;i++){
        fprintf(s,"%f",vals[i*m]);
        for(j=1;j<m;j++)
            fprintf(s,"%c%f",f,vals[i*m+j]);
        putc(rd,s);
    }

    return 0;
}

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
double *readCSV(FILE *stream, char* delim, int *m, int *n){

    // Input buffer
    char buffer[10240];

    // Set number of fields/rows to 0
    *m=0;*n=0;

    // Buffer strtok

    // Return value
    int currSize=10240;
    double *ret = malloc(sizeof(double)*currSize);

    int i;

    char *tok,*aux;

    // Read line by line 
    while(fgets(buffer,10240,stream)){

        // If number of fields have not been determined yet -> count
        if(*m==0){
            // Count number of elements in the vector (start in 1 since l = #delim +1)
            for(i=0,*m=1; buffer[i] ; i++) if(buffer[i]==*delim) (*m)++;
        }

        // Realloc if this line is over currentSize
        if( ( ((*n)+1) * (*m) )>currSize){
            currSize*=2;
            ret=realloc(ret,sizeof(double)*currSize);
            //TODO: ERROR check
        }

        // Strtok (parse line)
        tok=strtok(buffer,delim);
        i=0;
        while(tok){
            ret[((*n)*(*m))+(i++)] = strtod(tok,&aux);
            tok=strtok(NULL,delim);
        }

        // If i != *m -> omit line and send warn (via stderr)
        // FIXME: Error report / return instead of stderr
        if(i!=(*m)){
            fprintf(stderr,"ERROR: missing fields at line %d\n",(*n));
        }
        else
        {
            // Accept line
            (*n)++;
        }
    }

    return ret;
}
