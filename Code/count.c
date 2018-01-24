
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "dvbs2.h"
int main(void){
int i,j;
int no_nonzero_row[9000]={0};
for (j=0;j<9000;j++)
{	for (i=0;i<48599;i++)
	{
		if (rw_rowindex[i]==j)
		{
			no_nonzero_row[j]++;
		}
	}
}

FILE *stream;
stream = fopen( "no.h", "w" );
printf("int Nnonzero_row[] = ");
for (i=0;i<9000;i++)
 fprintf( stream, "%d, ", no_nonzero_row[i] ); //把数组值输出到文件fprintf.out中

}
