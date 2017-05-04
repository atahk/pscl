/* NAME: dmatTOdvec
 * DESCRIPTION: Takes a dmatrix and puts it in dvector form where 
 * each col is stored after each other. The dvector has elements 
 * vtr[0] to vtr[rows * columns -1]
 */




double*dmatTOdvec(double *vtr, double **dmtrx, int rows, int columns)
{
	int i, j, counter;
	counter = 0;
	for(j=0; j < columns; j++) {
	  for (i=0; i < rows; i++) {
	    vtr[counter] = dmtrx[i][j];
	    counter++;								
	  }
	}
	return vtr;
}
