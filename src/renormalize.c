void renormalizeVector(double *z, int p, double m)
{
  int j;
  
  for(j=0;j<p;j++){
    z[j] = z[j]/m;
  }
}

