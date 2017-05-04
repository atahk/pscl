void usage();

double check(double **data, int **ok, int n, int m);

double dtnorm(const double mu, const double sd, const double y);

void updatex(double **ystar, int **ok, double **beta, 
	     double **x, double **xp, double **xpv,
	     int n, int m, int d, 
	     int impute);
void makexreg(double **xreg, double **x, int n, int d, int q);

void updateb(double **ystar, int **ok, double **beta, double **xreg,
	     double **bp, double **bpv,
	     int n, int m, int d, int impute);

void updatey(double **ystar, double **y, double **x, double **beta,
	       int n, int m, int d, int iter);

void choldc(double **a, int n, double p[]);
void xchol(double **aorig, double **chol, int n, double *p, double **a);
void printmat(double **mat, int nr, int nc);

void rmvnorm(double *theta, double *mu, double **sigma, int k,
	     double *xprod, double **chol, double *z, double *p, double **a);
void rmvnorm_m(double **theta, double *mu, double **chol, int k);

void crossprod(double **x, int n, int p, double **xpx);
void crossprodslow(double **x, int n, int p, double **xpx);
void crossxy(double **x, double *y, int n, int k, double *xpy);
void crossxyd(double **x, double *y, int n, int k, double *xpy);
void crossxyi(double **beta, double **w, int m, int d, int p, double *bpw);
void crossxyj(double **x, double **y, int n, int k, int p, double *xpy);
void crossall(double **x, double **ystar, int n, int d, int j, 
	      double **xpx, double *xpy);
void crosscheck(double **x, double **ystar, int **ok,
		int n, int d, int j, 
		double **xpx, double *xpy);
void crosscheckx(double **beta, double **w, int **ok,
		 int m, int d, int i, 
		 double **bpb, double *bpw);

void gaussj(double **a, int n, double *b, int m);

void bayesreg(double **xpx, double *xpy, 
	      double *bp, double **priormat,
	      double *bpost, double **vpost,
	      int p);

void bayesregFull(double **xpx, double *xpy,
		  double sd,
		  double *bp, double **priormat,
		  double *bpost, double **vpost,
		  int p);

void renormalizeVector(double *z, int p, double m);

double r_sd(double s, double df);
