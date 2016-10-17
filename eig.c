#include "eig.h"
#define A(i,j)  ((i)*n + j)
#define eps 1.0e-8
#define eq(a,b) ((a-b<eps) && (b-a<eps))
#define jord_c 5

extern int debug;
extern int restr;

enum FUNC formula(char* str)
{
    if(strcmp(str, "symm") == 0)
        return  symm;    
    if(strcmp(str, "positive_symm") == 0)
        return positive_symm;
    if(strcmp(str, "hilbert") == 0)
        return hilbert;
	if(strcmp(str, "upper") == 0)
		return upper;
	if(strcmp(str,"disg")==0)
		return disg;
	if(strcmp(str,"jord")==0)
		return jord;
    return errf;
}

void fill(enum FUNC xin,int n, double* a)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            switch(xin) {
            case symm:
                a[A(i,j)] = abs(i-j);
                break;
            case positive_symm:
                a[A(i,j)] = 1+abs(i-j);
                break;
            case hilbert:
                a[A(i,j)] = 1./(double)(i+j+1);
                break;
			case upper:
				if(i==j) a[A(i,j)] = 1;
				if(i>j) a[A(i,j)] = 0;
				if(i<j) a[A(i,j)] = -1;
				break;
			case disg:
				if(i>j)
					a[A(i,j)] = n-i;
				else 
					a[A(i,j)] = n-j;
                break;
			case jord:
				if(i==j)
					a[A(i,j)] = 1;
				else if(i==j-1)
					a[A(i,j)] = jord_c;
				else
					a[A(i,j)] = 0;
			case errf:
				break;
              }
        }
    }
}
void printm(FILE* fout, int n, double* a)
{
  if(restr > n-2)
    {
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
               fprintf(fout, "%.4f ",a[i*n+j]);
            }
            fprintf(fout, "\n");
        }
    }
    else
    {
        for(int i = 0;i<restr;i++)
        {
            for(int j = 0; j<restr;j++)
            {
               fprintf(fout,"%.2f ",a[i*n+j]);
            }
           fprintf(fout,".. %.2f \n",a[i*n + n -1]);
        }
       fprintf(fout,"..\n");
        for(int j = 0; j<restr;j++)
        {
            fprintf(fout,"%.2f ",a[(n-1)*n+j]);
        }
       fprintf(fout,".. %.2f  \n",a[(n-1)*n + n -1]);
       
    }

}

int tridiag(int n, double* a)
{   

    double sphi, cphi, olx, oly, olbii, olbij, olbji, olbjj, sq;
    for(int j = 0;j<n-2;j++)
    {
        for(int i = j+2; i<n;i++)
        { 
		    sq = sqrt(a[A(j+1,j)]*a[A(j+1,j)] + a[A(i,j)]*a[A(i,j)]);
            if(eq(sq,0)) return 0;
            cphi = a[A(j+1,j)]/sq;
            sphi = - a[A(i,j)]/sq;
            for(int k = j; k<n;k++)
            {
                olx = a[A(j+1,k)];
                oly = a[A(i,k)];
                a[A(j+1,k)] = olx*cphi - oly*sphi;
                a[A(i,k)] = olx*sphi + oly*cphi;
            }
            olbii = a[A(j+1,j+1)];
            olbij = a[A(j+1,i)];
            olbji = a[A(i,j+1)];
            olbjj = a[A(i,i)];
            a[A(j+1,j+1)] = olbii*cphi - olbij*sphi;
            a[A(i,j+1)] = olbji*cphi - olbjj*sphi;
            a[A(j+1,i)] = olbii*sphi + olbij*cphi;
            a[A(i,i)] = olbji*sphi + olbjj*cphi;
            for(int k=j;k<n;k++)
            {
                a[A(k,j+1)] = a[A(j+1,k)];
                a[A(k,i)] = a[A(i,k)];
            }

        }
    }


    if(debug)
    {
        fprintf(stderr,"\n\nTridiagonal form:\n");
	    printm(stderr,n,a);
    }
    return 1;
}
