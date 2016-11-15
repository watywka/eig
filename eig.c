#include "eig.h"
#define A(i,j)  ((i)*n + j)
#define eps 1.0e-15
#define eq(a,b) ((a-b<eps) && (b-a<eps))
#define jord_c 5

extern int debug;
extern int restr;

enum FUNC formula(char* str)
{
    if(strcmp(str,"test")==0)
        return test;
    if(strcmp(str,"test2")==0)
        return test2;
    return errf;
}

void fill(enum FUNC xin,int n, double* a)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            switch(xin) {
                case test:
                    if(i == n -1)
                    {
                        a[A(i,j)] = j;
                    }
                    else if(j == n-1)
                    {
                        a[A(i,j)] = i;
                    }
                    else if(i==j)
                    {
                        a[A(i,j)] = 1;
                    }
                    else 
                    {
                        a[A(i,j)] = 0;
                case test2:
                    if(i == n -1)
                    {
                        a[A(i,j)] = j+1;
                    }
                    else if(j == n-1)
                    {
                        a[A(i,j)] = i+1;
                    }
                    else if(i==j)
                    {
                        a[A(i,j)] = 1;
                    }
                    else 
                    {
                        a[A(i,j)] = 0;
                    }
                   }
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
    fprintf(fout,"\n");

}

double norm(int n, double* a)
{
    double max = 0;
    double tmp = 0;
    for(int i = 0; i < n; i++)
    {
        max+=abs(a[A(0,i)]);
    }
    for(int j = 1; j < n; j++)
    {
        tmp = 0;
        for(int i=0;i<n;i++)
        {
            tmp+=abs(a[A(j,i)]);
        }
        if(tmp>=max) max = tmp;
    }
    return max;
}

int tridiag(int n, double* a, double e)
{

    double sphi, cphi, olx, oly, olbii, olbij, olbji, olbjj, sq, tmp1, tmp2, nor, v1, v2, norma;
    int re=0;
    int r = n;
    double D = 0, k = 0, p = 0;
    norma = norm(n,a);
    
    //
    //Privedeniye k tridiag vrascheneim
    //

    for(int j = 0;j<n-2;j++)
    {
        for(int i = j+2; i<n;i++)
        {
            sq = sqrt(a[A(j+1,j)]*a[A(j+1,j)] + a[A(i,j)]*a[A(i,j)]);
            if(eq(sq,0)) continue;
            else
            {
                cphi = a[A(j+1,j)]/sq;
                sphi = - a[A(i,j)]/sq;
            }
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

    //
    //Применение к подматрице 0,..,r алгоритма
    //

    while(r>2)
    {
        while(fabs(a[A(r-2,r-1)])>e*norma)
        {
            re++;
            for(int i = 0; i < r-1; i++)
            {
                if(eq(a[A(i+1,i)],0))
                {
                    tmp1 = 0;
                    tmp2 = 1;
                }
                else
                {
                    sq = sqrt(a[A(i,i)]*a[A(i,i)] + a[A(i+1,i)]*a[A(i+1,i)]);
                    nor = sqrt((a[A(i,i)] - sq)*(a[A(i,i)] - sq) + a[A(i+1,i)]*a[A(i+1,i)]);
                    tmp1 =(a[A(i,i)] - sq) /nor;
                    tmp2 = a[A(i+1,i)]/nor;
                }
                a[A(i,i)] = a[A(i,i)] -2*tmp1*(tmp1*a[A(i,i)]+tmp2*a[A(i+1,i)]);;
                a[A(i+1,i)] = 0;
                olx = tmp1 * a[A(i,i+1)] + tmp2 * a[A(i+1,i+1)];
                a[A(i,i+1)] = a[A(i,i+1)] - 2*tmp1*olx;
                a[A(i+1,i+1)] = a[A(i+1,i+1)] - 2*tmp2*olx;
                if (i<r-2) a[A(i+1,i+2)] = a[A(i+1,i+2)] - 2*tmp2*(tmp1 * a[A(i,i+2)] + tmp2* a[A(i+1,i+2)]);
                if (i>0)
                {
                    olx = v1*a[A(i-1,i-1)] + v2 * a[A(i-1,i)];
                    a[A(i-1,i-1)] = a[A(i-1,i-1)] - 2*v1*olx;
                    a[A(i-1,i)] = a[A(i-1,i)] - 2*v2*olx;
                    olx = v2*a[A(i,i)];
                    a[A(i,i-1)] = - 2*v1*olx;
                    a[A(i,i)] = a[A(i,i)] - 2*v2*olx;
                }
                v1 = tmp1;
                v2 = tmp2;
            }
            olx = v1*a[A(r-2,r-2)] + v2*a[A(r-2,r-1)];
            a[A(r-2,r-2)] = a[A(r-2,r-2)] - 2*v1*olx;
            olx = v2*a[A(r-1,r-1)];
            a[A(r-1,r-2)] = - 2*v1*olx;
            a[A(r-1,r-1)] = a[A(r-1,r-1)] - 2*v2*olx;
            for(int i = 0; i < r-1; i++)
            {
                a[A(i,i+1)] = a[A(i+1,i)];
            }
        }
        r--;
    }
    k = a[0]+a[n+1];
    p = a[0]*a[n+1] - a[1]*a[n];
    D = sqrt( k*k - 4 * p);
    a[1] = a[n] = 0;
    a[0] = (k + D)/2;
    a[n+1] = (k - D)/2; 
    printf("%d\n",re);
    return 1;
}
