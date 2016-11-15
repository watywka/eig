#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "eig.h"

int debug;
int restr;

int main(int argc, char* argv[])
{
    int n = 10;
    double trdif=0;
    double atr = 0, tatr=0;
    double alength = 0, talength = 0, ldiff =0;
    double  epsillon=0.0000001;
    double  *a, *ta;
    int opt;
    int fflag = 0,xflag = 0;
    struct timespec begin, end;
    FILE* fin = NULL;
    enum FUNC xin;
    opterr = 0;
    restr = 5;
    debug = 0;
    if(argc == 1)
    {
        fprintf(stderr,"Usage: eig [options] \n   -f   file with the matrix\n OR\n   -x   name of formula\n   -d   debug mode\n   -r   output restrictions\n   -n   size of matrix\n   -e   precision\n");
        return 1;
    }
    while((opt = getopt(argc,argv,"dx:f:r:n:e:")) != -1)
    {
        switch(opt) {
            case 'n':
                if(sscanf(optarg,"%d",&n) != 1)
                {
                    fprintf(stderr,"Cannot read matrix size\n");
                    return -1;
                }
                if(n<1)
                {
                    fprintf(stderr,"Wrong matrix size\n");
                    return -1;
                }
                break;
            case 'd':
                debug=1;
                break;
            case 'e':
                if(sscanf(optarg,"%lf",&epsillon) != 1)
                {
                    fprintf(stderr,"Cannot read precision value\n");
                    return -1;
                }
                if(epsillon<0)
                {
                    fprintf(stderr,"Negative precision value\n");
                    return -1;
                }
                break;
            case 'x':
                if(xflag == 1)
                {
                    fprintf(stderr,"Multiple usage of -x option is not allowed\n");
                    return -1;
                }
                if(fflag == 1)
                {
                    fprintf(stderr,"-x and -f options cannot be used simultaneously\n");
                    fclose(fin);
                    return -1;
                }
                xflag = 1;
                xin = formula(optarg);
                if(xin == errf)
                {
                    fprintf(stderr,"Unknown formula\n");
                    return -1;
                }
                break;
            case 'f':
                if(fflag == 1)
                {
                    fprintf(stderr,"Multiple usage of -f option is not allowed\n");
                    fclose(fin);
                    return -1;
                }
                if(xflag == 1)
                {
                    fprintf(stderr,"-x and -f options cannot be used simultaneously\n");
                    return -1;
                }
                fflag = 1;
                fin = fopen(optarg,"r");
                if(!fin) {
                    fprintf(stderr,"File not found\n");
                    return -1;
                }
                break;
            case 'r':
                if( sscanf(optarg,"%d",&restr) != 1) {
                    fprintf(stderr,"Wrong option usage\n");
                    return -1;
                }
                if(restr<1)
                {
                    fprintf(stderr,"Restriction must be positive\n");
                    return -1;
                }
                break;
            case '?':
                if((optopt == 'r') || (optopt == 'f') || (optopt == 'x') || (optopt == 'n') || (optopt == 'e')) {
                    fprintf(stderr,"-%c requires an argument\n",optopt);
                    return -1;
                }
                fprintf(stderr, "%c: unknown option\n", optopt);
                return -1;
                break;
            default:
                return -1;
        }
    }
    if(optind<argc)
    {
        fprintf(stderr,"Wrong options format\n");
        if(fflag)
            fclose(fin);
        return -1;
    }
    if( (!fflag) && (!xflag))
    {
        fprintf(stderr,"No file or formula provided\n");
        return -1;
    }
    if(xflag)
    {
        if(!(a=(double*) malloc(sizeof(double)*n*n))  ||(!(ta =(double*) malloc(sizeof(double)*n*n))))
        {
            if(a) free(a);
            if(ta) free(ta);
            fprintf(stderr,"Can't allocate memory\n");
            return -1;
        }
        fill(xin, n, a);
        fill(xin, n, ta);
    }
    if(fflag)
    {
        if(fscanf(fin,"%d", &n)!=1)
        {
            fprintf(stderr,"Corrupted file\n");
            fclose(fin);
            return -1;
        }
        if(n<1)
        {
            fprintf(stderr,"Wrong matrix size\n");
            fclose(fin);
            return -1;
        }
        if(!(a=(double*) malloc(sizeof(double)*n*n)) ||(!(ta =(double*) malloc(sizeof(double)*n*n))))
        {
            if(a) free(a);
            if(ta) free(ta);
            fprintf(stderr,"Can't allocate memory\n");
            fclose(fin);
            return -1;
        }
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                if(fscanf(fin,"%lf", &a[i*n+j])!=1)
                {
                    fprintf(stderr,"Corrupted file\n");
                    fclose(fin);
                    free(ta);
                    free(a);
                    return -1;
                }
                ta[i*n+j]=a[i*n+j];
            }
        }
        fclose(fin);
    }
    printf("Matrix:\n");
    printm(stdout,n,a);
    clock_gettime(CLOCK_MONOTONIC,&begin);
    fprintf(stdout,"\n");
    if(! tridiag(n,a,epsillon)) fprintf(stderr,"Unable\n");
    printm(stdout,n,a);
    clock_gettime(CLOCK_MONOTONIC, &end);
    printf("\nTime spent:%fs\n",((end.tv_sec-begin.tv_sec)+(double)(end.tv_nsec-begin.tv_nsec)/1000000000));
    
    for (int i=0;i<n;i++)
    {
        atr +=a[i*n+i];
        tatr += ta[i*n+i];
        for(int j=0;j<n;j++)
        {
            alength += a[i*n+j]*a[i*n+j];
            talength += ta[i*n+j]*ta[i*n+j];
        }
    }
    ldiff = alength - talength;
    trdif = atr - tatr;
    printf("trace A: %lf\n eigenvalues sum: %lf \n diff: %lf \n", atr, tatr, trdif);
    printf("length A: %lf\n new length: %lf \n diff: %lf \n", alength, talength, ldiff);
    free(ta);
    free(a);
    return 1;

}
