#ifndef EIG_H
#define EIG_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

enum FUNC
{
    symm,
    positive_symm,
    hilbert,
	upper,
	disg,
	jord,
    errf
};
extern int debug;
enum FUNC formula(char* str);
void fill(enum FUNC xin,int n, double* a);
double norm(int n, double* a);
int tridiag(int n, double* a, double e);
void printm(FILE* fout, int n, double* a);

#endif
