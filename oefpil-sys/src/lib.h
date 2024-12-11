#ifndef OEFPIL_SYS
#define OEFPIL_SYS

#include "oefpil/oefpil.h"
#include <stdio.h>

void oefpil_tcm_blockdiag_set_tile_half(
	int n,
	int bn,
	double* A,
	int* tilemap,
	int ii,
	const double* B);

void oefpil_tcm_full_set_tile_half(
	int n,
	int bn,
	double* A,
	int* tilemap,
	int ii,
	int jj,
	const double* B);

FILE* stdin_file();
FILE* stdout_file();
FILE* stderr_file();

#endif
