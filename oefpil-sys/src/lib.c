#include "lib.h"

void oefpil_tcm_blockdiag_set_tile_half(
	int n,
	int bn,
	double* A,
	int* tilemap,
	int ii,
	const double* B)
{
	(void)bn;
	A += ii * n * n;
	for (int i = 0; i < n; i++)
		for (int j = 0; j <= i; j++)
			A[j * n + i] = A[i * n + j] = *B++;
	tilemap[ii] = TILE_FULL;
}

void oefpil_tcm_full_set_tile_half(
	int n,
	int bn,
	double* A,
	int* tilemap,
	int ii,
	int jj,
	const double* B)
{
	int o = bn * n;
	double* AN = A + (ii * o + jj) * n;
	double* AT = A + (jj * o + ii) * n;
	for (int i = 0; i < n; i++)
		for (int j = 0; j <= i; j++)
			AT[j * o + i] = AN[i * o + j] = AT[i * o + j] = AN[j * o + i] = *B++;
	tilemap[jj * bn + ii] = tilemap[ii * bn + jj] = TILE_FULL;
}

FILE* stdin_file() {
	return stdin;
}
FILE* stdout_file() {
	return stdout;
}
FILE* stderr_file() {
	return stderr;
}
