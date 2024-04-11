#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

void print_matrix(double** T, int rows, int cols);

int main(int agrc, char* agrv[]) {
    double* a0;  // auxiliary 1D for 2D matrix a
    double** a;  // 2D matrix for sequential computation

    double* a10;
    double** a1;

    int n;  // input size
    int T;  // thread number
    int i, j, k;
    int indk;
    double c, amax;

    struct timeval start_time, end_time;
    long seconds, microseconds;
    double elapsed;

    const int BLOCK_SIZE = 4;

    if (agrc == 3) {
        n = atoi(agrv[1]);
        T = atoi(agrv[2]);
        printf("The matrix size:  %d * %d \n", n, n);
        printf("The thread number:  %d\n", T);
    } else {
        printf(
            "Usage: %s n\n\n \
             n: the matrix size\n \
             T: the thread number\n\n",
            agrv[0]);
        return 1;
    }

    omp_set_num_threads(T);

    printf("Creating and initializing matrices...\n\n");
    /*** Allocate contiguous memory for 2D matrices ***/
    a0 = (double*)malloc(n * n * sizeof(double));
    a = (double**)malloc(n * sizeof(double*));

    a10 = (double*)malloc(n * n * sizeof(double));
    a1 = (double**)malloc(n * sizeof(double*));

    for (i = 0; i < n; i++) {
        a[i] = a0 + i * n;
        a1[i] = a10 + i * n;
    }

    srand(time(0));
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            a1[i][j] = a[i][j] = (double)rand() / RAND_MAX;
        }
    }

    printf("Starting sequential computation...\n\n");
    /**** Sequential computation *****/
    gettimeofday(&start_time, 0);
    for (i = 0; i < n - 1; i++) {
        // find and record k where |a(k,i)|=max|a(j,i)|
        amax = a[i][i];
        indk = i;
        for (k = i + 1; k < n; k++) {
            if (fabs(a[k][i]) > fabs(amax)) {
                amax = a[k][i];
                indk = k;
            }
        }

        // exit with a warning that a is singular
        if (amax == 0) {
            printf("matrix is singular!\n");
            exit(1);
        } else if (indk != i) {
            // swap row i and row k
            for (j = 0; j < n; j++) {
                c = a[i][j];
                a[i][j] = a[indk][j];
                a[indk][j] = c;
            }
        }

        // store multiplier in place of A(k,i)
        for (k = i + 1; k < n; k++) {
            a[k][i] = a[k][i] / a[i][i];
        }

        // subtract multiple of row a(i,:) to zero out a(j,i)
        for (k = i + 1; k < n; k++) {
            c = a[k][i];
            for (j = i + 1; j < n; j++) {
                a[k][j] -= c * a[i][j];
            }
        }
    }
    gettimeofday(&end_time, 0);

    // print the running time
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    printf("sequential calculation time: %f\n\n", elapsed);

    printf("Starting unrolling and blocking computation...\n\n");
    /**** Unrolling and block computation *****/

    int END, END_R, END_C;
    double c0, c1, c2, c3;
    double r0, r1, r2, r3;
    double* cp;
    double atop;
    int m, tmp, t;

    // clang-format off

    gettimeofday(&start_time, 0);
    for (i = 0; i < n - BLOCK_SIZE + 1; i += BLOCK_SIZE) {
        END = i + BLOCK_SIZE;
        // apply BLAS2 version of GEPP to get A(i:n , i:END)
        // update row(i:i+BLOCK_SIZE)
        for (j = i; j < END; j++) {
            // find and record k where |a(k,i)|=max|a(j,i)|
            amax = a1[j][j];
            indk = j;
            for (k = j + 1; k < n; k++) {
                if (fabs(a1[k][j]) > fabs(amax)) {
                    amax = a1[k][j];
                    indk = k;
                }
            }

            // exit with a warning that a is singular
            if (amax == 0) {
                printf("matrix is singular!\n");
                exit(1);
            } else if (indk != i) {
                // swap row i and row k with pointer
                cp = a1[j];
                a1[j] = a1[indk];
                a1[indk] = cp;
            }

            // store multiplier in place of A(k,i), update col
            atop = a1[j][j];
			#pragma omp parallel for shared(atop)
            for (k = j + 1; k < n; k++) {
                a1[k][j] = a1[k][j] / atop;
            }

            // update A(j+1:n, j:END) for next swap
			#pragma omp parallel for
            for (k = j + 1; k < n; ++k) {
                c = a1[k][j];
                for (m = j + 1; m < END; m++) {
                    a1[k][m] -= a1[j][m] * c;
                }
            }

            // print_matrix(a1, n, n);
        }

        // calculate A(i:END, END:n)
        for (j = 0; j < BLOCK_SIZE - 1; j++) {
			#pragma omp parallel for
            for (k = i + j + 1; k < END; k++) {
                tmp = i + j;
                c0 = a[k][tmp];
                for (m = END; m < n; m++) {
                    a1[k][m] -= c0 * a1[tmp][m];
                }
            }
        }

        // calculate trailing matrix in BLAS 3
		#pragma omp parallel for shared(BLOCK_SIZE, END_R, END_C) lastprivate(j)
        for (j = END; j < n - BLOCK_SIZE + 1; j += BLOCK_SIZE) {
            END_R = j + BLOCK_SIZE;
            for (k = END; k < n - BLOCK_SIZE + 1; k += BLOCK_SIZE) {
                END_C = k + BLOCK_SIZE;
                for (m = j; m < END_R; m++) {
                    c0 = a1[m][i];
                    c1 = a1[m][i + 1];
                    c2 = a1[m][i + 2];
                    c3 = a1[m][i + 3];
                    for (t = k; t < END_C; t++) {
                        a1[m][t] = a1[m][t] - c0 * a1[i][t] - c1 * a1[i + 1][t] - c2 * a1[i + 2][t] - c3 * a1[i + 3][t];
                    }
                }
            }
        }

		k = j;

        // rest col and row(< 3) BLAS 2
        for (; k < n; k++) {
            r0 = a1[i][k];
            r1 = a1[i + 1][k];
            r2 = a1[i + 2][k];
            r3 = a1[i + 3][k];

			#pragma omp parallel for shared(r0, r1, r2, r3)
            for (m = END; m < n - (n % BLOCK_SIZE ? n % BLOCK_SIZE : BLOCK_SIZE); m++) {
                a1[m][k] = a1[m][k] - r0 * a1[m][i] - r1 * a1[m][i + 1] - r2 * a1[m][i + 2] - r3 * a1[m][i + 3];
            }
        }

        for (; j < n; j++) {
            c0 = a1[j][i];
            c1 = a1[j][i + 1];
            c2 = a1[j][i + 2];
            c3 = a1[j][i + 3];

			#pragma omp parallel for shared(c0, c1, c2, c3)
            for (t = END; t < n; t++) {
                a1[j][t] = a1[j][t] - c0 * a1[i][t] - c1 * a1[i + 1][t] - c2 * a1[i + 2][t] - c3 * a1[i + 3][t];
            }
        }

    }

    for (; i < n - 1; i++) {
        // find and record k where |a(k,i)|=max|a(j,i)|
        amax = a1[i][i];
        indk = i;
        for (k = i + 1; k < n; k++) {
            if (fabs(a1[k][i]) > fabs(amax)) {
                amax = a1[k][i];
                indk = k;
            }
        }

        // exit with a warning that a is singular
        if (amax == 0) {
            printf("matrix is singular!\n");
            exit(1);
        } else if (indk != i) {
            // swap row i and row k with pointer
            cp = a1[i];
            a1[i] = a1[indk];
            a1[indk] = cp;
        }

		// #pragma omp parallel for
        for (k = i + 1; k < n; k++) {
            a1[k][i] /= a1[i][i];
        }

        // subtract multiple of row a(i,:) to zero out a(j,i)
        for (k = i + 1; k < n; k++) {
            c = a1[k][i];
			#pragma omp parallel for shared(c)
            for (j = i + 1; j < n; j++) {
                a1[k][j] -= c * a1[i][j];
            }
        }
    }
    gettimeofday(&end_time, 0);

    // print the running time
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    printf("unrolling and blocking calculation time: %f\n\n", elapsed);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fabs(a[i][j] - a1[i][j]) > 1.0E-5) {
                printf("check fail\n");
                exit(1);
            }
        }
    }
}

void print_matrix(double** T, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        // printf("{");
        for (int j = 0; j < cols; j++) {
            printf("%.2f\t", T[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}
