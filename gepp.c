#include <math.h>
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
    int i, j, k;
    int indk;
    double c, amax;
    double* cp;

    struct timeval start_time, end_time;
    long seconds, microseconds;
    double elapsed;

    if (agrc == 2) {
        n = atoi(agrv[1]);
        printf("The matrix size:  %d * %d \n", n, n);
    } else {
        printf(
            "Usage: %s n\n\n"
            " n: the matrix size\n\n",
            agrv[0]);
        return 1;
    }

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

    //    printf("matrix a: \n");
    //    print_matrix(a, n, n);

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

    const int BLOCK_SIZE = 4;
    register double c0, c1, c2, c3;
    register double a00, a01, a02, a03;

    gettimeofday(&start_time, 0);
    for (i = 0; i < n - 1; i++) {
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

        // store multiplier in place of A(k,i)
        double aii = a1[i][i];
        for (k = i + 1; k < n - 4; k += 4) {
            a1[k][i] = a1[k][i] / aii;
            a1[k + 1][i] = a1[k + 1][i] / aii;
            a1[k + 2][i] = a1[k + 2][i] / aii;
            a1[k + 3][i] = a1[k + 3][i] / aii;
        }

        for (; k < n; k++) {
            a1[k][i] = a1[k][i] / aii;
        }

        // print_matrix(a1, n, n);



        // subtract multiple of row a(i,:) to zero out a(j,i)
        for (k = i + 1; k < n - 4; k += 4) {
            c0 = a1[k][i];
            c1 = a1[k + 1][i];
            c2 = a1[k + 2][i];
            c3 = a1[k + 3][i];

            for (j = i + 1; j < n - 4; j += 4) {
                a00 = a1[i][j];
                a01 = a1[i][j + 1];
                a02 = a1[i][j + 2];
                a03 = a1[i][j + 3];

                a1[k][j] -= c0 * a00;
                a1[k][j + 1] -= c0 * a01;
                a1[k][j + 2] -= c0 * a02;
                a1[k][j + 3] -= c0 * a03;

                a1[k + 1][j] -= c1 * a00;
                a1[k + 1][j + 1] -= c1 * a01;
                a1[k + 1][j + 2] -= c1 * a02;
                a1[k + 1][j + 3] -= c1 * a03;

                a1[k + 2][j] -= c2 * a00;
                a1[k + 2][j + 1] -= c2 * a01;
                a1[k + 2][j + 2] -= c2 * a02;
                a1[k + 2][j + 3] -= c2 * a03;

                a1[k + 3][j] -= c3 * a00;
                a1[k + 3][j + 1] -= c3 * a01;
                a1[k + 3][j + 2] -= c3 * a02;
                a1[k + 3][j + 3] -= c3 * a03;
            }

            for (; j < n; j++) {
                a00 = a1[i][j];
                a1[k][j] -= c0 * a00;
                a1[k + 1][j] -= c1 * a00;
                a1[k + 2][j] -= c2 * a00;
                a1[k + 3][j] -= c3 * a00;
            }
        }

        for (; k < n; k++) {
            c = a1[k][i];

            for (j = i + 1; j < n - 4; j += 4) {
                a1[k][j] -= c * a1[i][j];
                a1[k][j + 1] -= c * a1[i][j + 1];
                a1[k][j + 2] -= c * a1[i][j + 2];
                a1[k][j + 3] -= c * a1[i][j + 3];
            }

            for (; j < n; j++) {
                a1[k][j] -= c * a1[i][j];
            }
        }
    }
    gettimeofday(&end_time, 0);

    // print the running time
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    printf("sequential calculation time: %f\n\n", elapsed);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fabs(a[i][j] - a1[i][j]) > 1.0E-10) {
                printf("check fail");
                exit(1);
            }
        }
    }
}

void print_matrix(double** T, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%.2f   ", T[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}
