#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

void writeTempToFile(double **T, int nx, int ny, int step, int rank) {
    char filename[64];
    sprintf(filename, "T_x_y_%06d_rank%d.dat", step, rank);
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("File opening failed");
        exit(EXIT_FAILURE);
    }

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            fprintf(fp, "%lf ", T[j][i]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 2) {
        if (rank == 0) printf("Usage: %s input_file\n", argv[0]);
        MPI_Finalize();
        return -1;
    }

    FILE *input = fopen(argv[1], "r");
    if (!input) {
        perror("Input file error");
        MPI_Finalize();
        return -1;
    }

    int nx, ny;
    double x_start, x_end, y_start, y_end;
    double t_start, dt_out, dt, t_stop;
    double alpha;
    int left, right;

    fscanf(input, "%d %d", &nx, &ny);
    fscanf(input, "%lf %lf %lf %lf", &x_start, &x_end, &y_start, &y_end);
    fscanf(input, "%lf %lf %lf %lf", &t_start, &dt_out, &dt, &t_stop);
    fscanf(input, "%lf", &alpha);
    fscanf(input, "%d %d", &left, &right);
    fclose(input);

    int local_ny = ny / size;
    int remainder = ny % size;
    int start_j = rank * local_ny + (rank < remainder ? rank : remainder);
    int my_ny = local_ny + (rank < remainder ? 1 : 0);

    double dx = (x_end - x_start) / (nx - 1);
    double dy = (y_end - y_start) / (ny - 1);
    int nsteps = (int)((t_stop - t_start) / dt);

    double **T = malloc(my_ny * sizeof(double *));
    double **T_new = malloc(my_ny * sizeof(double *));
    for (int j = 0; j < my_ny; j++) {
        T[j] = malloc(nx * sizeof(double));
        T_new[j] = malloc(nx * sizeof(double));
    }

    for (int j = 0; j < my_ny; j++) {
        for (int i = 0; i < nx; i++) {
            T[j][i] = 0.0;
            T_new[j][i] = 0.0;
        }
    }

    for (int j = 0; j < my_ny; j++) {
        T[j][0] = left;
        T[j][nx - 1] = right;
        T_new[j][0] = left;
        T_new[j][nx - 1] = right;
    }

    writeTempToFile(T, nx, my_ny, 0, rank);

    double *top = malloc(nx * sizeof(double));
    double *bottom = malloc(nx * sizeof(double));

    for (int step = 1; step <= nsteps; step++) {
        if (rank > 0) MPI_Sendrecv(T[0], nx, MPI_DOUBLE, rank - 1, 0,
                                   top, nx, MPI_DOUBLE, rank - 1, 0,
                                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (rank < size - 1) MPI_Sendrecv(T[my_ny - 1], nx, MPI_DOUBLE, rank + 1, 0,
                                          bottom, nx, MPI_DOUBLE, rank + 1, 0,
                                          MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int j = 0; j < my_ny; j++) {
            for (int i = 1; i < nx - 1; i++) {
                double d2Tdx2 = (T[j][i - 1] - 2 * T[j][i] + T[j][i + 1]) / (dx * dx);
                double d2Tdy2;
                if (j == 0) {
                    if (rank == 0)
                        d2Tdy2 = (T[1][i] - 2 * T[0][i]) / (dy * dy);
                    else
                        d2Tdy2 = (top[i] - 2 * T[0][i] + T[1][i]) / (dy * dy);
                } else if (j == my_ny - 1) {
                    if (rank == size - 1)
                        d2Tdy2 = (T[j - 1][i] - 2 * T[j][i]) / (dy * dy);
                    else
                        d2Tdy2 = (T[j - 1][i] - 2 * T[j][i] + bottom[i]) / (dy * dy);
                } else {
                    d2Tdy2 = (T[j - 1][i] - 2 * T[j][i] + T[j + 1][i]) / (dy * dy);
                }
                T_new[j][i] = T[j][i] + alpha * dt * (d2Tdx2 + d2Tdy2);
            }
        }

        for (int j = 0; j < my_ny; j++) {
            for (int i = 1; i < nx - 1; i++) {
                T[j][i] = T_new[j][i];
            }
        }

        if (step % 1 == 0 || step == nsteps) {
            writeTempToFile(T, nx, my_ny, step, rank);
        }
    }

    for (int j = 0; j < my_ny; j++) {
        free(T[j]);
        free(T_new[j]);
    }
    free(T);
    free(T_new);
    free(top);
    free(bottom);

    MPI_Finalize();
    return 0;
}
