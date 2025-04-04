#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void writeTempToFile(double **T, int nx, int ny, int step) {
    char filename[64];
    sprintf(filename, "T_x_y_%06d.dat", step);
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
    printf("Saved %s\n", filename);
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Usage: %s input_file\n", argv[0]);
        return -1;
    }

    FILE *input = fopen(argv[1], "r");
    if (!input) {
        perror("Input file error");
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

    printf("Inputs are: %d %lf %lf %lf %lf %lf\n", nx, x_start, x_end, y_start, dt, t_stop);
    printf("Inputs are: %d %lf %lf\n", ny, y_start, y_end);

    double dx = (x_end - x_start) / (nx - 1);
    double dy = (y_end - y_start) / (ny - 1);
    int nsteps = (int)((t_stop - t_start) / dt);
    int output_interval = 1; // write every step

    printf("dx: %f, dy: %f\n", dx, dy);
    printf("dt: %f, nsteps: %d\n", dt, nsteps);

    double **T = malloc(ny * sizeof(double *));
    double **T_new = malloc(ny * sizeof(double *));
    for (int j = 0; j < ny; j++) {
        T[j] = malloc(nx * sizeof(double));
        T_new[j] = malloc(nx * sizeof(double));
    }

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            T[j][i] = 0.0;
            T_new[j][i] = 0.0;
        }
    }

    for (int j = 0; j < ny; j++) {
        T[j][0] = left;
        T[j][nx - 1] = right;
        T_new[j][0] = left;
        T_new[j][nx - 1] = right;
    }

    writeTempToFile(T, nx, ny, 0);

    for (int step = 1; step <= nsteps; step++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                double d2Tdx2 = (T[j][i - 1] - 2 * T[j][i] + T[j][i + 1]) / (dx * dx);
                double d2Tdy2 = (T[j - 1][i] - 2 * T[j][i] + T[j + 1][i]) / (dy * dy);
                T_new[j][i] = T[j][i] + alpha * dt * (d2Tdx2 + d2Tdy2);
            }
        }

        // Copy T_new to T
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                T[j][i] = T_new[j][i];
            }
        }

        if (step % output_interval == 0 || step == nsteps) {
            writeTempToFile(T, nx, ny, step);
        }
    }

    for (int j = 0; j < ny; j++) {
        free(T[j]);
        free(T_new[j]);
    }
    free(T);
    free(T_new);

    return 0;
}