#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned int uint;

double *field_mem;

double *E_x, *E_y, *E_z;
double *B_x, *B_y, *B_z;
double *inv_eps;  // precomputed inverse permittivity (1/ε)
double *inv_mu;   // precomputed inverse permeability (1/μ)

double   dt = -1;
double   dx = -1, dy = -1, dz = -1;
unsigned grid_dim_x = 0, grid_dim_y = 0, grid_dim_z = 0;
double   time = 0;  // current simulation time

void create_field(unsigned dim_x, unsigned dim_y, unsigned dim_z) {
    time       = 0.0;  // Initialize simulation time
    grid_dim_x = dim_x + 1;
    grid_dim_y = dim_y + 1;
    grid_dim_z = dim_z + 1;

    unsigned total_elements = grid_dim_x * grid_dim_y * grid_dim_z;

    // single contiguous memory block to be partitioned for field arrays
    size_t total_size = total_elements * 8 * sizeof(double);
    field_mem         = (double *) malloc(total_size);
    if (!field_mem) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    memset(field_mem, 0, total_size);

    // block partitioning
    E_x     = field_mem;
    E_y     = E_x + total_elements;
    E_z     = E_y + total_elements;
    B_x     = E_z + total_elements;
    B_y     = B_x + total_elements;
    B_z     = B_y + total_elements;
    inv_eps = B_z + total_elements;
    inv_mu  = inv_eps + total_elements;
}

void free_field() {
    if (field_mem) {
        free(field_mem);
        field_mem = NULL;
    }
}

void update_E_field() {
    const unsigned stride_x = grid_dim_y * grid_dim_z;
    const unsigned stride_y = grid_dim_z;
    const unsigned stride_z = 1;

    // first and last E field layers updated later as PEC boundary, and idx starts at E_*[1, 1, 1]
    uint idx = stride_z + stride_x + stride_y;
    for (int i = 1; i < grid_dim_x - 1; i++) {
        for (int j = 1; j < grid_dim_y - 1; j++) {
            for (int k = 1; k < grid_dim_z - 1; k++) {
                E_x[idx] += dt * ((B_z[idx] - B_z[idx - stride_y]) / dy - (B_y[idx] - B_y[idx - stride_z]) / dz) * inv_eps[idx];
                E_y[idx] += dt * ((B_x[idx] - B_x[idx - stride_z]) / dz - (B_z[idx] - B_z[idx - stride_x]) / dx) * inv_eps[idx];
                E_z[idx] += dt * ((B_y[idx] - B_y[idx - stride_x]) / dx - (B_x[idx] - B_x[idx - stride_y]) / dy) * inv_eps[idx];
                idx++;
            }
            idx += 2;         // skip PEC boundaries @ z = 0 and z = grid_dim_z - 1
        }
        idx += 2 * stride_y;  // skip PEC boundaries @ y = 0 and y = grid_dim_y - 1
    }
}

void update_B_field() {
    const unsigned stride_x = grid_dim_y * grid_dim_z;
    const unsigned stride_y = grid_dim_z;
    const unsigned stride_z = 1;

    uint idx = 0;
    // rear B-field layers are buffers (-> for grid_dim_* - 1)
    for (int i = 0; i < grid_dim_x - 1; i++) {
        for (int j = 0; j < grid_dim_y - 1; j++) {
            for (int k = 0; k < grid_dim_z - 1; k++) {
                B_x[idx] -= dt * ((E_z[idx + stride_y] - E_z[idx]) / dy - (E_y[idx + stride_z] - E_y[idx]) / dz) * inv_mu[idx];
                B_y[idx] -= dt * ((E_x[idx + stride_z] - E_x[idx]) / dz - (E_z[idx + stride_x] - E_z[idx]) / dx) * inv_mu[idx];
                B_z[idx] -= dt * ((E_y[idx + stride_x] - E_y[idx]) / dx - (E_x[idx + stride_y] - E_x[idx]) / dy) * inv_mu[idx];
                idx++;
            }
            idx++;        // skip B-field-buffer
        }
        idx += stride_y;  // skip B-field-buffer
    }
}

// Initialize B field to t = -dt/2 for proper leapfrog starting point
void initialize_leapfrog() {
    const double   half_dt  = dt * 0.5;
    const unsigned stride_x = grid_dim_y * grid_dim_z;
    const unsigned stride_y = grid_dim_z;
    const unsigned stride_z = 1;

    uint idx = 0;
    for (int i = 0; i < grid_dim_x - 1; i++) {
        for (int j = 0; j < grid_dim_y - 1; j++) {
            for (int k = 0; k < grid_dim_z - 1; k++) {
                B_x[idx] -= half_dt * ((E_z[idx + stride_y] - E_z[idx]) / dy - (E_y[idx + stride_z] - E_y[idx]) / dz) * inv_mu[idx];
                B_y[idx] -= half_dt * ((E_x[idx + stride_z] - E_x[idx]) / dz - (E_z[idx + stride_x] - E_z[idx]) / dx) * inv_mu[idx];
                B_z[idx] -= half_dt * ((E_y[idx + stride_x] - E_y[idx]) / dx - (E_x[idx + stride_y] - E_x[idx]) / dy) * inv_mu[idx];
                idx++;
            }
            idx++;        // skip B-field-buffer
        }
        idx += stride_y;  // skip B-field-buffer
    }
}

void advance_EM_field() {
    update_E_field();
    time += dt;
    update_B_field();
}

int main() {
    create_field(50, 50, 50);

    initialize_leapfrog();

    const double end_time = 100.0;
    while (time < end_time) {
        advance_EM_field();
    }

    free_field();
    return 0;
}
