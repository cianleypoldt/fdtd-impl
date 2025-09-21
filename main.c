#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned int uint;

typedef struct {
    unsigned x_min, x_max;
    unsigned y_min, y_max;
    unsigned z_min, z_max;
} grid_cutout;

double *field_mem = NULL;

double *E_x, *E_y, *E_z;
double *B_x, *B_y, *B_z;
double *inv_eps;  // (1/ε)
double *inv_mu;   // (1/μ)

enum boundry_condition {
    PML = 0,
    PEC,
    PMC,
    WRAPPED
};

enum boundry_condition boundry = 0;

grid_cutout regular_E_field_layers = (grid_cutout) { 0, 0, 0, 0, 0, 0 };
grid_cutout regular_B_field_layers = (grid_cutout) { 0, 0, 0, 0, 0, 0 };

double   dt = 0;
double   dx = 0, dy = 0, dz = 0;
unsigned grid_dim_x = 0, grid_dim_y = 0, grid_dim_z = 0;
double   time = 0;

void alloc_field() {
    assert(grid_dim_x > 0 && grid_dim_y > 0 && grid_dim_z > 0);

    unsigned total_elements = grid_dim_x * grid_dim_y * grid_dim_z;
    size_t   total_size     = total_elements * 8 * sizeof(double);
    field_mem               = (double *) malloc(total_size);
    if (!field_mem) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    memset(field_mem, 0, total_size);

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

void create_field(double size_x, double size_y, double size_z, char *boundry_by_name) {
    grid_dim_x = 10, grid_dim_y = 10, grid_dim_z = 10;

    if (strcmp(boundry_by_name, "PML") == 0) {
        regular_B_field_layers = (grid_cutout) { 1, grid_dim_x - 2, 1, grid_dim_y - 2, 1, grid_dim_z - 2 };
        regular_E_field_layers = (grid_cutout) { 1, grid_dim_x - 2, 1, grid_dim_y - 2, 1, grid_dim_z - 2 };
        boundry                = PML;
    }
    if (strcmp(boundry_by_name, "PEC") == 0) {
        regular_B_field_layers = (grid_cutout) { 1, grid_dim_x - 2, 1, grid_dim_y - 2, 1, grid_dim_z - 2 };
        regular_E_field_layers = (grid_cutout) { 1, grid_dim_x - 2, 1, grid_dim_y - 2, 1, grid_dim_z - 2 };
        boundry                = PEC;
    }
    if (strcmp(boundry_by_name, "PMC") == 0) {
        regular_B_field_layers = (grid_cutout) { 1, grid_dim_x - 2, 1, grid_dim_y - 2, 1, grid_dim_z - 2 };
        regular_E_field_layers = (grid_cutout) { 1, grid_dim_x - 2, 1, grid_dim_y - 2, 1, grid_dim_z - 2 };
        boundry                = PMC;
    }
    if (strcmp(boundry_by_name, "WRAPPED") == 0) {
        regular_B_field_layers = (grid_cutout) { 0, grid_dim_x - 1, 0, grid_dim_y - 1, 0, grid_dim_z - 1 };
        regular_E_field_layers = (grid_cutout) { 0, grid_dim_x - 1, 0, grid_dim_y - 1, 0, grid_dim_z - 1 };
        boundry                = WRAPPED;
    }
}

void update_E_field() {
    const unsigned stride_x = grid_dim_y * grid_dim_z;
    const unsigned stride_y = grid_dim_z;
    const unsigned stride_z = 1;

    // first and last E field layers updated later on as PEC boundary, index initialized at E_*[1, 1, 1]
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
    uint           idx      = 0;

    // rear B field layers are buffers (-> for grid_dim_* - 1)
    for (int i = 0; i < grid_dim_x - 1; i++) {
        for (int j = 0; j < grid_dim_y - 1; j++) {
            for (int k = 0; k < grid_dim_z - 1; k++) {
                B_x[idx] -= dt * ((E_z[idx + stride_y] - E_z[idx]) / dy - (E_y[idx + stride_z] - E_y[idx]) / dz) * inv_mu[idx];
                B_y[idx] -= dt * ((E_x[idx + stride_z] - E_x[idx]) / dz - (E_z[idx + stride_x] - E_z[idx]) / dx) * inv_mu[idx];
                B_z[idx] -= dt * ((E_y[idx + stride_x] - E_y[idx]) / dx - (E_x[idx + stride_y] - E_x[idx]) / dy) * inv_mu[idx];
                idx++;
            }
            idx++;        // skip B field-buffer
        }
        idx += stride_y;  // skip B field-buffer
    }
}

// Initialize B field to t = -dt/2 for leapfrog starting point
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
            idx++;        // skip B field-buffer
        }
        idx += stride_y;  // skip B field-buffer
    }
}

void advance_EM_field() {
    update_E_field();
    time += dt;
    update_B_field();
}

int main() {
    create_field(50, 50, 50, "PEC");

    initialize_leapfrog();

    const double end_time = 0.001;
    while (time < end_time) {
        advance_EM_field();
    }

    free_field();
    return 0;
}
