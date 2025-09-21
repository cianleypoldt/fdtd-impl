#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const double vacuum_permeability = 4.0 * M_PI * 1e-7;  // μ₀ in H/m
const double vacuum_permittivity = 8.85418782e-12;     // ε₀ in F/m

double *field_mem  = NULL;
int     cell_count = 0;

int nx = 0, ny = 0, nz = 0;

double *Ex, *Ey, *Ez;
double *Bx, *By, *Bz;
double *eps;
double *mu;

double dt = 0;
double dx = 0, dy = 0, dz = 0;
double simulation_time = 0;

int stride_x = 0;
int stride_y = 0;
int stride_z = 0;

void alloc_field(int x_dim, int y_dim, int z_dim) {
    nx = x_dim + 1;
    ny = y_dim + 1;
    nz = z_dim + 1;

    cell_count        = nx * ny * nz;
    size_t total_size = cell_count * 8 * sizeof(double);
    field_mem         = (double *) malloc(total_size);
    if (!field_mem) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    memset(field_mem, 0, total_size);

    Ex  = field_mem;
    Ey  = Ex + cell_count;
    Ez  = Ey + cell_count;
    Bx  = Ez + cell_count;
    By  = Bx + cell_count;
    Bz  = By + cell_count;
    eps = Bz + cell_count;
    mu  = eps + cell_count;

    // initialize eps and mu
    for (int i = 0; i < cell_count; i++) {
        *(eps + i) = vacuum_permittivity;
        *(mu + i)  = vacuum_permeability;
    }

    stride_x = ny * nz;
    stride_y = nz;
    stride_z = 1;
}

void free_field() {
    if (field_mem) {
        free(field_mem);
        field_mem = NULL;
    }
}

void update_E(double time_step) {
    // first and last E field layers updated later on as PEC boundary, index initialized at E_*[1, 1, 1]
    int idx = stride_z + stride_x + stride_y;
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                Ex[idx] += time_step * ((Bz[idx] - Bz[idx - stride_y]) / dy - (By[idx] - By[idx - stride_z]) / dz) / eps[idx];
                Ey[idx] += time_step * ((Bx[idx] - Bx[idx - stride_z]) / dz - (Bz[idx] - Bz[idx - stride_x]) / dx) / eps[idx];
                Ez[idx] += time_step * ((By[idx] - By[idx - stride_x]) / dx - (Bx[idx] - Bx[idx - stride_y]) / dy) / eps[idx];
                idx++;
            }
            idx += 2;         // skip PEC boundaries @ z = 0 and z = grid_dim_z - 1
        }
        idx += 2 * stride_y;  // skip PEC boundaries @ y = 0 and y = grid_dim_y - 1
    }
}

void update_B(double time_step) {
    // rear B field layers are buffers (-> for grid_dim_* - 1)
    int idx = 0;
    for (int i = 0; i < nx - 1; i++) {
        for (int j = 0; j < ny - 1; j++) {
            for (int k = 0; k < nz - 1; k++) {
                Bx[idx] -= time_step * ((Ez[idx + stride_y] - Ez[idx]) / dy - (Ey[idx + stride_z] - Ey[idx]) / dz) / mu[idx];
                By[idx] -= time_step * ((Ex[idx + stride_z] - Ex[idx]) / dz - (Ez[idx + stride_x] - Ez[idx]) / dx) / mu[idx];
                Bz[idx] -= time_step * ((Ey[idx + stride_x] - Ey[idx]) / dx - (Ex[idx + stride_y] - Ex[idx]) / dy) / mu[idx];
                idx++;
            }
            idx++;        // skip B field-buffer
        }
        idx += stride_y;  // skip B field-buffer
    }
}

void apply_PEC_border_condition() {
    // Apply PEC boundary conditions for x-planes
    int idx_begin = 0;
    int idx_end   = (nx - 1) * stride_x;
    for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
            Ex[idx_begin] = -Ex[idx_begin + stride_x];
            Ey[idx_begin] = -Ey[idx_begin + stride_x];
            Ez[idx_begin] = -Ez[idx_begin + stride_x];

            Ex[idx_end] = -Ex[idx_end - stride_x];
            Ey[idx_end] = -Ey[idx_end - stride_x];
            Ez[idx_end] = -Ez[idx_end - stride_x];

            idx_begin++;
            idx_end++;
        }
        idx_begin += stride_y - nz;
        idx_end += stride_y - nz;
    }

    // Apply PEC boundary conditions for y-planes
    for (int i = 0; i < nx; i++) {
        idx_begin = i * stride_x;
        idx_end   = i * stride_x + (ny - 1) * stride_y;
        for (int k = 0; k < nz; k++) {
            Ex[idx_begin] = -Ex[idx_begin + stride_y];
            Ey[idx_begin] = -Ey[idx_begin + stride_y];
            Ez[idx_begin] = -Ez[idx_begin + stride_y];

            Ex[idx_end] = -Ex[idx_end - stride_y];
            Ey[idx_end] = -Ey[idx_end - stride_y];
            Ez[idx_end] = -Ez[idx_end - stride_y];

            idx_begin++;
            idx_end++;
        }
    }

    // Apply PEC boundary conditions for z-planes
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            idx_begin = i * stride_x + j * stride_y;
            idx_end   = idx_begin + (nz - 1);

            Ex[idx_begin] = -Ex[idx_begin + stride_z];
            Ey[idx_begin] = -Ey[idx_begin + stride_z];
            Ez[idx_begin] = -Ez[idx_begin + stride_z];

            Ex[idx_end] = -Ex[idx_end - stride_z];
            Ey[idx_end] = -Ey[idx_end - stride_z];
            Ez[idx_end] = -Ez[idx_end - stride_z];
        }
    }
}

int main() {
}
