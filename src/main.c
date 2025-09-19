#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct em_field {
  double *data;             // Contiguous memory block for all arrays
  double *e_x, *e_y, *e_z;  // Pointers to respective portions of data
  double *h_x, *h_y, *h_z;
  double *eps;
  double *mu;
  int dim_x, dim_y, dim_z;  // Dimensions of the field
} em_field;

em_field *alloc_field(int dim_x, int dim_y, int dim_z) {
  int total_elements = dim_x * dim_y * dim_z;
  em_field *field = (em_field *)malloc(sizeof(em_field));
  if (!field) return NULL;

  field->dim_x = dim_x;
  field->dim_y = dim_y;
  field->dim_z = dim_z;

  // Allocate a single contiguous block for all arrays
  size_t total_size = total_elements * 8 * sizeof(double);  // 8 arrays total
  field->data = (double *)malloc(total_size);
  if (!field->data) {
    free(field);
    return NULL;
  }

  memset(field->data, 0, total_size);

  // Set pointers to respective portions of the data block
  field->e_x = field->data;
  field->e_y = field->e_x + total_elements;
  field->e_z = field->e_y + total_elements;
  field->h_x = field->e_z + total_elements;
  field->h_y = field->h_x + total_elements;
  field->h_z = field->h_y + total_elements;
  field->eps = field->h_z + total_elements;
  field->mu = field->eps + total_elements;

  return field;
}

void update_e_field(em_field *field, double dt) {
  int nx = field->dim_x;
  int ny = field->dim_y;
  int nz = field->dim_z;
  int idx, idx_xp, idx_yp, idx_zp;
  double curl_h_x, curl_h_y, curl_h_z;

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        idx = i + j * nx + k * nx * ny;
        idx_yp = i + ((j + 1) % ny) * nx + k * nx * ny;
        idx_zp = i + j * nx + ((k + 1) % nz) * nx * ny;
        idx_xp = ((i + 1) % nx) + j * nx + k * nx * ny;

        // Calculate curl of H
        curl_h_x = (field->h_z[idx_yp] - field->h_z[idx]) - (field->h_y[idx_zp] - field->h_y[idx]);

        curl_h_y = (field->h_x[idx_zp] - field->h_x[idx]) - (field->h_z[idx_xp] - field->h_z[idx]);

        curl_h_z = (field->h_y[idx_xp] - field->h_y[idx]) - (field->h_x[idx_yp] - field->h_x[idx]);

        // Update E field components
        field->e_x[idx] += dt * curl_h_x / field->eps[idx];
        field->e_y[idx] += dt * curl_h_y / field->eps[idx];
        field->e_z[idx] += dt * curl_h_z / field->eps[idx];
      }
    }
  }
}

void update_h_field(em_field *field, double dt) {
  int nx = field->dim_x;
  int ny = field->dim_y;
  int nz = field->dim_z;
  int idx, idx_xm, idx_ym, idx_zm;
  double curl_e_x, curl_e_y, curl_e_z;

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        idx = i + j * nx + k * nx * ny;
        idx_ym = i + ((j - 1 + ny) % ny) * nx + k * nx * ny;
        idx_zm = i + j * nx + ((k - 1 + nz) % nz) * nx * ny;
        idx_xm = ((i - 1 + nx) % nx) + j * nx + k * nx * ny;

        // Calculate curl of E
        curl_e_x = (field->e_z[idx] - field->e_z[idx_ym]) - (field->e_y[idx] - field->e_y[idx_zm]);

        curl_e_y = (field->e_x[idx] - field->e_x[idx_zm]) - (field->e_z[idx] - field->e_z[idx_xm]);

        curl_e_z = (field->e_y[idx] - field->e_y[idx_xm]) - (field->e_x[idx] - field->e_x[idx_ym]);

        // Update H field components
        field->h_x[idx] -= dt * curl_e_x / field->mu[idx];
        field->h_y[idx] -= dt * curl_e_y / field->mu[idx];
        field->h_z[idx] -= dt * curl_e_z / field->mu[idx];
      }
    }
  }
}
