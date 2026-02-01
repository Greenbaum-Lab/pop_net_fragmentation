#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>

int sum_of_n_k(int n, int p) {
    int sum = 0;
    for (int k = 0; k < p; k++) {
        sum += (n - k);
    }
    return sum;
}

int min(int a, int b) {
    return a < b ? a : b;
}

int max(int a, int b) {
    return a > b ? a : b;
}

double calculate_first_coefficients(double* matrix, int n, int j, int i, int same_pop, int lower_bound, int upper_bound, int* counter) {
    if (j == same_pop) {
        double sum = 0;
        for (int k = 0; k < n; k++) {
            sum += matrix[i * n + k];
        }
        return 1 + sum;
    }
    if (lower_bound <= j && j <= upper_bound) {
        (*counter)++;
        return -1 * matrix[i * n + (i + *counter - 1)];
    }
    for (int p = 0; p < i; p++) {
        if (j == (i - p) + sum_of_n_k(n, p)) {
            return -1 * matrix[i * n + p];
        }
    }
    return 0;
}

double calculate_last_coefficients(double* matrix, int n, int j, int cur_pop, int other_pop) {
    int p, t_index, t, not_t, min_t_p, max_t_p, sum;
    int t_values[2] = {other_pop, cur_pop};
    if (j == sum_of_n_k(n, other_pop) + (cur_pop - other_pop)) {
        double sum_cur_pop = 0, sum_other_pop = 0;
        for (int k = 0; k < n; k++) {
            sum_cur_pop += matrix[cur_pop * n + k];
            sum_other_pop += matrix[other_pop * n + k];
        }
        return sum_cur_pop + sum_other_pop;
    }
    for (p = 0; p < n; p++) {
        for (t_index = 0; t_index < 2; t_index++) {
            t = t_values[t_index];
            if (t == other_pop) {
                not_t = cur_pop;
            } else {
                not_t = other_pop;
            }
            if (p != not_t) {
                min_t_p = min(t, p);
                max_t_p = max(t, p);
                if (j == sum_of_n_k(n, min_t_p) + max_t_p - min_t_p) {
                    return -1 * matrix[not_t * n + p];
                }
            }
        }
    }
    return 0;
}


double* coefficient_matrix_from_migration(double* migration_matrix, int n) {
    int n_last_equations = (n * (n - 1)) / 2;
    int n_first_equations = n;
    int mat_size = n_first_equations + n_last_equations;

    double* coefficient_mat = (double*) malloc(mat_size * mat_size * sizeof(double));

    int* counter = (int*) malloc(sizeof(int));
    for (int i = 0; i < n_first_equations; i++) {
        *counter = 1;
        int same_population = sum_of_n_k(n, i);
        int lower_bound = same_population + 1;
        int upper_bound = sum_of_n_k(n, i + 1) - 1;
        for (int j = 0; j < mat_size; j++) {
            coefficient_mat[i * mat_size + j] = calculate_first_coefficients(migration_matrix, n, j, i, same_population, lower_bound, upper_bound, counter);
        }
    }

    free(counter);

    int cur_population = 1;
    int other_population = 0;
    for (int i = n_first_equations; i < mat_size; i++) {
        if (other_population == cur_population) {
            other_population = 0;
            cur_population += 1;
        }
        for (int j = 0; j < mat_size; j++) {
            coefficient_mat[i * mat_size + j] = calculate_last_coefficients(migration_matrix, n, j, cur_population, other_population);
        }
        other_population += 1;
    }

    return coefficient_mat;
}

double* produce_solution_vector(int n) {
    int n_first = n;
    int n_last = (n * (n - 1)) / 2;
    double* b = (double*) malloc((n_first + n_last) * sizeof(double));
    for (int i = 0; i < n_first; i++) {
        b[i] = 1;
    }
    for (int i = n_first; i < n_first + n_last; i++) {
        b[i] = 2;
    }
    return b;
}

double* coalescence_from_migration(double* migration_matrix, int n) {
    double* A = coefficient_matrix_from_migration(migration_matrix, n);
    double* b = produce_solution_vector(n);

    gsl_vector *x = gsl_vector_alloc(n * (n + 1) / 2);
    gsl_permutation * p = gsl_permutation_alloc(n * (n + 1) / 2);

    gsl_matrix_view m = gsl_matrix_view_array(A, n * (n + 1) / 2, n * (n + 1) / 2);
    gsl_vector_view bv = gsl_vector_view_array(b, n * (n + 1) / 2);

    int s;

    gsl_linalg_LU_decomp(&m.matrix, p, &s);
    gsl_linalg_LU_solve(&m.matrix, p, &bv.vector, x);

    double* T_mat = (double*) malloc(n * n * sizeof(double));
    int cur_ind = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            T_mat[i * n + j] = gsl_vector_get(x, cur_ind);
            T_mat[j * n + i] = gsl_vector_get(x, cur_ind);
            cur_ind++;
        }
    }

    gsl_permutation_free(p);
    gsl_vector_free(x);
    free(A);
    free(b);

    return T_mat;
}

