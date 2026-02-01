import ctypes
import numpy as np
import math
from collections import deque
import scipy.linalg as la

lib = ctypes.cdll.LoadLibrary('./libmigration.so')
lib.coefficient_matrix_from_migration.restype = ctypes.POINTER(ctypes.c_double)
lib.coefficient_matrix_from_migration.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int]
lib.coalescence_from_migration.restype = ctypes.POINTER(ctypes.c_double)
lib.coalescence_from_migration.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int]

def coefficient_matrix_from_migration_wrapper(migration_matrix):
    n = migration_matrix.shape[0]
    mat_size = n + (n * (n - 1)) // 2  # size of the coefficient matrix
    migration_matrix_c = migration_matrix.flatten().ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    result_c = lib.coefficient_matrix_from_migration(migration_matrix_c, n)
    result = np.ctypeslib.as_array(result_c, shape=(mat_size*mat_size,)).reshape((mat_size, mat_size))

    return result


class Coalescence:
    def __init__(self, matrix: np.ndarray) -> None:
        """
        Initialize a coalescence times matrix object
        :param matrix: input Coalescence time matrix
        """
        self.matrix = matrix
        self.shape = matrix.shape[0]

    def produce_fst(self) -> np.ndarray:
        """
        produces and returns the corresponding Fst matrix
        :return: The corresponding Fst matrix
        """
        F_mat = np.zeros((self.shape, self.shape))
        for i in range(self.shape):
            for j in range(i + 1, self.shape):
                t_S = (self.matrix[i, i] + self.matrix[j, j]) / 2
                t_T = (self.matrix[i, j] + t_S) / 2
                if np.isinf(t_T):
                    F_i_j = 1
                else:
                    F_i_j = (t_T - t_S) / t_T
                F_mat[i, j], F_mat[j, i] = F_i_j, F_i_j
        return F_mat


class Migration:
    def __init__(self, matrix: np.ndarray) -> None:
        """
        Initialize a migration matrix object
        :param matrix: input migration matrix
        """
        self.matrix = matrix
        self.shape = matrix.shape[0]


    def produce_coalescence(self) -> np.ndarray:
        n = self.matrix.shape[0]
        migration_matrix_c = self.matrix.flatten().ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        result_c = lib.coalescence_from_migration(migration_matrix_c, n)
        T_mat = np.ctypeslib.as_array(result_c, shape=(n * n,)).reshape((n, n))
        return T_mat
    def calculate_first_coefficients(self, j: int, i: int, same_pop: int, lower_bound: int, upper_bound: int,
                                     p_list: list, counter: list) -> float:
        """
        calculates the coefficients for the first n equations
        :param j: column of coefficient matrx
        :param i: row of coefficient matrix
        :param same_pop: The column corresponding to T(i,i)
        :param lower_bound: j values in range [lower_bound, upper_bound] correspond to coefficient -M(i,i+counter)
        :param upper_bound: j values in range [lower_bound, upper_bound] correspond to coefficient -M(i,i+counter)
        :param p_list: all values that are smaller than i
        :param counter: counts the number of time j was in the interval [lower_bound,upper_bound]
        :return: The coefficient in place [i,j], for i in [n-1]
        """
        n = self.matrix.shape[0]
        if j == same_pop:
            return 1 + np.sum(self.matrix[i, :])
        if lower_bound <= j <= upper_bound:
            counter[0] += 1
            return -1 * self.matrix[i, i + counter[0] - 1]
        for p in p_list:
            if j == (i - p) + np.sum([n - k for k in range(p)]):
                return -1 * self.matrix[i, p]
        return 0


    def calculate_last_coefficients(self, j, cur_pop, other_pop) -> float:
        """
        calculates the coefficients for the last (n choose 2) rows in the coefficient matrx
        :param j: the column in the coefficient matrix
        :param cur_pop: the index of the population that corresponds to the current value
        :param other_pop: the index of the other population that corresponds to the current value
        :return: The coefficient in the coefficient matrix according to certain conditions deduced from
        Wilkinson-Herbots' equations.
        """
        n = self.matrix.shape[0]
        if j == np.sum([n - k for k in range(other_pop)]) + (cur_pop - other_pop):
            return float(np.sum(self.matrix[[cur_pop, other_pop], :]))
        for p in range(n):
            for t in [other_pop, cur_pop]:
                if t == other_pop:
                    not_t = cur_pop
                else:
                    not_t = other_pop
                if p != not_t:
                    min_t_p = min(t, p)
                    max_t_p = max(t, p)
                    if j == np.sum([n - k for k in range(min_t_p)]) + max_t_p - min_t_p:
                        return -1 * self.matrix[not_t, p]
        return 0


    def produce_solution_vector(self):
        """
        produce the solution vector(b), according to Wilkinson-Herbot's equations
        :return: solution vector b
        """
        n = self.shape
        n_first = np.repeat(1, n)
        n_last = np.repeat(2, comb(n, 2))
        return np.hstack((n_first, n_last))


def comb(n: int, k: int) -> int:
    """
    calculate and return n Choose k
    :param n: number of objects
    :param k: number of selected objects
    :return: n Choose k
    """
    return int(math.factorial(n) / (math.factorial(k) * math.factorial(n - k)))


def find_fst(m: np.ndarray) -> np.ndarray:
    """
    Receives a migration matrix with one connected component(a squared, positive matrix with zeroes on the diagonal),
    and returns it's corresponding Fst matrix according to Wilkinson-Herbot's equations.
    :param m: Migration matrix- squared, positive, with zeroes on the diagonal.
    :return: Corresponding Fst matrix according to Wilkinson-Herbot's equations. If there is no solution, an error will
    occur.
    """
    if m.shape[0] == 1:
        return np.array([[0]])
    M = Migration(m)
    t = M.produce_coalescence()
    T = Coalescence(t)
    return T.produce_fst()


def find_coalescence(m: np.ndarray) -> np.ndarray:
    """
       Receives a migration matrix with one connected component
       (a squared, positive matrix with zeroes on the diagonal), and returns it's corresponding Coalescent times
       (T) matrix according to Wilkinson-Herbot's equations.
       :param m: Migration matrix- squared, positive, with zeroes on the diagonal.
       :return: Corresponding T matrix according to Wilkinson-Herbot's equations. If there is no solution,
       an error will occur.
       """
    if m.shape[0] == 1:
        return np.array([[1]])
    M = Migration(m)
    return M.produce_coalescence()


def find_components(matrix: np.ndarray) -> dict:
    """
    Find connected components in a connected graph represented by adjacency matrix.
    :param matrix: adjacency matrix representing a directed graph
    :return:something
    """
    components = 1
    n = matrix.shape[0]
    queue = deque()
    visited = set()
    not_visited = set([i for i in range(1, n)])
    visited.add(0)
    comp_dict = {components: [0]}
    queue.append(0)
    while len(not_visited) != 0:
        while len(queue) != 0:
            cur_vertex = queue.popleft()
            for i in range(n):
                if i not in visited and (matrix[cur_vertex, i] != 0 or matrix[i, cur_vertex] != 0):
                    queue.append(i)
                    visited.add(i)
                    not_visited.remove(i)
                    comp_dict[components].append(i)
        for vertex in not_visited:
            components += 1
            queue.append(vertex)
            visited.add(vertex)
            not_visited.remove(vertex)
            comp_dict[components] = [vertex]
            break
    return comp_dict


def split_migration_matrix(migration_matrix: np.ndarray, connected_components: list) -> list:
    """
    Splits a migration matrix to sub-matrices according to it's connected components.
    :param migration_matrix: A valid migration matrix.
    :param connected_components: list of lists, where each list represents a connected component's vertices
    (populations).
    :return: A list of sub-matrices, where each sub-matrix is the migration matrix of a connected component. Note that
    in order to interpret which populations are described in each sub matrix the connected components list is needed.
    """
    sub_matrices = []
    for component in connected_components:
        sub_matrix = migration_matrix[np.ix_(component, component)]
        sub_matrices.append(sub_matrix)

    return sub_matrices


def split_migration(migration_matrix: np.ndarray) -> tuple:
    """
    Finds a migration matrix connected components, and splits the matrix to its connected components.
    :param migration_matrix: A valid migration matrix.
    :return: A tuple (sub_matrices, components). Sub matrices is a list of numpy arrays, where each array is a
    component's migration matrix. components is a list of lists, where each list represents a component vertices
    (populations). The order of the components corresponds to the order of the sub-matrices.
    """
    components = list(find_components(migration_matrix).values())
    sub_matrices = split_migration_matrix(migration_matrix, components)
    return sub_matrices, components


def reassemble_matrix(sub_matrices: list, connected_components: list, which: str) -> np.ndarray:
    num_nodes = sum(len(component) for component in connected_components)
    if which == "fst":
        adjacency_matrix = np.ones((num_nodes, num_nodes), dtype=float)
    else:
        adjacency_matrix = np.full((num_nodes, num_nodes), np.inf)

    for component, sub_matrix in zip(connected_components, sub_matrices):
        indices = np.array(component)
        adjacency_matrix[np.ix_(indices, indices)] = sub_matrix

    return adjacency_matrix


def m_to_f(m: np.ndarray) -> np.ndarray:
    """
    Receives a migration matrix(a squared, positive matrix with zeroes on the diagonal) with any number
    of connected components, and returns it's corresponding Fst matrix according to Wilkinson-Herbot's equations and
    Slatkin equations.
    :param m: Migration matrix- squared, positive, with zeroes on the diagonal.
    :return: Corresponding Fst matrix according to Wilkinson-Herbot's equations. If there is no solution, an error will
    occur.
    """
    split = split_migration(m)
    sub_matrices, components = split[0], split[1]
    f_matrices = []
    for matrix in sub_matrices:
        f_matrices.append(find_fst(matrix))
    return reassemble_matrix(f_matrices, components, "fst")


def m_to_t(m: np.ndarray) -> np.ndarray:
    """
       Receives a migration matrix(a squared, positive matrix with zeroes on the diagonal) with any number
       of connected components, and returns its corresponding Coalescent times (T) matrix according to
       Wilkinson-Herbot's equations.
       :param m: Migration matrix- squared, positive, with zeroes on the diagonal.
       :return: Corresponding T matrix according to Wilkinson-Herbot's equations. If there is no solution,
       an error will occur.
       """
    split = split_migration(m)
    sub_matrices, components = split[0], split[1]
    t_matrices = []
    for matrix in sub_matrices:
        t_matrices.append(find_coalescence(matrix))
    return reassemble_matrix(t_matrices, components, "coalescence")

def transform_matrix(m: np.ndarray) -> tuple:
    """
       Receives a migration matrix (a squared, positive matrix with zeroes on the diagonal) with any number
       of connected components, and returns its corresponding Coalescent times (T) matrix according to
       Wilkinson-Herbot's equations, and it's corresponding Fst matrix(F) according to Slatkin equations.
       :param m: Migration matrix- squared, positive, with zeroes on the diagonal.
       :return:  A tuple (T,F). Corresponding T matrix according to Wilkinson-Herbot's equations,
       Corresponding F matrix according to Slatkin equations.
       If there is no solution, an error will occur.
       """
    split = split_migration(m)
    sub_matrices, components = split[0], split[1]
    t_matrices = []
    f_matrices = []
    for matrix in sub_matrices:
        t_matrix = find_coalescence(matrix)
        t_matrices.append(t_matrix)
        T = Coalescence(t_matrix)
        f_matrices.append(T.produce_fst())
    return reassemble_matrix(t_matrices, components, "coalescence"), reassemble_matrix(f_matrices, components, "fst")


def check_conservative(m: np.ndarray):
    """
    checks if a given migration matrix is conservative.
    :param m: A migration matrix
    :return: True if m is conservative, False otherwise.
    """
    for i in range(m.shape[0]):
        if np.sum(m[i, :]).round(8) != np.sum(m[:, i]).round(8):
            print("Found non-conservative migration matrix.")
            return False
    print("Migration matrix is conservative.")
    return True


def conservative_from_normal(binary_m,      # 0/1 mask of allowed directed edges
                             mu=1.0, sigma=0.3,
                             lower=0.2, upper=4,
                             seed=None):
    """
    1. Draw each allowed edge weight ~ N(mu, sigma^2), truncated to (lower, upper).
    2. Orthogonally project onto the conservative sub-space (row sum = col sum).
    3. If any entry leaves (lower, upper) by more than 1e-12, clip & re-project once.

    Returns
    -------
    M : (n,n) ndarray  -- conservative migration matrix.
    """
    rng = np.random.default_rng(seed)

    # ---------- step 1 : raw draw -----------------------------------
    nz   = np.argwhere(binary_m == 1)          # list of directed edges
    n, _ = binary_m.shape
    E    = len(nz)

    x = rng.normal(mu, sigma, size=E)
    # truncate to (lower, upper)
    x = np.clip(x, lower+1e-12, upper-1e-12)

    # ---------- build balance matrix  A  ----------------------------
    A = np.zeros((n, E))
    for j, (r, c) in enumerate(nz):
        A[r, j] =  1         # outflow
        A[c, j] = -1         # inflow

    # pre-compute projector  P = I – Aᵀ (AAᵀ)⁻¹ A
    # (AAᵀ) is n×n and of rank n-1 when graph is strongly connected
    ATA_inv = la.pinv(A @ A.T, rtol=1e-10)
    P = np.eye(E) - A.T @ ATA_inv @ A

    # ---------- step 2 : first projection ---------------------------
    m = P @ x

    # ---------- step 3 : clip & one Dykstra correction --------------
    m = np.clip(m, lower, upper)
    # Project again to restore perfect balance (cheap; AAᵀ has same inverse)
    m = P @ m

    # ---------- pack back into square matrix ------------------------
    M = np.zeros_like(binary_m, dtype=float)
    M[nz[:, 0], nz[:, 1]] = m
    return M