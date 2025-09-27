use ndarray::{Array, Array2, ArrayBase, Axis, Dim, OwnedRepr};
use ndarray_linalg::Solve;

fn create_matrix(x_vec: Vec<f64>) -> ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>> {
    let x_len = x_vec.len();

    let mut matrix = Array2::<f64>::default((x_len, x_len));
    // Заполняем матрицу СЛАУ (левая часть матрицы)
    for (index, mut row) in matrix.axis_iter_mut(Axis(0)).enumerate() {
        let mut x_degree = 1.0;
        let x_value = x_vec[index];
        for col in row.iter_mut() {
            *col = x_degree;
            x_degree *= x_value;
        }
    }

    matrix
}
