use ndarray::{Array, Array2, ArrayBase, Axis, Dim, OwnedRepr};
use ndarray_linalg::Solve;

fn create_matrix(x_vec: &Vec<f64>) -> ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>> {
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

fn solve_slay(x_vec: &Vec<f64>, f_vec: Vec<f64>) -> Vec<f64> {
    let matrix = create_matrix(x_vec);

    // Заполняем матрицу СЛАУ (правая часть матрицы)
    let f_column = Array::from_vec(f_vec);

    // Решение системы
    matrix
        .solve_into(f_column)
        .expect("Решение не найдено")
        .into_iter()
        .collect()
}

pub fn vandermonde_interpolation(x_vec: Vec<f64>, f_vec: Vec<f64>) -> (Vec<f64>, Vec<f64>) {
    let slay_solved = solve_slay(&x_vec, f_vec.clone());

    let mut prev_x = 0.0;
    let mut res_x_vec: Vec<f64> = Vec::new();
    let mut res_f_vec: Vec<f64> = Vec::new();

    // Находим промежуточные решения
    for (index, current_x) in x_vec.iter().enumerate() {
        if index == 0 {
            res_x_vec.push(*current_x);
            res_f_vec.push(f_vec[index]);
            prev_x = *current_x;
            continue;
        }

        let mid_x = (prev_x + current_x) / 2.0;
        let mut new_mid_x = 1.0;
        let mut new_mid_f = 0.0;
        for coefficients in slay_solved.iter().rev() {
            new_mid_f += coefficients * new_mid_x;
            new_mid_x *= mid_x;
        }

        res_x_vec.push(mid_x);
        res_f_vec.push(new_mid_f);

        res_x_vec.push(*current_x);
        res_f_vec.push(f_vec[index]);
    }

    (res_x_vec, res_f_vec)
}
