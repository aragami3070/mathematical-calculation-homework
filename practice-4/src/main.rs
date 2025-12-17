use ndarray::{Array1, Array2, ArrayBase, Dim, OwnedRepr, array};
use ndarray_linalg::Solve;

/// создает матрицу для нахождения коэффициентов кубических сплайнов
fn create_spline_matrix(x: &Array1<f64>, f: &Array1<f64>) -> (Array2<f64>, Array1<f64>) {
    let n = x.len() - 1; // количество интервалов
    let size = 4 * n;
    let h = x[1] - x[0];

    let mut a = Array2::<f64>::zeros((size, size));
    let mut b = Array1::<f64>::zeros(size);

    // s_i(x_i) = f_i
    for i in 0..n {
        let row = i;
        a[[row, 4 * i]] = 1.0; // a_i
        b[row] = f[i];
    }

    // s_i(x_{i+1}) = f_{i+1}
    for i in 0..n {
        let row = n + i;
        a[[row, 4 * i]] = 1.0;
        a[[row, 4 * i + 1]] = h;
        a[[row, 4 * i + 2]] = h.powi(2);
        a[[row, 4 * i + 3]] = h.powi(3);
        b[row] = f[i + 1];
    }

    // s'_i(x_{i+1}) = s'_{i+1}(x_{i+1})
    for i in 0..(n - 1) {
        let row = 2 * n + i;

        a[[row, 4 * i + 1]] = 1.0;
        a[[row, 4 * i + 2]] = 2.0 * h;
        a[[row, 4 * i + 3]] = 3.0 * h.powi(2);
        a[[row, 4 * (i + 1) + 1]] = -1.0;
        b[row] = 0.0;
    }

    // s''_i(x_{i+1}) = s''_{i+1}(x_{i+1})
    for i in 0..(n - 1) {
        let row = 3 * n - 1 + i;

        a[[row, 4 * i + 2]] = 2.0;
        a[[row, 4 * i + 3]] = 6.0 * h;
        a[[row, 4 * (i + 1) + 2]] = -2.0;
        b[row] = 0.0;
    }

    // s''_0(x_0) = 0 и s''_{n-1}(x_n) = 0
    let row1 = 4 * n - 2;
    a[[row1, 2]] = 2.0; // 2*c_0

    let row2 = 4 * n - 1;
    a[[row2, 4 * (n - 1) + 2]] = 2.0;
    a[[row2, 4 * (n - 1) + 3]] = 6.0 * h;

    (a, b)
}

/// решает систему уравнений для нахождения коэффициентов сплайнов
fn solve_spline_coefficients(
    x: &Array1<f64>,
    f: &Array1<f64>,
) -> ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>, f64> {
    let (a, b) = create_spline_matrix(x, f);
    let coeffs = a.solve_into(b).expect("Failed to solve system");
    let n = x.len() - 1;

    coeffs.into_shape_clone((n, 4)).expect("Reshape error")
}

fn main() {
    let x = array![0.0, 1.0, 2.0, 3.0];
    let f = array![1.0, 2.0, 9.0, 28.0];

    let (a, b) = create_spline_matrix(&x, &f);
    let coefficients = solve_spline_coefficients(&x, &f);

    let size = 4 * (x.len() - 1);

    println!("Матрица системы {}x{}:", size, size);
    for (i, row) in a.axis_iter(ndarray::Axis(0)).enumerate() {
        print!("eq_{:<2}: [", i);
        for val in row.iter() {
            print!("{:8.3} ", val);
        }
        println!("]");
    }

    println!("\nВектор правой части:");
    for (i, val) in b.iter().enumerate() {
        println!("eq_{:<2}: {:.3}", i, val);
    }

    println!("\nКоэффициенты сплайнов (по строкам: [a_i, b_i, c_i, d_i] для каждого интервала):");
    for (i, coeff_row) in coefficients.axis_iter(ndarray::Axis(0)).enumerate() {
        println!(
            "interval {}: a_i = {:.3}, b_i = {:.3}, c_i = {:.3}, d_i = {:.3}",
            i, coeff_row[0], coeff_row[1], coeff_row[2], coeff_row[3]
        );
    }

    println!("\nРазмер матрицы: {:?}", a.dim());
    println!("Количество интервалов: {}", x.len() - 1);

    println!("Сплайны: ");
    println!("S_{}: {}", x[0], x[0]);
    for (i, coeff_row) in coefficients.axis_iter(ndarray::Axis(0)).enumerate() {
        let xi = x[i];
        let xip1 = x[i + 1];
        let dx = xip1 - xi;

        // midpoint local coordinate
        let t_mid = dx / 2.0;
        let s_mid = coeff_row[0]
            + coeff_row[1] * t_mid
            + coeff_row[2] * t_mid.powi(2)
            + coeff_row[3] * t_mid.powi(3);
        println!("S_{:.1}: {:.3}", xi + t_mid, s_mid);

        // right end local coordinate
        let t_end = dx;
        let s_end = coeff_row[0]
            + coeff_row[1] * t_end
            + coeff_row[2] * t_end.powi(2)
            + coeff_row[3] * t_end.powi(3);
        println!("S_{:.1}: {:.3}", xip1, s_end);
    }
}
