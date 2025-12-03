use ndarray::{Array1, arr2};
use ndarray_linalg::Determinant;
use std::f64;

fn norm_stop(xk: &Array1<f64>, xkp1: &Array1<f64>, epsilon: f64) -> bool {
    xk.iter()
        .zip(xkp1.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0. / 0., f64::max)
        < epsilon
}

fn main() {
    // Define matrix A
    let a = arr2(&[
        [18.0, 0.18, 0.18, 0.18, 0.18],
        [0.19, 19.0, 0.19, 0.19, 0.19],
        [0.20, 0.20, 20.0, 0.20, 0.20],
        [0.21, 0.21, 0.21, 21.0, 0.21],
        [0.22, 0.22, 0.22, 0.22, 22.0],
    ]);

    println!("Матрица A:");
    println!("{:?}", a);

    let det = a.det().expect("Failed to compute determinant");
    println!("Определитель матрицы A: {}", det);

    // Compute vector b with diagonal elements of A
    let b_diag: Array1<f64> = a.diag().to_owned();
    // b = A * b_diag (reshaped as column vector)
    let b = a.dot(&b_diag);

    println!("Столбец b:");
    println!("{:?}", b);

    // Create matrix alpha
    let n = a.nrows();
    let mut alpha = a.clone();
    for i in 0..n {
        for j in 0..n {
            if i == j {
                alpha[(i, j)] = 0.0;
            } else {
                alpha[(i, j)] = -a[(i, j)] / a[(i, i)];
            }
        }
    }

    println!("Матрица alpha:");
    println!("{:?}", alpha);

    // Create vector beta
    let mut beta = b.clone();
    for i in 0..n {
        beta[i] = b[i] / a[(i, i)];
    }

    println!("Столбец beta:");
    println!("{:?}", beta);

    let epsilon = 1e-9;
    println!("Считаем до точности epsilon= {}", epsilon);

    let mut xk = Array1::<f64>::zeros(n);

    println!("x^(0) = {:?}", xk);

    for i in 0..17 {
        let xkp1 = alpha.dot(&xk) + &beta;
        println!("x^({})= {:?}", i + 1, xkp1);
        if norm_stop(&xk, &xkp1, epsilon) {
            break;
        }
        xk = xkp1;
    }
}
