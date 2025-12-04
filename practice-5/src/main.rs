use ndarray::{Array1, Array2, s};
use ndarray_linalg::Determinant;

fn forward_elimination(mut a: Array2<f64>, mut b: Array1<f64>) -> (Array2<f64>, Array1<f64>) {
    let n = b.len();
    for k in 0..(n - 1) {
        for i in (k + 1)..n {
            if a[[i, k]] != 0.0 {
                let factor = a[[i, k]] / a[[k, k]];
                // Update row i from k+1 to end: A[i, k+1..n] -= factor * A[k, k+1..n]
                for col in (k + 1)..n {
                    let a_ik = a[[i, k]];
                    let a_kk = a[[k, k]];
                    // factor = a_ik / a_kk;
                    let factor = a_ik / a_kk;
                    a[[i, col]] -= factor * a[[k, col]];
                }
                a[[i, k]] = 0.0; // Zero out explicitly for clarity
                b[i] -= factor * b[k];
            }
        }
    }
    (a, b)
}

fn backward_substitution(u: &Array2<f64>, b: &Array1<f64>) -> Array1<f64> {
    let n = b.len();
    let mut x = Array1::<f64>::zeros(n);
    for i in (0..n).rev() {
        let dot = if i + 1 < n {
            u.slice(s![i, (i + 1)..n]).dot(&x.slice(s![(i + 1)..n]))
        } else {
            0.0
        };
        x[i] = (b[i] - dot) / u[[i, i]];
    }
    x
}

fn main() {
    let a = Array2::<f64>::from_shape_vec(
        (5, 5),
        vec![
            18.0, 0.18, 0.18, 0.18, 0.18,
            0.19, 19.0, 0.19, 0.19, 0.19,
            0.20, 0.20, 20.0, 0.20, 0.20,
            0.21, 0.21, 0.21, 21.0, 0.21,
            0.22, 0.22, 0.22, 0.22, 22.0,
        ],
    )
    .unwrap();

    println!("1) Матрица A:\n{:?}", a);
    let det = a.det().unwrap();
    println!("Определитель матрицы A = {}", det);

    // vector b = diagonal of A
    let n = a.shape()[0];
    let mut b = Array1::<f64>::ones(n);
    println!("\nКолонка b:");
    for i in 0..n {
        b[i] = a[[i, i]];
        println!("[{}]", b[i]);
    }

    let (a_tmp, b_tmp) = forward_elimination(a, b);
    println!("Матрица после прямого прохода:");
    println!("{:?}", a_tmp);

    let solution = backward_substitution(&a_tmp, &b_tmp);

    println!("Решение (A|b) методом Гаусса");
    println!("Вектор решения после обратного прохода:");
    for (i, val) in solution.iter().enumerate() {
        println!("x{} = {}", i, val);
    }

    println!("x =");
    for val in solution.iter() {
        println!("[{}]", val);
    }
}
