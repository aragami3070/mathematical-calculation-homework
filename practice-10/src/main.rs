use std::error::Error;

use ndarray::{Array1, Array2, s};
use polars::prelude::*;

const V: f64 = 18.0;
const H: f64 = 0.1;
const N: usize = 180;

fn ytoch(x: f64) -> f64 {
    V * x * x * (x - V)
}

fn f(x: f64) -> f64 {
    4.0 * V * x.powi(4) - 3.0 * V * V * x.powi(3) + 6.0 * V * x - 2.0 * V * V
}

fn p(x: f64) -> f64 {
    x * x
}

fn q(x: f64) -> f64 {
    x
}

fn phi_k(x: f64, k: usize) -> f64 {
    x.powi(k as i32) * (x - V)
}

fn dphi_k(x: f64, k: usize) -> f64 {
    let k_f64 = k as f64;
    (k_f64 + 1.0) * x.powi((k.saturating_sub(1)) as i32)
        - V * k_f64 * x.powi((k.saturating_sub(2)) as i32)
}

fn ddphi_k(x: f64, k: usize) -> f64 {
    let k_f64 = k as f64;
    k_f64 * (k_f64 + 1.0) * x.powi((k.saturating_sub(2)) as i32)
        - V * k_f64 * (k_f64 - 1.0) * x.powi((k.saturating_sub(3)) as i32)
}

fn forward_elimination(mut a: Array2<f64>, mut b: Array1<f64>) -> (Array2<f64>, Array1<f64>) {
    let n = b.len();

    for k in 0..n {
        // Поиск максимального элемента в столбце k
        let mut max_row = k;
        let mut max_val = a[[k, k]].abs();
        for i in (k + 1)..n {
            if a[[i, k]].abs() > max_val {
                max_val = a[[i, k]].abs();
                max_row = i;
            }
        }

        // Перестановка строк
        if max_row != k {
            // Перестановка строк в матрице A
            for col in 0..n {
                let temp = a[[k, col]];
                a[[k, col]] = a[[max_row, col]];
                a[[max_row, col]] = temp;
            }
            // Перестановка элементов в векторе b
            let temp_b = b[k];
            b[k] = b[max_row];
            b[max_row] = temp_b;
        }

        // Проверка на нулевой ведущий элемент
        if a[[k, k]].abs() < 1e-12 {
            println!(
                "Предупреждение: малый диагональный элемент A[{},{}] = {:.2e}",
                k,
                k,
                a[[k, k]]
            );
            a[[k, k]] = 1e-12;
        }

        // Исключение
        for i in (k + 1)..n {
            if a[[i, k]].abs() > 1e-15 {
                let factor = a[[i, k]] / a[[k, k]];
                a[[i, k]] = 0.0;
                for col in (k + 1)..n {
                    a[[i, col]] -= factor * a[[k, col]];
                }
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

fn main() -> Result<(), Box<dyn Error>> {
    // Создаем сетку
    let xk: Vec<f64> = (0..=N).map(|i| i as f64 * H).collect();

    println!("Проверка точного решения в ключевых точках:");
    let test_points = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    for &x in &test_points {
        println!("x = {:.0}: y_toch = {:.2}", x, ytoch(x));
    }

    // Инициализация матрицы A и вектора b
    let mut a_matrix = Array2::zeros((N, N));
    let mut b_vec = Array1::zeros(N);

    // Заполнение матрицы A и вектора b для ВНУТРЕННИХ точек (i=1..n)
    for i in 1..=N {
        let i_idx = i - 1;
        b_vec[i_idx] = f(xk[i]);
        for k in 1..=N {
            let k_idx = k - 1;
            let x = xk[i];
            let val = ddphi_k(x, k) + p(x) * dphi_k(x, k) + q(x) * phi_k(x, k);
            a_matrix[[i_idx, k_idx]] = val;
        }
    }

    println!("\nРазмерность матрицы A: {}x{}", N, N);
    println!("Размерность вектора b: {}", N);

    // Решение системы A * c = b методом Гаусса
    println!("\nРешение СЛАУ методом Гаусса");
    let (u, b_transformed) = forward_elimination(a_matrix, b_vec);
    let c_vec = backward_substitution(&u, &b_transformed);

    println!("Первые 10 коэффициентов a_k:");
    for i in 0..10.min(N) {
        println!("a_{} = {:e}", i + 1, c_vec[i]);
    }

    // Сравнение точного и приближенного решений через polars
    println!("\nСравнение решений в некоторых точках:");
    let mut comparison_data = Vec::new();

    for x in 0..=V as usize {
        let x_f64 = x as f64;
        let y_exact = ytoch(x_f64);
        let y_appr = y_approx(&c_vec, x_f64);
        let rel_error = if y_exact.abs() > 1e-12 {
            ((y_appr - y_exact) / y_exact).abs()
        } else {
            y_appr.abs()
        };

        comparison_data.push((x_f64, y_exact, y_appr, rel_error));
    }

    let df = df![
        "x" => (0..=V as usize).map(|x| x as f64).collect::<Vec<_>>(),
        "Точное y" => (0..=V as usize).map(|x| ytoch(x as f64)).collect::<Vec<_>>(),
        "Приближенное y" => (0..=V as usize).map(|x| y_approx(&c_vec, x as f64)).collect::<Vec<_>>(),
        "Относительная погрешность" => (0..=V as usize)
            .map(|x| {
                let x_f64 = x as f64;
                let y_exact = ytoch(x_f64);
                let y_appr = y_approx(&c_vec, x_f64);
                if y_exact.abs() > 1e-12 {
                    ((y_appr - y_exact) / y_exact).abs()
                } else {
                    y_appr.abs()
                }
            })
            .collect::<Vec<_>>(),
    ]?;

    println!("{}", df);
    polars::env::set_config(polars::env::ConfigOptions {
        table_row_count: Some(100),
        ..Default::default()
    });
    Ok(())
}

fn y_approx(c: &Array1<f64>, x: f64) -> f64 {
    let mut result = 0.0;
    for k in 1..=c.len() {
        result += c[k - 1] * phi_k(x, k);
    }
    result
}
