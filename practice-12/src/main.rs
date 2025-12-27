use ndarray::{Array1, Array2};
use std::error::Error;

fn solve_fredholm_degenerate_correct(v: f64) -> Result<(), Box<dyn Error>> {
    let a = 0.0;
    let b = 1.0;
    let n = 3;
    let lam = 1.0;

    let mut alpha = Array2::zeros((n as usize, n as usize));
    let mut gamma = Array1::zeros(n as usize);

    for i in 0..n {
        for k in 0..n {
            let deg = i + k;
            alpha[[i as usize, k as usize]] = 1.0 / (deg as f64 + 1.0);
        }
        let deg1 = 1 + i;
        let deg2 = 2 + i;
        let deg3 = 3 + i;
        gamma[i as usize] = v
            * ((4.0 / 3.0) / (deg1 as f64 + 1.0)
                + 0.25 / (deg2 as f64 + 1.0)
                + 0.2 / (deg3 as f64 + 1.0));
    }

    let eye = Array2::eye(n as usize);
    let mut alpha_t = alpha.clone();
    alpha_t.swap_axes(0, 1);
    let a_matrix = &eye + &(lam * &alpha_t);
    let q = gauss_method(&a_matrix, &gamma)?;

    // Тестовые точки
    let num_points = 18usize;
    let mut x_test = Vec::with_capacity(num_points);
    for i in 0..num_points {
        let x = a + (b - a) * (i as f64) / ((num_points - 1) as f64);
        x_test.push(x);
    }

    let mut y_num_vals = Vec::with_capacity(num_points);
    let mut y_ex_vals = Vec::with_capacity(num_points);
    let mut errors = Vec::with_capacity(num_points);

    // Функции базиса
    let a0 = |x: f64| x;
    let a1 = |x: f64| x * x;
    let a2 = |x: f64| x.powi(3);

    for &x in &x_test {
        // f(x)
        let mut result = v * (4.0 / 3.0 * x + 0.25 * x * x + 0.2 * x.powi(3));
        // Вычитаем поправку
        result -= lam * q[0] * a0(x);
        result -= lam * q[1] * a1(x);
        result -= lam * q[2] * a2(x);

        let y_num = result;
        let y_ex = v * x;

        y_num_vals.push(y_num);
        y_ex_vals.push(y_ex);
        errors.push((y_num - y_ex).abs());
    }

    let max_len = x_test.len();
    println!("┌{0:─<10}┬{0:─<14}┬{0:─<14}┬{0:─<12}┐", "─".repeat(8));
    println!(
        "│ {:<8} │ {:<12} │ {:<12} │ {:<10} │",
        "x", "y_метода", "y_точн", "погреш"
    );
    println!("├{0:─<10}┼{0:─<14}┼{0:─<14}┼{0:─<12}┤", "─".repeat(8));

    for i in 0..max_len {
        println!(
            "│ {:<8.1} │ {:<12.4} │ {:<12.4} │ {:<10.4} │",
            x_test[i], y_num_vals[i], y_ex_vals[i], errors[i]
        );
    }

    println!("└{0:─<10}┴{0:─<14}┴{0:─<14}┴{0:─<12}┘", "─".repeat(8));
    Ok(())
}

fn gauss_method(a: &Array2<f64>, b: &Array1<f64>) -> Result<Array1<f64>, Box<dyn Error>> {
    let n = a.nrows();
    let mut mat = a.to_owned();
    let mut vecb = b.to_owned();

    for i in 0..n {
        // Поиск главного элемента
        let mut max_row = i;
        for k in (i + 1)..n {
            if mat[[k, i]].abs() > mat[[max_row, i]].abs() {
                max_row = k;
            }
        }

        // Перестановка строк
        if max_row != i {
            for col in 0..n {
                let temp = mat[[i, col]];
                mat[[i, col]] = mat[[max_row, col]];
                mat[[max_row, col]] = temp;
            }
            let temp = vecb[i];
            vecb[i] = vecb[max_row];
            vecb[max_row] = temp;
        }

        // Прямой ход
        for k in (i + 1)..n {
            if mat[[i, i]].abs() < 1e-12 {
                return Err("Singular matrix".into());
            }
            let c = -mat[[k, i]] / mat[[i, i]];
            for j in i..n {
                if i == j {
                    mat[[k, j]] = 0.0;
                } else {
                    mat[[k, j]] += c * mat[[i, j]];
                }
            }
            vecb[k] += c * vecb[i];
        }
    }

    // Обратный ход
    let mut x = Array1::zeros(n);
    for i in (0..n).rev() {
        let mut sum = 0.0;
        for j in (i + 1)..n {
            sum += mat[[i, j]] * x[j];
        }
        x[i] = (vecb[i] - sum) / mat[[i, i]];
    }

    Ok(x)
}

fn main() -> Result<(), Box<dyn Error>> {
    solve_fredholm_degenerate_correct(18.0)
}
