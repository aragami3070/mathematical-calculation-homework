use ndarray::{Array1, Array2};
use libm::powf;

fn rhs_func(x: f32, variant: f32) -> f32 {
    variant * (4.0 / 3.0 * x + 0.25 * powf(x, 2.0) + 0.2 * powf(x, 3.0))
}

fn simpson_integrate<F>(mut f: F, a: f32, b: f32, n: usize) -> f32
where
    F: FnMut(f32) -> f32,
{
    let h = (b - a) / n as f32;
    let mut sum = f(a) + f(b);
    let mut x = a + h;

    for i in 1..n {
        let weight = if i % 2 == 0 { 2.0 } else { 4.0 };
        sum += weight * f(x);
        x += h;
    }
    sum * h / 3.0
}

fn build_alpha(size: usize) -> Array2<f32> {
    let mut alpha = Array2::zeros((size, size));
    for i in 0..size {
        for j in 0..size {
            let integrand = |t: f32| -> f32 {
                let ai = match i {
                    0 => t,
                    1 => powf(t, 2.0),
                    _ => powf(t, 3.0),
                };
                let bj = match j {
                    0 => t,
                    1 => powf(t, 2.0),
                    _ => powf(t, 3.0),
                };
                ai * bj
            };
            alpha[[i, j]] = simpson_integrate(integrand, 0.0, 1.0, 1000);
        }
    }
    alpha
}

fn build_gamma(size: usize, variant: f32) -> Array1<f32> {
    let mut gamma = Array1::zeros(size);
    for i in 0..size {
        let integrand = |t: f32| -> f32 {
            let bi = match i {
                0 => t,
                1 => powf(t, 2.0),
                _ => powf(t, 3.0),
            };
            rhs_func(t, variant) * bi
        };
        gamma[i] = simpson_integrate(integrand, 0.0, 1.0, 1000);
    }
    gamma
}

fn gauss_method(mut mat: Array2<f32>, mut vecb: Array1<f32>) -> Array1<f32>{
    let n = mat.nrows();
    for i in 0..n {
        let mut max_row = i;
        for k in (i + 1)..n {
            if mat[[k, i]].abs() > mat[[max_row, i]].abs() {
                max_row = k;
            }
        }

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
        for k in (i + 1)..n {
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

    let mut x = Array1::zeros(n);
    for i in (0..n).rev() {
        let mut sum = 0.0;
        for j in (i + 1)..n {
            sum += mat[[i, j]] * x[j];
        }
        x[i] = (vecb[i] - sum) / mat[[i, i]];
    }
    x
}

fn fredholm_solver(variant: f32, rank: usize) -> (Array1<f32>, Array1<f32>, Array1<f32>, Array1<f32>) {
    let alpha = build_alpha(rank);
    let gamma = build_gamma(rank, variant);
    let mut system_matrix = Array2::eye(rank);
    system_matrix += &alpha;

    let coeffs = gauss_method(system_matrix, gamma);

    let step = 0.1;
    let x_vals: Vec<f32> = (0..=10).map(|i| i as f32 * step).collect();
    let x_vals_arr = Array1::from_vec(x_vals.clone());

    let mut y_num = Array1::zeros(x_vals.len());
    for (i, &x) in x_vals.iter().enumerate() {
        let mut y_val = rhs_func(x, variant);
        for j in 0..rank {
            let aj = match j {
                0 => x,
                1 => powf(x, 2.0),
                _ => powf(x, 3.0),
            };
            y_val -= coeffs[j] * aj;
        }
        y_num[i] = y_val;
    }

    let y_true = &x_vals_arr * variant;
    let err = (&y_num - &y_true).mapv(f32::abs);

    (x_vals_arr, y_num, y_true, err)
}

fn main() {
    let v = 18.0;
    println!("Вычисление для варианта {}", v as i32);

    let (x, y_calc, y_exact, error) = fredholm_solver(v, 3);

    println!("\nРешение интегрального уравнения Фредгольма (вырожденное ядро)");
    println!("x:      {}",
        x.iter()
         .map(|xi| format!("{:.6}", xi))
         .collect::<Vec<_>>()
         .join(" "));
    println!("y_мет:  {}",
        y_calc.iter()
              .map(|yi| format!("{:.6}", yi))
              .collect::<Vec<_>>()
              .join(" "));
    println!("y_точн: {}",
        y_exact.iter()
               .map(|yi| format!("{:.6}", yi))
               .collect::<Vec<_>>()
               .join(" "));
    println!("погреш: {}",
        error.iter()
             .map(|ei| format!("{:.6}", ei))
             .collect::<Vec<_>>()
             .join(" "));
}
