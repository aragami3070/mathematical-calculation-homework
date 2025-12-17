use ndarray::Array1;

// Function f(x, y, V)
fn f(x: f64, y: f64, v: f64) -> f64 {
    2.0 * v * x + v * x.powi(2) - y
}

// Euler method
fn euler_method(x0: f64, y0: f64, h: f64, n: usize, v: f64) -> (Array1<f64>, Array1<f64>) {
    let mut x = Array1::<f64>::zeros(n + 1);
    let mut y = Array1::<f64>::zeros(n + 1);

    x[0] = x0;
    y[0] = y0;

    for i in 0..n {
        x[i + 1] = x[i] + h;
        y[i + 1] = y[i] + h * f(x[i], y[i], v);
    }

    (x, y)
}

// Improved Euler method
fn improved_euler_method(x0: f64, y0: f64, h: f64, n: usize, v: f64) -> (Array1<f64>, Array1<f64>) {
    let mut x = Array1::<f64>::zeros(n + 1);
    let mut y = Array1::<f64>::zeros(n + 1);

    x[0] = x0;
    y[0] = y0;

    for i in 0..n {
        x[i + 1] = x[i] + h;
        let y_half = y[i] + (h / 2.0) * f(x[i], y[i], v);
        let x_half = x[i] + h / 2.0;
        y[i + 1] = y[i] + h * f(x_half, y_half, v);
    }

    (x, y)
}

// Exact solution
fn exact_solution(x: &Array1<f64>, v: f64) -> Array1<f64> {
    x.mapv(|xi| v * xi.powi(2))
}

fn main() {
    let x0 = 1.0;
    let v = 18.0;
    let y0 = v;
    let h = 0.001;
    let n = 10;

    let (x_euler, y_euler) = euler_method(x0, y0, h, n, v);
    let (x_improved, y_improved) = improved_euler_method(x0, y0, h, n, v);
    let y_exact = exact_solution(&x_euler, v);

    let error_euler = &y_euler - &y_exact;
    let error_euler = error_euler.mapv(f64::abs);
    let error_improved = &y_improved - &y_exact;
    let error_improved = error_improved.mapv(f64::abs);

    println!("\nМетод Эйлера:");
    println!("{}", "=".repeat(130));
    print!("x:      ");
    for xi in x_euler.iter() {
        print!("{:10.7} ", xi);
    }
    println!();

    print!("y_M:    ");
    for yi in y_euler.iter() {
        print!("{:10.7} ", yi);
    }
    println!();

    print!("y_T:    ");
    for yi in y_exact.iter() {
        print!("{:10.7} ", yi);
    }
    println!();

    print!("Погрешн:");
    for e in error_euler.iter() {
        print!("{:10.7} ", e);
    }
    println!();
    println!("{}", "=".repeat(130));

    println!("\nУсовершенствованный метод Эйлера:");
    println!("{}", "=".repeat(130));
    print!("x:      ");
    for xi in x_improved.iter() {
        print!("{:10.7} ", xi);
    }
    println!();

    print!("y_M:    ");
    for yi in y_improved.iter() {
        print!("{:10.7} ", yi);
    }
    println!();

    print!("y_T:    ");
    for yi in y_exact.iter() {
        print!("{:10.7} ", yi);
    }
    println!();

    print!("Погрешн:");
    for e in error_improved.iter() {
        print!("{:10.7} ", e);
    }
    println!();
    println!("{}", "=".repeat(130));
}
