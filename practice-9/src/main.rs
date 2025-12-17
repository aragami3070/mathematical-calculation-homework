fn main() {
    let y = move |x: f64, v: f64| v * x * x * (x - v);
    let derivative = move |x: f64, v: f64| {
        let x3 = x * x * x;
        let v2 = v * v;
        -(4.0 * v * x3 * x - 3.0 * v2 * x3 + 6.0 * v * x - 2.0 * v2)
    };
    let p = move |x: f64| -x * x;
    let q = move |x: f64| -x;

    let v = 18.0;

    let n = 10;
    let x_0 = 0.0;
    let h = v / n as f64;

    let x: Vec<f64> = (0..=n).map(|i| x_0 + i as f64 * h).collect();
    let exact: Vec<f64> = x.iter().map(|x_i| y(*x_i, v)).collect();

    let mut f = vec![0.0; n + 1];
    let mut s = vec![0.0; n + 1];
    let mut t = vec![0.0; n + 1];
    let mut r = vec![0.0; n + 1];
    let mut f1 = vec![0.0; n + 1];
    let mut s1 = vec![0.0; n + 1];
    let mut y_vec = vec![0.0; n + 1];
    let mut e = vec![0.0; n + 1];

    for i in 1..n {
        f[i] = 0.5 * (1.0 + 0.5 * h * p(x[i]));
        s[i] = 0.5 * (1.0 - 0.5 * h * p(x[i]));
        t[i] = 1.0 + 0.5 * h * h * q(x[i]);
        r[i] = 0.5 * h * h * derivative(x[i], v);
    }

    f1[1] = 0.0;
    s1[1] = 0.0;

    for i in 1..n {
        let denom = t[i] - f[i] * f1[i];
        f1[i + 1] = s[i] / denom;
        s1[i + 1] = (r[i] + f[i] * s1[i]) / denom;
    }

    y_vec[n] = 0.0;
    for i in (1..=n - 1).rev() {
        y_vec[i] = f1[i + 1] * y_vec[i + 1] + s1[i + 1];
    }

    let mut max_e = 0.0;
    let mut max_e_index = 0;
    for i in 0..=n {
        e[i] = (y_vec[i] - exact[i]).abs();
        if e[i] > max_e {
            max_e = e[i];
            max_e_index = i
        }
    }

    println!("{:<6} {:>12} {:>12} {:>13}", "x", "y", "exact", "e");
    for i in 0..=n {
        println!(
            "{:<6.2} {:>12.2} {:>12.2} {:>13.8}",
            x[i], y_vec[i], exact[i], e[i]
        );
    }
    println!("Максимальный e: {}", max_e);
    println!("Номер максимального е: {}", max_e_index)
}
