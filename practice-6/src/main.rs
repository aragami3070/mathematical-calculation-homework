use ndarray::{Array1, Array2, array};

fn main() {
    // Исходные данные
    let a: Array2<f64> = array![
        [18.0, 0.18, 0.0, 0.0, 0.0],
        [0.19, 19.0, 0.19, 0.0, 0.0],
        [0.0, 0.0, 20.0, 0.20, 0.20],
        [0.0, 0.0, 0.21, 21.0, 0.21],
        [0.0, 0.0, 0.0, 0.22, 22.0],
    ];

    // Вычисляем вектор b
    let n = a.shape()[0];
    let mut b = Array1::<f64>::zeros(n);
    for i in 0..n {
        b[i] = a[[i, i]];
    }

    // b = dot(A, b) in column vector form
    let b = a.dot(&b);

    println!("Матрица A:");
    println!("{:?}", a);
    println!("Столбец b:");
    println!("{:?}", b);

    let mut q = Array1::<f64>::zeros(n);
    let mut p = Array1::<f64>::zeros(n - 1);

    p[0] = -a[[0, 1]] / a[[0, 0]];
    q[0] = b[0] / a[[0, 0]];

    println!("Прямая прогонка");
    println!("Список Pi и Qi");

    for i in 1..p.len() {
        p[i] = a[[i, i + 1]] / (-a[[i, i]] - a[[i, i - 1]] * p[i - 1]);
    }
    for i in 1..q.len() {
        q[i] = (a[[i, i - 1]] * q[i - 1] - b[i]) / (-a[[i, i]] - a[[i, i - 1]] * p[i - 1]);
    }
    println!("{:?},\n{:?}", p, q);

    println!("Обратная прогонка");
    let mut x = Array1::<f64>::zeros(n);
    x[n - 1] = q[n - 1];
    println!("x 5 = {}", x[n - 1]);

    for i in (0..n - 1).rev() {
        x[i] = p[i] * x[i + 1] + q[i];
        println!("x {} = {}", i + 1, x[i]);
    }

    println!("{:?}", x);
}
