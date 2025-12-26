#import "conf.typ" : conf
#show: conf.with(title: [Отчёт\ по практической подготовке], type: "pract", info: (author: (
  name: [Смирнова Егора Ильича],
  faculty: [компьютерных наук и информационных технологий],
  group: "351",
  sex: "male",
), inspector: (degree: "", name: "")), settings: (title_page: (enabled: true), contents_page: (enabled: true)))

= Построение интерполяционного многочлена в общем виде

*Условие*

Необходимо найти интерполяционный многочлен в общем виде.

#table(rows: 2, columns: (1fr, 1fr, 1fr, 1fr, 1fr))[*$x$*][
  0
][1][2][3][*$f(x)$*][1][2][9][28]

*Код*

```rust
use std::{error::Error, fs, path::Path};
type X = f64;
type F = f64;
pub fn read(path: &str) -> Result<(Vec<X>, Vec<F>), Box<dyn Error>> {
    let file_content = fs::read_to_string(Path::new(path))?;
    let split_content = file_content.split('\n').collect::<Vec<&str>>();
    let x_vec: Vec<X> = split_content[0]
        .split_whitespace()
        .map(|s| s.parse::<X>())
        .collect::<Result<_, _>>()?;
    let f_vec: Vec<X> = split_content[1]
        .split_whitespace()
        .map(|s| s.parse::<X>())
        .collect::<Result<_, _>>()?;
    Ok((x_vec, f_vec))
}
pub fn write(path: &str, x_vec: &[f64], f_vec: &[f64]) -> Result<(), Box<dyn Error>>{
  let mut x_str = "X: ".to_string();
  for x in x_vec {
    x_str.push_str(&x.to_string());
  }
  fs::write(path, x_str)?;
  let mut f_str = "F: ".to_string();
  for f in f_vec {
    f_str.push_str(&f.to_string());
    f_str.push(' ');
  }
  fs::write(path, f_str)?;
  Ok(())
}
use ndarray::{Array, Array2, ArrayBase, Axis, Dim, OwnedRepr};
use ndarray_linalg::Solve;
fn create_matrix(x_vec: &[f64]) -> ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>> {
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
fn solve_slay(x_vec: &[f64], f_vec: Vec<f64>) -> Vec<f64> {
    let matrix = create_matrix(x_vec);
    // Заполняем матрицу СЛАУ (правая часть матрицы)
    let f_column = Array::from_vec(f_vec);
```
```rust
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
        prev_x = *current_x;
        res_x_vec.push(mid_x);
        res_f_vec.push(new_mid_f);
        res_x_vec.push(*current_x);
        res_f_vec.push(f_vec[index]);
    }
    (res_x_vec, res_f_vec)
}
use std::process;
use crate::file_works::write;
mod file_works;
mod solve;
fn main() {
    // Считываем таблицу из файла
    let (x_vec, f_vec) = match file_works::read("assets/input.txt") {
        Ok(pair) => pair,
        Err(err) => {
            eprintln!("Error: {err}");
            process::exit(1);
        }
    };
    let (res_x, res_f) = solve::vandermonde_interpolation(x_vec, f_vec);
    if let Err(err) = write("assets/output.txt", &res_x, &res_f) {
        eprintln!("Error: {err}");
        process::exit(1);
    }
    println!("Результат:");
    println!("{res_x:?}");
    println!("{res_f:?}")
}
```
#pagebreak()
*Результат*

Входные
#image("images/01.png")
#image("images/02.png")

= Интерполяционный многочлен в форме Лагранжа

*Условие*

По данным интерполяции из предыдущего задания построить интерполяционный многочлен
в форме Лагранжа.

*Код*

```rust
use itertools::Itertools;
use std::io;
fn main() {
    println!("Введите x:");
    let mut input = String::new();
    io::stdin()
        .read_line(&mut input)
        .expect("Failed input x_list");
    let x_list: Vec<f64> = input
        .split_whitespace()
        .map(|x| x.trim().parse::<f64>().expect("Failed parse x_list"))
        .collect();
    println!("Введите f:");
    input.clear();
    io::stdin()
        .read_line(&mut input)
        .expect("Failed input f_list");
    let f_list: Vec<f64> = input
        .split_whitespace()
        .map(|f| f.trim().parse::<f64>().expect("Failed parse f_list"))
        .collect();
    let x_extended_list: Vec<f64> = x_list
        .iter()
        .tuple_windows()
        .flat_map(|(x_0, x_1)| [*x_0, (x_0 + x_1) / 2.0])
        .chain(x_list.last().copied())
        .collect();
    let mut numerator = 1.0;
    let mut denominator = 1.0;
    let mut result: Vec<f64> = Vec::new();
    for j in 0..x_extended_list.len() {
        let mut sum = 0.0;
        for k in 0..f_list.len() {
            let f_k = f_list[k];
```
```rust
            for i in 0..x_list.len() {
                if k != i {
                    numerator *= x_extended_list[j] - x_list[i];
                    denominator *= x_list[k] - x_list[i];
                }
            }
            sum += f_k * numerator / denominator;
            numerator = 1.0;
            denominator = 1.0;
        }
        result.push(sum);
    }
    println!("Результат:");
    println!("x: {x_extended_list:?}");
    println!("f: {result:?}")
}
```

*Результат*
#image("images/03.png", width: 60%)

= Интерполяционный многочлен в форме Ньютона

*Условие*

По данным интерполяции из предыдущего задания построить интерполяционный многочлен
в форме Ньютона.

*Код*

```rust
use itertools::Itertools;
use std::io;
fn main() {
    let mut input = String::new();
    println!("Введите x:");
    io::stdin().read_line(&mut input).expect("Failed input x_list");
    let x_list: Vec<f64> = input
        .split_whitespace()
        .map(|x| x.trim().parse::<f64>().expect("Failed parse x_list"))
        .collect();
    println!("Введите f:");
    input.clear();
    io::stdin() .read_line(&mut input) .expect("Failed input f_list");
    let f_list: Vec<f64> = input
        .split_whitespace()
        .map(|f| f.trim().parse::<f64>().expect("Failed parse f_list"))
        .collect();
    let x_extended_list: Vec<f64> = x_list
        .iter()
        .tuple_windows()
        .flat_map(|(x_0, x_1)| [*x_0, (x_0 + x_1) / 2.0])
        .chain(x_list.last().copied())
        .collect();
```
```rust
    let mut result: Vec<f64> = Vec::new();
    let mut coefficients: Vec<Vec<f64>> = Vec::new();
    let f_list_len = f_list.len();
    for ind in 0..f_list_len {
        if ind == 0 {
            coefficients.push(f_list.clone());
            continue;
        }
        let mut new_vec: Vec<f64> = Vec::new();
        for i in 0..f_list_len - ind {
            let numerator = coefficients[ind - 1][i + 1] - coefficients[ind - 1][i];
            let denominator = x_list[i + ind] - x_list[i];
            new_vec.push(numerator / denominator);
        }
        coefficients.push(new_vec);
    }
    for x in x_extended_list.clone().iter() {
        let mut end_ind = 0;
        let mut sum = 0.0;
        for sub_list in coefficients.iter() {
            let mut mult_x = 1.0;
            for (i, x_i) in x_list.iter().enumerate() {
                if i == end_ind {
                    break;
                }
                mult_x *= x - x_i
            }
            sum += sub_list[0] * mult_x;
            end_ind += 1;
        }
        result.push(sum);
    }
    println!("Таблица входных:");
    println!("x: {x_extended_list:?}");
    println!("f: {result:?}");
    println!("Результат:");
    println!("x: {x_extended_list:?}");
    println!("N_(x): {result:?}")
}
```

*Результат*
#image("images/04.png", width: 55%)

= Интерполяция кубическими сплайнами

*Условие*

По данным интерполяции из предыдущего задания построить кусочно-непрерывную
склейку кубических сплайнов

#pagebreak()
*Код*

```rust
use ndarray_linalg::Solve;
/// создает матрицу для нахождения коэффициентов кубических сплайнов
fn create_spline_matrix(x: &Array1<f64>, f: &Array1<f64>) -> (Array2<f64>, Array1<f64>) {
    let n = x.len() - 1; // количество интервалов
    let size = 4 * n;
    let h = x[1] - x[0];
    let mut a = Array2::<f64>::zeros((size, size));
    let mut b = Array1::<f64>::zeros(size);
    for i in 0..n {
        let row = i;
        a[[row, 4 * i]] = 1.0; // a_i
        b[row] = f[i];
    }
    for i in 0..n {
        let row = n + i;
        a[[row, 4 * i]] = 1.0;
        a[[row, 4 * i + 1]] = h;
        a[[row, 4 * i + 2]] = h.powi(2);
        a[[row, 4 * i + 3]] = h.powi(3);
        b[row] = f[i + 1];
    }
    for i in 0..(n - 1) {
        let row = 2 * n + i;
        a[[row, 4 * i + 1]] = 1.0;
        a[[row, 4 * i + 2]] = 2.0 * h;
        a[[row, 4 * i + 3]] = 3.0 * h.powi(2);
        a[[row, 4 * (i + 1) + 1]] = -1.0;
        b[row] = 0.0;
    }
    for i in 0..(n - 1) {
        let row = 3 * n - 1 + i;
        a[[row, 4 * i + 2]] = 2.0;
        a[[row, 4 * i + 3]] = 6.0 * h;
        a[[row, 4 * (i + 1) + 2]] = -2.0;
        b[row] = 0.0;
    }
    let row1 = 4 * n - 2;
    a[[row1, 2]] = 2.0; // 2*c_0
    let row2 = 4 * n - 1;
    a[[row2, 4 * (n - 1) + 2]] = 2.0;
    a[[row2, 4 * (n - 1) + 3]] = 6.0 * h;
    (a, b)
}
/// решает систему уравнений
fn solve_spline_coefficients(
    x: &Array1<f64>,
    f: &Array1<f64>,
) -> ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>, f64> {
    let (a, b) = create_spline_matrix(x, f);
    let coeffs = a.solve_into(b).expect("Failed to solve system");
    let n = x.len() - 1;
    coeffs.into_shape_clone((n, 4)).expect("Reshape error")
}
```
```rust
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
        let t_mid = dx / 2.0;
        let s_mid = coeff_row[0]
            + coeff_row[1] * t_mid
            + coeff_row[2] * t_mid.powi(2)
            + coeff_row[3] * t_mid.powi(3);
        println!("S_{:.1}: {:.3}", xi + t_mid, s_mid);
        let t_end = dx;
        let s_end = coeff_row[0]
            + coeff_row[1] * t_end
            + coeff_row[2] * t_end.powi(2)
            + coeff_row[3] * t_end.powi(3);
        println!("S_{:.1}: {:.3}", xip1, s_end);
    }
}
```

#pagebreak()
*Результат*
#image("images/05.png", width: 80%)
= Метод Гаусса решения СЛАУ

*Условие*

Решить следующую СЛАУ методом Гаусса:

$
  A = mat(
    18, 0.18, 0.18, 0.18, 0.18;0.19, 19, 0.19, 0.19, 0.19;0.20, 0.20, 20, 0.20, 0.20;0.21, 0.21, 0.21, 21, 0.21;0.22, 0.22, 0.22, 0.22, 22
  ) quad
  b = mat(18;19;20;21;22).
$

*Код*
```rust
use ndarray::{Array1, Array2, s};
use ndarray_linalg::Determinant;
fn forward_elimination(mut a: Array2<f64>, mut b: Array1<f64>) -> (Array2<f64>, Array1<f64>) {
    let n = b.len();
    for k in 0..(n - 1) {
        for i in (k + 1)..n {
            if a[[i, k]] != 0.0 {
                let factor = a[[i, k]] / a[[k, k]];
                for col in (k + 1)..n {
                    let a_ik = a[[i, k]];
```
```rust
                    let a_kk = a[[k, k]];
                    let factor = a_ik / a_kk;
                    a[[i, col]] -= factor * a[[k, col]];
                }
                a[[i, k]] = 0.0;
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
    println!("Матрица A:\n{:?}", a);
    let det = a.det().unwrap();
    println!("Определитель матрицы A = {}", det);
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
```

#pagebreak()
*Результат*
```
Матрица A:
[[18.0, 0.18, 0.18, 0.18, 0.18],
 [0.19, 19.0, 0.19, 0.19, 0.19],
 [0.2, 0.2, 20.0, 0.2, 0.2],
 [0.21, 0.21, 0.21, 21.0, 0.21],
 [0.22, 0.22, 0.22, 0.22, 22.0]]

Определитель матрицы A = 3156982.6488520275

Колонка b:
[18]
[19]
[20]
[21]
[22]

Матрица после прямого прохода:
[[18.0, 0.18, 0.18, 0.18, 0.18],
 [0.0, 18.9981, 0.1881, 0.1881, 0.1881],
 [0.0, 0.0, 19.996039603960398, 0.19603960396039605, 0.19603960396039605],
 [0.0, 0.0, 0.0, 20.993823529411767, 0.20382352941176468],
 [0.0, 0.0, 0.0, 0.0, 21.991456310679613]]

Решение (A|b) методом Гаусса
Вектор решения после обратного прохода:
x0 = 0.9615384615384615
x1 = 0.9615384615384613
x2 = 0.9615384615384615
x3 = 0.9615384615384612
x4 = 0.9615384615384615
x =
[0.9615384615384615]
[0.9615384615384613]
[0.9615384615384615]
[0.9615384615384612]
[0.9615384615384615]
```

= Метод прогонки решения СЛАУ (трехдиагональных)

*Условие*

Решить следующую СЛАУ методом прогонки
$
  mat(-18, 0.18, 0, 0, 0;0.18, -19, 0.19, 0, 0;0, 0.19, -20, 0.20, 0;0, 0, 0.20, -21, 0.21;0, 0, 0, 0.21, -22) x = mat(18;19;20;21;22).
$

*Код*

```rust
use ndarray::{Array1, Array2, array};
fn main() {
    // Исходные данные
    let a: Array2<f64> = array![
        [18.0, 0.18, 0.0, 0.0, 0.0],
        [0.19, 19.0, 0.19, 0.0, 0.0],
```
```rust
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
```

*Результат*
```
Матрица A:
[[18.0, 0.18, 0.0, 0.0, 0.0],
 [0.19, 19.0, 0.19, 0.0, 0.0],
 [0.0, 0.0, 20.0, 0.2, 0.2],
 [0.0, 0.0, 0.21, 21.0, 0.21],
 [0.0, 0.0, 0.0, 0.22, 22.0]]

Столбец b:
[327.42, 368.22, 408.59999999999997, 449.82, 488.62]
Прямая прогонка
Список Pi и Qi
[-0.01, -0.010001000100010001, -0.01, -0.01000100010001]
[18.19, 19.2000200020002, 20.43, 21.217821782178216, 22.00002200440088]
Обратная прогонка
x 5 = 22.00002200440088
x 4 = 20.997799559911982
x 3 = 20.22002200440088
x 2 = 18.997799559911982
x 1 = 18.000022004400883
[18.000022004400883, 18.997799559911982, 20.22002200440088, 20.997799559911982, 22.00002200440088],
```

= Метод простой итерации

*Условие*

Решить СЛАУ методом простой итерации

При решении СЛАУ $A x = b$, где $A$ --- квадратная матрица, которую можно
преобразовать к виду:
$
  mat(
    0, -a_12 / a_11, dots, -a_(1 n) / a_11;-a_21 / a_22, 0, dots, -a_(2 n) / a_22;dots.v, dots.v, dots.down, dots.v;-a_(n 1)/a_(n n), -a_(n 2)/a_(n n), dots, 0
  )
  x = mat(b_1 / a_11;b_2 / a_22;dots.v;b_n / a_(n n)).
$

*Код*
```rust
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
    let b_diag: Array1<f64> = a.diag().to_owned();
    let b = a.dot(&b_diag);
    println!("Столбец b:");
    println!("{:?}", b);
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
    let mut beta = b.clone();
```
```rust
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
        println!("x^({}) = {:?}", i + 1, xkp1);
        if norm_stop(&xk, &xkp1, epsilon) {
            break;
        }
        xk = xkp1;
    }
}
```

*Результат*
```
Матрица A:
[[18.0, 0.18, 0.18, 0.18, 0.18],
 [0.19, 19.0, 0.19, 0.19, 0.19],
 [0.2, 0.2, 20.0, 0.2, 0.2],
 [0.21, 0.21, 0.21, 21.0, 0.21],
 [0.22, 0.22, 0.22, 0.22, 22.0]]
Определитель матрицы A: 3156982.6488520275
Столбец b:
[338.76, 376.39000000000004, 416.0, 457.59000000000003, 501.16]

Матрица alpha:
[[0.0, -0.01, -0.01, -0.01, -0.01],
 [-0.01, 0.0, -0.01, -0.01, -0.01],
 [-0.01, -0.01, 0.0, -0.01, -0.01],
 [-0.01, -0.01, -0.01, 0.0, -0.01],
 [-0.01, -0.01, -0.01, -0.01, 0.0]]
Столбец beta:
[18.82, 19.810000000000002, 20.8, 21.790000000000003, 22.78]
Считаем до точности epsilon= 0.000000001
x^(0) = [0.0, 0.0, 0.0, 0.0, 0.0]
x^(1) = [18.82, 19.810000000000002, 20.8, 21.790000000000003, 22.78]
x^(2) = [17.9682, 18.968100000000003, 19.968, 20.967900000000004, 21.9678]
x^(3) = [18.001282, 19.001281000000002, 20.00128, 21.001279000000004, 22.001278000000003]
x^(4) = [17.99994882, 18.999948810000003, 19.999948800000002, 20.99994879, 21.99994878]
x^(5) = [18.0000020482, 19.0000020481, 20.000002048, 21.000002047900004, 22.000002047800002]
x^(6) = [17.999999918082, 18.999999918081002, 19.99999991808, 20.999999918079002, 21.999999918078]
x^(7) = [18.00000000327682, 19.000000003276813, 20.0000000032768, 21.000000003276792, 22.00000000327678]
x^(8) = [17.99999999986893, 18.99999999986893, 19.99999999986893, 20.99999999986893, 21.99999999986893]
x^(9) = [18.000000000005244, 19.000000000005244, 20.000000000005244, 21.000000000005244, 22.000000000005244]
```

= Метод Эйлера и усовершенствованный метод Эйлера
#pagebreak()
*Условие:*

Решить задачу Коши методом Эйлера и усовершенствованным методом эйлера
$
  cases(
    y'(x) = 2 V x + V x^2 - y(x),
    y(1) = V
  )
$

где $ y_"точн"(x) = V x^2$, $V$ — номер варианта.

*Код:*
```rust
use ndarray::Array1;
fn f(x: f64, y: f64, v: f64) -> f64 {
    2.0 * v * x + v * x.powi(2) - y
}
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
```
```rust
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
```

*Результат:*
```
Метод Эйлера:
==================================================================================================================================
x:       1.0000000  1.0010000  1.0020000  1.0030000  1.0040000  1.0050000  1.0060000  1.0070000  1.0080000  1.0090000  1.0100000
y_M:    18.0000000 18.0360000 18.0720360 18.1081081 18.1442161 18.1803602 18.2165403 18.2527564 18.2890085 18.3252966 18.3616208
y_T:    18.0000000 18.0360180 18.0720720 18.1081620 18.1442880 18.1804500 18.2166480 18.2528820 18.2891520 18.3254580 18.3618000
Погрешн: 0.0000000  0.0000180  0.0000360  0.0000539  0.0000719  0.0000898  0.0001077  0.0001256  0.0001435  0.0001614  0.0001792
==================================================================================================================================

Усовершенствованный метод Эйлера:
==================================================================================================================================
x:       1.0000000  1.0010000  1.0020000  1.0030000  1.0040000  1.0050000  1.0060000  1.0070000  1.0080000  1.0090000  1.0100000
y_M:    18.0000000 18.0360180 18.0720720 18.1081620 18.1442880 18.1804500 18.2166480 18.2528820 18.2891520 18.3254580 18.3618000
```
```
y_T:    18.0000000 18.0360180 18.0720720 18.1081620 18.1442880 18.1804500 18.2166480 18.2528820 18.2891520 18.3254580 18.3618000
Погрешн: 0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
==================================================================================================================================
```

= Разностный метод
*Условие:*

Решить краевую задачу разностным метдом
$
  cases(
    y'' + x^2y' + x y = 4 V x^4 - 3 V T x^3 + 6 V x - 2 V T,
    y'(0) = y(T) = 0,
    V = T = 18
  )
$

*Код:*
```rust
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
```
```rust
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
```

*Результат:*
```
x                 y        exact             e
0.00           0.00        -0.00    0.00000000
1.80       -1121.72      -944.78  176.93764366
3.60       -3342.43     -3359.23   16.80414160
5.40       -6968.68     -6613.49  355.19096528
7.20      -10062.34    -10077.70   15.35246882
9.00      -13674.70    -13122.00  552.70179365
10.80     -15105.59    -15116.54   10.95700167
12.60     -16185.13    -15431.47  753.65305586
14.40     -13431.26    -13436.93    5.66892704
16.20      -9458.86     -8503.06  955.80164848
18.00          0.00         0.00    0.00000000

Максимальный e: 955.801648478242
Номер максимального е: 9
```
= Краевая задача методом неопределенных коэффициентов
*Условие:*

Решить краевую задачу методом неопределенных коэффициентов.
$
  cases(
    y'' + x^2y' + x y = 4 V x^4 - 3 V T x^3 + 6 V x - 2 V T,
    y'(0) = y(T) = 0,
    y_#[точн] (x) = V x^2 (x - T),
    T = V
  )
$
*Код:*
```rust
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
            for col in 0..n {
                let temp = a[[k, col]];
                a[[k, col]] = a[[max_row, col]];
                a[[max_row, col]] = temp;
            }
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
    let xk: Vec<f64> = (0..=N).map(|i| i as f64 * H).collect();
    println!("Проверка точного решения в ключевых точках:");
    let test_points = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    for &x in &test_points {
        println!("x = {:.0}: y_toch = {:.2}", x, ytoch(x));
    }
    let mut a_matrix = Array2::zeros((N, N));
    let mut b_vec = Array1::zeros(N);
    // Заполнение матрицы A и вектора b для внутренних точек (i=1..n)
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
    println!("\nРешение СЛАУ методом Гаусса");
    let (u, b_transformed) = forward_elimination(a_matrix, b_vec);
    let c_vec = backward_substitution(&u, &b_transformed);
    println!("Первые 10 коэффициентов a_k:");
    for i in 0..10.min(N) {
        println!("a_{} = {:e}", i + 1, c_vec[i]);
    }
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
```
*Результат:*
```
```

= Метод неопределенных коэффициентов
Решите следующее интегральное уравнение:
$
  y(x) + 1 * integral_0^1 (x t + x^2 t^2 + x^3 t^3) y(t) d t
    = V (4/3 x + 1/4 x^2 + 1/5 x^3)
$

*Код:*
```rust
```
*Результат:*

= Метод квадратур
Решите следующее интегральное уравнение:
$
  y(x) + 1 * integral_0^1 (x t + x^2 t^2 + x^3 t^3) y(t) d t
    = V (4/3 x + 1/4 x^2 + 1/5 x^3)
$

*Код:*
```rust
```
*Результат:*

