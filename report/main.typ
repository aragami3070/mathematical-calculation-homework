#import "conf.typ" : conf
#show: conf.with(
  title: [Отчёт\ по практической подготовке],
  type: "pract",
  info: (
      author: (
        name: [Смирнова Егора Ильича],
        faculty: [компьютерных наук и информационных технологий],
        group: "351",
        sex: "male"
      ),
      inspector: (
        degree: "",
        name: ""
      )
  ),
  settings: (
    title_page: (
      enabled: true
    ),
    contents_page: (
      enabled: true
    )
  )
)

= Построение интерполяционного многочлена в общем виде

*Условие*

Необходимо найти интерполяционный многочлен в общем виде.

#table(rows: 2, columns: (1fr, 1fr, 1fr, 1fr, 1fr))[*$x$*][
  0][1][2][3][*$f(x)$*][1][2][9][28]

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

```
```rust
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
```
```rust
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
            for i in 0..x_list.len() {
                if k != i {
                    numerator *= x_extended_list[j] - x_list[i];
                    denominator *= x_list[k] - x_list[i];
                }
            }
```
```rust
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
#image("images/03.png")

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
```
```rust
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
#image("images/04.png")


= Интерполяция кубическими сплайнами

*Условие*

Необходимо построить интерполяционный многочлен с помощью кубических сплайнов
(алгебраических многочленов третьей степени, где сплайн --- фрагмент, отрезок чего-либо).


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
```
```rust
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
```
```rust
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

*Результат*
#image("images/05.png")
#image("images/06.png")

= Метод Гаусса решения СЛАУ

*Условие*

Решить следующую СЛАУ методом Гаусса
Метод Гаусса должен решать уравнения вида $A x = b$, где $A$ - матрица.
Для упрощения тестирования матрица $А$ примет вид:

$
  A = mat(18, 0.18, 0.18, 0.18, 0.18;
      0.19, 19, 0.19, 0.19, 0.19;
      0.20, 0.20, 20, 0.20, 0.20;
      0.21, 0.21, 0.21, 21, 0.21;
      0.22, 0.22, 0.22, 0.22, 22) quad
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
```
```rust
    println!("1) Матрица A:\n{:?}", a);
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

*Результат*
#image("images/07.png", width: 70%)
#image("images/08.png", width: 70%)

=	Метод прогонки решения СЛАУ (трехдиагональных)

*Условие*
В данном случае решается система линейных уравнений вида $A x = b$, где $A$ --- матрица вида:

$
    mat(
      -18, 0.18, 0, 0, 0;
      0.18, -19, 0.19, 0, 0;
      0, 0.19, -20, 0.20, 0;
      0, 0, 0.20, -21, 0.21;
      0, 0, 0, 0.21, -22) x = mat(18; 19; 20; 21; 22).
$

*Код*

```rust
use ndarray::{Array1, Array2, array};

fn main() {
    // Исходные данные
    let a: Array2<f64> = array![
        [18.0, 0.18, 0.18, 0.18, 0.18],
        [0.19, 19.0, 0.19, 0.19, 0.19],
        [0.20, 0.20, 20.0, 0.20, 0.20],
        [0.21, 0.21, 0.21, 21.0, 0.21],
        [0.22, 0.22, 0.22, 0.22, 22.0],
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


```
```rust
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
#image("images/09.png")

= Метод простой итерации

*Условие*

При решении СЛАУ вида $A x = b$, где $A$ --- квадратная матрица, мы можем преобразовать ее к эквивалентному виду:

$
  mat(0, - a_12/a_11, ..., -a_(1 n) / a_11;
      -a_21/a_22, 0, ..., -a_(2 n)/a_22;
      dots.v, dots.v, dots.down, dots.v;
    -a_(n 1)/a_(n n), -a_(n 2)/a_(n n), ..., 0
    )
  x = mat(b_1 / a_11; b_2 / a_22; dots.v; b_n / a_(n n)).
$

Таким образом исходная система допускает представление в виде:

$
 alpha x + beta = x,
$

а критерий остановки вычислений:

$
  ||x^(k) - x^(k-1)|| < e.
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
```
```rust
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
```

*Результат*

```
```
