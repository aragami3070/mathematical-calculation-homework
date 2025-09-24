use ndarray_interp::interp1d::*;
use ndarray::array;

fn main() {
    // Исходные точки по X
    let x = array![0.0, 1.0, 2.0, 3.0];

    // Значения функции в этих точках
    let y = array![1.0, 2.0, 3.0, 28.0];

    // Создаём интерполятор (линейная интерполяция с экстраполяцией)
    let interp = Interp1DBuilder::new(y)
        .x(x)
        .strategy(Linear::new().extrapolate(true))
        .build()
        .unwrap();

    // Пример интерполяции в точке 1.5
    let value = interp.interp_scalar(1.5).unwrap();

    println!("Интерполированное значение в 1.5: {}", value);
}

