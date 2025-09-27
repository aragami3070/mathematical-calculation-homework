use std::process;


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

	println!("{res_x:?}");
	println!("{res_f:?}");
}
