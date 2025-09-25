use core::f64;
use std::io;
fn main() {
    let mut input = String::new();
    io::stdin()
        .read_line(&mut input)
        .expect("Failed input x_list");
    let x_list: Vec<f64> = input
        .split_whitespace()
        .map(|x| x.trim().parse::<f64>().expect("Failed parse x_list"))
        .collect();

    io::stdin()
        .read_line(&mut input)
        .expect("Failed input f_list");
    let f_list: Vec<f64> = input
        .split_whitespace()
        .map(|f| f.trim().parse::<f64>().expect("Failed parse f_list"))
        .collect();

	let mut numerator = 1.0;
	let mut denominator = 1.0;

	let mut result: Vec<f64> = Vec::new();

    for j in 0..x_list.len() {
		let mut sum = 0.0;
        for k in 0..f_list.len() {
            let f_k = f_list[k];
            for i in 0..x_list.len() {
                if k != i {
					numerator *= f_k * (x_list[j] - x_list[i]);
					denominator *= x_list[k] - x_list[i];
                }
			}
			sum += numerator / denominator;
        }
		result.push(sum);
    }

	println!("{result:?}")
}
