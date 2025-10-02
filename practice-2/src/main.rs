use itertools::Itertools;
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
            sum += f_k * numerator / denominator;
            numerator = 1.0;
            denominator = 1.0;
        }
        result.push(sum);
    }

    println!("{result:?}")
}
