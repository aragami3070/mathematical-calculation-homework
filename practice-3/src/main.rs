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

    for x in x_extended_list {
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

    for res in result {
        print!("{res} ")
    }
}
