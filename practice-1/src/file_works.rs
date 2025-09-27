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
