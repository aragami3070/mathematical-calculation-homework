use std::{error::Error, fs, path::Path};

type X = f64;
type F = f64;

pub fn read(name: &str) -> Result<(Vec<X>, Vec<F>), Box<dyn Error>> {
    let file_content = fs::read_to_string(Path::new(name))?;

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
