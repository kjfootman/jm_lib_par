use std::fmt;
use std::convert::From;
use std::ops::Mul;
use std::path::PathBuf;
// use std::sync::{Mutex, Arc};
use rayon::prelude::*;
use regex::Regex;
// use crate::vector::Vector;
use crate::linear_algebra::vector::Vector;

pub struct Matrix {
    m: usize,
    n: usize,
    AA: Vec<f64>,
    JA: Vec<usize>,
    IA: Vec<usize>
}

/***********************************************************************************************************/
impl Matrix {
    pub fn new() -> Matrix {
        Matrix { 
            m: 0, 
            n: 0, 
            AA: Vec::new(), 
            JA: Vec::new(), 
            IA: Vec::new() 
        }
    }

//-----------------------------------------------------------------------------------------------------------//
    pub fn import_file(path: &str) -> Matrix {
        let path = PathBuf::from(path);
        
        if path.exists() {
            println!("import matrix from {}", path.display());
            let mut m = usize::MAX; 
            let mut n = usize::MAX; 
            let mut z = usize::MAX;
            let mut re;
            let text = std::fs::read_to_string(path)
                .expect("can not read file");
            let mut data = vec![];
            
            // header
            re = Regex::new(r"%%MatrixMarket matrix coordinate\s(\w+)\s(\w+)").unwrap(); 
            for cap in re.captures_iter(&text) {
                println!("{}, {}", &cap[1], &cap[2]);
            }

            // parse dimension
            re = Regex::new(r"(\d+)\s(\d+)\s(\d+)\n").unwrap(); 
            for cap in re.captures_iter(&text) {
                m = cap[1].trim().parse::<usize>().unwrap();
                n = cap[2].trim().parse::<usize>().unwrap();
                z = cap[3].trim().parse::<usize>().unwrap();
            }
            if m != n {
                println!("caution: not square a matrix");
            }

            // parse each row to the tuple of (i, j, value)
            re = Regex::new(r"(\d+)\s(\d+)\s([-|\s].+)").unwrap(); 
            for cap in re.captures_iter(&text) {
                data.push((
                    cap[1].parse::<usize>().unwrap(),
                    cap[2].parse::<usize>().unwrap(),
                    cap[3].trim().parse::<f64>().unwrap(),
                ));
            }

            // sort by row
            data.sort_by(|a, b| {
                a.0.cmp(&b.0).then(a.1.cmp(&b.1))
            });

            // assemble matrix
            let mut AA = Vec::with_capacity(z);
            let mut JA = Vec::with_capacity(z);
            let mut IA = Vec::with_capacity(m);
            let mut row = usize::MAX;

            for (i, j, value) in &data {
                if i - 1 != row {
                    row = i - 1;
                    IA.push(AA.len());
                }
                AA.push(*value);
                JA.push(j - 1);
            }
            IA.push(data.len());

            println!("done.");
            Matrix::from(AA, JA, IA)
        } else {            
            println!("can not find file {}", path.display());
            Matrix::new()
        }
    }

//-----------------------------------------------------------------------------------------------------------//
    pub fn from<
        T: IntoParallelIterator, 
        U: IntoParallelIterator, 
        V: IntoParallelIterator
    >(AA: T, JA: U, IA: V) 
    -> Self where T::Item: Into<f64>, U::Item: Into<usize>, V::Item: Into<usize> { 
        let AA = AA.into_par_iter()
            .map(|v| v.into())
            .collect::<Vec<f64>>();
        let JA = JA.into_par_iter()
            .map(|v| v.into())
            .collect::<Vec<usize>>();
        let IA = IA.into_par_iter()
            .map(|v| v.into())
            .collect::<Vec<usize>>();

        let m = IA.len() - 1;
        let n = *JA.par_iter().max().unwrap() + 1;

        Matrix {
            m,
            n,
            AA,
            JA,
            IA
        }
    }

//-----------------------------------------------------------------------------------------------------------//
    pub fn num_rows(&self) -> usize {
        self.m
    }

    pub fn num_cols(&self) -> usize {
        self.n
    }

    pub fn AA(&self) -> &Vec<f64> {
        &self.AA
    }

    pub fn JA(&self) -> &Vec<usize> {
        &self.JA
    }
    
    pub fn IA(&self) -> &Vec<usize> {
        &self.IA
    }
}

/***********************************************************************************************************/
impl Mul<&Vector> for &Matrix {
    type Output = Vector;

    fn mul(self, rhs: &Vector) -> Self::Output {
        assert!(self.m == rhs.num_rows());
        
        let m = self.m;
        let mut AA = vec![0f64; m];

        AA.par_iter_mut().enumerate()
            .for_each(|(i, v)| {
                for j in self.IA[i]..self.IA[i+1] {
                    *v += self.AA[j] * rhs[self.JA[j]];
                }
            });

        Vector::from(AA)
    }
}

/***********************************************************************************************************/
impl<I: IntoIterator> From<I> for Matrix 
    where I::Item: IntoIterator, <<I as IntoIterator>::Item as IntoIterator>::Item: Into<f64> {
    fn from(value: I) -> Self {

        for v in value {
            for v2 in v {
                let tmp: f64 = v2.into();
                print!("{} ", tmp);
            }
            println!();
        }

        Matrix::new()
    }
}

/***********************************************************************************************************/
impl fmt::Display for Matrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let m = self.m;
        let n = self.n;

        writeln!(f, "Matrix")?;
        writeln!(f, "dim: [{} x {}]", m, n)?;

        let width = f.width().unwrap_or_else(|| 0);
        let precision = f.precision().unwrap_or_else(|| 0);

        for i in 0..m {
            let mut row = vec![0.0f64; n];
            for j in self.IA[i]..self.IA[i+1] {
                row[self.JA[j]] = self.AA[j];
            }

            for val in row {
                write!(f, "{val:>width$.pre$} ", width = width, pre = precision)?;
            }

            writeln!(f)?;
        }

        Ok(())
    } 
}