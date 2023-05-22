use std::fmt;
use std::convert::From;
use std::ops::Mul;
use std::sync::{Mutex, Arc};
use rayon::prelude::*;
use crate::vector::Vector;

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
        let m = self.m;
        let mid = m / 2;
        let mut AA = vec![0f64; m];
        let data = Arc::new(Mutex::new(&mut AA));

        rayon::join(|| {
            let data = Arc::clone(&data);
            let mut v = data.lock().unwrap();

            for i in 0..mid {
                for j in self.IA[i]..self.IA[i+1] {
                    v[i] += self.AA[j] * rhs[self.JA[j]];
                }
            }
        }, || {
            let data = Arc::clone(&data);
            let mut v = data.lock().unwrap();

            for i in mid..m {
                for j in self.IA[i]..self.IA[i+1] {
                    v[i] += self.AA[j] * rhs[self.JA[j]];
                }
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