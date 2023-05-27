use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use crate::linear_algebra::matrix::Matrix;
use crate::linear_algebra::vector::Vector;

pub fn SGS(A: &Matrix, v: &Vector) -> Vector {
    // Symetric Gauss Seidel
    let m = A.num_rows();
    let y = Vector::from(vec![0.0; m]);
    let mut t = Vector::from(vec![0f64; m]);
    
    t
}

//-----------------------------------------------------------------------------------------------------------//
pub fn Jacobi(A: &Matrix, b: &Vector) -> Vector {
    let m = A.num_rows();
    // let mut x = Vector::from(vec![1.0f64; m]);
    let mut x = Vector::from(vec![0f64; m]);

    x.par_iter_mut().enumerate()
        .for_each(|(i, x)| {
            let j1 = A.IA()[i];
            let j2 = A.IA()[i+1];
            let idx = A.JA()[j1..j2].iter().position(|&j| i == j);

            match idx {
                Some(idx) => {
                    *x = b[i];
                    *x /= A.AA()[j1..j2][idx];
                },
                None => panic!("can't find the diagonal index")
            }
        });

    x
}