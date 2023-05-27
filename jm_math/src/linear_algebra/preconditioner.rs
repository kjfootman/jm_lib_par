use rayon::prelude::*;
use crate::linear_algebra::matrix::Matrix;
use crate::linear_algebra::vector::Vector;

pub fn SGS(A: &Matrix, v: &Vector) -> Vector {
    let m = A.num_rows();
    let y = Vector::from(vec![0.0; m]);
    let mut t = Vector::from(vec![0f64; m]);

    // A.IA().par_iter().for_each(|&i| {
    //     let j1 = A.IA()[i];
    //     let j2 = A.IA()[i+1];
    //     let mut j = j1;

    //     t[i] = v[i];

    //     while j < i {
    //         t[i] -= A.AA()[j] * y.AA()[A.JA()[j]];
    //         j += 1;
    //     }
    // });
    t.par_iter_mut().enumerate()
        .for_each(|(i, t)| {
            let j1 = A.IA()[i];
            let j2 = A.IA()[i+1];
            let mut j = j1;

            *t = v[i];

            while j < i {
                *t -= A.AA()[j] * y[A.IA()[j]];
                j += 1;
            }

            *t /= A.AA()[j];
        });
    
    t
}