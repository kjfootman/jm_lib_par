use rayon::prelude::*;
// use std::sync::{Arc, Mutex};
use crate::linear_algebra::matrix::Matrix;
use crate::linear_algebra::vector::Vector;

//-----------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------//
pub fn none<'a>(_A: &Matrix, b: &'a Vector) -> &'a Vector {

    b
}

//-----------------------------------------------------------------------------------------------------------//
pub fn Jacobi(A: &Matrix, b: &Vector) -> Vector {
    let m = A.num_rows();
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

            // let idx = j1;
        });

    x
}

//-----------------------------------------------------------------------------------------------------------//
pub fn Gauss_Seidel(A: &Matrix, b: &Vector) -> Vector {
    let m = A.num_rows();
    let mut x = Vector::from(vec![0f64; m]);

    for i in 0..m {
        let j1 = A.IA()[i];
        // let j2 = A.IA()[i+1];
        let mut j = j1;
        let mut jrow: usize;

        x[i] = b[i];

        loop {
            jrow = A.JA()[j];

            if jrow >= i {
                break;
            }

            x[i] -= A.AA()[j] * x[jrow]; 

            j += 1;
        }

        if jrow != i || A.AA()[j] == 0.0 {
            panic!("diagonal element error");
        }

        x[i] /= A.AA()[j];
    }

    x
}

//-----------------------------------------------------------------------------------------------------------//
pub fn level_schduling(A: &Matrix, b: &Vector) -> Vector {
    // lower part of A
    let m = A.num_rows();
    let mut depth = vec![0usize; m];
    let mut uptr = vec![0usize; m];
    let mut q = (0..m).collect::<Vec<usize>>();
    let mut x = Vector::from(vec![0f64; m]);

    for i in 0..m {
        let j1 = A.IA()[i];
        let mut j = j1;
        let mut jrow: usize; 

        loop {
            jrow = A.JA()[j];

            if jrow >= i {
                break;
            }

            j += 1;
        }

        uptr[i] = j;
        if jrow != i {
            panic!("{}", "diagonal element error");
        }

        // depth for vertices
        depth[i] = 1 + A.JA()[j1..=j].iter().map(|&i| depth[i]).max().unwrap();
    }

    q.par_sort_by_key(|&i| depth[i]);

    let mut tmp = depth[q[0]];
    let mut level = Vec::with_capacity(*q.last().unwrap());

    level.push(0);
    for i in 0..q.len() {
        if depth[q[i]] != tmp {
            tmp = depth[q[i]];
            level.push(i); // 비효율 발생
        }
    }
    level.push(q.len());

    // println!("num_levels: {}", level.len() - 1);
    //

    for lev in 0..level.len() - 1 {
        // todo: parallelize following code,
        // for k in level[lev]..level[lev+1] {
        //     let i = q[k];
        //     x[i] = b[i];
        //     for j in A.IA()[i]..uptr[i] {
        //         x[i] -= A.AA()[j] * x[A.JA()[j]];
        //     }

        //     x[i] /= A.AA()[uptr[i]];
        // }

        let sub = (level[lev]..level[lev+1]).into_par_iter().map(|k| {
            let i = q[k];
            let mut value = b[i];

            for j in A.IA()[i]..uptr[i] {
                value -= A.AA()[j] * x[A.JA()[j]];
            }

            value /= A.AA()[uptr[i]];

            (i, value)
        }).collect::<Vec<_>>();

        sub.iter().for_each(|(i, value)| {
            x[*i] = *value;
        });
    }

// test for rayon
//----------------------------------------//

//----------------------------------------//
    x
}