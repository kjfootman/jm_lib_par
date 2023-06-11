use core::panic;

use rayon::prelude::*;
// use std::sync::{Arc, Mutex};
use crate::linear_algebra::matrix::Matrix;
use crate::linear_algebra::vector::Vector;

pub enum Preconditioner {
    Jacobi,
    GS,
    SGS,
    SOR(f64),
    SSOR(f64),
    ILU,
    None
}

impl Preconditioner {
    pub fn from(&self, A: &Matrix) -> Option<Matrix> {
        match self {
            Preconditioner::GS => {
                println!("Gauss-Seidel preconditioner");
                Some(GS(A))
            },
            Preconditioner::None => {
                println!("preconditioner not used");
                None
            },
            _ => panic!("not available preconditioner")
        }
    }
}

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

    println!("num_levels: {}", level.len() - 1);
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

    x
}

//-----------------------------------------------------------------------------------------------------------//
pub fn GS(A: &Matrix) -> Matrix {
    // Gauss Seidel preconditioner
    // return lower part of A as a preconditioner

    let m = A.num_rows();
    let mut AA = Vec::with_capacity(A.AA().len());
    let mut JA = Vec::with_capacity(A.JA().len());
    let mut IA = vec![0usize; m+1];
    let mut UPTR = vec![usize::MAX; m];

    for i in 0..m {
        let j1 = A.IA()[i];
        let mut j = j1;
        let mut jrow;

        // for lower elements
        loop {
            jrow = A.JA()[j];

            if jrow >= i {
                break;
            }

            AA.push(A.AA()[j] * AA[UPTR[jrow]]);
            JA.push(jrow);

            j += 1;
        }

        // for upper elements
        if jrow != i || A.AA()[j] == 0f64 {
            panic!("diagonal element error");
        }

        UPTR[i] = AA.len();
        AA.push(1f64 / A.AA()[j]);
        JA.push(jrow);
        IA[i+1] = AA.len();
    }

    AA.shrink_to_fit();
    JA.shrink_to_fit();

    let mut M = Matrix::from(AA, JA, IA);
    M.set_dia_ptr(UPTR);

    M
    // Matrix::new()
}

//-----------------------------------------------------------------------------------------------------------//
pub fn LU_solve(M: &Matrix, v: &Vector) -> Vector {
    // LU solver
    // L: unit lower matrix
    // U: upper matrix

    let m = M.num_rows();
    let UPTR = match M.UPTR() {
        Some(uptr) => {
            uptr
        },
        None => {
            panic!("can not find diagonal pointer");
        }
    };
    let mut x = Vector::from(vec![0f64; m]);

    // foward sweep
    for i in 0..m {
        x[i] = v[i];
        for j in M.IA()[i]..UPTR[i] {
            x[i] -= M.AA()[j] * x[M.JA()[j]];
        }
    }

    // backward sweep
    for i in (0..m).rev() {
        for j in UPTR[i]+1..M.IA()[i+1] {
            x[i] -= M.AA()[j] * x[M.JA()[j]];
        }
        x[i] *= M.AA()[UPTR[i]];
    }

    x
}