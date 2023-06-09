use rayon::prelude::*;
use crate::linear_algebra::{
    vector::Vector,
    matrix::Matrix,
    preconditioner as precon,
    preconditioner::Preconditioner
};

pub fn GMRES(iMax: usize, tol: f64, restart: usize, A: &Matrix, b: &Vector, preconditioner: Preconditioner) -> Vector {
    assert!(A.num_cols() == b.num_rows());

    let m = b.num_rows();
    let bl = b.l2_norm();
    let mut iter = 0;
    let mut residual = f64::MAX;
    let mut x = Vector::from(vec![0.0; m]);
    let P = preconditioner.from(A);

    while iter < iMax && residual > tol {
        assert!(A.num_cols() == b.num_rows());

        let mut V: Vec<Vector> = Vec::with_capacity(restart);
        let mut H = Vec::with_capacity(restart);
        let mut g = Vec::from(vec![0.0; restart + 1]);
        let r = b - &(A * &x);

        g[0] = r.l2_norm();
        V.push(&r / g[0]);
        
        // * Arnoldi's process - Modified Grame-Schmidt
        for j in 0..restart {
            let mut h = vec![0.0; j + 2];

            //* right preconditioning w = A * M_inv * w
            let mut w = match &P {
                Some(M) => A * &precon::LU_solve(M, &V[j]),
                None => A * &V[j]
            };
            
            for i in 0..=j {
                h[i] = &w * &V[i];
                w -= &(h[i] * &V[i]);
            }

            h[j+1] = w.l2_norm();
            H.push(h);

            if H[j][j+1].abs() < tol {
                println!("lucky breakdown");
                break;
            }

            w = &w / H[j][j+1];
            V.push(w);
        }

        // * Given's rotation
        Givens_rotation(&mut H, &mut g, &tol);

        // * Upper triangular matrix solve
        let y = upper_triangular_solve(&H, &g);

        let mut z = Vector::from(vec![0.0; m]);
        for i in 0..H.len() {
            z += &(y[i] * &V[i]);
        }
        
        //* right preconditioning x = x + M_inv * z
        x = match &P {
            Some(M) => &x + &precon::LU_solve(M, &z),
            None => &x + &z
        };

        iter += 1;
        residual = g[H.len()].abs() / bl;
    }

    let log = Log {
        solver: "GMRES",
        precon: preconditioner,
        restart: Some(restart),
        iMax,
        iter,
        residual
    };

    display(log);
    
    x
}

//-----------------------------------------------------------------------------------------------------------//
pub fn HGMRES(iMax: usize, tol: f64, restart: usize, A: &Matrix, b: &Vector, preconditioner: Preconditioner) -> Vector {
    assert!(A.num_cols() == b.num_rows());

    let m = b.num_rows();
    let bl = b.l2_norm();
    let mut iter = 0;
    let mut residual = f64::MAX;
    let mut x = Vector::from(vec![0.0; m]); 
    let P = preconditioner.from(A);

    while iter < iMax && residual > tol {
        // let mut degree = 0;
        let mut W: Vec<Vector> = Vec::with_capacity(restart);
        let mut H = Vec::with_capacity(restart);
        let mut g = Vec::from(vec![0.0; restart + 1]);
        let mut z = b - &(A * &x);

        for j in 0..restart + 1 {
            //* calculate Householder vector
            // match Householder_vec(j, &z) {
            //     Some(w) => {
            //         // degree = j;
            //         W.push(w);
            //     },
            //     None => {
            //         println!("householder tmp");
            //         let mut tmp = z.par_iter().map(|&v| v).collect::<Vec<f64>>();
            //         tmp.push(0.0);
            //         H.push(tmp);
            //         break;
            //     } 
            // }
            W.push(Householder_vec(j, &z).unwrap());

            //* calculate h[j-1] = P[j]z[j] 
            if j != 0 {
                let mut h = z[0..=j].iter()
                    .map(|&v| v)
                    .collect::<Vec<f64>>();
                h[j] = z[j];
                
                let sigma = (&z.AA()[j..m], &W[j][0..m-j]).into_par_iter()
                    .map(|(v1, v2)| v1 * v2).sum::<f64>();
                h[j] = z[j] - 2.0 * sigma * W[j][0];
                
                H.push(h);

                if H[j-1][j].abs() < tol {
                    println!("lucky breakdown");
                    break;
                }
            } else {
                g[0] = -z[0].signum() * z.l2_norm();
            }

            //* calculate basis vector v[j] = P(0) P(1) .. P(j) e(j)
            let mut v = Vector::from(vec![0.0; m]);
            v[j] = 1.0;

            // test
            for n in (0..=j).rev() {
                let sigma = (&v.AA()[n..m], &W[n][0..m-n]).into_par_iter()
                    .map(|(v1, v2)| v1 * v2).sum::<f64>();

                v[n..m].par_iter_mut().enumerate()
                    .for_each(|(i, v)|{
                        // println!("n: {}, i: {}", n, i);
                        *v -= 2.0 * sigma * W[n][i]
                    });
            }

            //* calculate z = P(j) .. P(1) P(0) A v(j)
            //* right preconditioning z = A * M_inv * v(j)
            z = match &P {
                Some(M) => A * &precon::LU_solve(M, &v),
                None => A * &v
            };
            // z = A * &v;

            for n in 0..=j {
                let sigma = (&z.AA()[n..m], &W[n][0..m-n]).into_par_iter()
                    .map(|(v1, v2)| v1 * v2).sum::<f64>();

                z[n..m].par_iter_mut().enumerate()
                    .for_each(|(i, v)|{
                        // println!("n: {}, i: {}", n, i);
                        *v -= 2.0 * sigma * W[n][i]
                    });
            }
        }
        
        // * Girven's rotation
        Givens_rotation(&mut H, &mut g, &tol);

        // * upper triangular matrix solve
        let y = upper_triangular_solve(&H, &g);

        // * update solution z = P(j)(y(j)e(j) - z)
        z = Vector::from(vec![0.0; m]);
        for n in (0..H.len()).rev() {
            z[n] += y[n];

            let sigma = (&z.AA()[n..m], &W[n][0..m-n]).into_par_iter()
                .map(|(v1, v2)| v1 * v2).sum::<f64>();

            z[n..m].par_iter_mut().enumerate()
                .for_each(|(i, v)|{
                    *v -= 2.0 * sigma * W[n][i]
                });
        }

        //* right preconditioning x = x + M_inv * z
        x = match &P {
            Some(M) => &x + &precon::LU_solve(M, &z),
            None => &x + &z
        };
        // x += &z;

        // println!("iteration: {}, residual: {:.4E}", self.iter, g[H.len()].abs() / bl);
        iter += 1;
        residual = g[H.len()].abs() / bl;
    } 

    let log = Log {
        solver: "HGMRES",
        precon: preconditioner,
        restart: Some(restart),
        iMax,
        iter,
        residual
    };

    display(log);

    // let mut print = format!(
    //     "MSolver: HGMRES({}) {} iteration: {:5}, residual: {:.4E}",
    //     restart,
    //     "-".repeat(10),
    //     iter,
    //     residual
    // );

    // if iter == iMax {
    //     print.push_str(", ***** maximum iteration exceeded!");
    // }

    // println!("{print}");

    x
}

fn Householder_vec(i: usize, v: &Vector) -> Option<Vector> {
    let n = v.len();
    let AA = v[i..n].par_iter().map(|&v| v).collect::<Vec<f64>>();
    let mut w = Vector::from(AA);

    w[0] = v.get(i)? + v.get(i)?.signum() * w.l2_norm(); 
    Some(&w / w.l2_norm())
}

//-----------------------------------------------------------------------------------------------------------//
fn Givens_rotation(H: &mut Vec<Vec<f64>>, g: &mut Vec<f64>, tol: &f64) {
    let v = H.last().expect("error in last colum of H")
                    .last().expect("error in las element of H[j]"); 

    let dim = if v.abs() > *tol {
        H.len()
    } else {
        H.len() - 1
    };

    for j in 0..dim {
        let hCol = &H[j];
        let r = hCol[j];
        let h = hCol[j+1];
        let l = (r.powf(2.0) + h.powf(2.0)).sqrt();
        let cos = r / l;
        let sin = -h / l;
        let g0 = g[j];
        // let g1 = g[j+1];
        
        for k in j..H.len() {
            let r = H[k][j];
            let h = H[k][j+1];
            H[k][j] = r * cos - h * sin;
            H[k][j+1] = r * sin + h * cos;
        }

        g[j] = g0 * cos /*- g1 * sin*/;
        g[j+1] = g0 * sin /*+ g1 * cos*/;
    }
}

//-----------------------------------------------------------------------------------------------------------//
fn upper_triangular_solve(H: &Vec<Vec<f64>>, g: &Vec<f64>) -> Vec<f64> {
    let mut y = vec![0.0; H.len()];

    for i in (0..H.len()).rev() {
        y[i] = g[i];
        for j in i+1..H.len() {
            y[i] -= H[j][i] * y[j];
        }
        y[i] /= H[i][i];
    }
    
    y
}

//-----------------------------------------------------------------------------------------------------------//
pub fn CG(iMax: usize, tol: f64, A: &Matrix, b:&Vector) -> Vector {
    // conjugate gradient solver
    assert!(A.num_cols() == b.num_rows());        
    
    let m = b.num_rows();
    let bl =  b.l2_norm();
    let mut iter = 0;
    let mut residual = f64::MAX;
    let mut x = Vector::from(vec![0.0; m]);
    let mut r = b - &(A * &x);
    let mut p = Vector::from(r.clone());
    let mut rsold = &r * &r;

    while iter < iMax && residual > tol {
        let Ap = A * &p;
        let alpha = rsold / (&p * &Ap);
        
        x += &(alpha * &p);
        r -= &(alpha * &Ap);
        // * rsnew = &r * &z;
        let rsnew = &r * &r;

        residual = rsnew.sqrt() / bl;
        // residual = (b - &(A * &x)).l2_norm() / bl;

        p = &r + &((rsnew / rsold) * &p);
        rsold = rsnew;
        
        iter += 1;
    }
    
    let mut print = format!(
        "MSolver: CG {} iteration: {:5}, residual: {:.4E}",
        "-".repeat(10),
        iter,
        residual
    );

    if iter == iMax {
        print.push_str(", ***** maximum iteration exceeded!");
    }
    println!("{print}");

    x
}

//-----------------------------------------------------------------------------------------------------------//
pub fn Gauss_Seidel(iMax: usize, tol: f64, A: &Matrix, b: &Vector) -> Vector {
    assert!(A.num_cols() == b.num_rows());
    
    let m = b.num_rows();
    let bl = b.l2_norm();
    let AA = A.AA();
    let JA = A.JA();
    let IA = A.IA();
    let mut iter = 0;
    let mut residual = f64::MAX;
    let mut x = Vector::from(vec![0.0; m]);

    while iter < iMax && residual > tol {
        for i in 0..m {
            let mut denominator = 0.0;
            x[i] = b[i];
            for j in IA[i]..IA[i+1]{
                if i != JA[j] {
                    x[i] -= AA[j] * x[JA[j]];
                } else {
                    denominator = AA[j];
                }
            }
            x[i] /= denominator;
        }

        residual = (b - &(A * &x)).l2_norm() / bl;
        iter += 1;
    }

    let mut print = format!(
        "MSolver: Gauss-Seidel {} iteration: {:5}, residual: {:.4E}",
        "-".repeat(10),
        iter,
        residual
    );

    if iter == iMax {
        print.push_str(", ***** maximum iteration exceeded!");
    }

    println!("{print}");

    x
}

//-----------------------------------------------------------------------------------------------------------//
struct Log {
    solver: &'static str,
    precon: Preconditioner,
    restart: Option<usize>,
    iMax: usize,
    iter: usize,
    residual: f64
}

// fn display(solver: &str, preconditioner: Preconditioner, restart: Option<usize>, iter: usize, residual: f64) -> String {
fn display(log: Log) {
    let solver = log.solver;
    let precondition = match log.precon {
        Preconditioner::GS => " with GS precondition",
        Preconditioner::SGS => " with SGS precondition",
        Preconditioner::ILU => " with ILU precondition",
        _ => ""
    };
    let restart = match log.restart {
        Some(restart) => {
            format!("({restart})")
        },
        None => String::new()
    };
    let iter = log.iter;
    let residual = log.residual;
    let mut format = format!(
        "MSolver: {solver}{restart}{precondition} {sep} iteration: {iter:5}  residual: {residual:.4E}",
        sep="-".repeat(10)
    );

    if iter == log.iMax {
        format += " ***** warning: maximum iteration exeeded!";
    }

    println!("{format}");
}