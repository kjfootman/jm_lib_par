use crate::vector::Vector;
use crate::matrix::Matrix;

pub fn GMRES(iMax: usize, tol: f64, restart: usize, A: &Matrix, b: &Vector) -> Vector {
    let m = b.num_rows();
    let bl = b.l2_norm();
    let mut iter = 0;
    let mut residual = f64::MAX;
    let mut x = Vector::from(vec![0.0; m]);

    while iter < iMax && residual > tol {
        let mut V: Vec<Vector> = Vec::with_capacity(restart);
        let mut H = Vec::with_capacity(restart);
        let mut g = Vec::from(vec![0.0; restart + 1]);
        let r = b - &(A * &x);

        g[0] = r.l2_norm();
        V.push(&r / g[0]);

        // * Arnoli's method
        for j in 0..restart {
            let mut v = A * &V[j];
            let mut h = vec![0.0; j + 2];

            for i in 0..=j {
                h[i] = &v * &V[i];
                v -= &(h[i] * &V[i]);
            }

            h[j+1] = v.l2_norm();
            H.push(h);

            // degree = j + 1;

            if H[j][j+1].abs() < tol {
                println!("lucky breakdown");
                // degree = j;
                break;
            }

            v = &v / H[j][j+1];
            V.push(v);
        }

        // * Given's rotation
        Givens_rotation(&mut H, &mut g, &tol);

        // * Upper triangular matrix solve
        let y = upper_triangular_solve(&H, &g);

        for i in 0..H.len() {
            x += &(y[i] * &V[i]);
        }

        iter += 1;
        residual = g[H.len()].abs() / bl;
    }

    let mut print = format!(
        "MSolver: GMRES({}) {} iteration: {:5}, residual: {:.2E}",
        restart,
        "-".repeat(10),
        iter,
        residual
    );

    if iter == iMax {
        print.push_str(", ***** maximum iteration exceeded!");
    }
    println!("{print}");
    
    // println!("GMRES({})-------iter:{}", restart, iter);

    x
}

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
        let l = (r.powi(2) + h.powi(2)).sqrt();
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