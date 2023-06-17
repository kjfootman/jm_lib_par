use std::sync::Mutex;

use jm_math::prelude::*;
use devtimer::{run_benchmark, DevTime};
use rayon::prelude::*;

const DIM: usize = 300_000_000;

pub fn test1() {
    let n = 10;
    let v1 = Vector::from(vec![1.0; DIM]);
    let v2 = Vector::from(vec![1.0; DIM]);
    
    //* Vector struct method */
    let bench_result = run_benchmark(n, |_| {
        let v = &v1 - &v2;
        v.par_iter().sum::<f64>();
    });
    let time0 = bench_result.get_average() as f64 * 1.0E-9;
    println!();

    //* direct operation method */
    let bench_result = run_benchmark(n, |_| {
        let v = (v1.AA(), v2.AA()).into_par_iter().map(|(v1, v2)| {
            v1 - v2
        }).collect::<Vec<f64>>();

        v.par_iter().sum::<f64>();
    });
    let time1 = bench_result.get_average() as f64 * 1.0E-9;

    println!();
    println!("Vector struct operation: {:>10.4} sec", time0);
    println!("direct operation: {:>10.4} sec", time1);
}

//-----------------------------------------------------------------------------------------------------------//
pub fn test2() {
    let n = 10;
    let v1 = Vector::from(vec![1.0; DIM]);
    let v2 = Vector::from(vec![1.0; DIM]);

    //* single thread operation */
    let bench_result = run_benchmark(n, |_| {
        let sum = v1.iter().zip(v2.iter()).map(|(v1, v2)| {
            v1 * v2
        }).sum::<f64>();
        
        println!("{}", sum);
    });
    let time0 = bench_result.get_average() as f64 * 1.0E-9;
    println!();

    //* multi thread operation */
    let bench_result = run_benchmark(n, |_| {
        let r = &v1 * &v2;
        println!("{}", r); 
    });
    let time1 = bench_result.get_average() as f64 * 1.0E-9;

    println!();
    println!("single thread operation: {:>10.4} sec", time0);
    println!("multi thread operation: {:>10.4} sec", time1);
}

//-----------------------------------------------------------------------------------------------------------//
#[allow(non_snake_case)]
pub fn test3() {
    let n = 10;
    let m = 10_000_000;
    // let m = 10;
    let mut AA = vec![3.0, 1.0];
    let mut JA = vec![0usize, 1];
    let mut IA = vec![0];

    for i in 1..m - 1 {
        IA.push(AA.len());

        AA.push(1.0);
        AA.push(3.0);
        AA.push(1.0);
        
        JA.push(i - 1);
        JA.push(i);
        JA.push(i + 1);
    }

    IA.push(AA.len());
    AA.append(&mut vec![1.0, 3.0]);
    JA.append(&mut vec![m - 2, m - 1]);
    IA.push(AA.len());

    let v = Mutex::new(Vector::from(vec![0.0; m]));

    //* single thread operation */
    let bench_result = run_benchmark(n, |_| {
        let mut v = v.lock().unwrap();

        for i in 0..m {
            for j in IA[i]..IA[i+1] {
                v[i] += AA[j] * v[JA[j]];
            }
        }

        v.AA().iter().sum::<f64>();
    });
    let time0 = bench_result.get_average() as f64 * 1.0E-9;
    
    let M = Matrix::from(AA, JA, IA);
    
    //* multi thread operation */
    let bench_result = run_benchmark(n, |_| {
        let v = v.lock().unwrap();
        let v = &M * &v;

        v.par_iter().sum::<f64>();
    });
    let time1 = bench_result.get_average() as f64 * 1.0E-9;

    println!();
    println!("single thread operation: {:>10.4} sec", time0);
    println!("multi thread operation: {:>10.4} sec", time1);

    // let v = Vector::from(vec![1f64; m]);
    // println!("{:.2}", &M * &v);
}

//-----------------------------------------------------------------------------------------------------------//
#[allow(non_snake_case)]
pub fn test4() {
    // let AA = [3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3];
    // let JA = [0usize, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4];
    // let IA = [0usize, 2, 5, 8, 11 , 13];

    // let v = Vector::from([4, 5, 5, 5, 4]); 
    // let A = Matrix::from(AA, JA, IA);

    // let x = msolver::GMRES(1000, 1E-13, 2, &A, &v);
    // let x = msolver::HGMRES(1000, 1E-13, 2, &A, &v);
    // println!("{}", x);

    let m = 3_000_000;
    let (A, b) = tri_diagonal(m);
    // println!("{:.2}", A);
    // println!("{:.2}", b);
    let bench_result = run_benchmark(10, |_| {
        let x = msolver::GMRES(1000, 1.0E-13, 5, &A, &b, Preconditioner::GS);
        println!("{:.2}", x.AA().par_iter().sum::<f64>());
    });
    let time0 = bench_result.get_average() as f64 * 1.0E-9;

    let bench_result = run_benchmark(10, |_| {
        let x = msolver::CG(1000, 1.0E-13, &A, &b);
        println!("{:.2}", x.AA().par_iter().sum::<f64>());
    });
    let time1 = bench_result.get_average() as f64 * 1.0E-9;


    println!();
    println!("GMRES: {:>10.4} sec", time0);
    println!("CG: {:>10.4} sec", time1);
}

//-----------------------------------------------------------------------------------------------------------//
#[allow(non_snake_case)]
pub fn test5() {
    let m = 7_000_000;
    let tol = 1.0E-10;
    let (A, b) = tri_diagonal(m);

    // let x = msolver::HGMRES(1000, 1.0E-13, 3, &A, &b);
    // println!("{}", x.par_iter().sum::<f64>());

    let bench_result = run_benchmark(1, |_| {
        let x = msolver::CG(1000, tol, &A, &b);
        println!("{:.6}", x.par_iter().sum::<f64>());
    });
    let time1 = bench_result.get_average() as f64 * 1.0E-9;

    let bench_result = run_benchmark(1, |_| {
        let x = msolver::Gauss_Seidel(1000, tol, &A, &b);
        println!("{:.6}", x.par_iter().sum::<f64>());
    });
    let time2 = bench_result.get_average() as f64 * 1.0E-9;

    let bench_result = run_benchmark(1, |_| {
        let x = msolver::GMRES(1000, tol, 5, &A, &b, Preconditioner::GS);
        println!("{:.6}", x.par_iter().sum::<f64>());
    });
    let time0 = bench_result.get_average() as f64 * 1.0E-9;

    let bench_result = run_benchmark(1, |_| {
        let x = msolver::HGMRES(1000, tol, 5, &A, &b, Preconditioner::GS);
        println!("{:.6}", x.par_iter().sum::<f64>());
    });
    let time3 = bench_result.get_average() as f64 * 1.0E-9;

    println!();
    println!("CG: {:>10.4} sec", time1);
    println!("GS: {:>10.4} sec", time2);
    println!("GMRES: {:>10.4} sec", time0);
    println!("HGMRES: {:>10.4} sec", time3);
}

//-----------------------------------------------------------------------------------------------------------//
#[allow(non_snake_case)]
pub fn test6() {
    let n = 1000;
    let m = 100_000_000;
    let v1 = vec![1.0; m];
    let v2 = vec![1.0; m];
    let mut v3 = vec![0.0; m];

    v3[n..m].par_iter_mut().enumerate()
        .for_each(|(i, v)| *v -= v1[i] * v2[i]);

    println!("sum: {}", v3.par_iter().sum::<f64>());

    let AA = vec![5, 1, 2, 1, 1, 1, 4];
    let JA = vec![0usize, 2, 1, 2, 0, 1, 2];
    let IA = vec![0usize, 2, 4, 7];
    let b = Vector::from(vec![8, 7, 15]);

    let A = Matrix::from(AA, JA, IA);
    let x = msolver::HGMRES(1000, 1.0E-13, 3, &A, &b, Preconditioner::GS);
    println!("{:.2}", x);
}

//-----------------------------------------------------------------------------------------------------------//
#[allow(non_snake_case)]
pub fn test7() {
    //* preconditioner verification
    let n = 1;
    let AA = vec![1f64; 29];
    let JA = vec![
        0usize, 0, 1, 1, 2, 2, 3, 0, 4, 1, 4, 5, 2, 5, 6,
        3, 6, 7, 4, 8, 5, 8, 9, 6, 9, 10, 7, 10, 11
    ];
    let IA = vec![0usize, 1, 3, 5, 7, 9, 12, 15, 18, 20, 23, 26, 29];
    let M = Matrix::from(AA, JA, IA);

    let v = (0..M.num_rows()).map(|i| M.AA()[M.IA()[i]..M.IA()[i+1]].iter().sum()).collect::<Vec<f64>>();
    let v = Vector::from(v);

    // let x = preconditioner::level_schduling(&M, &v);
    // println!("{:.2}", x);

    // let M = Matrix::import_mtx("./res/nos4.mtx");
    // let M = Matrix::import_mtx("res/bcsstm12.mtx");
    // let (M, v) = tri_diagonal(7);
    // let v = (0..M.num_rows()).map(|i| M.AA()[M.IA()[i]..M.IA()[i+1]].iter().sum()).collect::<Vec<f64>>();
    // let v = Vector::from(v);
    // println!("{:>13.2}", M);

    // let bench_result = run_benchmark(n, |_| {
    //     let x = msolver::GMRES(1000, 1.0E-13, 10, &M, &v, Preconditioner::GS);
    //     println!("{:.6}", x.iter().sum::<f64>());
    // });
    // let time0 = bench_result.get_average() as f64 * 1.0E-9;

    // let bench_result = run_benchmark(n, |_| {
    //     let x = msolver::HGMRES(1000, 1.0E-13, 10, &M, &v, Preconditioner::GS);
    //     println!("{:.6}", x.iter().sum::<f64>());
    // });
    // let time1 = bench_result.get_average() as f64 * 1.0E-9;

    // let bench_result = run_benchmark(n, |_| {
    //     let x = msolver::CG(1000, 1.0E-13, &M, &v);
    //     println!("{:.6}", x.iter().sum::<f64>());
    // });
    // let time2 = bench_result.get_average() as f64 * 1.0E-9;
    
    // println!();
    // println!("GMRES: {:>10.4} sec", time0);
    // println!("HGMRES: {:>10.4} sec", time1);
    // println!("CG: {:>10.4} sec", time2);


    let M = Matrix::import_mtx("res/bcsstk14.mtx");
    // let M = Matrix::import_mtx("res/nos4.mtx");
    // let M = Matrix::import_mtx("res/jpwh_991.mtx");
    // let M = Matrix::import_mtx("res/bcsstm12.mtx");
    // let M = Matrix::import_mtx("res/orsirr_1.mtx");
    let m = M.num_cols();
    let v = Vector::from(vec![1f64; m]);
    let v = &M * &v;

    let x = msolver::CG(1000, 1.0E-6, &M, &v);
    println!("{:.4}\n", x.iter().sum::<f64>());
    let x = msolver::GMRES(1000, 1.0E-6, 10, &M, &v, Preconditioner::SGS);
    println!("{:.4}\n", x.iter().sum::<f64>());
    let x = msolver::GMRES(1000, 1.0E-6, 10, &M, &v, Preconditioner::ILU);
    println!("{:.4}\n", x.iter().sum::<f64>());
    let x = msolver::HGMRES(1000, 1.0E-6, 10, &M, &v, Preconditioner::SGS);
    println!("{:.4}", x.iter().sum::<f64>());
}

//-----------------------------------------------------------------------------------------------------------//
#[allow(non_snake_case)]
pub fn test8() {
    let AA = vec![1, 2, 3, 4, 5, 6, 7, 8, 9];
    let JA = vec![0usize, 3, 1, 2, 3, 0, 2, 3, 3];
    let IA = vec![0usize, 2, 5, 8, 9];
    // let M = Matrix::import_mtx("./res/bcsstk01.mtx");
    let A = Matrix::from(AA, JA, IA);
    let v = Vector::from(vec![1, 3, 13, 9]);

    println!("{:.2}", A);

    let M = preconditioner::SGS(&A);
    println!("{:.2}", M);
    
    let M = preconditioner::ILU(&A);
    println!("{:.2}", M);
}

//-----------------------------------------------------------------------------------------------------------//
#[allow(non_snake_case)]
fn tri_diagonal(m: usize) -> (Matrix, Vector) {
    let a = 2.0;
    let mut AA = vec![a, 1.0];
    let mut JA = vec![0usize, 1];
    let mut IA = vec![0];
    let mut v = vec![a + 2.0; m - 2];

    for i in 1..m - 1 {
        IA.push(AA.len());

        AA.push(1.0);
        AA.push(a);
        AA.push(1.0);
        
        JA.push(i - 1);
        JA.push(i);
        JA.push(i + 1);
    }

    IA.push(AA.len());
    AA.append(&mut vec![1.0, a]);
    JA.append(&mut vec![m - 2, m - 1]);
    IA.push(AA.len());

    let M = Matrix::from(AA, JA, IA);

    v.insert(0, 1.0 + a);
    v.push(1.0 + a);
    let v = Vector::from(v);

    (M, v)
}