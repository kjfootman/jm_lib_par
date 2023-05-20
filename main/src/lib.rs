use std::sync::Mutex;

use jm_math::prelude::*;
use devtimer::run_benchmark;
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

    let m = 300;
    let (A, b) = tri_diagonal(m);
    // println!("{:.2}", A);
    // println!("{:.2}", b);
    let bench_result = run_benchmark(10, |_| {
        let x = msolver::HGMRES(1000, 1.0E-13, 5, &A, &b);
        println!("{:.2}", x.AA().par_iter().sum::<f64>());
    });
    let time0 = bench_result.get_average() as f64 * 1.0E-9;

    let bench_result = run_benchmark(10, |_| {
        let x = msolver::Conjugate_gradient(1000, 1.0E-13, &A, &b);
        println!("{:.2}", x.AA().par_iter().sum::<f64>());
    });
    let time1 = bench_result.get_average() as f64 * 1.0E-9;


    println!();
    println!("GMRES: {:>10.4} sec", time0);
    println!("CG: {:>10.4} sec", time1);
}

//-----------------------------------------------------------------------------------------------------------//
#[allow(non_snake_case)]
fn tri_diagonal(m: usize) -> (Matrix, Vector) {
    // let m = 10_000_000;
    // let m = 10;
    let a = 2.3;
    let mut AA = vec![a, 1.0];
    let mut JA = vec![0usize, 1];
    let mut IA = vec![0];
    let mut v = vec![5f64; m - 2];

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