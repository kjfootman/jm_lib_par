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

    let v = Vector::from([true, false]);
    println!("{}", v);
}