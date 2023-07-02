// use core::num;
// use rayon::ThreadPool;
use rayon::prelude::*;
use std::time::Duration;
use std::thread;
use std::sync::Mutex;

fn main() {
    // test1();
    // test2();
    // test3();
    test4();
}

fn test4() {
    // par_iter() collect
    let v = (0..50).into_par_iter().map(|i| i).collect::<Vec<_>>();
    for (i, v) in v.iter().enumerate() {
        println!("v[{}] = {}", i, v);
    }
}

fn test3() {
    // rayon threadpool 
    let pool = rayon::ThreadPoolBuilder::new().num_threads(2).build().unwrap();

    pool.install(|| rayon::scope(|s| {
        s.spawn(|_| {
            let sec = 5;
            let id = pool.current_thread_index().unwrap();
            println!("thread{} started: sleep {sec}", id);
            thread::sleep(Duration::from_secs(sec));
            println!("thread{} ended", id);
        });

        s.spawn(|_| {
            let sec = 1;
            let id = pool.current_thread_index().unwrap();
            println!("thread{} started: sleep {sec}", id);
            thread::sleep(Duration::from_secs(sec));
            println!("thread{} ended", id);
        });

        s.spawn(|_| {
            let sec = 2;
            let id = pool.current_thread_index().unwrap();
            println!("thread{} started: sleep {sec}", id);
            thread::sleep(Duration::from_secs(sec));
            println!("thread{} ended", id);
        });
    }));
}


fn test2() {
    // remainder
    let numerator = 55;
    let denominator = 14; 
    let quotient = numerator / denominator;
    let remainder = numerator % denominator;
    let mut v = vec![usize::MAX; denominator];

    for i in 0..remainder {
        v[i] = quotient + 1;
    }

    for i in remainder..denominator {
        v[i] = quotient;
    }

    println!("numerator: {}, denominator: {}", numerator, denominator);
    println!("quotient: {}, remainder: {}", quotient, remainder);
    println!("{:?}", v);
}

fn test1() {
    // rayon: data sharing test
    let perm = [5, 3, 4, 0, 1, 2];
    let inverse = Mutex::new(vec![usize::MAX; perm.len()]);
    let mid = perm.len() / 2;

    rayon::scope(|s| {
        s.spawn(|_| {
            // some heaby process
            thread::sleep(Duration::from_secs(10));

            let mut inverse = inverse.lock().unwrap();
            println!("thread1 update inverse");
            for i in 0..mid {
                inverse[perm[i]] = i;
            }
            println!("{:?}", inverse);
        });

        s.spawn(|_| {
            // some heaby process
            thread::sleep(Duration::from_secs(5));

            let mut inverse = inverse.lock().unwrap();
            println!("thread2 update inverse");
            for i in mid..inverse.len() {
                inverse[perm[i]] = i;
            }
            println!("{:?}", inverse);
        });
    });

    let tmp = inverse.into_inner().unwrap();
    println!("main thread: {:?}", tmp);
}