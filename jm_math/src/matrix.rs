use std::fmt;
use crate::vector::Vector;
use std::ops::{Mul, MulAssign};

pub struct Matrix {
    m: usize,
    n: usize,
    AA: Vec<f64>,
    JA: Vec<usize>,
    IA: Vec<usize>
}

/***********************************************************************************************************/
impl Matrix {
    pub fn new() {

    }
}