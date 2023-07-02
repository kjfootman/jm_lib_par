use std::ops::{
    Add,
    AddAssign,
    Sub,
    SubAssign,
    Deref,
    DerefMut
};

#[derive(PartialEq, Debug)]
pub struct Point {
    x: f64,
    y: f64,
    z: f64,
}

/***********************************************************************************************************/
impl Point {
    pub fn new(x: f64, y: f64, z: f64) -> Point {
        Point { 
            x,
            y,
            z,
        }
    }
}