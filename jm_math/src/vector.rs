use std::fmt;
use std::convert::From;
use std::ops::{Add, AddAssign, Sub, SubAssign, Mul, Div};
use std::ops::{Deref, DerefMut};
use rayon::prelude::*;

pub struct Vector {
    m: usize,
    AA: Vec<f64>
}

/***********************************************************************************************************/
impl Vector {
    pub fn new() -> Self {
        Self { 
            m: 0, 
            AA: Vec::new() 
        }
    }

    pub fn from_iter<I: IntoIterator>(v: I) -> Vector where I::Item: Into<f64> {
        let AA = v.into_iter().map(|v| {
            v.into()
        }).collect::<Vec<f64>>();

        Vector::from(AA)
    }

//-----------------------------------------------------------------------------------------------------------//
    pub fn dim(&self) -> usize {
        self.m
    }

    pub fn AA(&self) -> &Vec<f64> {
        &self.AA
    }
}

/***********************************************************************************************************/
impl Add<&Vector> for &Vector {
    type Output = Vector;

    fn add(self, rhs: &Vector) -> Self::Output {
        let AA = (self.AA(), rhs.AA()).into_par_iter()
            .map(|(v1, v2)| {
                v1 + v2
            }).collect::<Vec<f64>>();

        Vector::from(AA)
    }
}

/***********************************************************************************************************/
impl AddAssign<&Vector> for Vector {
    fn add_assign(&mut self, rhs: &Vector) {
        let AA = (self.AA(), rhs.AA()).into_par_iter()
            .map(|(v1, v2)| {
                v1 + v2
            }).collect::<Vec<f64>>();

        *self = Vector::from(AA);
    }
}

/***********************************************************************************************************/
impl Sub<&Vector> for &Vector {
    type Output = Vector;

    fn sub(self, rhs: &Vector) -> Self::Output {
        let AA = (self.AA(), rhs.AA()).into_par_iter()
            .map(|(v1, v2)| {
                v1 - v2
            }).collect::<Vec<f64>>();

        Vector::from(AA)
    }    
}

/***********************************************************************************************************/
impl SubAssign<&Vector> for Vector {
    fn sub_assign(&mut self, rhs: &Vector) {
        let AA = (self.AA(), rhs.AA()).into_par_iter()
            .map(|(v1, v2)| {
                v1 - v2
            }).collect::<Vec<f64>>();

        *self = Vector::from(AA);
    }
}

/***********************************************************************************************************/
impl Mul<&Vector> for &Vector {
    type Output = f64;

    fn mul(self, rhs: &Vector) -> Self::Output {
        let sum = (self.AA(), rhs.AA()).into_par_iter()
            .map(|(v1, v2)| {
                v1 * v2
            }).sum::<f64>();
        
        sum
    }
}

/***********************************************************************************************************/
impl Mul<&Vector> for f64 {
    type Output = Vector;

    fn mul(self, rhs: &Vector) -> Self::Output {
        let AA = rhs.AA().par_iter()
            .map(|v| {
                self * v
            }).collect::<Vec<f64>>();
        
        Vector::from(AA)
    }
}

/***********************************************************************************************************/
impl Div<f64> for Vector {
    type Output = Vector;

    fn div(self, rhs: f64) -> Self::Output {
        let AA = self.AA().par_iter()
            .map(|v| {
                v / rhs
            }).collect::<Vec<f64>>();

        Vector::from(AA)
    }
}

/***********************************************************************************************************/
impl Deref for Vector {
    type Target = Vec<f64>;

    fn deref(&self) -> &Self::Target {
        &self.AA
    }
}

/***********************************************************************************************************/
impl DerefMut for Vector {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.AA
    }
}

/***********************************************************************************************************/
impl<I: IntoIterator> From<I> for Vector where I::Item: Into<f64> {
    fn from(value: I) -> Self {
        // let AA: Vec<f64> = value.into();
        // let m = AA.len();
        let AA = value.into_iter()
            .map(|v| v.into())
            .collect::<Vec<f64>>();
        let m = AA.len();
        Self { 
            m, 
            AA 
        }
    }
}

/***********************************************************************************************************/
impl fmt::Display for Vector {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let m = self.m;

        writeln!(f, "Vector")?;
        writeln!(f, "dim: [{} x 1]", m)?;

        let width = f.width().unwrap_or_else(|| 0);
        let precision = f.precision().unwrap_or_else(|| 2);

        for (i, val) in self.AA.iter().enumerate() {
            writeln!(f, "v[{i:}] = {val:>width$.pre$}", width = width, pre = precision)?
        }

        Ok(())
    }
}