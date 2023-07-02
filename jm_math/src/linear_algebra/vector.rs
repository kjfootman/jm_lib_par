use std::{
    fmt,
    convert::From
};
use std::ops::{
    Add, 
    AddAssign, 
    Sub, 
    SubAssign, 
    Mul, 
    Div, 
    Neg,
    Deref,
    DerefMut,
};
use rayon::prelude::*;

#[derive(PartialEq, Debug, Clone)]
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

//-----------------------------------------------------------------------------------------------------------//
    pub fn num_rows(&self) -> usize {
        self.m
    }

    pub fn AA(&self) -> &Vec<f64> {
        &self.AA
    }

    pub fn l2_norm(&self) -> f64 {
        self.AA.par_iter()
            .map(|v| v * v)
            .sum::<f64>().sqrt()
    }

//-----------------------------------------------------------------------------------------------------------//
    pub fn permutate(self, perm: &Vec<usize>) -> Vector {
        let AA = perm.par_iter().map(|&i| self.AA[i]).collect::<Vec<_>>();

        Vector::from(AA)
    }
}

/***********************************************************************************************************/
impl Add<&Vector> for &Vector {
    type Output = Vector;

    fn add(self, rhs: &Vector) -> Self::Output {
        let AA = (self.AA(), rhs.AA()).into_par_iter()
            .map(|(v1, v2)| v1 + v2)
            .collect::<Vec<f64>>();

        Vector::from(AA)
    }
}

/***********************************************************************************************************/
impl AddAssign<&Vector> for Vector {
    fn add_assign(&mut self, rhs: &Vector) {
        let AA = (self.AA(), rhs.AA()).into_par_iter()
            .map(|(v1, v2)| v1 + v2)
            .collect::<Vec<f64>>();

        *self = Vector::from(AA);
    }
}

/***********************************************************************************************************/
impl Sub<&Vector> for &Vector {
    type Output = Vector;

    fn sub(self, rhs: &Vector) -> Self::Output {
        let AA = (self.AA(), rhs.AA()).into_par_iter()
            .map(|(v1, v2)| v1 - v2)
            .collect::<Vec<f64>>();

        Vector::from(AA)
    }    
}

/***********************************************************************************************************/
impl SubAssign<&Vector> for Vector {
    fn sub_assign(&mut self, rhs: &Vector) {
        let AA = (self.AA(), rhs.AA()).into_par_iter()
            .map(|(v1, v2)| v1 - v2 )
            .collect::<Vec<f64>>();

        *self = Vector::from(AA);
    }
}

/***********************************************************************************************************/
impl Mul<&Vector> for &Vector {
    type Output = f64;

    fn mul(self, rhs: &Vector) -> Self::Output {
        let sum = (self.AA(), rhs.AA()).into_par_iter()
            .map(|(v1, v2)| v1 * v2)
            .sum::<f64>();
        
        sum
    }
}

/***********************************************************************************************************/
impl Mul<&Vector> for f64 {
    type Output = Vector;

    fn mul(self, rhs: &Vector) -> Self::Output {
        let AA = rhs.AA().par_iter()
            .map(|v| self * v)
            .collect::<Vec<f64>>();
        
        Vector::from(AA)
    }
}

/***********************************************************************************************************/
impl Div<f64> for &Vector {
    type Output = Vector;

    fn div(self, rhs: f64) -> Self::Output {
        let AA = self.AA().par_iter()
            .map(|v| v / rhs)
            .collect::<Vec<f64>>();

        Vector::from(AA)
    }
}

/***********************************************************************************************************/
impl Neg for Vector {
    type Output = Vector;

    fn neg(self) -> Self::Output {
        let AA = self.AA.par_iter()
            .map(|v| -v)
            .collect::<Vec<f64>>();

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

// #[cfg(test)]
// mod tests {
//     use super::*;
    
//     #[test]
//     fn permuate_() {
//         let perm = vec![1usize, 0, 2];    
//         let v = Vector::from(vec![0f64, 1f64, 2f64]);
//         let v = v.permutate(&perm);
//         assert_eq!(Vector::from(vec![1f64, 0f64, 2f64]), v);
//     }
// }