#[allow(non_snake_case)]
mod vector;
#[allow(non_snake_case)]
mod matrix;
#[allow(non_snake_case)]
pub mod msolver;

pub mod prelude {
    pub use crate::vector::Vector;
    pub use crate::matrix::Matrix;
    pub use crate::msolver;
}
