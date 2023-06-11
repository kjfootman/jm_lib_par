mod linear_algebra;

pub mod prelude {
    pub use crate::linear_algebra::vector::Vector;
    pub use crate::linear_algebra::matrix::Matrix;
    pub use crate::linear_algebra::msolver;
    pub use crate::linear_algebra::preconditioner::{self, Preconditioner};
}
