#[allow(non_snake_case)]
mod mesh;

pub mod prelude {
    pub use crate::mesh::point::Point;
}

pub fn add(left: usize, right: usize) -> usize {
    left + right
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::mesh::point::Point;

//     #[test]
//     fn point_define() {
//         let p = Point::new(1.0, 2.0, 3.0);
//         assert_eq!(3, p);
//     }
// }
