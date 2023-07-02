use std::{
    fmt,
    collections::HashSet,
    convert::From,
    ops::Mul,
    path::PathBuf,
};
use rayon::prelude::*;
use regex::Regex;
use crate::linear_algebra::vector::Vector;

#[derive(Debug, Clone)]
pub struct Matrix {
    m: usize,
    n: usize,
    AA: Vec<f64>,
    JA: Vec<usize>,
    IA: Vec<usize>,
    UPTR: Option<Vec<usize>>
}

impl Matrix {
    pub fn new() -> Matrix {
        Matrix { 
            m: 0, 
            n: 0, 
            AA: Vec::new(), 
            JA: Vec::new(), 
            IA: Vec::new(),
            UPTR: None
        }
    }

//-----------------------------------------------------------------------------------------------------------//
    fn degree(&self, i: usize) -> usize {
        // degree of a node i
        let j1 = self.IA[i];
        let j2 = self.IA[i+1];
        let rows = &self.JA[j1..j2];

        if rows.contains(&i) {
            rows.len() - 1
        } else {
            rows.len()
        }
    }
    
//-----------------------------------------------------------------------------------------------------------//
    pub fn neighbor(&self, i: usize) -> Vec<usize> {
        // todo: return neighbor nodes
        (&self.JA[self.IA[i]..self.IA[i+1]]).to_vec()
    }

//-----------------------------------------------------------------------------------------------------------//
    pub fn bandwidth(&self) -> usize {
        // return bandwidth of a matrix, max|i - j|
        let m = self.m;

        let bw = (0..m).into_par_iter().map(|i| {
            let j1 = self.IA[i];
            let j2 = self.IA[i+1];

            self.JA[j1..j2].iter().map(|&j| i.abs_diff(j))
                .max_by(|x, y| x.cmp(y)).unwrap()

        }).max().expect("can not find a max value");

        bw
    }

//-----------------------------------------------------------------------------------------------------------//
    pub fn RCM(&self) -> Vec<usize> {
        // reverse Cuthill-McKee
        assert!(self.num_cols() == self.num_rows());

        let m = self.num_rows();
        let mut next: usize;
        let mut levset = HashSet::new();
        let mut marker = vec![0usize; m];
        let mut iperm = vec![usize::MAX; m];
        let mut iter = 0;

        //* find a initial node which has the minimum degree
        let i = (0..m).into_par_iter().min_by_key(|&i| self.degree(i))
            .expect("can not find a initial node");

        levset.insert(i);
        next = 1;
        marker[i] = 1;
        iperm[0] = i;

        while next < m  && iter < m{
            // println!("iter: {iter}, next: {next}, m: {m}");

            let mut next_levset = HashSet::new(); 
            let mut nodes = levset.into_iter().collect::<Vec<usize>>();

            nodes.sort_by(|&i, &j| 
                self.degree(i).cmp(&self.degree(j)).then(i.cmp(&j))
            );
            
            for j in nodes {
                //* neighbors oredered by its degree
                let mut neighbor = self.neighbor(j);
                neighbor.sort_by(|&i, &j|
                    self.degree(i).cmp(&self.degree(j)).then(i.cmp(&j))
                );

                //* visiting neighbour i of a node j in degree increasing order
                for i in neighbor {
                    if marker[i] == 0 {
                        next_levset.insert(i);
                        marker[i] = 1;
                        iperm[next] = i;
                        next += 1;
                    }
                }
            }
            
            levset = next_levset;

            iter += 1;
        }

        iperm.reverse();
        iperm
    }

//-----------------------------------------------------------------------------------------------------------//
    pub fn permutate_par(self, perm: &Vec<usize>) -> Matrix {
        // permutate matrix with permutation vector
        if self.m != perm.len() {
            panic!("theh length of permuation arrays is not equal to the number of rows");
        };

        let m = self.m;
        let mut inverse = vec![usize::MAX; m];
        
        // set inverse array of perm
        for i in 0..m {
            inverse[perm[i]] = i;
        }

        let argsort = (0..m).into_par_iter()
            .map(|i| {
                let j1 = self.IA[perm[i]];
                let j2 = self.IA[perm[i]+1];
                let ja = &self.JA[j1..j2];

                // sort cols by inverse
                let mut argsort = (0..j2-j1).collect::<Vec<_>>();
                argsort.sort_by_key(|&j| inverse[ja[j]]);

                argsort
            }).collect::<Vec<Vec<_>>>();

        let AA = argsort.par_iter().enumerate()
            .map(|(i, indices)| {
                let j1 = self.IA[perm[i]];
                let j2 = self.IA[perm[i]+1];
                let aa = &self.AA[j1..j2];
                let aa = indices.iter().map(|&j| aa[j]).collect::<Vec<_>>();

                aa
            }).flatten().collect::<Vec<_>>();
        
        let JA = argsort.par_iter().enumerate()
            .map(|(i, indices)| {
                let j1 = self.IA[perm[i]];
                let j2 = self.IA[perm[i]+1];
                let ja = &self.JA[j1..j2];
                let ja = indices.iter().map(|&j| inverse[ja[j]]).collect::<Vec<_>>();

                ja
            }).flatten().collect::<Vec<_>>();
        
        // todo: to be improved by parallelism
        let mut IA = Vec::with_capacity(m + 1);
        IA.push(0);
        for i in 0..m {
            IA.push(IA[i] + argsort[i].len())
        }

        Matrix::from(AA, JA, IA)
    }

//-----------------------------------------------------------------------------------------------------------//
    // pub fn permutate(self, perm: &Vec<usize>) -> Matrix {
    //     // permutate matrix with permutation vector
    //     // assert!(self.m == perm.len());

    //     let m = self.m;
    //     let mut inverse = vec![usize::MAX; m];
    //     let mut AA = Vec::with_capacity(self.AA.len());
    //     let mut JA = Vec::with_capacity(self.JA.len());
    //     let mut IA = Vec::with_capacity(self.IA.len());

    //     // set inverse array of perm
    //     for i in 0..m {
    //         inverse[perm[i]] = i;
    //     }
        
    //     IA.push(0);
    //     for i in 0..m {
    //         let j1 = self.IA[perm[i]];
    //         let j2 = self.IA[perm[i]+1];
    //         let aa = &self.AA[j1..j2];
    //         let ja = &self.JA[j1..j2];

    //         // sort cols by inverse
    //         let mut argsort = (0..j2-j1).collect::<Vec<_>>();
    //         argsort.sort_by_key(|&j| inverse[ja[j]]);

    //         // argsort();
    //         argsort.iter().for_each(|&j| {
    //             AA.push(aa[j]);
    //             JA.push(inverse[ja[j]]);
    //         });

    //         IA.push(AA.len());
    //         // println!("{:?}", argsort);
    //     }

    //     Matrix::from(AA, JA, IA)
    // }

//-----------------------------------------------------------------------------------------------------------//
    pub fn import_mtx(path: &str) -> Matrix {
        let path = PathBuf::from(path);
        
        if path.exists() {
            println!("import matrix from {}", path.display());
            let mut m = usize::MAX; 
            let mut n = usize::MAX; 
            let mut z = usize::MAX;
            let mut symmetry = false;
            let mut re;
            let text = std::fs::read_to_string(path)
                .expect("can not read file");
            let mut data = vec![];
            
            // header
            re = Regex::new(r"%%MatrixMarket matrix coordinate\s(\w+)\s(\w+)").unwrap(); 
            for cap in re.captures_iter(&text) {
                symmetry = cap[2].to_string().eq("symmetric");
                // println!("{}, {}", &cap[1], &cap[2]);
            }

            // parse dimension
            re = Regex::new(r"(\d+)\s(\d+)\s(\d+)\n").unwrap(); 
            for cap in re.captures_iter(&text) {
                m = cap[1].trim().parse::<usize>().unwrap();
                n = cap[2].trim().parse::<usize>().unwrap();
                z = cap[3].trim().parse::<usize>().unwrap();
            }
            if m != n {
                println!("caution: not square a matrix");
            }

            // parse each row to the tuple of (i, j, value)
            re = Regex::new(r"(\d+)\s(\d+)\s([-|\s].+)").unwrap(); 
            for cap in re.captures_iter(&text) {
                data.push((
                    cap[1].trim().parse::<usize>().unwrap(),
                    cap[2].trim().parse::<usize>().unwrap(),
                    cap[3].trim().parse::<f64>().unwrap(),
                ));
            }

            //* if symmetric matrix */
            if symmetry {
                let upper = data.par_iter()
                    .filter(|(i, j, _)| i != j)
                    .map(|(i, j, value)| {
                        (*j, *i, *value)
                    }).collect::<Vec<_>>();
                data.par_extend(upper.par_iter());
            }

            // sort by row
            data.sort_by(|a, b| {
                a.0.cmp(&b.0).then(a.1.cmp(&b.1))
            });

            // assemble matrix
            let mut AA = Vec::with_capacity(z);
            let mut JA = Vec::with_capacity(z);
            let mut IA = Vec::with_capacity(m);
            let mut row = usize::MAX;

            for (i, j, value) in &data {
                if i - 1 != row {
                    row = i - 1;
                    IA.push(AA.len());
                }
                AA.push(*value);
                JA.push(j - 1);
            }
            IA.push(data.len());

            println!("done.");
            Matrix::from(AA, JA, IA)
        } else {            
            println!("can not find file {}", path.display());
            Matrix::new()
        }
    }

//-----------------------------------------------------------------------------------------------------------//
    pub fn from<
        T: IntoParallelIterator, 
        U: IntoParallelIterator, 
        V: IntoParallelIterator
    >(AA: T, JA: U, IA: V) 
    -> Self where T::Item: Into<f64>, U::Item: Into<usize>, V::Item: Into<usize> { 
        let AA = AA.into_par_iter()
            .map(|v| v.into())
            .collect::<Vec<f64>>();
        let JA = JA.into_par_iter()
            .map(|v| v.into())
            .collect::<Vec<usize>>();
        let IA = IA.into_par_iter()
            .map(|v| v.into())
            .collect::<Vec<usize>>();

        let m = IA.len() - 1;
        let n = *JA.par_iter().max().unwrap() + 1;

        Matrix {
            m,
            n,
            AA,
            JA,
            IA,
            UPTR: None
        }
    }

//-----------------------------------------------------------------------------------------------------------//
    pub fn num_rows(&self) -> usize {
        self.m
    }

    pub fn num_cols(&self) -> usize {
        self.n
    }

    pub fn AA(&self) -> &Vec<f64> {
        &self.AA
    }

    pub fn JA(&self) -> &Vec<usize> {
        &self.JA
    }
    
    pub fn IA(&self) -> &Vec<usize> {
        &self.IA
    }

    pub fn UPTR(&self) -> Option<&Vec<usize>> {
        self.UPTR.as_ref()
    }

    pub fn set_dia_ptr(&mut self, UPTR: Vec<usize>) {
        self.UPTR = Some(UPTR);
    }
}

impl Mul<&Vector> for &Matrix {
    type Output = Vector;

    fn mul(self, rhs: &Vector) -> Self::Output {
        assert!(self.m == rhs.num_rows());
        
        let m = self.m;
        let mut AA = vec![0f64; m];

        AA.par_iter_mut().enumerate()
            .for_each(|(i, v)| {
                for j in self.IA[i]..self.IA[i+1] {
                    *v += self.AA[j] * rhs[self.JA[j]];
                }
            });

        Vector::from(AA)
    }
}

impl<I: IntoIterator> From<I> for Matrix 
    where I::Item: IntoIterator, <<I as IntoIterator>::Item as IntoIterator>::Item: Into<f64> {
    fn from(value: I) -> Self {

        for v in value {
            for v2 in v {
                let tmp: f64 = v2.into();
                print!("{} ", tmp);
            }
            println!();
        }

        Matrix::new()
    }
}

impl fmt::Display for Matrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let m = self.m;
        let n = self.n;

        writeln!(f, "Matrix")?;
        writeln!(f, "dim: [{} x {}]", m, n)?;

        let width = f.width().unwrap_or_else(|| 0);
        let precision = f.precision().unwrap_or_else(|| 0);

        for i in 0..m {
            let mut row = vec![0.0f64; n];
            for j in self.IA[i]..self.IA[i+1] {
                row[self.JA[j]] = self.AA[j];
            }

            for val in row {
                write!(f, "{val:>width$.pre$} ", width = width, pre = precision)?;
            }

            writeln!(f)?;
        }

        Ok(())
    } 
}

// #[cfg(test)]
// mod tests {
//     #[test]
//     fn RCM_test() {
//         assert!(false);
//     }
// }