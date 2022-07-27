//all the methods of calculating averages and variances etc. 

use {
    rand_pcg::Pcg64,
    net_ensembles::sampling::bootstrap::bootstrap_copyable
};


//Credit:: Yannick

pub const BOOTSTRAP_SAMPLES: usize = 400;
pub struct MyVarianceBootstrap
{
    pub mean: f64,
    pub mean_err: f64,
    pub var: f64,
    pub var_err: f64
}

impl MyVarianceBootstrap{

    pub fn mean(&self) -> f64
    {
        self.mean
    }

    pub fn variance_of_mean(&self) -> f64
    {
        self.var
    }

    pub fn from_slice(
        slice: &[u32], 
        frac: Option<f64>, 
        bootstrap_rng1: &mut Pcg64,
        bootstrap_rng2: &mut Pcg64
    ) -> Self
    {
        let (mean, mean_err) = bootstrap_copyable(
            bootstrap_rng1,
            BOOTSTRAP_SAMPLES,
            slice,
            |slice| calc_average(slice, frac)
        );
        let (var, var_err) = bootstrap_copyable(
            bootstrap_rng2,
            BOOTSTRAP_SAMPLES,
            slice,
            |slice| calc_variance(slice, mean, frac)
        );
        Self{
            mean,
            var,
            mean_err,
            var_err
        }
    }
}
#[derive(Clone)]
pub struct MyVariance
{
    pub mean: f64,
    pub var: f64
}

impl MyVariance{

    pub fn mean(&self) -> f64
    {
        self.mean
    }

    pub fn variance_of_mean(&self) -> f64
    {
        self.var
    }

    pub fn from_slice(slice: &[u32], frac: Option<f64>) -> Self
    {
        let mean = calc_average(slice, frac);
        let var = calc_variance(slice, mean, frac);
        Self{
            mean,
            var
        }
    }
}


pub struct Measured
{
    pub var_m: MyVariance,
    pub var_c: MyVariance,
}

pub fn calc_average(slice: &[u32], frac: Option<f64>) -> f64
{
    let mut sum = 0_u64;
    for val in slice
    {
        sum += *val as u64;
    }

    let len = slice.len() as u64;
    let rest = sum % len;
    let div = sum / len;

    let res = div as f64 + (rest as f64) / (len as f64);
    match frac{
        None => res,
        Some(f) => res / f
    }
}

pub fn calc_variance(slice: &[u32], average: f64, frac: Option<f64>) -> f64
{
    let mut var_sum = 0.0;

    match frac{
        None => {
            for &val in slice{
                let dif = average - val as f64;
                var_sum += dif * dif;
            }
        },
        Some(v) => {
            for &val in slice{
                let dif = average - val as f64 / v;
                var_sum += dif * dif;
            }
        }
    }


    var_sum / slice.len() as f64
}

pub fn vector_average(slice:&Vec<u32>,frac:Option<f64>) -> f64{
    let mut sum = 0_u64;
    for val in slice
    {
        sum += *val as u64;
    }

    let len = slice.len() as u64;
    let rest = sum % len;
    let div = sum / len;

    let res = div as f64 + (rest as f64) / (len as f64);
    match frac{
        None => res,
        Some(f) => res / f
    }

}
pub fn vector_variance(slice: &Vec<u32>, average: f64, frac: Option<f64>) -> f64
{
    let mut var_sum = 0.0;

    match frac{
        None => {
            for &val in slice{
                let dif = average - val as f64;
                var_sum += dif * dif;
            }
        },
        Some(v) => {
            for &val in slice{
                let dif = average - val as f64 / v;
                var_sum += dif * dif;
            }
        }
    }


    var_sum / slice.len() as f64
}



pub fn chi_squared_no_errors(yvals:Vec<f64>,fitvals:Vec<f64>) -> f64{
    let mut sum = 0.;
    for j in 0..yvals.len(){
        let diff = yvals[j]-fitvals[j];
        sum += diff*diff;
    }
    sum
}


