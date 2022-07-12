use{
    std::{
        time::Instant
    },
    structopt::StructOpt,
    indicatif::*
};

pub mod sir_model;
pub mod life_span_data_collection;
pub mod grid;
pub mod life_span_size_fitting;
//use crate::life_span_size_fitting::*;
pub mod misc_types;
use crate::misc_types::*;
pub mod scan_lambda;
pub mod scan_lambda_gamma;
pub mod stats_methods;
pub mod critical_lambda;
pub mod json_parsing;
pub mod prop_test;
pub mod lockdown_methods;
pub mod scan_lambda_lock_thresh;
//pub mod time_graph;
pub mod lifespanhist;
pub mod large_deviations;
pub mod simple_sampling;


pub const VERSION: &str = env!("CARGO_PKG_VERSION");

fn main() {
    let start_time = Instant::now();
    let opt = CmdOption::from_args();
    match opt{
        CmdOption::ScanLambda(o) => o.execute(),
        CmdOption::ScanLambdaGamma(o) => o.execute(), 
        CmdOption::LifeSpan(o) => o.execute(),
        CmdOption::LifespanSizeFitting(o) => o.execute(),
        CmdOption::CriticalLambda(o) => o.execute(),
        CmdOption::PropTest(o) => o.execute(),
        CmdOption::ScanLambdaThresh(o) => o.execute(),
        CmdOption::LargeDevSimpleSample(o) => o.execute(),
        CmdOption::SimpleSample(o) => o.execute(),
        CmdOption::LargeDeviationsLD(o) => o.execute(start_time),
        CmdOption::LargeDeviationsLDContinue(o) => o.execute(start_time),

    }
    println!("Execution took {}",humantime::format_duration(start_time.elapsed()))
    
}

pub fn indication_bar(len: u64) -> ProgressBar
{
        // for indication on when it is finished
        let bar = ProgressBar::new(len);
        bar.set_style(ProgressStyle::default_bar()
            .template("{msg} [{elapsed_precise} - {eta_precise}] {wide_bar}"));
        bar
}


#[derive(Debug, StructOpt, Clone)]
#[structopt(about = "Simulations for the SIR Model with lockdowns!")]
pub enum CmdOption 
{
    ScanLambda(scan_lambda::ScanLambda),
    ScanLambdaThresh(scan_lambda_lock_thresh::ScanLambdaThresh),
    ScanLambdaGamma(scan_lambda_gamma::ScanLambdaGamma),
    LifeSpan(lifespanhist::LifeSpan),
    LifespanSizeFitting(life_span_size_fitting::LifespanSizeFitting),
    CriticalLambda(critical_lambda::CriticalLambda),
    PropTest(prop_test::PropTest),
    LargeDevSimpleSample(large_deviations::SimpleSampleCMDopts),
    SimpleSample(simple_sampling::SimpleSample),
    LargeDeviationsLD(large_deviations::BALDOptsLD),
    LargeDeviationsLDContinue(large_deviations::BALDContinueCmdOpts)

    
}