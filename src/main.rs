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
pub mod life_span_size_fitting_lambdascan;
pub mod life_span_size_fitting_threshscan;
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
pub mod scan_gamma;
pub mod scan_lock_params;
pub mod critical_threshold;
pub mod connectedcomponent;
pub mod simplecurves;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

fn main() {
    let start_time = Instant::now();
    let opt = CmdOption::from_args();
    match opt{
        CmdOption::ScanLambda(o) => o.execute(),
        CmdOption::ScanGamma(o) => o.execute(),
        CmdOption::ScanLambdaGamma(o) => o.execute(), 
        CmdOption::ScanThreshParams(o) => o.execute(),
        CmdOption::LifeSpanHist(o) => o.execute(),
        CmdOption::LifeSpanSizeFittingLambda(o) => o.execute(),
        CmdOption::LifeSpanSizeFittingThresh(o) => o.execute(),
        CmdOption::CriticalLambda(o) => o.execute(),
        CmdOption::CriticalThresh(o) => o.execute(),
        CmdOption::PropTest(o) => o.execute(),
        CmdOption::ScanLambdaThresh(o) => o.execute(),
        CmdOption::LargeDevSimpleSample(o) => o.execute(),
        CmdOption::SimpleSample(o) => o.execute(),
        CmdOption::LargeDeviationsLD(o) => o.execute(start_time),
        CmdOption::LargeDeviationsLDContinue(o) => o.execute(start_time),
        CmdOption::REESLd(o) => o.execute(start_time),
        CmdOption::ConnectedComponent(o) => o.execute(),
        CmdOption::SimpleCurves(o) => o.execute(),

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
    ScanGamma(scan_gamma::ScanGamma),
    ScanLambdaThresh(scan_lambda_lock_thresh::ScanLambdaThresh),
    ScanLambdaGamma(scan_lambda_gamma::ScanLambdaGamma),
    ScanThreshParams(scan_lock_params::ScanLock),
    LifeSpanHist(lifespanhist::LifeSpan),
    LifeSpanSizeFittingLambda(life_span_size_fitting_lambdascan::LifespanSizeFitting),
    LifeSpanSizeFittingThresh(life_span_size_fitting_threshscan::LifespanSizeFittingThresh),
    CriticalLambda(critical_lambda::CriticalLambda),
    CriticalThresh(critical_threshold::CriticalThresh),
    PropTest(prop_test::PropTest),
    LargeDevSimpleSample(large_deviations::SimpleSampleCMDopts),
    SimpleSample(simple_sampling::SimpleSampleScan),
    LargeDeviationsLD(large_deviations::LDOptsLD),
    LargeDeviationsLDContinue(large_deviations::LDContinueCmdOpts),
    REESLd(large_deviations::ReesOpts),
    ConnectedComponent(connectedcomponent::ConnectedComponent),
    SimpleCurves(simplecurves::SimpleCurves),

    
}