pub mod sir_model;
use crate::sir_model::sir_states::InfectionState;
use crate::sir_model::base_model::BaseSWModel;
use crate::sir_model::*;
use crate::sir_model::base_model_options::BaseSwOptions;
use{
    std::{
        time::Instant
    },
    structopt::StructOpt,
    indicatif::*
};
use net_ensembles::{SwEnsemble};

use net_ensembles::rand::SeedableRng;
pub mod grid;
use crate::grid::*;
pub mod misc_types;
use crate::misc_types::*;

pub mod scan_lambda;
use crate::scan_lambda::*;
pub mod stats_methods;
use crate::stats_methods::*;



pub mod json_parsing;
use crate::json_parsing::*;

use rand_pcg::Pcg64;
//use rand::SeedableRng;
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

fn main() {
    let num_nodes = 200;
    let rng = Pcg64::seed_from_u64(45);
    let e = SwEnsemble::<InfectionState,Pcg64>::new(num_nodes,0.1,rng);

    let base = BaseSWModel{
        ensemble: e,
        gamma: 0.14,
        lambda:0.1763,
        n:num_nodes,};
    let mut model = SimpleSample::from_base(base, 45);
    let x = model.propagate_until_completion_max();
    println!("{}",x)
}




pub fn indication_bar(len: u64) -> ProgressBar
{
        // for indication on when it is finished
        let bar = ProgressBar::new(len);
        bar.set_style(ProgressStyle::default_bar()
            .template("{msg} [{elapsed_precise} - {eta_precise}] {wide_bar}"));
        bar
}