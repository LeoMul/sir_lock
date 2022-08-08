use std::fmt::Display;

use{
    super::*,
    structopt::StructOpt,
    std::num::*,
    crate::json_parsing::*,
    serde::{Serialize, Deserialize},
    serde_json::Value,
    crate::lockdown_methods::*,
    crate::misc_types::*,

};
#[derive(Debug, StructOpt, Clone)]
///3D Plot, Lambda, Gamma
pub struct ScanLambdaGamma{
    #[structopt(long)]
    json: Option<String>,

    #[structopt(long)]
    num_threads:Option<NonZeroUsize>
}

impl ScanLambdaGamma{
    pub fn parse(&self) -> (ScanLambdaGammaParams, Value){
        parse(self.json.as_ref())
    }
    pub fn execute(&self){
        let (opt, json) = self.parse();
        run_simulation(opt,json,self.num_threads)
    }
}


#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ScanLambdaGammaParams{
    pub system_size: NonZeroUsize,
    pub lambda_range: F64RangeBuilder,
    pub gamma_range: F64RangeBuilder,
    //pub recovery_prob: f64,
    pub graph_type: GraphType,
    pub samples_per_step: u64,
    //pub measure: MeasureType,
    pub fraction: bool,
    pub graph_seed: u64,
    pub sir_seed: u64,
    pub lockdown: LockdownType,
    pub lock_thresh: f64,
    pub rel_thresh: f64,
    pub compare_nolock: bool,
    pub initial_infected: usize
}

impl Default for ScanLambdaGammaParams{
    fn default() -> Self{
        let lambda_range_def = F64RangeBuilder{
            start:0.0,
            end: 1.0,
            steps: NonZeroUsize::new(20).unwrap()
        };
        let gamma_range_def = F64RangeBuilder{
            start:0.01,
            end: 1.0,
            steps: NonZeroUsize::new(20).unwrap()
        };
        Self{
            lambda_range:lambda_range_def,
            gamma_range:gamma_range_def,
            system_size: DEFAULT_SYSTEM_SIZE,
            //recovery_prob:DEFAULT_RECOVERY_PROB,
            graph_type: GraphType::SmallWorld(0.1),
            samples_per_step: DEFAULT_SAMPLES_PER_STEP,
            //measure: MeasureType::C,
            fraction: true,
            graph_seed: DEFAULT_GRAPH_SEED,
            sir_seed: DEFAULT_SIR_SEED,
            lockdown:LockdownType::Random(0.5541),
            lock_thresh:0.165,
            rel_thresh:0.0206,
            compare_nolock: false,
            initial_infected: DEFAULT_INITIAL_INFECTED

        }
    }
}

impl ScanLambdaGammaParams{
    pub fn name<E>(&self, measure: MeasureType, file_ending:E , num_threads:Option<NonZeroUsize>) -> String where E:Display{
        let k = match num_threads{
            None => "".to_owned(),
            Some(v) => format!("k{}",v)
        };
        format!(
            "ver{}LamGamScan_{}_N{}t{}-{}_{}r{}-{}_{}InInf{}SamStep{}_Graph{}_GSeed{}_SS{}_THR{}_LOCK{}_LT{}_RT{}_COMP{}.{}",
            crate::VERSION,
            measure.name(),
            self.system_size,
            //self.recovery_prob,
            self.lambda_range.start,
            self.lambda_range.end,
            self.lambda_range.steps,
            self.gamma_range.start,
            self.gamma_range.end,
            self.gamma_range.steps,
            self.initial_infected,
            self.samples_per_step,
            self.graph_type.name(),
            self.graph_seed,
            self.sir_seed,
            k,
            lockdown_naming_string(self.lockdown),
            self.lock_thresh,
            self.rel_thresh,
            self.compare_nolock,
            file_ending


        )
    }
}