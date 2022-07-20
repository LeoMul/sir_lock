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
///Perform a scan over Lambda and the Lockdown threshold. One network.
pub struct ScanLambdaThresh{
    #[structopt(long)]
    json: Option<String>,

    #[structopt(long)]
    num_threads:Option<NonZeroUsize>
}

impl ScanLambdaThresh{
    pub fn parse(&self) -> (ScanLambdaThreshParams, Value){
        parse(self.json.as_ref())
    }
    pub fn execute(&self){
        let (opt, json) = self.parse();
        run_simulation(opt,json,self.num_threads)
    }
}


#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ScanLambdaThreshParams{
    pub system_size: NonZeroUsize,
    pub lambda_range: F64RangeBuilder,
    pub lockt_range: F64RangeBuilder,
    pub recovery_prob: f64,
    pub graph_type: GraphType,
    pub samples_per_step: u64,
    pub fraction: bool,
    pub graph_seed: u64,
    pub sir_seed: u64,
    pub lockdown: LockdownParameters,
}

impl Default for ScanLambdaThreshParams{
    fn default() -> Self{
        let lambda_range_def = F64RangeBuilder{
            start:0.0,
            end: 1.0,
            steps: NonZeroUsize::new(20).unwrap()
        };
        let lock_thresh_range_def = F64RangeBuilder{
            start:0.06,
            end: 1.0,
            steps: NonZeroUsize::new(20).unwrap()
        };
        Self{
            lambda_range:lambda_range_def,
            lockt_range:lock_thresh_range_def,
            system_size: DEFAULT_SYSTEM_SIZE,
            recovery_prob:DEFAULT_RECOVERY_PROB,
            graph_type: GraphType::Barabasi(2,10),
            samples_per_step: DEFAULT_SAMPLES_PER_STEP,
            fraction: true,
            graph_seed: DEFAULT_GRAPH_SEED,
            sir_seed: DEFAULT_SIR_SEED,
            lockdown:LockdownParameters{
                lock_style: LockdownType::Random(123131315,0.6),
                lock_threshold: 0.1,
                release_threshold: 0.05,
            },
        }
    }
}

impl ScanLambdaThreshParams{
    pub fn name<E>(&self, measure: MeasureType, file_ending:E , num_threads:Option<NonZeroUsize>) -> String where E:Display{
        let k = match num_threads{
            None => "".to_owned(),
            Some(v) => format!("k{}",v)
        };
        format!(
            "ver{}LamThreshScan_{}_N{}t{}-{}_{}r{}LockThresh{}-{}_{}SamStep{}_Graph{}_GSeed{}_SS{}_THR{}_LOCK{}.{}",
            crate::VERSION,
            measure.name(),
            self.system_size,
            //self.recovery_prob,
            self.lambda_range.start,
            self.lambda_range.end,
            self.lambda_range.steps,
            self.recovery_prob,
            self.lockt_range.start,
            self.lockt_range.end,
            self.lockt_range.steps,
            self.samples_per_step,
            self.graph_type.name(),
            self.graph_seed,
            self.sir_seed,
            k,
            lockdown_naming_string(self.lockdown.lock_style),
                    file_ending


        )
    }
}