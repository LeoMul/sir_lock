use std::fmt::Display;

use{
    super::*,
    structopt::StructOpt,
    std::num::*,
    crate::json_parsing::*,
    serde::{Serialize, Deserialize},
    serde_json::Value,
    
    crate::misc_types::*,
    crate::lockdown_methods::*
};



#[derive(Debug, StructOpt, Clone)]
pub struct ScanLambda{
    #[structopt(long)]
    json: Option<String>,

    #[structopt(long)]
    num_threads:Option<NonZeroUsize>
}

impl ScanLambda{
    pub fn parse(&self) -> (ScanLambdaParams, Value){
        parse(self.json.as_ref())
    }
    pub fn execute(&self){
        let (opt, json) = self.parse();
        run_simulation(opt,json,self.num_threads)
    }
}


#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ScanLambdaParams{
    pub system_size: NonZeroUsize,
    pub lambda_range: F64RangeBuilder,
    pub recovery_prob: f64,
    pub graph_type: GraphType,
    pub samples_per_step: u64,
    pub measure: MeasureType,
    pub fraction: bool,
    pub graph_seed: u64,
    pub sir_seed: u64,
    pub lockdown: LockdownParameters,
    pub complock: LockdownParameters,
    pub compare: bool,
}

impl Default for ScanLambdaParams{
    fn default() -> Self{
        let lambda_range_def = F64RangeBuilder{
            start:0.0,
            end: 1.0,
            steps: NonZeroUsize::new(20).unwrap()
        };
        Self{
            lambda_range:lambda_range_def,
            system_size: DEFAULT_SYSTEM_SIZE,
            recovery_prob:DEFAULT_RECOVERY_PROB,
            graph_type: GraphType::Barabasi(2,10),
            samples_per_step: DEFAULT_SAMPLES_PER_STEP,
            measure: MeasureType::C,
            fraction: true,
            graph_seed: DEFAULT_GRAPH_SEED,
            sir_seed: DEFAULT_SIR_SEED,
            lockdown:LockdownParameters{
                lock_style: LockdownType::Targeted,
                dynamic_bool: true,
                lock_threshold: 0.1,
                release_threshold: 0.05,
            },
            complock:LockdownParameters{
                lock_style: LockdownType::None,
                dynamic_bool: false,
                lock_threshold: 0.1,
                release_threshold: 0.05,
            },
            compare: true,
            

        }
    }
}

impl ScanLambdaParams{
    pub fn name<E>(&self, file_ending:E , num_threads:Option<NonZeroUsize>) -> String where E:Display{
        let k = match num_threads{
            None => "".to_owned(),
            Some(v) => format!("k{}",v)
        };
        format!(
            "ver{}LamScan_{}_N{}r{}t{}-{}_{}SamStep{}_Graph{}_GSeed{}_SS{}_THR{}_LOCK{}_LT{}_RT{}_COMP{}{}.{}",
            crate::VERSION,
            self.measure.name(),
            self.system_size,
            self.recovery_prob,
            self.lambda_range.start,
            self.lambda_range.end,
            self.lambda_range.steps,
            self.samples_per_step,
            self.graph_type.name(),
            self.graph_seed,
            self.sir_seed,
            k,
            lockdown_naming_string(self.lockdown.lock_style),
            self.lockdown.lock_threshold,
            self.lockdown.release_threshold,
            lockdown_naming_string(self.complock.lock_style),
            self.compare,
            file_ending


        )
    }
}