use std::fmt::Display;

use{
    super::*,
    structopt::StructOpt,
    std::num::*,
    crate::json_parsing::*,
    serde::{Serialize, Deserialize},
    serde_json::Value,
   
    crate::misc_types::*,

};
#[derive(Debug, StructOpt, Clone)]
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
            sir_seed: DEFAULT_SIR_SEED

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
            "ver{}LS_{}_{}_N{}t{}-{}_{}r{}-{}_{}v{}S_G{}GS{}SS{}{}{}.{}",
            crate::VERSION,
            measure.name(),
            "none", //NO VACCINES
            self.system_size,
            //self.recovery_prob,
            self.lambda_range.start,
            self.lambda_range.end,
            self.lambda_range.steps,
            self.gamma_range.start,
            self.gamma_range.end,
            self.gamma_range.steps,
            0, //VACCINES
            self.samples_per_step,
            self.graph_type.name(),
            self.graph_seed,
            self.sir_seed,
            k,
            file_ending


        )
    }
}