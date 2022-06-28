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
pub struct CriticalLambda{
    #[structopt(long)]
    json: Option<String>,

    #[structopt(long)]
    num_threads:Option<NonZeroUsize>
}

impl CriticalLambda{
    pub fn parse(&self) -> (CriticalLambdaParams, Value){
        parse(self.json.as_ref())
    }
    pub fn execute(&self){
        let (opt, json) = self.parse();
        run_simulation(opt,json,self.num_threads)
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct CriticalLambdaParams{
    pub system_size_range: UsizeRangeBuilder,
    pub recovery_prob: f64,
    pub lambda_range: F64RangeBuilder,
    pub graph_type: GraphType,
    pub num_networks: u64,
    pub samples_per_lambda: u64,
    pub fraction: bool,
    pub graph_seed: u64,
    pub sir_seed: u64,
}
impl Default for CriticalLambdaParams{
    fn default() -> Self{
        let system_size_range_def = UsizeRangeBuilder{
            start: 2100,
            end: 3600,
            steps: NonZeroUsize::new(3).unwrap()
        };
        let trans_prob_range = F64RangeBuilder{
            start: 0.05,
            end:0.25,
            steps: NonZeroUsize::new(20).unwrap() 
        };
        Self{
            system_size_range: system_size_range_def,
            recovery_prob: DEFAULT_RECOVERY_PROB,
            lambda_range: trans_prob_range,
            graph_type:GraphType::SmallWorld(0.1),
            num_networks: 10000,
            samples_per_lambda: 1000,
            fraction: true,
            graph_seed:DEFAULT_GRAPH_SEED,
            sir_seed: DEFAULT_SIR_SEED
        }
    }
}
impl CriticalLambdaParams{
    pub fn name<E>(&self, file_ending:E , num_threads:Option<NonZeroUsize>) -> String where E:Display{
        let k = match num_threads{
            None => "".to_owned(),
            Some(v) => format!("k{}",v)
        };
        format!(
            "ver{}CriticalLam_Size{}to{}_{}Lam{}to{}_{}_Sam{}NumNet{}_GT{}_GS{}_SIRS{}_THR{}.{}",
            crate::VERSION,
            self.system_size_range.start,
            self.system_size_range.end,
            self.system_size_range.steps,
            //self.recovery_prob,
            self.lambda_range.start,
            self.lambda_range.end,
            self.lambda_range.steps,
            self.samples_per_lambda,
            self.num_networks,
            
            self.graph_type.name(),
            self.graph_seed,
            self.sir_seed,
            k,
            file_ending


        )
    }
}