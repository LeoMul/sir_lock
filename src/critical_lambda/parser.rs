use std::fmt::Display;

use{
    super::*,
    structopt::StructOpt,
    std::num::*,
    crate::json_parsing::*,
    serde::{Serialize, Deserialize},
    serde_json::Value,
   
    crate::misc_types::*,
    crate::lockdown_methods::*,
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
    pub system_size_range: Vec<usize>,
    pub recovery_prob: f64,
    pub lambda_range: F64RangeBuilder,
    pub graph_type: GraphType,
    pub num_networks: u64,
    pub fraction: bool,
    pub graph_seed: u64,
    pub sir_seed: u64,
    pub lockdown: LockdownParameters,
}
impl Default for CriticalLambdaParams{
    fn default() -> Self{
        let system_size_range_def =vec![200,600,1200,4800];
        let trans_prob_range = F64RangeBuilder{
            start: 0.05,
            end:0.25,
            steps: NonZeroUsize::new(200).unwrap() 
        };
        Self{
            system_size_range: system_size_range_def,
            recovery_prob: DEFAULT_RECOVERY_PROB,
            lambda_range: trans_prob_range,
            graph_type:GraphType::Barabasi(2,10),
            num_networks: 100000,
            fraction: true,
            graph_seed:DEFAULT_GRAPH_SEED,
            sir_seed: DEFAULT_SIR_SEED,
            lockdown: LockdownParameters{
                lock_style: LockdownType::LimitContacts(2),
                dynamic_bool: false,
                lock_threshold: 0.1,
                release_threshold: 0.05,
            }
            
        }
    }
}
impl CriticalLambdaParams{
    pub fn name<E>(&self, file_ending:E , num_threads:Option<NonZeroUsize>,particular_n:usize) -> String where E:Display{
        let k = match num_threads{
            None => "".to_owned(),
            Some(v) => format!("k{}",v)
        };
        format!(
            "ver{}CriticalLam_ThisFileN{}Size{}to{}Lam{}to{}_{}_NumNet{}_GT{}_GS{}_SIRS{}_THR{}_LOCK{}.{}",
            crate::VERSION,
            particular_n,
            self.system_size_range[0],
            self.system_size_range[self.system_size_range.len()-1],
            //self.recovery_prob,
            self.lambda_range.start,
            self.lambda_range.end,
            self.lambda_range.steps,
            
            self.num_networks,
            
            self.graph_type.name(),
            self.graph_seed,
            self.sir_seed,
            k,
            lockdown_naming_string(self.lockdown.lock_style),
            file_ending


        )
    }
}