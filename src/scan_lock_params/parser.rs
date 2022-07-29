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
/// Do a scan over thresholds and severity. 
pub struct ScanLock{
    #[structopt(long)]
    json: Option<String>,

    #[structopt(long)]
    num_threads:Option<NonZeroUsize>
}

impl ScanLock{
    pub fn parse(&self) -> (ScanLockParams, Value){
        parse(self.json.as_ref())
    }
    pub fn execute(&self){
        let (opt, json) = self.parse();
        run_simulation(opt,json,self.num_threads)
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ScanLockParams{
    pub lockdown_severities: Vec<f64>,
    pub recovery_prob: f64,
    pub trans_prob: f64,
    pub thresholdrange: F64RangeBuilder,
    pub graph_type: GraphType,
    pub system_size:NonZeroUsize,
    pub num_samples: u64,
    pub fraction: bool,
    pub graph_seed: u64,
    pub sir_seed: u64,
    pub lockdown: LockdownParameters,
    pub initial_infected: usize,
}
impl Default for ScanLockParams{
    fn default() -> Self{
        let lockdown_severities =vec![0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
        let thresholdrange = F64RangeBuilder{
            start: 0.02,
            end:0.2,
            steps: NonZeroUsize::new(50).unwrap() 
        };
        Self{
            lockdown_severities,
            recovery_prob: DEFAULT_RECOVERY_PROB,
            trans_prob: 0.2,
            thresholdrange,            
            graph_type:GraphType::SmallWorld(0.1),
            num_samples: 100000,
            system_size: NonZeroUsize::new(2048).unwrap() ,
            fraction: true,
            graph_seed:DEFAULT_GRAPH_SEED,
            sir_seed: DEFAULT_SIR_SEED,
            lockdown: LockdownParameters{
                lock_style: LockdownType::Random(0.6),
                lock_threshold: 0.1,
                release_threshold: 0.05,
            },
            initial_infected: DEFAULT_INITIAL_INFECTED,
    
            
        }
    }
}
impl ScanLockParams{
    pub fn name<E>(&self, file_ending:E , num_threads:Option<NonZeroUsize>,particular_sev:f64) -> String where E:Display{
        let k = match num_threads{
            None => "".to_owned(),
            Some(v) => format!("k{}",v)
        };
        
        format!(
            "ver{}ScanLock_ThisFileSev{}SevList{}to{}Thresh{}to{}_{}Lam{}_Gam{}_InInf{}NumSamples{}_GT{}_GS{}_SIRS{}_THR{}_LOCK{}.{}",
            crate::VERSION,
            particular_sev*100.0,
            self.lockdown_severities[0],
            self.lockdown_severities[self.lockdown_severities.len()-1],
            //
            self.thresholdrange.start,
            self.thresholdrange.end,
            self.thresholdrange.steps,
            self.trans_prob,
            self.recovery_prob,
            self.initial_infected,
            self.num_samples,
            
            self.graph_type.name(),
            self.graph_seed,
            self.sir_seed,
            k,
            lockdown_naming_string(self.lockdown.lock_style),
            file_ending


        )
    }
}