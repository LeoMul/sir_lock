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
///Scan over lambda &system size to find max lifespans!
pub struct LifespanSizeFitting{
    #[structopt(long)]
    json: Option<String>,

    #[structopt(long)]
    num_threads:Option<NonZeroUsize>
}

impl LifespanSizeFitting{
    pub fn parse(&self) -> (LifespanSizeFittingParams, Value){
        parse(self.json.as_ref())
    }
    pub fn execute(&self){
        let (opt, json) = self.parse();
        run_simulation(opt,json,self.num_threads)
    }
}


#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct LifespanSizeFittingParams{
    pub system_size_range: Vec<usize>,
    pub recovery_prob: f64,
    pub trans_prob_range: F64RangeBuilder,
    pub graph_type: GraphType,
    pub num_networks: u64,
    pub fraction: bool,
    pub graph_seed: u64,
    pub sir_seed: u64,
    pub lockdown: LockdownParameters,
    pub lifespanpercent: f64,
    pub initial_infected: usize
}

impl Default for LifespanSizeFittingParams{
    fn default() -> Self{
        //let system_size_range_def = UsizeRangeBuilder{
          //  start: 2100,
            //end: 3600,
            //steps: NonZeroUsize::new(3).unwrap()
        //};
        let system_size_range_def =vec![200,400,600,1000,1200,1600,2000,2400,2800,3200];
        let trans_prob_range = F64RangeBuilder{
            start: 0.05,
            end:0.25,
            steps: NonZeroUsize::new(20).unwrap() 
        };
        Self{
            system_size_range: system_size_range_def,
            recovery_prob: DEFAULT_RECOVERY_PROB,
            trans_prob_range,
            graph_type:GraphType::Barabasi(2,10),
            num_networks: 10000,
            fraction: true,
            graph_seed:DEFAULT_GRAPH_SEED,
            sir_seed: DEFAULT_SIR_SEED,
            lockdown: LockdownParameters{
                lock_style: LockdownType::Random(2151515153,0.6),
                lock_threshold: 0.1,
                release_threshold: 0.05,
            },
            lifespanpercent: 0.98,
            initial_infected: DEFAULT_INITIAL_INFECTED
        }
    }
}

impl LifespanSizeFittingParams{
    pub fn name<E>(&self, file_ending:E , num_threads:Option<NonZeroUsize>,system_size:usize) -> String where E:Display{
        let k = match num_threads{
            None => "".to_owned(),
            Some(v) => format!("k{}",v)
        };
        format!(
            "ver{}LifeSpanFitting_Size{}_r{}_t{}{}-{}_InInf{}_NumNet{}_gr{}_gs{}_sir{}_thr{}_LSPER{}_LOCK{}.{}",
            crate::VERSION,
            system_size,
            self.recovery_prob,
            self.trans_prob_range.start,
            self.trans_prob_range.end,
            self.trans_prob_range.steps,
            self.initial_infected,
            self.num_networks,
            self.graph_type.name(),
            self.graph_seed,
            self.sir_seed,
            k,
            self.lifespanpercent,
            lockdown_naming_string(self.lockdown.lock_style),
            file_ending


        )
    }
}