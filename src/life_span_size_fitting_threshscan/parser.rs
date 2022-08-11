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
///Scan over thresh &system size to find max lifespans!
pub struct LifespanSizeFittingThresh{
    #[structopt(long)]
    json: Option<String>,

    #[structopt(long)]
    num_threads:Option<NonZeroUsize>
}

impl LifespanSizeFittingThresh{
    pub fn parse(&self) -> (LifespanSizeFittingThreshParams, Value){
        parse(self.json.as_ref())
    }
    pub fn execute(&self){
        let (opt, json) = self.parse();
        run_simulation(opt,json,self.num_threads)
    }
}


#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct LifespanSizeFittingThreshParams{
    pub system_size_range: Vec<usize>,
    pub recovery_prob: f64,
    pub trans_prob: f64,
    pub lockdownthreshrange: F64RangeBuilder,
    pub graph_type: GraphType,
    pub num_networks: u64,
    pub fraction: bool,
    pub graph_seed: u64,
    pub sir_seed: u64,
    pub lockdowntype: LockdownType,
    pub releasetype: ReleaseType,
    pub lifespanpercent: f64,
    pub initial_infected: usize
}

impl Default for LifespanSizeFittingThreshParams{
    fn default() -> Self{
        //let system_size_range_def = UsizeRangeBuilder{
          //  start: 2100,
            //end: 3600,
            //steps: NonZeroUsize::new(3).unwrap()
        //};
        let system_size_range_def =vec![200,400,600,1000,1200,1600,2000,2400,2800,3200];
        let lockdownthreshrange = F64RangeBuilder{
            start: 0.05,
            end:0.25,
            steps: NonZeroUsize::new(20).unwrap() 
        };
        Self{
            system_size_range: system_size_range_def,
            recovery_prob: DEFAULT_RECOVERY_PROB,
            trans_prob: 0.2,
            lockdownthreshrange,
            graph_type:GraphType::SmallWorld(0.1),
            num_networks: 10000,
            fraction: true,
            graph_seed:DEFAULT_GRAPH_SEED,
            sir_seed: DEFAULT_SIR_SEED,
            lockdowntype:LockdownType::Random(0.5491),
            releasetype:ReleaseType::FracOfLock(0.125),
            lifespanpercent: 0.98,
            initial_infected: DEFAULT_INITIAL_INFECTED
        }
    }
}

impl LifespanSizeFittingThreshParams{
    pub fn name<E>(&self, file_ending:E , num_threads:Option<NonZeroUsize>,system_size:usize) -> String where E:Display{
        let k = match num_threads{
            None => "".to_owned(),
            Some(v) => format!("k{}",v)
        };
        let s = if let LockdownType::Random(n) = self.lockdowntype{
            format!("percent{n}")
        }else{
            "".to_owned()
        };
        format!(
            "ver{}LifeSpanFittingThreshScan_Size{}_r{}_t{}lt{}{}-{}_InInf{}_NumNet{}_gr{}_gs{}_sir{}_thr{}_LSPER{}_LOCK{}{s}reltype{}.{}",
            crate::VERSION,
            system_size,
            self.recovery_prob,
            self.trans_prob,
            self.lockdownthreshrange.start,
            self.lockdownthreshrange.end,
            self.lockdownthreshrange.steps,
            self.initial_infected,
            self.num_networks,
            self.graph_type.name(),
            self.graph_seed,
            self.sir_seed,
            k,
            self.lifespanpercent,
            lockdown_naming_string(self.lockdowntype),
            self.releasetype.name(),
    
            file_ending


        )
    }
}