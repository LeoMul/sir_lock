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
/// Do a scan over system sizes and thresh's. Then fit for crit thresh!
pub struct CriticalThresh{
    #[structopt(long)]
    json: Option<String>,

    #[structopt(long)]
    num_threads:Option<NonZeroUsize>
}

impl CriticalThresh{
    pub fn parse(&self) -> (CriticalThreshParams, Value){
        parse(self.json.as_ref())
    }
    pub fn execute(&self){
        let (opt, json) = self.parse();
        run_simulation(opt,json,self.num_threads)
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct CriticalThreshParams{
    pub system_size_range: Vec<usize>,
    pub recovery_prob: f64,
    pub lambda:f64,
    pub thresh_range: F64RangeBuilder,
    pub graph_type: GraphType,
    pub num_networks: u64,
    pub fraction: bool,
    pub graph_seed: u64,
    pub sir_seed: u64,
    pub lockdown: LockdownParameters,
    pub initial_infected: usize,
    pub bootbool: bool,
    pub bootsamples: usize,
    //pub energy: MeasureType
}
impl Default for CriticalThreshParams{
    fn default() -> Self{
        let system_size_range_def =vec![3000,400,600,800,1000,1200,1600,2000,2500,200];
        let thresh_range = F64RangeBuilder{
            start: 0.01,
            end:0.07,
            steps: NonZeroUsize::new(150).unwrap() 
        };
        Self{
            system_size_range: system_size_range_def,
            recovery_prob: DEFAULT_RECOVERY_PROB,
            lambda: 0.2,
            graph_type:GraphType::SmallWorld(0.1),
            num_networks: 100000,
            fraction: true,
            graph_seed:DEFAULT_GRAPH_SEED,
            sir_seed: DEFAULT_SIR_SEED,
            lockdown: LockdownParameters{
                lock_style: LockdownType::Random(0.35),
                lock_threshold: 0.1,
                release_threshold: 0.05,
            },
            thresh_range,
            initial_infected: DEFAULT_INITIAL_INFECTED,
            bootbool: false,
            bootsamples: crate::stats_methods::stats::BOOTSTRAP_SAMPLES,
            //energy: MeasureType::C
            
        }
    }
}
impl CriticalThreshParams{
    pub fn name<E>(&self, file_ending:E , num_threads:Option<NonZeroUsize>,particular_n:usize,energy:MeasureType) -> String where E:Display{
        let k = match num_threads{
            None => "".to_owned(),
            Some(v) => format!("k{}",v)
        };
        let string = if self.bootbool{
            format!("BootSamples{}",self.bootsamples)
        }
        else{
            "".to_owned()
        };
        let s = if let LockdownType::Random(n) = self.lockdown.lock_style{
            format!("percent{n}")
        }else{
            "".to_owned()
        };
        format!(
            "ver{}CriticalThresh{}_ThisFileN{}Size{}to{}Thresh{}to{}_{}_Lam{}_Gam{}_InInf{}{string}NumNet{}_GT{}_GS{}_SIRS{}_THR{}_LOCK{}{s}.{}",
            crate::VERSION,
            energy.name(),
            particular_n,
            self.system_size_range[0],
            self.system_size_range[self.system_size_range.len()-1],
            //
            self.thresh_range.start,
            self.thresh_range.end,
            self.thresh_range.steps,
            self.lambda,
            self.recovery_prob,
            self.initial_infected,
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