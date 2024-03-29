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
///Perform a scan over Gamma. One network used.
pub struct ScanGamma{
    #[structopt(long)]
    json: Option<String>,

    #[structopt(long)]
    num_threads:Option<NonZeroUsize>
}

impl ScanGamma{
    pub fn parse(&self) -> (ScanGammaParams, Value){
        parse(self.json.as_ref())
    }
    pub fn execute(&self){
        let (opt, json) = self.parse();
        run_simulation(opt,json,self.num_threads)
    }
}


#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ScanGammaParams{
    pub system_size: NonZeroUsize,
    pub gamma_range: F64RangeBuilder,
    pub transmission_prob: f64,
    pub graph_type: GraphType,
    pub samples_per_step: u64,
    pub measure: MeasureType,
    pub fraction: bool,
    pub graph_seed: u64,
    pub sir_seed: u64,
    pub lockdown: LockdownParameters,
    pub complock: LockdownParameters,
    pub compare: bool,
    pub initial_infected:usize
}

impl Default for ScanGammaParams{
    fn default() -> Self{
        let gamma_range_def = F64RangeBuilder{
            start:0.0,
            end: 1.0,
            steps: NonZeroUsize::new(20).unwrap()
        };
        Self{
            gamma_range:gamma_range_def,
            system_size: DEFAULT_SYSTEM_SIZE,
            transmission_prob:DEFAULT_LAMBDA,
            graph_type: GraphType::SmallWorld(0.1),
            samples_per_step: DEFAULT_SAMPLES_PER_STEP,
            measure: MeasureType::C,
            fraction: true,
            graph_seed: DEFAULT_GRAPH_SEED,
            sir_seed: DEFAULT_SIR_SEED,
            lockdown:LockdownParameters{
                lock_style: LockdownType::Random(DEFAULT_RANDOM_LOCKDOWN_FRAC),
                lock_threshold: DEFAULT_LOCKDOWN_THRESHOLD,
                release_threshold: DEFAULT_RELEASE_THRESHOLD,
            },
            complock:LockdownParameters{
                lock_style: LockdownType::None,
                lock_threshold: 0.1,
                release_threshold: 0.05,
            },
            compare: false,
            initial_infected: DEFAULT_INITIAL_INFECTED
            

        }
    }
}

impl ScanGammaParams{
    pub fn name<E>(&self, file_ending:E , num_threads:Option<NonZeroUsize>) -> String where E:Display{
        let k = match num_threads{
            None => "".to_owned(),
            Some(v) => format!("k{}",v)
        };
        format!(
            "ver{}LamScan_{}_N{}t{}r{}-{}_{}InInf{}SamStep{}_Graph{}_GSeed{}_SS{}_THR{}_LOCK{}_LT{}_RT{}_COMP{}{}.{}",
            crate::VERSION,
            self.measure.name(),
            self.system_size,
            self.transmission_prob,
            self.gamma_range.start,
            self.gamma_range.end,
            self.gamma_range.steps,
            self.initial_infected,
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