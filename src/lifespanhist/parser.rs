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
///Find the histogram and cumulative dist of disease lifespans
pub struct LifeSpan{
    #[structopt(long)]
    json: Option<String>,

    #[structopt(long)]
    num_threads:Option<NonZeroUsize>
}
impl LifeSpan{
    pub fn parse(&self) -> (LifeSpanParams, Value){
        parse(self.json.as_ref())
    }
    pub fn execute(&self){
        let (opt, json) = self.parse();
        run_simulation(opt,json,self.num_threads)
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct LifeSpanParams{
    pub system_size: NonZeroUsize,
    pub recovery_prob: f64,
    pub trans_prob:f64,
    pub graph_type: GraphType,
    pub samples: u64,
    pub fraction: bool,
    pub graph_seed: u64,
    pub sir_seed: u64,
    pub initial_infected:usize
}

impl Default for LifeSpanParams{
    fn default() -> Self{
        Self{
            system_size:DEFAULT_SYSTEM_SIZE,
            recovery_prob: 0.14,
            trans_prob: 0.25,
            graph_type: GraphType::SmallWorld(0.1),
            samples: DEFAULT_SAMPLES_PER_STEP,
            fraction: true,
            graph_seed: DEFAULT_GRAPH_SEED,
            sir_seed: DEFAULT_SIR_SEED,
            initial_infected: DEFAULT_INITIAL_INFECTED

        }
    }
}


impl LifeSpanParams{
    pub fn name<E>(&self,something_else:E, file_ending:E , num_threads:Option<NonZeroUsize>) -> String where E:Display{
        let k = match num_threads{
            None => "".to_owned(),
            Some(v) => format!("k{}",v)
        };
        format!(
            //"ver{}LS_{}_{}_N{}t{}-{}_{}r{}-{}_{}v{}S_G{}GS{}SS{}{}{}.{}",
            "ver{}TIMEGRAPH_{}_S{}_R{}_T{}_InInf{}_G{}_SAM{}_GS{}_SS{}_THREADS{}.{}",
            crate::VERSION,
            something_else,
            self.system_size,
            self.recovery_prob,
            self.trans_prob,
            self.initial_infected,
            self.graph_type.name(),
            self.samples,
            self.graph_seed,
            self.sir_seed,
            k,
            file_ending)
    }
}