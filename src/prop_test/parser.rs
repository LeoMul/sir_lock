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
///Check a single propagation is working
pub struct PropTest{
    #[structopt(long)]
    json: Option<String>,

    #[structopt(long)]
    num_threads:Option<NonZeroUsize>
}

impl PropTest{
    pub fn parse(&self) -> (PropTestParams, Value){
        parse(self.json.as_ref())
    }
    pub fn execute(&self){
        let (opt, json) = self.parse();
        run_simulation(opt,json,self.num_threads)
    }
}


#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct PropTestParams{
    pub system_size: NonZeroUsize,
    pub lambda: f64,
    pub recovery_prob: f64,
    pub graph_type: GraphType,
    pub samples_per_step: u64,
    pub measure: MeasureType,
    pub fraction: bool,
    pub graph_seed: u64,
    pub sir_seed: u64,
}

impl Default for PropTestParams{
    fn default() -> Self{
        
        Self{
            lambda:0.1763,
            system_size: DEFAULT_SYSTEM_SIZE,
            recovery_prob:DEFAULT_RECOVERY_PROB,
            graph_type: GraphType::SmallWorld(0.1),
            samples_per_step: DEFAULT_SAMPLES_PER_STEP,
            measure: MeasureType::C,
            fraction: true,
            graph_seed: DEFAULT_GRAPH_SEED,
            sir_seed: DEFAULT_SIR_SEED

        }
    }
}