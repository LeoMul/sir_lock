use{
    structopt::StructOpt,
    serde::{Serialize, Deserialize},
    serde_json::Value,
    std::num::*,
    crate::lockdown_methods::*,
    crate::misc_types::*,
    crate::json_parsing::*,
};

#[derive(Debug, StructOpt, Clone)]
/// Do a simple sampling simulation and get curves
pub struct SimpleCurves
{
    /// Specify the json file with the options
    /// If not given, an example json will be printed
    #[structopt(long)]
    json: Option<String>,

    /// Number of threads to use
    #[structopt(long)]
    num_threads: Option<NonZeroUsize>
}

impl SimpleCurves {
    pub fn parse(&self) -> (SimpleCurvesParam, Value)
    {
        parse(self.json.as_ref())
    }

    pub fn execute(&self) 
    {
        let (param, json) = self.parse();
        super::execute::execute_sir(param, json, self.num_threads)
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct SimpleCurvesParam
{
    pub system_size: NonZeroUsize,
    pub lambda: f64,
    pub recovery_prob: f64,
    pub graph_type: GraphType,

    pub lockdown: LockdownParameters,

    //pub vaccine_doses: usize,
    pub graph_seed: u64,
    pub sir_seed: u64,

    //pub vaccine_seed: u64,
    pub initial_infected: usize,
    pub samples: usize
}
impl SimpleCurvesParam
{
    pub fn quick_name(
        &self) -> String
    {
        //let j = match num_threads
        //{
        //    None => "".to_owned(),
        //    Some(v) => format!("j{}", v)
        //};
        format!(
            "v{}SimpleCurves_Lock{}R{}L{}N{}Rec{}Trans{}InInf{}Sam{}Graph{}GS{}SS{}",
            crate::VERSION,
            lockdown_naming_string(self.lockdown.lock_style),
            self.lockdown.release_threshold,
            self.lockdown.lock_threshold,
            self.system_size,
            self.recovery_prob,
            self.lambda,
            self.initial_infected,
            self.samples,
            self.graph_type.name(),
            self.graph_seed,
            self.sir_seed
        )
    
    
    }
}
impl Default for SimpleCurvesParam
{
    fn default() -> Self {
        
        Self{
            lambda: DEFAULT_LAMBDA,
            recovery_prob: DEFAULT_RECOVERY_PROB,
            samples: DEFAULT_SAMPLES_SIMPLE_SAMPLE,
            graph_seed: DEFAULT_GRAPH_SEED,
            graph_type: GraphType::SmallWorld(0.1),
            system_size: DEFAULT_SYSTEM_SIZE,
            sir_seed: DEFAULT_SIR_SEED,
            lockdown: LockdownParameters{
                lock_style: LockdownType::Random(0.6),
                lock_threshold: 0.1,
                release_threshold: 0.05,
            },
            initial_infected: DEFAULT_INITIAL_INFECTED
        }
    }
}
        