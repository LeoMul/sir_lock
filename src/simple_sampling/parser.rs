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
/// Do a simple sampling simulation and get P(M) and P(C)
pub struct SimpleSampleScan
{
    /// Specify the json file with the options
    /// If not given, an example json will be printed
    #[structopt(long)]
    json: Option<String>,

    /// Number of threads to use
    #[structopt(long)]
    num_threads: Option<NonZeroUsize>
}

impl SimpleSampleScan {
    pub fn parse(&self) -> (SimpleSampleParam, Value)
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
pub struct SimpleSampleParam
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
impl SimpleSampleParam
{
    pub fn quick_name(
        &self,  
        num_threads: Option<NonZeroUsize>,
        measure_type: MeasureType
    ) -> String
    {
        let j = match num_threads
        {
            None => "".to_owned(),
            Some(v) => format!("j{}", v)
        };
        let s = if let LockdownType::Random(n) = self.lockdown.lock_style{
            format!("percent{n}")
        }else{
            "".to_owned()
        };
        format!(
            "v{}SimpleSampling_Lock{}{s}R{}L{}_Measure{}N{}Rec{}Trans{}InInf{}Sam{}Graph{}GS{}SS{}THR{}.dat",
            crate::VERSION,
            lockdown_naming_string(self.lockdown.lock_style),
            self.lockdown.release_threshold,
            self.lockdown.lock_threshold,
            measure_type.name(),
            self.system_size,
            self.recovery_prob,
            self.lambda,
            self.initial_infected,
            self.samples,
            self.graph_type.name(),
            self.graph_seed,
            self.sir_seed,
            j
        )
    
    
    }
}
impl Default for SimpleSampleParam
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
                lock_style: LockdownType::Random(DEFAULT_RANDOM_LOCKDOWN_FRAC),
                lock_threshold: DEFAULT_LOCKDOWN_THRESHOLD,
                release_threshold: DEFAULT_RELEASE_THRESHOLD,
            },
            initial_infected: DEFAULT_INITIAL_INFECTED
        }
    }
}
        