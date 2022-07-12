//use crate::simple_sample::parser::DEFAULT_SAMPLES_SIMPLE_SAMPLE;
use crate::misc_types::*;
use{
    structopt::StructOpt,
    std::num::*,
    super::*,
    crate::sir_model::*,
    crate::{GraphType, MeasureType},
    serde::{Serialize, Deserialize},
    crate::lockdown_methods::*,
};

#[derive(StructOpt, Debug, Clone)]
pub struct SimpleSampleCMDopts
{
    /// Specify the json file with the options
    /// If not given, an example json will be printed
    #[structopt(long)]
    pub json: Option<String>,

    /// Number of threads to use
    #[structopt(long)]
    pub num_threads: Option<NonZeroUsize>,
}

impl SimpleSampleCMDopts
{
    pub fn execute(self)
    {
        typical_event_sampling(self)
    }
}


#[derive(Serialize, Deserialize)]
pub struct SimpleSampleldParam
{
    pub system_size: NonZeroUsize,
    pub recovery_prob: f64,
    pub graph_type: GraphType,
    pub graph_seed: u64,
    pub sir_seed: u64,
    pub lambda: f64,
    pub large_deviation_param: LargeDeviationParam,
    pub lockdown: LockdownParameters,
    pub energy: MeasureType,
    pub samples: usize,
    pub randomize: Randomize
}

impl SimpleSampleldParam
{
    pub fn quick_name(&self) -> String 
    {
        let rand = match self.randomize{
            Randomize::Random => "Rand".to_owned(),
            Randomize::Markov(m) => {
                format!("msteps{}every{}", m.step_size, m.every)
            }
        };

        format!(
            "{}_LargeDevSimSamN{}_{rand}_trans{}_gamma{}_time{}S{}_order_rand_1.dat",
            crate::VERSION,
            self.system_size,
            self.lambda,
            self.recovery_prob,
            self.large_deviation_param.time_steps,
            self.samples
        )
    }
}

#[derive(Serialize, Debug, Deserialize, Clone, Copy)]
pub enum Randomize {
    Random,
    Markov(MarkovParam)
}

impl Default for Randomize
{
    fn default() -> Self
    {
        Randomize::Markov(
            MarkovParam{
                step_size: NonZeroUsize::new(DEFAULT_MARKOV_STEP_SIZE).unwrap(),
                every: ONE
            }
        )
    }
}

#[derive(Serialize, Debug, Deserialize, Clone, Copy)]
pub struct MarkovParam
{
    pub step_size: NonZeroUsize,
    pub every: NonZeroUsize,
}

impl Default for SimpleSampleldParam
{
    fn default() -> Self
    {
        Self{
            system_size: DEFAULT_SYSTEM_SIZE,
            recovery_prob: DEFAULT_RECOVERY_PROB,
            graph_type: GraphType::Barabasi(2,10),
            graph_seed: DEFAULT_GRAPH_SEED,
            sir_seed: DEFAULT_SIR_SEED,
            lambda: DEFAULT_LAMBDA,
            large_deviation_param: LargeDeviationParam
            {
                time_steps: ONE,
                markov_seed: DEFAULT_MARKOV_SEED
            },
            lockdown: LockdownParameters{
                lock_style: LockdownType::Random(12151515,0.6),
                dynamic_bool: true,
                lock_threshold: 0.1,
                release_threshold: 0.05,
            },
            energy: MeasureType::C,
            samples: DEFAULT_SAMPLES_SIMPLE_SAMPLE,
            randomize: Randomize::default()
        }
    }
}