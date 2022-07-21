use{
    structopt::StructOpt,
    std::{num::*, fs::File, io::BufReader},
    serde::{Serialize, Deserialize},
    crate::sir_model::*,
    net_ensembles::sampling::{HistU32Fast, IntervalOrder, HistogramPartition},
    crate::misc_types::*,
    crate::lockdown_methods::*,
};




#[derive(StructOpt, Debug, Clone)]
/// Continue a Large deviation simulation!
pub struct LDContinueCmdOpts
{
    /// Json file containing options.
    /// If not given it will print an example file instead
    #[structopt(long)]
    json: Option<String>,

    /// Number of threads to use
    #[structopt(long)]
    pub num_threads: Option<NonZeroUsize>,

    /// Do not create a savestate at the end
    #[structopt(long)]
    pub no_save: bool,
}


impl LDContinueCmdOpts
{

    pub fn execute(&self, instant: std::time::Instant)
    {
        if let Some(num) = self.num_threads
        {
            rayon::ThreadPoolBuilder::new()
                .num_threads(num.get())
                .build_global()
                .unwrap();
        }

        let opts: LDContinueOpts = match self.json.as_ref()
        {
            None => {
                let def = LDContinueOpts::default();
                serde_json::to_writer_pretty(
                    std::io::stdout(), 
                    &def
                ).unwrap();
                std::process::exit(0)
            },
            Some(filename) => {
                let file = File::open(filename).expect("Unable to open json file");
                let reader = BufReader::new(file);
                serde_json::from_reader(reader).expect("Parsing error for json")
            }
        };
        
        
        crate::large_deviations::load_high_degree_rewl(
            opts, 
            instant, 
            self.no_save, 
        )

    }

}

#[derive(Clone, Serialize, Deserialize)]
pub struct LDContinueOpts
{
    pub file_name: String,
    pub seconds: Option<NonZeroU64>,
    pub minutes: Option<NonZeroU64>,
    pub days: Option<NonZeroU64>,
    pub hours: Option<NonZeroU64>,
    pub change_step_size: Option<NonZeroUsize>,
    pub rewltype: RewlType
}
impl Default for LDContinueOpts
{
    fn default() -> Self
    {
        Self{
            file_name: "insert file with binary save".to_owned(),
            seconds: None,
            minutes: None,
            hours: None,
            days: None,
            change_step_size: None,
            rewltype: RewlType::None
        }
    }
}




impl LDContinueOpts
{
    pub fn seconds(&self) -> u64
    {
        let seconds = ((get_or_0(self.days) * 24
            + get_or_0(self.hours) )* 60
            + get_or_0(self.minutes)) * 60
            + get_or_0(self.seconds);

        if seconds == 0 {
            eprintln!("Error, you have to specify a time for the simulation. Time has to be at least 1 second");
            panic!("Insufficient time!");
        }
        seconds
    }
}
#[inline]
fn get_or_0(val: Option<NonZeroU64>) -> u64
{
    match val {
        None => 0,
        Some(v) => v.get()
    }
}

pub struct LDLdOpts
{
    pub energy: MeasureType,
    pub allowed_seconds: u64,
    pub quick_name:  Box<dyn Fn (Option<usize>, &str, LargeDeviationMode) -> String>,
    pub value: Vec<serde_json::Value>,
    pub no_save: bool
}
#[derive(StructOpt, Debug, Clone)]
/// Large deviation REWL using lockdowns!
pub struct LDOptsLD
{
    /// Specify the json file with the options
    /// If not given, an example json will be printed
    #[structopt(long)]
    pub json: Option<String>,

    /// Number of threads to use
    #[structopt(long)]
    pub num_threads: Option<NonZeroUsize>,

    /// Print other example
    #[structopt(long)]
    pub example2: bool,

    /// Do not create a savestate at the end
    #[structopt(long)]
    pub no_save: bool,
}

impl LDOptsLD{
    pub fn execute(self, instant: std::time::Instant)
    {
        if self.example2{
            let o = LDLDparam::example2();
            let _ = serde_json::to_writer_pretty(
                std::io::stdout(),
                &o
            );
            std::process::exit(0);
        }
        super::execute_large_dev(self, instant)
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct LDLDparam
{
    pub system_size: NonZeroUsize,
    pub recovery_prob: f64,
    pub graph_type: GraphType,
    pub graph_seed: u64,
    pub sir_seed: u64,
    pub f_threshold: f64,
    pub lambda: f64,
    pub large_deviation_param: LargeDeviationParam,
    pub lockdownparams:LockdownParameters,
    pub histograms: HistogramCreator,
    pub sweep_size: NonZeroUsize,
    pub step_size: usize,
    pub walkers_per_interval: NonZeroUsize,
    pub energy: MeasureType,
    pub seconds: u64,
    pub minutes: u64,
    pub hours: u64,
    pub days: u64,
    pub times_continued: Option<NonZeroUsize>
}

pub enum LargeDeviationMode {
    Rewl,
    Rees
}

impl LargeDeviationMode
{
    fn name(&self) -> &'static str
    {
        match self{
            Self::Rewl => "Rewl",
            Self::Rees => "Rees"
        }
    }
}

impl LDLDparam
{
    pub fn quick_name(&self, interval: Option<usize>, end: &str, mode: LargeDeviationMode) -> String
    {
        let interval = match interval{
            None => "".to_owned(),
            Some(i) => format!("_{i}")
        };

        let num_intervals = self.histograms.num_intervals();
        let interval_info = self.histograms.start_end_inlcusive();

        let times = match self.times_continued
        {
            None => "".to_owned(),
            Some(val) => format!("_x{val}")
        };

        format!(
            "v{}_{}_LDHD_N{}r{}t{}LockType{}rel{}lockth{}G{}GS{}SS{}_NI{num_intervals}_{interval_info}{interval}{times}MarkovStepSize{}.{end}",
            crate::VERSION,
            mode.name(),
            self.system_size,
            self.recovery_prob,
            self.lambda,
            lockdown_naming_string(self.lockdownparams.lock_style),
            self.lockdownparams.release_threshold,
            self.lockdownparams.lock_threshold,
            self.graph_type.name(),
            self.graph_seed,
            self.step_size,
            self.sir_seed
        )
    }

    pub fn times_plus_1(&mut self)
    {
        self.times_continued = match self.times_continued{
            None => Some(ONE),
            Some(val) => NonZeroUsize::new(val.get() + 1)
        };
    }

    pub fn example2() -> Self
    {
        let mut opt = Self::default();
        let i1 = Interval{
            start: 1,
            end_inlcusive: 120
        };
        let i2 = Interval{
            start: 100,
            end_inlcusive: 200
        };
        opt.histograms = HistogramCreator::Manual(
            vec![i1, i2]
        );
        opt
    }

    pub fn allowed_seconds(&self) -> u64
    {
        self.seconds + 60 * self.minutes + 60*60 * self.hours + 60*60*24*self.days
    }
}

impl Default for LDLDparam
{
    fn default() -> Self
    {
        Self{
            system_size: DEFAULT_SYSTEM_SIZE,
            recovery_prob: DEFAULT_RECOVERY_PROB,
            graph_type: GraphType::Barabasi(2,10),
            graph_seed: DEFAULT_GRAPH_SEED,
            sir_seed: DEFAULT_SIR_SEED,
            f_threshold: DEFAULT_F_THRESHOLD,
            lambda: DEFAULT_LAMBDA,
            large_deviation_param: LargeDeviationParam
            {
                time_steps: ONE,
                markov_seed: DEFAULT_MARKOV_SEED
            },
            lockdownparams: LockdownParameters{
                lock_style: LockdownType::Random(191905810985091580,0.6),
                lock_threshold: 0.1,
                release_threshold: 0.05
            },
            histograms: HistogramCreator::Automatic(
                Intervals{
                    start: None,
                    end_inlcusive: None,
                    num_intervals: 2,
                    overlap: None,
                    greedy_search_steps: None
                }
            ),
            walkers_per_interval: ONE,
            step_size: DEFAULT_MARKOV_STEP_SIZE,
            sweep_size: DEFAULT_SWEEP_SIZE,
            energy: MeasureType::C,
            seconds: 0,
            minutes: 5,
            hours: 1,
            days: 0,
            times_continued: None
        }
    }
}

#[derive(Clone, Copy, Deserialize, Serialize)]
pub struct Interval
{
    pub start: u32,
    pub end_inlcusive: u32
}

impl Interval{
    pub fn is_valid(&self) -> bool
    {
        self.start < self.end_inlcusive
    }

    pub fn to_hist(&self) -> HistU32Fast
    {
        HistU32Fast::new_inclusive(self.start, self.end_inlcusive)
            .expect("unable to create hist")
    }
}

#[derive(Clone, Copy, Deserialize, Serialize)]
pub struct Intervals
{
    pub start: Option<u32>,
    pub end_inlcusive: Option<u32>,
    pub num_intervals: usize,
    pub overlap: Option<NonZeroUsize>,
    pub greedy_search_steps: Option<NonZeroUsize>
}

#[derive(Clone, Deserialize, Serialize)]
pub enum HistogramCreator
{
    Manual(Vec<Interval>),
    Automatic(Intervals)
}

impl HistogramCreator
{
    pub fn greedy_search_steps(&self) -> Option<NonZeroUsize>
    {
        match self{
            Self::Automatic(i) => {
                i.greedy_search_steps
            },
            _ => None
        }
    }

    pub fn num_intervals(&self) -> usize
    {
        match self{
            Self::Manual(v) => v.len(),
            Self::Automatic(i) => i.num_intervals
        }
    }

    pub fn start_end_inlcusive(&self) -> String
    {
        match self{
            Self::Automatic(a) => {
                format!(
                    "{}-{}",
                    a.start.unwrap_or(1),
                    if let Some(v) = a.end_inlcusive {
                        v.to_string()
                    } else {
                        "end".to_string()
                    }
                )
            },
            Self::Manual(v) => {
                let first = v.first().unwrap();
                let mut start = first.start;
                let mut end = first.end_inlcusive;
                for i in v.iter()
                {
                    if i.start < start{
                        start = i.start;
                    }
                    if i.end_inlcusive > end {
                        end = i.end_inlcusive;
                    }
                }
                format!("{start}-{end}")
            }
        }
    }

    pub fn create(&self, min_possible: u32, max_possible: u32, walkers_per_interval: NonZeroUsize) -> Vec<HistU32Fast>
    {
        let mut hists: Vec<_> = match self
        {
            Self::Manual(intervals) =>
            {
                assert!(
                    intervals.iter().all(Interval::is_valid),
                    "An interval is invalid!"
                );
                intervals.iter()
                    .map(Interval::to_hist)
                    .collect()
            },
            Self::Automatic(auto) => {
                let left = auto.start.unwrap_or(min_possible);
                let right = auto.end_inlcusive.unwrap_or(max_possible);
                let encapsulating = HistU32Fast::new_inclusive(left, right)
                    .expect("Histogram Creation Error");
                
                encapsulating.overlapping_partition(auto.num_intervals, auto.overlap.unwrap_or(ONE).get())
                    .expect("Partition error")
            }
        };
        if walkers_per_interval.get() > 1
        {
            unimplemented!()
        }

        hists.sort_unstable_by(HistU32Fast::left_compare);

        for (left, right) in hists.iter().zip(hists[1..].iter())
        {
            let range_left = left.range_inclusive();
            assert!(
                range_left.contains(&right.left()),
                "Overlap Missing"
            );
        }

        hists
    }
}

pub fn calc_m_ba(model: &mut BALargeDeviationWithLocks) -> Option<u32>
{
    Some(
        if !model.ld_model.markov_changed{
            model.ld_model.energy
        }
        else{
            model.ld_model.energy = model.ld_energy_m();
            model.ld_model.energy 
        })
}

pub fn calc_c_ba(model: &mut BALargeDeviationWithLocks) -> Option<u32>
{
    Some(
        if !model.ld_model.markov_changed{
            model.ld_model.energy
        }
        else{
            model.ld_energy_m();
            model.ld_model.energy = model.calculate_ever_infected() as u32;
            model.ld_model.energy 
        })
}
pub fn calc_m_sw(model: &mut SWLargeDeviationWithLocks) -> Option<u32>
{
    Some(
        if !model.ld_model.markov_changed{
            model.ld_model.energy
        }
        else{
            model.ld_model.energy = model.ld_energy_m();
            model.ld_model.energy 
        })
}

pub fn calc_c_sw(model: &mut SWLargeDeviationWithLocks) -> Option<u32>
{
    Some(
        if !model.ld_model.markov_changed{
            model.ld_model.energy
        }
        else{
            model.ld_energy_m();
            model.ld_model.energy = model.calculate_ever_infected() as u32;
            model.ld_model.energy 
        })
}

pub fn energy_function_returner_ba(measure_type:MeasureType) -> impl Fn(&mut BALargeDeviationWithLocks) -> Option<u32>  + Sync + Send + Copy{



    match measure_type{
        MeasureType::C => calc_c_ba,
        MeasureType::M => {
            calc_m_ba
        }
    }
    

}

pub fn energy_function_returner_sw(measure_type:MeasureType) -> impl Fn(&mut SWLargeDeviationWithLocks) -> Option<u32>  + Sync + Send + Copy{



    match measure_type{
        MeasureType::C => calc_c_sw,
        MeasureType::M => {
            calc_m_sw
        }
    }
    

}