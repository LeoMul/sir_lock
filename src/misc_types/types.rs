use{
    serde::{Serialize, Deserialize},
    
    std::{
        num::*, 
        
    
        ops::RangeInclusive,
       
    },
    crate::stats_methods::*,
    crate::grid::*,
    std::io::Write,
};


pub const DEFAULT_SYSTEM_SIZE: NonZeroUsize = unsafe{NonZeroUsize::new_unchecked(200)};
pub const DEFAULT_RECOVERY_PROB: f64 = 0.14;
pub const DEFAULT_GRAPH_SEED: u64 = 875629289;
pub const DEFAULT_SIR_SEED: u64 = 1489264107025;
pub const DEFAULT_INITIAL_INFECTED:usize = 5;
//pub const DEFAULT_VACCINE_SEED: u64 = 4896709264107025;
pub const DEFAULT_SAMPLES_PER_STEP: u64 = 5000;
pub const ONE: NonZeroUsize = unsafe{NonZeroUsize::new_unchecked(1)};
pub const DEFAULT_GRAPH:GraphType = GraphType::SmallWorld(0.1);
pub const DEFAULT_F_THRESHOLD: f64 = 0.0000001;
pub const DEFAULT_LAMBDA: f64 = 0.2;
pub const DEFAULT_MARKOV_SEED: u64 = 782063498562509862;
pub const DEFAULT_SWEEP_SIZE: NonZeroUsize = unsafe{NonZeroUsize::new_unchecked(2222)};
pub const DEFAULT_MARKOV_STEP_SIZE: usize = 100;
pub const DEFAULT_SAMPLES_SIMPLE_SAMPLE: usize = 10000;
pub const DEFAULT_RANDOM_LOCKDOWN_SEED: u64 = 123131315;
pub const DEFAULT_RANDOM_LOCKDOWN_FRAC: f64 = 0.6;

//pub type GenGraphSIR = net_ensembles::GenericGraph<crate::sir_model::sir_states::InfectionState, net_ensembles::graph::NodeContainer<crate::sir_model::sir_states::InfectionState>>;
#[derive(Clone, Serialize, Deserialize, Copy)]
pub struct LargeDeviationParam
{
    pub time_steps: NonZeroUsize,
    pub markov_seed: u64,
    pub initial_infected: usize
}



#[derive(Serialize, Deserialize, Clone, Debug)]
pub enum GraphType{
    SmallWorld(f64),
    // place holder for later graph types I might want to use
    // This makes sure I have to implement i
    Barabasi(usize,usize),
    Invalid
}

impl GraphType{
    pub fn name(&self) -> String
    {
        match self
        {
            Self::SmallWorld(p) => format!("SmallWorld{}", p),
            Self::Barabasi(q,r) => format!("Barabasi{}{}",q,r),
            Self::Invalid => unimplemented!()
        }
    }
}
#[derive(Serialize, Deserialize, Clone, Debug)]
pub enum RewlType{
    SmallWorld,
    Barabasi,
    None
}

#[derive(Clone, Copy, Serialize, Deserialize)]
pub struct ExchangeInfo
{
    pub index: usize,
    pub old_val: f64
}

#[derive(Clone, Copy, Serialize, Deserialize)]
pub enum MarkovStep
{
    RotateLeft,
    RotateRight,
    Transmission(ExchangeInfo),
    Recovery(ExchangeInfo),
    SwapTrans((usize, usize)),
    SwapRec((usize, usize)),
    MovePatientZero(usize,usize,bool),
    MovePatientZeroRandom(usize,usize)
}
#[derive(Serialize, Deserialize)]
pub struct LockdownMarkovMove{
    pub lockdown_index: usize,
    pub not_lockdown_index: usize
}
#[derive(Serialize, Deserialize)]
pub enum MarkovStepWithLocks{
    BaseMarkovStep(MarkovStep),
    LockdownStep(LockdownMarkovMove) 
}
impl From<MarkovStep> for MarkovStepWithLocks{
    fn from(other:MarkovStep) -> Self{
        Self::BaseMarkovStep(other)
    }
}
#[derive(Clone, Copy, Serialize, Deserialize, Default)]
pub struct StepTracker
{
    accepted: u64,
    rejected: u64
}

impl StepTracker
{
    pub fn new() -> Self 
    {
        Self::default()
    }

    pub fn accept(&mut self)
    {
        self.accepted += 1;
    }

    pub fn reject(&mut self)
    {
        self.rejected += 1;
        if self.rejected > self.total_steps(){
            panic!("{} {}", self.accepted, self.rejected);
        }
    }

    #[inline]
    pub fn total_steps(&self) -> u64
    {
        self.accepted + self.rejected
    }

    pub fn rejection_rate(self) -> f64
    {
        self.rejected as f64 / (self.total_steps()) as f64
    }

    pub fn acceptance_rate(self) -> f64
    {
        self.accepted as f64 / (self.total_steps()) as f64
    }

    pub fn write<W: Write>(&self, mut writer: W, name: &str)
    {
        let _ = writeln!(
            writer,
            "#{name} total_steps {} rejected {} acceptance_rate {}", 
            self.total_steps(),
            self.rejected,
            self.acceptance_rate()
        );
    }

    pub fn reset(&mut self)
    {
        self.accepted = 0;
        self.rejected = 0;
    }
}

#[derive(Clone, Serialize, Deserialize, Default)]
pub struct AcceptanceTracker
{
    lockdown: StepTracker,
    patient_move: StepTracker,
    patient_move_random: StepTracker,
    rotate_left: StepTracker,
    rotate_right: StepTracker,
    swap_trans_or_rec: StepTracker,
    base_move: StepTracker,
}

impl AcceptanceTracker{
    pub fn new() -> Self 
    {
        Self::default()
    }

    pub fn reset(&mut self)
    {
        self.lockdown.reset();
        self.patient_move.reset();
        self.patient_move_random.reset();
        self.rotate_left.reset();
        self.rotate_right.reset();
        self.swap_trans_or_rec.reset();
        self.base_move.reset(); 
    }

    pub fn write_stats<W: Write>(&self, mut writer: W)
    {
        self.lockdown.write(&mut writer, "LockdownMove");
        self.patient_move_random.write(&mut writer, "PatientEdgeMoveRandom");
        self.patient_move.write(&mut writer, "PatientMove");
        self.rotate_left.write(&mut writer, "RotateLeft");
        self.rotate_right.write(&mut writer, "RotateRight");
        self.swap_trans_or_rec.write(&mut writer, "SwapTransOrRec");
        self.base_move.write(writer, "BaseMove");
    }

    pub fn accept(&mut self, step: &MarkovStepWithLocks)
    {
        match step
        {
            MarkovStepWithLocks::LockdownStep(_) => {
                self.lockdown.accept()
            },
            MarkovStepWithLocks::BaseMarkovStep(base_step) => 
            {
                match base_step
                {
                    MarkovStep::RotateLeft => {
                        self.rotate_left.accept()
                    },
                    MarkovStep::RotateRight => {
                        self.rotate_right.accept()
                    },
                    MarkovStep::SwapRec(..) | MarkovStep::SwapTrans(..) => 
                    {
                        self.swap_trans_or_rec.accept()
                    },
                    MarkovStep::Transmission(..) | MarkovStep::Recovery(..) => {
                        self.base_move.accept()
                    },
                    MarkovStep::MovePatientZero(.., moved) => {
                        if *moved
                        {
                            self.patient_move.accept()
                        }
                    },
                    MarkovStep::MovePatientZeroRandom(..) => {
                        self.patient_move_random.accept()
                    }
                    
                }
            }
        }
    }

    pub fn reject(&mut self, step: &MarkovStepWithLocks)
    {
        match step
        {
            MarkovStepWithLocks::LockdownStep(_) => {
                self.lockdown.reject()
            },
            MarkovStepWithLocks::BaseMarkovStep(base_step) => 
            {
                match base_step
                {
                    MarkovStep::RotateLeft => {
                        self.rotate_left.reject()
                    },
                    MarkovStep::RotateRight => {
                        self.rotate_right.reject()
                    },
                    MarkovStep::SwapRec(..) | MarkovStep::SwapTrans(..) => 
                    {
                        self.swap_trans_or_rec.reject()
                    },
                    MarkovStep::Transmission(..) | MarkovStep::Recovery(..) => {
                        self.base_move.reject()
                    },
                    MarkovStep::MovePatientZero(.., moved) => {
                        if *moved{
                            self.patient_move.reject()
                        }
                    },
                    MarkovStep::MovePatientZeroRandom(..) => {
                        self.patient_move_random.reject()
                    }
                    
                }
            }
        }
    }
}





























#[derive(Serialize, Deserialize, Clone, Debug, Copy)]
pub enum MeasureType {
    // ever infected
    C,
    // max infected
    M,
}

impl MeasureType{
    pub fn name(self) -> &'static str
    {
        match self{
            Self::C => "C",
            Self::M => "M",
        }
    }

    pub fn is_c(self) -> bool 
    {
        matches!(self, Self::C)
    }
}
pub struct Measured
{
    pub var_m: MyVariance,
    pub var_c: MyVariance,
}



#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct F64RangeBuilder
{
    pub start: f64,
    pub end: f64,
    pub steps: NonZeroUsize
}

impl F64RangeBuilder{
    pub fn get_range(&self) -> GridRangeF64
    {
        GridRangeF64::new(
            self.start, 
            self.end, 
            self.steps.get()
        )
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct UsizeRangeBuilder
{
    pub start: usize,
    pub end: usize,
    pub steps: NonZeroUsize
}

impl UsizeRangeBuilder {
    pub fn range(&self) -> RangeInclusive<usize>
    {
        self.start..=self.end
    }

    pub fn iter(&self) -> impl Iterator<Item=usize>
    {
        self.range().step_by(self.steps.get())
    }

    pub fn step_size(&self) -> NonZeroUsize
    {
        self.steps
    }

    pub fn end(&self) -> usize
    {
        self.end
    }

    pub fn len(&self) -> usize
    {
        (self.end - self.start) / self.steps.get()
    }
    pub fn is_empty(&self) -> bool{
        self.end == 0
    }
    pub fn get_range(&self) -> GridRangeUsize{
        GridRangeUsize::new(
            self.start,self.end,self.steps.get())
    }
}

