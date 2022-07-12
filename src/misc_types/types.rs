use{
    serde::{Serialize, Deserialize},
    
    std::{
        num::*, 
        
    
        ops::RangeInclusive,
       
    },
    crate::stats_methods::*,
    crate::grid::*,
};
pub const DEFAULT_SYSTEM_SIZE: NonZeroUsize = unsafe{NonZeroUsize::new_unchecked(200)};
pub const DEFAULT_RECOVERY_PROB: f64 = 0.14;
pub const DEFAULT_GRAPH_SEED: u64 = 875629289;
pub const DEFAULT_SIR_SEED: u64 = 1489264107025;
//pub const DEFAULT_VACCINE_SEED: u64 = 4896709264107025;
pub const DEFAULT_SAMPLES_PER_STEP: u64 = 5000;
pub const ONE: NonZeroUsize = unsafe{NonZeroUsize::new_unchecked(1)};

pub const DEFAULT_F_THRESHOLD: f64 = 0.0000001;
pub const DEFAULT_LAMBDA: f64 = 0.1763;
pub const DEFAULT_MARKOV_SEED: u64 = 782063498562509862;
pub const DEFAULT_SWEEP_SIZE: NonZeroUsize = unsafe{NonZeroUsize::new_unchecked(2222)};
pub const DEFAULT_MARKOV_STEP_SIZE: usize = 100;
pub const DEFAULT_SAMPLES_SIMPLE_SAMPLE: usize = 10000;



//pub type GenGraphSIR = net_ensembles::GenericGraph<crate::sir_model::sir_states::InfectionState, net_ensembles::graph::NodeContainer<crate::sir_model::sir_states::InfectionState>>;

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
            Self::SmallWorld(p) => format!("sw{}", p),
            Self::Barabasi(q,r) => format!("ba{}{}",q,r),
            Self::Invalid => unimplemented!()
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

