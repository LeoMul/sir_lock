use serde::{Serialize, Deserialize};

use {
    net_ensembles::SwEnsemble,
    super::*,
    rand_pcg::Pcg64,
    net_ensembles::{
        WithGraph, 
        GraphIteratorsMut,
        SwGraph,
        MeasurableGraphQuantities,
        SimpleSample,
        rand::SeedableRng,
        GraphIterators,
    },
};


pub type SmallWorldEnsemble = SwEnsemble<InfectionState,Pcg64>;
pub type SmallWorldGraph = SwGraph<InfectionState>;


#[derive(Clone, Serialize, Deserialize)]
pub struct SWModel{
    pub ensemble: SmallWorldEnsemble,
    pub lambda: f64,
    pub gamma: f64,
    pub n: usize,
}

impl SWModel{
    pub fn ensemble(&self) -> &SmallWorldEnsemble
    {
        &self.ensemble
    }
    
    pub fn ensemble_mut(&mut self) -> &mut SmallWorldEnsemble
    {
        &mut self.ensemble
    }

    pub fn set_lambda(&mut self, lambda: f64)
    {
        self.lambda = lambda;
    }

    pub fn set_gamma(&mut self, gamma: f64)
    {
        self.gamma = gamma;
    }
    pub fn set_all_to_sus(&mut self){
        //resets all states to S
        self.ensemble
            .contained_iter_mut()
            .for_each(|s| *s = InfectionState::Suspectible);

        
        debug_assert_eq!(self.n,self.sus_count())
    
    }
    pub fn sus_count(&self) -> usize{
        self.ensemble.contained_iter().filter(|s| s.sus_check()).count()
    }
    pub fn infect_patient(&mut self,patient:usize){
        
        //println!("infecting");
        *self.ensemble_mut().at_mut(patient) = InfectionState::Infected
        //infects patient 0
    }
    pub fn infect_many_patients(&mut self, vec:&[usize]){

        self.ensemble
            .contained_iter_mut()
            .for_each(|s| *s = InfectionState::Suspectible);

        for patient in vec{
            *self.ensemble_mut().at_mut(*patient) = InfectionState::Infected
        }
    }
    /// Called C in the paper
    pub fn calculate_ever_infected(&self) -> usize
    {
        self.ensemble.contained_iter()
            .filter(|&v| v.is_or_was_infected())
            .count()
    }
}
impl From<SWOptions> for SWModel
{
    fn from(param: SWOptions) -> Self {
        let graph_rng = Pcg64::seed_from_u64(param.graph_seed);
        
        let mut ensemble = SmallWorldEnsemble::new(
            param.system_size.get(),
            param.rewire_prob,
            graph_rng
        );
        let mut counter = 0_u32;
        while !ensemble.is_connected().unwrap()
        {
            counter += 1;
            ensemble.randomize();
        }
        if counter > 0 
        {
            println!("Randomized the ensemble {} additional times to find a connected network", counter);
        }

        Self{
            ensemble,
            lambda: param.lambda,
            gamma: param.gamma,
            n: param.system_size.get()
        }
    }
}
