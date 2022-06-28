use serde::{Serialize, Deserialize};

use {
    net_ensembles::*,
    super::*,
    rand_pcg::Pcg64,
    net_ensembles::{
        WithGraph, 
        GraphIteratorsMut,
        MeasurableGraphQuantities,
        SimpleSample,
        rand::SeedableRng,
        GraphIterators,
    },
};

pub type BarabasiEnsemble = BAensemble<InfectionState,Pcg64>;
//pub type SmallWorldGraph = BAGraph<InfectionState>;





#[derive(Clone, Serialize, Deserialize)]
pub struct BarabasiModel{
    pub ensemble: BarabasiEnsemble,
    pub lambda: f64,
    pub gamma: f64,
    pub n: usize,
}

impl BarabasiModel{
    pub fn ensemble(&self) -> &BarabasiEnsemble
    {
        &self.ensemble
    }
    
    pub fn ensemble_mut(&mut self) -> &mut BarabasiEnsemble
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
    pub fn infect_patient(&mut self,patient:usize){
        
        self.ensemble
            .contained_iter_mut()
            .for_each(|s| *s = InfectionState::Suspectible);
            //resets all states to S
        //println!("infecting");
        *self.ensemble_mut().at_mut(patient) = InfectionState::Infected
        //infects patient 0
        }
    /// Called C in the paper
    pub fn calculate_ever_infected(&self) -> usize
    {
        self.ensemble.contained_iter()
            .filter(|&v| v.is_or_was_infected())
            .count()
    }
}

impl From<BarabasiOptions> for BarabasiModel{
    fn from(param:BarabasiOptions) -> Self{
        let graph_rng = Pcg64::seed_from_u64(param.graph_seed);
        let mut ensemble = BarabasiEnsemble::new(param.system_size.get(),graph_rng,param.m,param.source_n);
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
            lambda:param.lambda,
            gamma:param.gamma,
            n: param.system_size.get()
        }

    }

}