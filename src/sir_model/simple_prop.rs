use std::ops::{Deref, DerefMut};
use net_ensembles::rand::prelude::Distribution;
use rand_pcg::Pcg64;
use net_ensembles::{GraphIteratorsMut, rand};
use net_ensembles::rand::SeedableRng;
use rand::distributions::Uniform;
use crate::lockdown_methods::*;
use super::*;
use net_ensembles::WithGraph;
//pub type Generic = net_ensembles::GenericGraph<sir_model::sir_states::InfectionState,net_ensembles::graph::NodeContainer::sir_model::sir_states::InfectionState>;

pub type GenGraphSIR = net_ensembles::GenericGraph<crate::sir_model::sir_states::InfectionState, net_ensembles::graph::NodeContainer<crate::sir_model::sir_states::InfectionState>>;

#[derive(Clone)]
pub struct SimpleSample{
    base_model: BaseSWModel,
    infected_list: Vec<usize>,
    new_infected_list:Vec<usize>,
    rng_type: Pcg64,
    recovered_list:Vec<usize>,
    suspectible_list:Vec<usize>
}
//deref and deref mut solve some problems of things being references/mutable. 
impl Deref for SimpleSample
{
    type Target = BaseSWModel;
    fn deref(&self) -> &Self::Target {
        &self.base_model
    }
}

impl DerefMut for SimpleSample{
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.base_model
    }
}

impl SimpleSample{
    pub fn reseed_sir_rng(&mut self, rng: &mut Pcg64)
    { //reseeding rng
        self.rng_type = Pcg64::from_rng(rng).unwrap();
    }
    pub fn from_base(
        base_model: BaseSWModel,
        sir_sample_seed: u64
    ) -> Self
    {
        let rng_type = Pcg64::seed_from_u64(sir_sample_seed);
        let mut res = Self{
            base_model,
            infected_list:Vec::new(),
            new_infected_list:Vec::new(),
            rng_type,
            recovered_list:Vec::new(),
            suspectible_list:Vec::new()};
        res.reset_simple_sample_sir_simulation();
        res
    }

    /// Note: vaccine list should not contain patiend zero or any of 
    /// its neighbors. This is not checked here.
    pub fn reset_simple_sample_sir_simulation(&mut self)
    {  
        let un = Uniform::new(0,self.base_model.ensemble.graph().vertex_count());  
        let index = un.sample(&mut self.rng_type);
        self.infect_patient(index);
        
        // now only node 0 is infected!
        self.infected_list.clear();
        self.infected_list.push(index);
        
    }

    fn propagate_one_time_step(&mut self){
        debug_assert!(self.new_infected_list.is_empty());

        let prob_dist = Uniform::new_inclusive(0.0,1.0);
        let lambda = self.lambda; 
        //let dist = Uniform::new_inclusive(0.0, 1.0);
        //ok, there's a lot going on here.
        //firstly, we're lopping over the infected patients;
        //secondly, we're looping over the neighbours to the infected patient who are susceptible.
        //so, we're isolating the self.base_model parameter, then getting the ensemble. THEN we are isolating those neighbours.
        //from the neighbours, we're filtering out those are actually suspectible and can therefore be infected.

        for &index in self.infected_list.iter(){
            for (n_index,neighbour) in self.base_model
                .ensemble
                .contained_iter_neighbors_mut_with_index(index)
                .filter(|(_,neighbour)| neighbour.sus_check()){
                
                //I had an error for ages here, but this import fixed it: use net_ensembles::rand::prelude::Distribution;
                let prob = prob_dist.sample(&mut self.rng_type);
                
                if prob < lambda{
                    //* dereferences neighbour
                    *neighbour = InfectionState::Infected;
                    self.new_infected_list.push(n_index);

                }
                    
                }
        }
        for i in (0..self.infected_list.len()).rev(){
            if prob_dist.sample(&mut self.rng_type) < self.gamma{

                let removed_index = self.infected_list.swap_remove(i);
                *self.ensemble.at_mut(removed_index) = InfectionState::Recovered;
            }
        }
        self.infected_list.append(&mut self.new_infected_list);
    }

    pub fn propagate_until_completion_max(&mut self) -> usize{
        self.reset_simple_sample_sir_simulation();
        let mut max_infected = self.infected_list.len();
        debug_assert_eq!(max_infected,1);
        loop{
            self.propagate_one_time_step();
            max_infected = max_infected.max(self.infected_list.len());
            if self.infected_list.is_empty(){
                break;
            }

        }
        max_infected

    }
    

    pub fn propagate_until_completion_time(&mut self) -> u32{
        self.reset_simple_sample_sir_simulation();
        let max_infected = self.infected_list.len();
        debug_assert_eq!(max_infected,1);
        //
        for i in 1..{
            self.propagate_one_time_step();
            if self.infected_list.is_empty(){
                return i
            }

        }
        unreachable!()}

    pub fn produce_time_data(&mut self) -> (Vec<usize>,Vec<usize>,Vec<usize>,usize){
        self.reset_simple_sample_sir_simulation();
        let mut rec:Vec<usize> = Vec::new();
        let mut sus:Vec<usize> = Vec::new();
        let mut inf:Vec<usize> = Vec::new();
        //let time:Vec<usize> = Vec::new();
        rec.push(self.recovered_list.len());
        sus.push(self.suspectible_list.len());
        inf.push(self.infected_list.len());
        let mut time:usize = 0;
        
        loop{
            time += 1;
            self.propagate_one_time_step();
            rec.push(self.recovered_list.len());
            sus.push(self.suspectible_list.len());
            inf.push(self.infected_list.len());
            if self.infected_list.is_empty(){
                break;
            }
            
            
            

        }



        return (sus,inf,rec,time)

    }
}

#[derive(Clone)]
pub struct SimpleSampleBarabasi{
    base_model: BarabasiModel,
    infected_list: Vec<usize>,
    new_infected_list:Vec<usize>,
    rng_type: Pcg64,
    recovered_list:Vec<usize>,
    suspectible_list:Vec<usize>
}
impl Deref for SimpleSampleBarabasi
{
    type Target = BarabasiModel;
    fn deref(&self) -> &Self::Target {
        &self.base_model
    }
}

impl DerefMut for SimpleSampleBarabasi{
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.base_model
    }
}
impl SimpleSampleBarabasi{
    pub fn reseed_sir_rng(&mut self, rng: &mut Pcg64)
    { //reseeding rng
        self.rng_type = Pcg64::from_rng(rng).unwrap();
    }
    pub fn from_base(
        base_model: BarabasiModel,
        sir_sample_seed: u64
    ) -> Self
    {
        let rng_type = Pcg64::seed_from_u64(sir_sample_seed);
        let mut res = Self{
            base_model,
            infected_list:Vec::new(),
            new_infected_list:Vec::new(),
            rng_type,
            recovered_list:Vec::new(),
            suspectible_list:Vec::new()};
        res.reset_simple_sample_sir_simulation();
        res
    }

    /// Note: vaccine list should not contain patiend zero or any of 
    /// its neighbors. This is not checked here.
    pub fn reset_simple_sample_sir_simulation(&mut self)
    {  
        let un = Uniform::new(0,self.base_model.ensemble.graph().vertex_count());  
        let index = un.sample(&mut self.rng_type);
        self.infect_patient(index);
        
        // now only node 0 is infected!
        self.infected_list.clear();
        self.infected_list.push(index);
        self.recovered_list.clear();
        
    }

    
    
    


    fn propagate_one_time_step(&mut self){
        debug_assert!(self.new_infected_list.is_empty());

        let prob_dist = Uniform::new_inclusive(0.0,1.0);
        let lambda = self.lambda; 
        //let dist = Uniform::new_inclusive(0.0, 1.0);
        //ok, there's a lot going on here.
        //firstly, we're lopping over the infected patients;
        //secondly, we're looping over the neighbours to the infected patient who are susceptible.
        //so, we're isolating the self.base_model parameter, then getting the ensemble. THEN we are isolating those neighbours.
        //from the neighbours, we're filtering out those are actually suspectible and can therefore be infected.

        for &index in self.infected_list.iter(){
            for (n_index,neighbour) in self.base_model
                .ensemble
                .contained_iter_neighbors_mut_with_index(index)
                .filter(|(_,neighbour)| neighbour.sus_check()){
                
                //I had an error for ages here, but this import fixed it: use net_ensembles::rand::prelude::Distribution;
                let prob = prob_dist.sample(&mut self.rng_type);
                
                if prob < lambda{
                    //* dereferences neighbour
                    *neighbour = InfectionState::Infected;
                    self.new_infected_list.push(n_index);

                }
                    
                }
        }
        for i in (0..self.infected_list.len()).rev(){
            if prob_dist.sample(&mut self.rng_type) < self.gamma{

                let removed_index = self.infected_list.swap_remove(i);
                *self.ensemble.at_mut(removed_index) = InfectionState::Recovered;
            }
        }
        self.infected_list.append(&mut self.new_infected_list);
    }

    pub fn propagate_until_completion_max(&mut self) -> usize{
        self.reset_simple_sample_sir_simulation();


        let mut max_infected = self.infected_list.len();

        debug_assert_eq!(max_infected,1);
        loop{
            self.propagate_one_time_step();

            max_infected = max_infected.max(self.infected_list.len());
            if self.infected_list.is_empty(){
                break;
            }

        }
        max_infected

    }
    
    //maybe think about having the lockeddown graph calculated once and then passed to this function!

    pub fn create_locked_down_network(&mut self, lockdownparams:LockdownParameters)
    -> GenGraphSIR{

        //this will only work right now during a simulation, not at the beginning because of the current implementation of reset_simple_sample.
        let graph = self.base_model.ensemble.graph().clone();
        //This transfers the SIR information to the new network also. //doesn't work? need to fix maybe
        let locked_down_graph = lockdown(lockdownparams,graph,&self.infected_list);
        

        locked_down_graph


    }

    pub fn transfer_sir_information(&mut self, mut locked_down_graph:GenGraphSIR) -> GenGraphSIR{
        //resets and updates the sir informaton of the locked_down_graph
        for state in locked_down_graph.contained_iter_mut(){
            *state = InfectionState::Suspectible;

        }

        for index in &self.infected_list{
            *locked_down_graph.at_mut(*index) = InfectionState::Infected;
        }
        for index in &self.recovered_list{
            *locked_down_graph.at_mut(*index) = InfectionState::Recovered;
        }
        locked_down_graph
        
    }


    pub fn propagate_until_completion_max_with_lockdown(&mut self,mut post_locked_down_graph:GenGraphSIR,lockparams:LockdownParameters) -> usize{
        
        //THIS IS WHY create_locked_down_network doesnt work rn, as this makes a new position for patient zero each time. And leaves room for duplicates.
        self.reset_simple_sample_sir_simulation();

        //maybe have post_locked_down_graph be part of the sample struct? would allow for easier updating and function passing? who knows
        //let (mut post_locked_down_graph,stat_bool) = self.create_locked_down_network(locktype);

        //NOTE: The act of cloning the graph, sets the SIR information to default in all of the nodes. It is necessary to update them.

        //wont need after change made
        post_locked_down_graph = self.transfer_sir_information(post_locked_down_graph);
        //println!("Old edges {}" ,self.base_model.ensemble.graph().edge_count());

        let lockdown_threshold = lockparams.lock_threshold;
        let release_threshold = lockparams.release_threshold;
        let dynamic_bool = lockparams.dynamic_bool;

        let mut max_infected = self.infected_list.len();
        debug_assert_eq!(max_infected,1);
        let mut lockdown_indicator = false;
        loop{
            debug_assert_eq!(max_infected,post_locked_down_graph.contained_iter().filter(|&state| *state == InfectionState::Infected).count());



            let prob_dist = Uniform::new_inclusive(0.0,1.0);
            let lambda = self.lambda; 
            let inf = self.infected_list.len() as f64/self.n as f64;
            //println!("{}",inf);
            //println!("{:?}",self.infected_list);

            //transfer SIR information when lockdown comes in/leaves

            if inf > lockdown_threshold && lockdown_indicator == false{
                lockdown_indicator = true;
                //post_locked_down_graph = self.transfer_sir_information(post_locked_down_graph);

                if dynamic_bool{
                    //If the lockdown type is one which must be updated once the new lockdown is brought in. 
                    //Otherwise the lockdown is static and the structure is constant throughout the propagation
                    post_locked_down_graph = self.create_locked_down_network(lockparams);
                    post_locked_down_graph = self.transfer_sir_information(post_locked_down_graph);
                    //println!("Locking down: new edges {}",post_locked_down_graph.edge_count());
                }

            }
            if inf < release_threshold &&lockdown_indicator == true{
                //println!("Releasing lockdown: edges {}",self.base_model.ensemble.graph().edge_count());
                lockdown_indicator = false;
            }

            let bool_dup_checker = false;
            if !bool_dup_checker{
                if contains_duplicates(self.infected_list.clone()){
                println!("safety dance");
                }
            }

            if !lockdown_indicator{
                

                for &index in self.infected_list.iter(){
                    for (n_index,neighbour) in self.base_model.ensemble.contained_iter_neighbors_mut_with_index(index).filter(|(_,neighbour)| neighbour.sus_check()){
                
                        let prob = prob_dist.sample(&mut self.rng_type);
                
                        if prob < lambda{
                            //* dereferences neighbour
                            //need to update the corresponding neighbour in the other graph!
                            *neighbour = InfectionState::Infected;
                            
                            //This might only be necessary if lockdown is static, as dynamic lockdown will calculate this anyway...
                            *post_locked_down_graph.at_mut(n_index) = InfectionState::Infected;
                            self.new_infected_list.push(n_index);} }
                }

            }
            else{
                //println!("locking down");
                

                for &index in self.infected_list.iter(){

                    for (n_index,neighbour) in post_locked_down_graph.contained_iter_neighbors_mut_with_index(index).filter(|(_,neighbour)| neighbour.sus_check()){
                        let prob = prob_dist.sample(&mut self.rng_type);
                
                        if prob < lambda{
                            //* dereferences neighbour
                            //need to update the corresponding neighbour in the other graph!
                            *neighbour = InfectionState::Infected;
                            *self.base_model.ensemble.at_mut(n_index) = InfectionState::Infected;
                            self.new_infected_list.push(n_index);} }
                }
            }

            //Recoveries are independent of the topology.
            for i in (0..self.infected_list.len()).rev(){
                if prob_dist.sample(&mut self.rng_type) < self.gamma{
        
                    let removed_index = self.infected_list.swap_remove(i);
                    *self.ensemble.at_mut(removed_index) = InfectionState::Recovered;
                    *post_locked_down_graph.at_mut(removed_index) = InfectionState::Recovered;
                    self.recovered_list.push(removed_index);
                        //println!("recovering");
                }
            }
            self.infected_list.append(&mut self.new_infected_list);
            self.new_infected_list = Vec::new();


            
            max_infected = max_infected.max(self.infected_list.len());

            if self.infected_list.is_empty(){
                //let inf = self.infected_list.len() as f64/self.n as f64;
                //println!("{}",inf);

                break
            }
            
        }
        max_infected

    }

    

    pub fn propagate_until_completion_time(&mut self) -> u32{
        self.reset_simple_sample_sir_simulation();
        let max_infected = self.infected_list.len();
        debug_assert_eq!(max_infected,1);
        //
        for i in 1..{
            self.propagate_one_time_step();
            if self.infected_list.is_empty(){
                return i
            }

        }
        unreachable!()}

    pub fn produce_time_data(&mut self) -> (Vec<usize>,Vec<usize>,Vec<usize>,usize){
        self.reset_simple_sample_sir_simulation();
        let mut rec:Vec<usize> = Vec::new();
        let mut sus:Vec<usize> = Vec::new();
        let mut inf:Vec<usize> = Vec::new();
        //let time:Vec<usize> = Vec::new();
        rec.push(self.recovered_list.len());
        sus.push(self.suspectible_list.len());
        inf.push(self.infected_list.len());
        let mut time:usize = 0;
        
        loop{
            time += 1;
            self.propagate_one_time_step();
            rec.push(self.recovered_list.len());
            sus.push(self.suspectible_list.len());
            inf.push(self.infected_list.len());
            if self.infected_list.is_empty(){
                break;
            }
            
            
            

        }



        (sus,inf,rec,time)

    }

}

pub fn contains_duplicates(x:Vec<usize>)-> bool{
    for i in 0..x.len(){
        for j in 0..x.len(){
            if i != j{
                if x[i] == x[j]{
                    return true
                }
            }
        }
    }
    return false
}