use std::ops::{Deref, DerefMut};
use net_ensembles::rand::prelude::Distribution;
use net_ensembles::sw_graph::SwContainer;
use rand_pcg::Pcg64;
use net_ensembles::{GraphIteratorsMut, rand, GenericGraph};
use net_ensembles::rand::SeedableRng;
use rand::distributions::Uniform;
use crate::lockdown_methods::*;
use super::*;
use net_ensembles::WithGraph;
//pub type Generic = net_ensembles::GenericGraph<sir_model::sir_states::InfectionState,net_ensembles::graph::NodeContainer::sir_model::sir_states::InfectionState>;
//use crate::misc_types::*;
//pub type GenGraphSIR = net_ensembles::GenericGraph<crate::sir_model::sir_states::InfectionState, net_ensembles::graph::NodeContainer<crate::sir_model::sir_states::InfectionState>>;

#[derive(Clone)]
pub struct SimpleSampleSW{
    base_model: SWModel,
    infected_list: Vec<usize>,
    new_infected_list:Vec<usize>,
    rng_type: Pcg64,
    initial_infected:usize,
    recovered_list:Vec<usize>,
    suspectible_list:Vec<usize>,
    lockgraph:SwSIR
}
//deref and deref mut solve some problems of things being references/mutable. 
impl Deref for SimpleSampleSW
{
    type Target = SWModel;
    fn deref(&self) -> &Self::Target {
        &self.base_model
    }
}

impl DerefMut for SimpleSampleSW{
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.base_model
    }
}

impl SimpleSampleSW{
    pub fn reseed_sir_rng(&mut self, rng: &mut Pcg64)
    { //reseeding rng
        self.rng_type = Pcg64::from_rng(rng).unwrap();
    }
    pub fn from_base(
        base_model: SWModel,
        sir_sample_seed: u64, initial_infected:usize
    ) -> Self
    {
        let rng_type = Pcg64::seed_from_u64(sir_sample_seed);
        let graph = base_model.ensemble.graph().clone();

        //placeholder graph to be replaced by the lockdowngraph at each lockdown.
        //let graph = GenericGraph::<InfectionState,_>::new(1);

        //This transfers the SIR information to the new network also. //doesn't work? need to fix maybe
        
        let mut res = Self{
            base_model,
            infected_list:Vec::new(),
            new_infected_list:Vec::new(),
            rng_type,
            initial_infected,
            recovered_list:Vec::new(),
            suspectible_list:Vec::new(),
            lockgraph:graph};
        res.reset_simple_sample_sir_simulation();
        res
    }

    pub fn create_new_locked_down_network(&mut self,lockdownparams:LockdownParameters){
        let graph = self.base_model.ensemble.graph().clone();
        let mut lock_graph = lockdown(lockdownparams,graph,&mut self.rng_type);
        self.lockgraph = self.transfer_sir_information(&mut lock_graph).clone();

        //println!("{}",self.lockgraph.edge_count())


    }


    /// Note: vaccine list should not contain patiend zero or any of 
    /// its neighbors. This is not checked here.
    pub fn reset_simple_sample_sir_simulation(&mut self)
    {  
        self.infected_list.clear();
        self.set_all_to_sus();
        
        let un = Uniform::new(0,self.base_model.ensemble.graph().vertex_count());

        while self.infected_list.len() < self.initial_infected{
            let index = un.sample(&mut self.rng_type);
            if !self.infected_list.iter().any(|i| *i == index){
                self.infect_patient(index);
                self.infected_list.push(index);

            }
            

        } 
        
        
    }
    pub fn reset_simple_sample_sir_simulation_legacy(&mut self){  
        let un = Uniform::new(0,self.base_model.ensemble.graph().vertex_count());  
        let index = un.sample(&mut self.rng_type);
        self.set_all_to_sus();
        self.infect_patient(index);
        
        // now only node 0 is infected!
        self.infected_list.clear();
        self.infected_list.push(index);
        //self.recovered_list.clear();
        
        
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


        pub fn create_locked_down_network(&mut self, lockdownparams:LockdownParameters)
        -> SwSIR{

        //this will only work right now during a simulation, not at the beginning because of the current implementation of reset_simple_sample.
        let graph = self.base_model.ensemble.graph().clone();
        //This transfers the SIR information to the new network also. //doesn't work? need to fix maybe
        lockdown(lockdownparams,graph,&mut self.rng_type)

        }

    pub fn transfer_sir_information<'a>(&mut self, locked_down_graph: &'a mut SwSIR) ->  &'a mut SwSIR{
        //resets and updates the sir informaton of the locked_down_graph
        

        self.ensemble.graph().contained_iter().zip(locked_down_graph.contained_iter_mut()).for_each(
            |(old,new)|{ *new = *old}
        );

        locked_down_graph
        
    }
    

    pub fn iterate_once_with_locks(&mut self,lockdown_indicator:bool){
        let prob_dist = Uniform::new_inclusive(0.0,1.0);
        let lambda = self.lambda; 
        if !lockdown_indicator{
            //println!("test");

            for &index in self.infected_list.iter(){
                for (n_index,neighbour) in self.base_model.ensemble.contained_iter_neighbors_mut_with_index(index).filter(|(_,neighbour)| neighbour.sus_check()){
            
                    let prob = prob_dist.sample(&mut self.rng_type);
            
                    if prob < lambda{
                        //* dereferences neighbour

                        //need to update the corresponding neighbour in the other graph!
                        *neighbour = InfectionState::Infected;
                        
                        //This might only be necessary if lockdown is static, as dynamic lockdown will calculate this anyway...
                        
                        //why dont we just transfer all the information between the graphs each time we switch, as in the current implementation
                        //the generate lockdown graph self.function does precisely this, making the next line not useful and we can comment it. hopefully.
                        *self.lockgraph.at_mut(n_index) = InfectionState::Infected;
                        self.new_infected_list.push(n_index);} }
            }

        }
        else{
            //println!("locking down");
            

            for &index in self.infected_list.iter(){

                for (n_index,neighbour) in self.lockgraph.contained_iter_neighbors_mut_with_index(index).filter(|(_,neighbour)| neighbour.sus_check()){
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
                *self.lockgraph.at_mut(removed_index) = InfectionState::Recovered;
                //self.recovered_list.push(removed_index);
                    //println!("recovering");
            }
        }
        self.infected_list.append(&mut self.new_infected_list);
        self.new_infected_list = Vec::new();

    }

    

    pub fn propagate_until_completion_max_with_lockdown(&mut self,lockparams:LockdownParameters) -> usize{
        
        //this fn makes use of the iterate once function and is identical to the others. Its purpose is mostly so I can figure out how to do the large deviation one.

        //THIS IS WHY create_locked_down_network doesnt work rn, as this makes a new position for patient zero each time. And leaves room for duplicates.
        self.reset_simple_sample_sir_simulation();

        //maybe have post_locked_down_graph be part of the sample struct? would allow for easier updating and function passing? who knows
        //let (mut post_locked_down_graph,stat_bool) = self.create_locked_down_network(locktype);

        //NOTE: The act of cloning the graph, sets the SIR information to default in all of the nodes. It is necessary to update them.
        //wont need after change made
        //let mut post_locked_down_graph = &mut self.create_locked_down_network(lockparams);
        //let post_locked_down_graph = &mut self.transfer_sir_information(post_locked_down_graph);

        self.create_new_locked_down_network(lockparams);
        //println!("Old edges {}" ,self.base_model.ensemble.graph().edge_count());

        let lockdown_threshold = lockparams.lock_threshold;
        //println!("{}",lockdown_threshold);
        let release_threshold = lockparams.release_threshold;

        let mut max_infected = self.infected_list.len();
        debug_assert_eq!(max_infected,self.initial_infected);
        let mut lockdown_indicator = false;
        loop{
            //debug_assert_eq!(max_infected,post_locked_down_graph.contained_iter().filter(|&state| *state == InfectionState::Infected).count());
            //let sus = self.sus_count();
            //let infected_integer = self.infected_list.len();
            //println!("sus {sus} inf {infected_integer}");


            //let prob_dist = Uniform::new_inclusive(0.0,1.0);
            //let lambda = self.lambda; 
            let inf = self.infected_list.len() as f64/self.n as f64;
            //println!("{}",inf);
            //println!("{:?}",self.infected_list);

            //transfer SIR information when lockdown comes in/leaves
            //let mut new_locked_down_infected:Vec<usize> = Vec::new();
            if inf > lockdown_threshold && !lockdown_indicator{
                lockdown_indicator = true;
                //post_locked_down_graph = self.transfer_sir_information(post_locked_down_graph);
                //println!("Engaging lockdown: edges {}",post_locked_down_graph.edge_count());

            }
            if inf < release_threshold &&lockdown_indicator{
                //println!("Releasing lockdown: edges {}",self.base_model.ensemble.graph().edge_count());
                lockdown_indicator = false;
            }

            //let bool_dup_checker = false;
            //if bool_dup_checker &&contains_duplicates(self.infected_list.clone()){
            //    
            //    println!("safety dance");
            //    
            //}

            self.iterate_once_with_locks(lockdown_indicator);
            max_infected = max_infected.max(self.infected_list.len());

            if self.infected_list.is_empty(){
                //let inf = self.infected_list.len() as f64/self.n as f64;
                //println!("{}",inf);

                break
            }
            
        }
        max_infected

    }
    pub fn propagate_until_completion_max_with_lockdown_single(&mut self,lockparams:LockdownParameters,mut lockgraph:SwSIR) -> usize{
        
        //this fn makes use of the iterate once function and is identical to the others. Its purpose is mostly so I can figure out how to do the large deviation one.

        //THIS IS WHY create_locked_down_network doesnt work rn, as this makes a new position for patient zero each time. And leaves room for duplicates.
        self.reset_simple_sample_sir_simulation();

        //maybe have post_locked_down_graph be part of the sample struct? would allow for easier updating and function passing? who knows
        //let (mut post_locked_down_graph,stat_bool) = self.create_locked_down_network(locktype);
    
        //NOTE: The act of cloning the graph, sets the SIR information to default in all of the nodes. It is necessary to update them.
        //wont need after change made
        //let mut post_locked_down_graph = &mut self.create_locked_down_network(lockparams);
        self.lockgraph = self.transfer_sir_information(&mut lockgraph).clone();

        //self.create_new_locked_down_network(lockparams);
        //println!("Old edges {}" ,self.base_model.ensemble.graph().edge_count());

        let lockdown_threshold = lockparams.lock_threshold;
        //println!("{}",lockdown_threshold);
        let release_threshold = lockparams.release_threshold;

        let mut max_infected = self.infected_list.len();
        debug_assert_eq!(max_infected,self.initial_infected);
        let mut lockdown_indicator = false;
        loop{
            //debug_assert_eq!(max_infected,post_locked_down_graph.contained_iter().filter(|&state| *state == InfectionState::Infected).count());
            //let sus = self.sus_count();
            //let infected_integer = self.infected_list.len();
            //println!("sus {sus} inf {infected_integer}");


            //let prob_dist = Uniform::new_inclusive(0.0,1.0);
            //let lambda = self.lambda; 
            let inf = self.infected_list.len() as f64/self.n as f64;
            //println!("{}",inf);
            //println!("{:?}",self.infected_list);

            //transfer SIR information when lockdown comes in/leaves
            //let mut new_locked_down_infected:Vec<usize> = Vec::new();
            if inf > lockdown_threshold && !lockdown_indicator{
                lockdown_indicator = true;
                //post_locked_down_graph = self.transfer_sir_information(post_locked_down_graph);
                //println!("Engaging lockdown: edges {}",post_locked_down_graph.edge_count());

            }
            if inf < release_threshold &&lockdown_indicator{
                //println!("Releasing lockdown: edges {}",self.base_model.ensemble.graph().edge_count());
                lockdown_indicator = false;
            }

            //let bool_dup_checker = false;
            //if bool_dup_checker &&contains_duplicates(self.infected_list.clone()){
            //    
            //    println!("safety dance");
            //    
            //}

            self.iterate_once_with_locks(lockdown_indicator);
            max_infected = max_infected.max(self.infected_list.len());

            if self.infected_list.is_empty(){
                //let inf = self.infected_list.len() as f64/self.n as f64;
                //println!("{}",inf);

                break
            }
            
        }
        max_infected

    }
    pub fn propagate_until_completion_max_with_lockdown_printing(&mut self,post_locked_down_graph:&mut SwSIR,lockparams:LockdownParameters,writer:&mut SirWriter){
        
        //this fn makes use of the iterate once function and is identical to the others. Its purpose is mostly so I can figure out how to do the large deviation one.

        //THIS IS WHY create_locked_down_network doesnt work rn, as this makes a new position for patient zero each time. And leaves room for duplicates.
        self.reset_simple_sample_sir_simulation();
        writer.write_current(self.ensemble().graph())
            .unwrap();
        //maybe have post_locked_down_graph be part of the sample struct? would allow for easier updating and function passing? who knows
        //let (mut post_locked_down_graph,stat_bool) = self.create_locked_down_network(locktype);

        //NOTE: The act of cloning the graph, sets the SIR information to default in all of the nodes. It is necessary to update them.

        //wont need after change made
        let post_locked_down_graph = &mut self.transfer_sir_information(post_locked_down_graph);
        //println!("Old edges {}" ,self.base_model.ensemble.graph().edge_count());

        let lockdown_threshold = lockparams.lock_threshold;
        //println!("{}",lockdown_threshold);
        let release_threshold = lockparams.release_threshold;

        let mut max_infected = self.infected_list.len();
        debug_assert_eq!(max_infected,self.initial_infected);
        let mut lockdown_indicator = false;
        loop{
            //debug_assert_eq!(max_infected,post_locked_down_graph.contained_iter().filter(|&state| *state == InfectionState::Infected).count());
            //let sus = self.sus_count();
            //let infected_integer = self.infected_list.len();
            //println!("sus {sus} inf {infected_integer}");


            //let prob_dist = Uniform::new_inclusive(0.0,1.0);
            //let lambda = self.lambda; 
            let inf = self.infected_list.len() as f64/self.n as f64;
            //println!("{}",inf);
            //println!("{:?}",self.infected_list);

            //transfer SIR information when lockdown comes in/leaves
            //let mut new_locked_down_infected:Vec<usize> = Vec::new();
            if inf > lockdown_threshold && !lockdown_indicator{
                lockdown_indicator = true;
                //post_locked_down_graph = self.transfer_sir_information(post_locked_down_graph);
                //println!("Engaging lockdown: edges {}",post_locked_down_graph.edge_count());

            }
            if inf < release_threshold &&lockdown_indicator{
                //println!("Releasing lockdown: edges {}",self.base_model.ensemble.graph().edge_count());
                lockdown_indicator = false;
            }

            //let bool_dup_checker = false;
            //if bool_dup_checker &&contains_duplicates(self.infected_list.clone()){
            //    
            //    println!("safety dance");
            //    
            //}

            self.iterate_once_with_locks(lockdown_indicator);
            writer.write_current(self.ensemble().graph())
            .unwrap();
            max_infected = max_infected.max(self.infected_list.len());

            if self.infected_list.is_empty(){
                //let inf = self.infected_list.len() as f64/self.n as f64;
                //println!("{}",inf);

                break
            }
            
        }
        writer.write_line().unwrap();

    }

    pub fn propagate_until_completion_time_with_locks(&mut self,post_locked_down_graph:&mut SwSIR,lockparams:LockdownParameters) -> u32{
        self.reset_simple_sample_sir_simulation();
        let max_infected = self.infected_list.len();
        debug_assert_eq!(max_infected,1);
        let post_locked_down_graph = self.transfer_sir_information(post_locked_down_graph);
        let lockdown_threshold = lockparams.lock_threshold;
        //println!("{}",lockdown_threshold);
        let release_threshold = lockparams.release_threshold;
        let mut lockdown_indicator = false;
        //let prob_dist = Uniform::new_inclusive(0.0,1.0);
        //let lambda = self.lambda; 
        for i in 1..{
            debug_assert_eq!(max_infected,post_locked_down_graph.contained_iter().filter(|&state| *state == InfectionState::Infected).count());
            let inf = self.infected_list.len() as f64/self.n as f64;
            if inf > lockdown_threshold && !lockdown_indicator{
                lockdown_indicator = true;


            }
            if inf < release_threshold &&lockdown_indicator{
                //println!("Releasing lockdown: edges {}",self.base_model.ensemble.graph().edge_count());
                lockdown_indicator = false;
            }
            self.iterate_once_with_locks(lockdown_indicator);
            if self.infected_list.is_empty(){
                return i
            }

        }
        unreachable!()
    }
    
    pub fn propagate_until_completion_time_with_locks_new_lockgraph_for_each_lockdown(&mut self,lockparams:LockdownParameters) -> u32{
        self.reset_simple_sample_sir_simulation();

        let max_infected = self.infected_list.len();
        debug_assert_eq!(max_infected,self.initial_infected);

        //I feel as though this is perhaps not the correct way to do this... 

        //let mut post_locked_down_graph = self.create_locked_down_network(lockparams);
        //let post_locked_down_graph = self.transfer_sir_information(&mut post_locked_down_graph);
        //

        let lockdown_threshold = lockparams.lock_threshold;
        //println!("{}",lockdown_threshold);
        let release_threshold = lockparams.release_threshold;
        let mut lockdown_indicator = false;
        //let prob_dist = Uniform::new_inclusive(0.0,1.0);
        //let lambda = self.lambda; 
        for i in 1..{
            
            let inf = self.infected_list.len() as f64/self.n as f64;
            if inf > lockdown_threshold && !lockdown_indicator{
                lockdown_indicator = true;
                
                self.create_new_locked_down_network(lockparams);

                //checking that the infected were transferred correctly:
                debug_assert_eq!(self.infected_list.len(),self.lockgraph.contained_iter().filter(|&state| *state == InfectionState::Infected).count());

                let test = self.lockgraph.degree(0);
                println!("{}",test.unwrap());
            }
            if inf < release_threshold &&lockdown_indicator{
                //println!("Releasing lockdown: edges {}",self.base_model.ensemble.graph().edge_count());
                lockdown_indicator = false;
            }

            self.iterate_once_with_locks(lockdown_indicator);

            if self.infected_list.is_empty(){
                return i
            }

        }
        unreachable!()
    }

    pub fn propagate_until_completion_max_with_locks_new_lockgraph_for_each_lockdown(&mut self,lockparams:LockdownParameters) -> usize{
        
        //this fn makes use of the iterate once function and is identical to the others. Its purpose is mostly so I can figure out how to do the large deviation one.

        //THIS IS WHY create_locked_down_network doesnt work rn, as this makes a new position for patient zero each time. And leaves room for duplicates.
        self.reset_simple_sample_sir_simulation();

        //maybe have post_locked_down_graph be part of the sample struct? would allow for easier updating and function passing? who knows

        //past leo's calls for help have been answered. i am doing this now but i should have done it a year ago.
        // i hate myself.


        //let (mut post_locked_down_graph,stat_bool) = self.create_locked_down_network(locktype);


        let lockdown_threshold = lockparams.lock_threshold;
        //println!("{}",lockdown_threshold);
        let release_threshold = lockparams.release_threshold;

        let mut max_infected = self.infected_list.len();
        debug_assert_eq!(max_infected,self.initial_infected);
        let mut lockdown_indicator = false;

        loop{
            //debug_assert_eq!(max_infected,post_locked_down_graph.contained_iter().filter(|&state| *state == InfectionState::Infected).count());
        
            let inf = self.infected_list.len() as f64/self.n as f64;
            //println!("{}",inf);
            //println!("{:?}",self.infected_list);

            //transfer SIR information when lockdown comes in/leaves
            //let mut new_locked_down_infected:Vec<usize> = Vec::new();
            if inf > lockdown_threshold && !lockdown_indicator{
                lockdown_indicator = true;
                //creates new lokcdown network as required by the referee.
                self.create_new_locked_down_network(lockparams);

                //checking that the infected were transferred correctly:
                debug_assert_eq!(self.infected_list.len(),self.lockgraph.contained_iter().filter(|&state| *state == InfectionState::Infected).count());
                //let test = self.lockgraph.degree(1000);
                //let edges = self.lockgraph.edge_count();
                //println!("Degree of node zero {}, edge count {}",test.unwrap(),edges);
            }
            if inf < release_threshold &&lockdown_indicator{
                //println!("Releasing lockdown: edges {}",self.base_model.ensemble.graph().edge_count());
                lockdown_indicator = false;
            }

            //let bool_dup_checker = false;
            //if bool_dup_checker &&contains_duplicates(self.infected_list.clone()){
            //    
            //    println!("safety dance");
            //    
            //}

            self.iterate_once_with_locks(lockdown_indicator);
            max_infected = max_infected.max(self.infected_list.len());

            if self.infected_list.is_empty(){
                //let inf = self.infected_list.len() as f64/self.n as f64;
                //println!("{}",inf);

                break
            }
            
        }
        max_infected

    }

    
    
    
    
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

#[derive(Clone)]
pub struct SimpleSampleBarabasi{
    base_model: BarabasiModel,
    infected_list: Vec<usize>,
    new_infected_list:Vec<usize>,
    rng_type: Pcg64,
    initial_infected:usize
    
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
        sir_sample_seed: u64, initial_infected:usize
    ) -> Self
    {
        let rng_type = Pcg64::seed_from_u64(sir_sample_seed);
        let mut res = Self{
            base_model,
            infected_list:Vec::new(),
            new_infected_list:Vec::new(),
            rng_type,
            initial_infected
            };
        res.reset_simple_sample_sir_simulation();
        res
    }

    
    pub fn reset_simple_sample_sir_simulation(&mut self)
    {  
        self.infected_list.clear();
        self.set_all_to_sus();
        let un = Uniform::new(0,self.base_model.ensemble.graph().vertex_count());

        while self.infected_list.len() < self.initial_infected{
            let index = un.sample(&mut self.rng_type);
            if !self.infected_list.iter().any(|i| *i == index){
                self.infect_patient(index);
                self.infected_list.push(index);

            }
            

        } 
        
        
    }
    pub fn reset_simple_sample_sir_simulation_legacy(&mut self){  
        let un = Uniform::new(0,self.base_model.ensemble.graph().vertex_count());  
        let index = un.sample(&mut self.rng_type);
        self.set_all_to_sus();
        self.infect_patient(index);
        
        // now only node 0 is infected!
        self.infected_list.clear();
        self.infected_list.push(index);
        //self.recovered_list.clear();
        
        
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
        lockdown(lockdownparams,graph,&mut self.rng_type)
        

        


    }

    pub fn transfer_sir_information<'a>(&mut self, locked_down_graph: &'a mut GenGraphSIR) ->  &'a mut GenGraphSIR{
        //resets and updates the sir informaton of the locked_down_graph
        

        self.ensemble.graph().contained_iter().zip(locked_down_graph.contained_iter_mut()).for_each(
            |(old,new)|{ *new = *old}
        );

        locked_down_graph
        
    }
    

    pub fn iterate_once_with_locks(&mut self, post_locked_down_graph:&mut GenGraphSIR,lockdown_indicator:bool){
        let prob_dist = Uniform::new_inclusive(0.0,1.0);
        let lambda = self.lambda; 
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
                //self.recovered_list.push(removed_index);
                    //println!("recovering");
            }
        }
        self.infected_list.append(&mut self.new_infected_list);
        self.new_infected_list = Vec::new();

    }

    

    pub fn propagate_until_completion_max_with_lockdown(&mut self,post_locked_down_graph:&mut GenGraphSIR,lockparams:LockdownParameters) -> usize{
        
        //this fn makes use of the iterate once function and is identical to the others. Its purpose is mostly so I can figure out how to do the large deviation one.

        //THIS IS WHY create_locked_down_network doesnt work rn, as this makes a new position for patient zero each time. And leaves room for duplicates.
        self.reset_simple_sample_sir_simulation();

        //maybe have post_locked_down_graph be part of the sample struct? would allow for easier updating and function passing? who knows
        //let (mut post_locked_down_graph,stat_bool) = self.create_locked_down_network(locktype);

        //NOTE: The act of cloning the graph, sets the SIR information to default in all of the nodes. It is necessary to update them.

        //wont need after change made
        let post_locked_down_graph = &mut self.transfer_sir_information(post_locked_down_graph);
        //println!("Old edges {}" ,self.base_model.ensemble.graph().edge_count());

        let lockdown_threshold = lockparams.lock_threshold;
        //println!("{}",lockdown_threshold);
        let release_threshold = lockparams.release_threshold;
        

        let mut max_infected = self.infected_list.len();
        debug_assert_eq!(max_infected,1);
        let mut lockdown_indicator = false;
        loop{
            //debug_assert_eq!(max_infected,post_locked_down_graph.contained_iter().filter(|&state| *state == InfectionState::Infected).count());



            //let prob_dist = Uniform::new_inclusive(0.0,1.0);
            //let lambda = self.lambda; 
            let inf = self.infected_list.len() as f64/self.n as f64;
            //println!("{}",inf);
            //println!("{:?}",self.infected_list);
            //let sus = self.sus_count();
            //let infected_integer = self.infected_list.len();
            //println!("sus {sus} inf {infected_integer}");
            //transfer SIR information when lockdown comes in/leaves
            //let mut new_locked_down_infected:Vec<usize> = Vec::new();
            if inf > lockdown_threshold && !lockdown_indicator{
                lockdown_indicator = true;
                //post_locked_down_graph = self.transfer_sir_information(post_locked_down_graph);



            }
            if inf < release_threshold &&lockdown_indicator{
                //println!("Releasing lockdown: edges {}",self.base_model.ensemble.graph().edge_count());
                lockdown_indicator = false;
            }

            let bool_dup_checker = false;
            if bool_dup_checker &&contains_duplicates(self.infected_list.clone()){
                
                println!("safety dance");
                
            }

            self.iterate_once_with_locks(post_locked_down_graph, lockdown_indicator);
            max_infected = max_infected.max(self.infected_list.len());

            if self.infected_list.is_empty(){
                //let inf = self.infected_list.len() as f64/self.n as f64;
                //println!("{}",inf);

                break
            }
            
        }
        max_infected

    }
    pub fn propagate_until_completion_max_with_lockdown_printing(&mut self,post_locked_down_graph:&mut GenGraphSIR,lockparams:LockdownParameters,writer:&mut SirWriter){
        
        //this fn makes use of the iterate once function and is identical to the others. Its purpose is mostly so I can figure out how to do the large deviation one.

        //THIS IS WHY create_locked_down_network doesnt work rn, as this makes a new position for patient zero each time. And leaves room for duplicates.
        self.reset_simple_sample_sir_simulation();
        writer.write_current(self.ensemble().graph())
            .unwrap();
        //maybe have post_locked_down_graph be part of the sample struct? would allow for easier updating and function passing? who knows
        //let (mut post_locked_down_graph,stat_bool) = self.create_locked_down_network(locktype);

        //NOTE: The act of cloning the graph, sets the SIR information to default in all of the nodes. It is necessary to update them.

        //wont need after change made
        let post_locked_down_graph = &mut self.transfer_sir_information(post_locked_down_graph);
        //println!("Old edges {}" ,self.base_model.ensemble.graph().edge_count());

        let lockdown_threshold = lockparams.lock_threshold;
        //println!("{}",lockdown_threshold);
        let release_threshold = lockparams.release_threshold;

        let mut max_infected = self.infected_list.len();
        debug_assert_eq!(max_infected,self.initial_infected);
        let mut lockdown_indicator = false;
        loop{
            //debug_assert_eq!(max_infected,post_locked_down_graph.contained_iter().filter(|&state| *state == InfectionState::Infected).count());
            //let sus = self.sus_count();
            //let infected_integer = self.infected_list.len();
            //println!("sus {sus} inf {infected_integer}");


            //let prob_dist = Uniform::new_inclusive(0.0,1.0);
            //let lambda = self.lambda; 
            let inf = self.infected_list.len() as f64/self.n as f64;
            //println!("{}",inf);
            //println!("{:?}",self.infected_list);

            //transfer SIR information when lockdown comes in/leaves
            //let mut new_locked_down_infected:Vec<usize> = Vec::new();
            if inf > lockdown_threshold && !lockdown_indicator{
                lockdown_indicator = true;
                //post_locked_down_graph = self.transfer_sir_information(post_locked_down_graph);
                //println!("Engaging lockdown: edges {}",post_locked_down_graph.edge_count());

            }
            if inf < release_threshold &&lockdown_indicator{
                //println!("Releasing lockdown: edges {}",self.base_model.ensemble.graph().edge_count());
                lockdown_indicator = false;
            }

            //let bool_dup_checker = false;
            //if bool_dup_checker &&contains_duplicates(self.infected_list.clone()){
            //    
            //    println!("safety dance");
            //    
            //}

            self.iterate_once_with_locks(post_locked_down_graph, lockdown_indicator);
            writer.write_current(self.ensemble().graph())
            .unwrap();
            max_infected = max_infected.max(self.infected_list.len());

            if self.infected_list.is_empty(){
                //let inf = self.infected_list.len() as f64/self.n as f64;
                //println!("{}",inf);

                break
            }
            
        }
        writer.write_line().unwrap();

    }
    pub fn propagate_until_completion_time_with_locks(&mut self,mut post_locked_down_graph:GenGraphSIR,lockparams:LockdownParameters) -> u32{
        
        self.reset_simple_sample_sir_simulation();
        let max_infected = self.infected_list.len();
        debug_assert_eq!(max_infected,1);
        let post_locked_down_graph = self.transfer_sir_information(&mut post_locked_down_graph);
        let lockdown_threshold = lockparams.lock_threshold;
        //println!("{}",lockdown_threshold);
        let release_threshold = lockparams.release_threshold;
        let mut lockdown_indicator = false;
        //let prob_dist = Uniform::new_inclusive(0.0,1.0);
        //let lambda = self.lambda; 
        for i in 1..{
            debug_assert_eq!(max_infected,post_locked_down_graph.contained_iter().filter(|&state| *state == InfectionState::Infected).count());
            let inf = self.infected_list.len() as f64/self.n as f64;
            if inf > lockdown_threshold && !lockdown_indicator{
                lockdown_indicator = true;
                //post_locked_down_graph = self.transfer_sir_information(post_locked_down_graph);

            }
            if inf < release_threshold &&lockdown_indicator{
                //println!("Releasing lockdown: edges {}",self.base_model.ensemble.graph().edge_count());
                lockdown_indicator = false;
            }
            self.iterate_once_with_locks(post_locked_down_graph, lockdown_indicator);
            if self.infected_list.is_empty(){
                return i
            }

        }
        unreachable!()
    }

    
    
    

}

pub fn contains_duplicates(x:Vec<usize>)-> bool{
    for i in 0..x.len(){
        for j in 0..x.len(){
            if i !=j && x[i] == x[j]{
                return true
            }
        }
    }
    false
}