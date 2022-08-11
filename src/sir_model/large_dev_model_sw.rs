use {
    std::{num::*, ops::{Deref, DerefMut}},
    serde::{Serialize, Deserialize},
    net_ensembles::rand::{distributions::{Uniform, Distribution}, SeedableRng},
    net_ensembles::{GraphIterators, MarkovChain, MeasurableGraphQuantities, WithGraph},
    rand_pcg::Pcg64,
    net_ensembles::sampling::{HasRng, histogram::*},
    rand_distr::Binomial,
    crate::lockdown_methods::*,
    crate::sir_model::small_world::*,
    crate::sir_model::sir_states::*,
    crate::sir_model::offset::*,
    crate::sir_model::sir_writer::*,
    crate::misc_types::*,
};
use rand::Rng;
//might need but leaving it here.
//use net_ensembles::Contained,
const LOCKDOWN_CHANGE: f64 = 0.01;
const ROTATE_LEFT: f64 =  0.005;
const ROTATE_RIGHT: f64 =  0.01;
const PATIENT_MOVE: f64 = 0.05;
const PATIENT_MOVE_RANDOM: f64 = 0.08;

use net_ensembles::GraphIteratorsMut;
#[derive(Clone, Serialize, Deserialize)]
pub enum LockdownIndicator{
    Lockdown,
    NotLockdown
}

impl LockdownIndicator{

    pub fn is_in_lockdown(&self) -> bool{
        matches!(self,LockdownIndicator::Lockdown)
    }
    pub fn is_not_in_lockdown(&self) -> bool{
        matches!(self,LockdownIndicator::NotLockdown)
    }

}

#[derive(Clone, Serialize, Deserialize)]
pub struct SWLargeDeviation
{
    base_model: SWModel,
    lock_graph: SwSIR,
    #[cfg(feature = "ldprint")]   
    pub printingvector: Vec<usize>,
    not_lockdown_pairs: Vec<[usize;2]>,
    lockdown_pairs: Vec<[usize;2]>,
    lockdownparams:LockdownParameters, 
    lockdown_indicator: LockdownIndicator,
    transmission_rand_vec: Vec<f64>,
    pub markov_changed: bool,
    recovery_rand_vec: Vec<f64>,
    time_steps: NonZeroUsize,
    unfinished_simulations_counter: u64,
    total_simulations_counter: u64,
    markov_rng: Pcg64,
    offset: Offset,
    pub energy: u32,
    pub old_energy:u32,
    // temporary storage, allocated here so that 
    // I do not need to allocate it again and again
    new_infected: Vec<usize>,
    // Number of infected neighbors of a node
    // but only if the node is susceptible,
    // otherwise this is 0, as there is no danger
    dangerous_neighbor_count: Vec<i32>,
    one_minus_lambda: f64,
    pub system_size: NonZeroUsize,
    /// counting how many nodes are currently infected
    currently_infected_count: u32,
    zero_starting_conditions: bool,
    last_extinction_index: usize,
    pub patient_zero_vec: Vec<usize>,
    max_degree: usize,
    hist_patient_zero: HistUsizeFast
}
impl Deref for SWLargeDeviation
{
    type Target = SWModel;
    fn deref(&self) -> &Self::Target
    {
        &self.base_model
    }
}

impl DerefMut for SWLargeDeviation
{
    fn deref_mut(&mut self) -> &mut Self::Target
    {
        &mut self.base_model
    }
}




impl SWLargeDeviation
{
    pub fn new(base_model: SWModel, param: LargeDeviationParam,lockdown:LockdownParameters) -> Self
    {


        let mut markov_rng = Pcg64::seed_from_u64(param.markov_seed);
        let pairs_struct = create_lock_pairs_lists(lockdown, base_model.ensemble.graph(),&mut markov_rng);
        
        let zero_starting_conditions = param.zero_starting_conditions;

        #[cfg(feature = "ldprint")]   
        let printingvector = vec![5,6,7,8,9,10,55,56,57,58,59,60,61,62,63,64,65,66,95,96,97,98,99,100,195,196,197,198,199,200,295,296,297,298,299,300,395,396,397,398,399,400,495,496,497,498,499,500];
        
        let lock_graph =create_locked_down_network_from_pair_list(&pairs_struct, base_model.ensemble.graph());
        let not_lockdown_pairs = pairs_struct.to_be_removed;
        let lockdown_pairs = pairs_struct.to_be_kept;
        let lockdown_indicator = LockdownIndicator::NotLockdown;
        
        let system_size = NonZeroUsize::new(base_model.ensemble.graph().vertex_count())
            .unwrap();
        let rng_vec_len = system_size.get() * param.time_steps.get();

        let uniform = Uniform::new_inclusive(0.0_f64, 1.0);
        let mut transmission_rand_vec: Vec<_> = (0..rng_vec_len)
            .map(|_| uniform.sample(&mut markov_rng))
            .collect();

        let mut recovery_rand_vec: Vec<_> = (0..rng_vec_len)
            .map(|_| uniform.sample(&mut markov_rng))
            .collect();

        if zero_starting_conditions{
            transmission_rand_vec.iter_mut().for_each(|item| *item = 0.0);
            recovery_rand_vec.iter_mut().for_each(|item| *item = 0.0);

        }
        println!("{}",transmission_rand_vec[0]);


        let offset = Offset::new(param.time_steps.get(), system_size.get());
        let dangerous_neighbor_count = vec![0; system_size.get()];
        //let dangerous_neighbor_count_lock = vec![0; system_size.get()];
        let mut patient_zero_vec = Vec::new();
        let un = Uniform::new(0,base_model.ensemble.graph().vertex_count());
        let mut hist = HistUsizeFast::new(0, base_model.ensemble.vertex_count()).unwrap();

        while patient_zero_vec.len() < param.initial_infected{
            let index = un.sample(&mut markov_rng);
            if !patient_zero_vec.iter().any(|i| *i == index){
                //patient_zero_vec.infect_patient(index);
                patient_zero_vec.push(index);
                hist.increment_quiet(index);
            }
        }

        let max_degree = base_model.ensemble.graph().degree_iter().max().unwrap();

        Self{
            one_minus_lambda: 1.0 - base_model.lambda,
            base_model,
            lock_graph,
            not_lockdown_pairs,
            lockdown_pairs,
            markov_changed: true,
            lockdownparams:lockdown,
            markov_rng,
            time_steps: param.time_steps,
            offset,

            //writer,
            #[cfg(feature = "ldprint")]   
            printingvector,
            lockdown_indicator,
            energy: u32::MAX,
            old_energy: u32::MAX,
            transmission_rand_vec,
            zero_starting_conditions,
            recovery_rand_vec,
            unfinished_simulations_counter: 0,
            total_simulations_counter: 0,
            new_infected: Vec::new(),
            dangerous_neighbor_count,
            system_size,
            currently_infected_count: 0,
            last_extinction_index: usize::MAX,
            patient_zero_vec,
            max_degree,
            hist_patient_zero: hist
        }
    }

    pub fn hist_patient_zero(&self) -> &HistUsizeFast
    {
        &self.hist_patient_zero
    }

    //Change this to use lockdown indicator. Only one implementation.
    pub fn create_dangerous_neighbours_trans_sir(&mut self){

        //let n = self.system_size.get();
        self.dangerous_neighbor_count.iter_mut().for_each(|item| *item = 0);
        //let graph = self.base_model.ensemble.graph();

        match self.lockdown_indicator{
            LockdownIndicator::Lockdown => {
                self.base_model.ensemble.contained_iter()
                    .zip(
                        self.lock_graph.contained_iter_mut()
                    ).for_each(
                        |(old,new)|
                        {
                            *new= *old
                        }
                    );


                 for (index,contained) in self.lock_graph.contained_iter().enumerate(){
                     //let contained_mut = unsafe{self.lock_graph.get_contained_unchecked(index)};
                     
                     if contained.inf_check(){
                         let contained_index_iter = self.lock_graph.contained_iter_neighbors_with_index(index);
                         for (j, contained) in contained_index_iter
                         {
                             // only susceptible nodes are relevant later
                             if contained.sus_check(){
                                 // keep track of neighbors of infected nodes
                                 unsafe{
                                     *self.dangerous_neighbor_count.get_unchecked_mut(j) += 1;
                                 }
                             }
                         }
                     
                     }
                 
                 
                 
                 }

            },
            LockdownIndicator::NotLockdown => {
                    self.lock_graph.contained_iter()
                .zip(
                    self.base_model.ensemble.contained_iter_mut()
                ).for_each(
                    |(old,new)|
                    {
                        *new= *old
                    }
                );

                
                for (index,contained) in self.base_model.ensemble.graph().contained_iter().enumerate(){
                    if contained.inf_check(){
                        let contained_index_iter = self.base_model.ensemble.contained_iter_neighbors_with_index(index);
                        for (j, contained) in contained_index_iter
                        {
                            // only susceptible nodes are relevant later
                            if contained.sus_check(){
                                // keep track of neighbors of infected nodes
                                unsafe{
                                    *self.dangerous_neighbor_count.get_unchecked_mut(j) += 1;
                                }
                            }
                        }
    
                    }
    
                }
            }

        }
    }
   



    pub fn ld_iterate_once(&mut self){

        let new_infected = &mut self.new_infected;
        debug_assert!(new_infected.is_empty());


        //current implementation: If the system goes into lockdown, the dangerous neighbour list is recalculated before this function is called.
        //Hence, this iter contains the correct SIR information.
        let iter_danger = self.dangerous_neighbor_count.iter().enumerate();

        //New Infections. The choice of topology is handled here.

            
        for (index, &count) in iter_danger{
            if count > 0{
                #[cfg(feature = "debug_assertions")]
                match self.lockdown_indicator{
                    LockdownIndicator::Lockdown => {
                        debug_assert!(self.lock_graph.at(index).sus_check());

                    },
                    LockdownIndicator::NotLockdown => {
                        debug_assert!(self.base_model.ensemble.at(index).sus_check());

                    }
                }

                // probability of NOT infecting this node
                let prob = self.one_minus_lambda.powi(count);
                // infect if random number is bigger
                if self.transmission_rand_vec[self.offset.lookup_index(index)] >= prob {
                    new_infected.push(index);
            
                }
            }
        }

         
        
        let gamma = self.base_model.gamma;

        
        match self.lockdown_indicator{
            LockdownIndicator::Lockdown => {
                let graph = &mut self.lock_graph;

                for index in 0..self.system_size.get() {
                    let contained_mut = graph.at_mut(index);
                    if !contained_mut.inf_check(){
                        continue;
                    }
                    if self.recovery_rand_vec[self.offset.lookup_index(index)] < gamma {
                        *contained_mut = InfectionState::Recovered;
                        self.currently_infected_count -= 1;

                        let contained_index_iter = graph
                            .contained_iter_neighbors_with_index(index);
                        
                            
                        for (j, contained) in contained_index_iter
                        {
                            // only susceptible nodes are relevant later
                             if contained.sus_check(){
                             // keep track of neighbors of infected nodes
                                self.dangerous_neighbor_count[j] -= 1;
                             }
                        }
                    }

                }

            },
            LockdownIndicator::NotLockdown => {

                let graph = &mut self.base_model.ensemble;

                for index in 0..self.system_size.get() {
                    let contained_mut = graph.at_mut(index);
                    if !contained_mut.inf_check(){
                        continue;
                    }
                    if self.recovery_rand_vec[self.offset.lookup_index(index)] < gamma {
                        *contained_mut = InfectionState::Recovered;
                        self.currently_infected_count -= 1;

                        let contained_index_iter = graph
                            .contained_iter_neighbors_with_index(index);
                        

                        for (j, contained) in contained_index_iter
                        {
                            // only susceptible nodes are relevant later
                             if contained.sus_check(){
                             // keep track of neighbors of infected nodes
                                self.dangerous_neighbor_count[j] -= 1;
                             }
                        }
                    }

                }
            }
        };         

        
        //we have already dealt with the choice of topology, so we need only update the graphs's sir information.

        match self.lockdown_indicator{
            LockdownIndicator::Lockdown => {
                for &i in new_infected.iter()
                {   
                    let contained_mut = self.lock_graph.at_mut(i);
                    *contained_mut = InfectionState::Infected;
                    self.currently_infected_count += 1;
                    self.dangerous_neighbor_count[i] = 0;
                    let contained_index_iter = self.lock_graph
                    .contained_iter_neighbors_with_index(i);

                    for (j, contained) in contained_index_iter
                    {
                        if contained.sus_check(){
                            self.dangerous_neighbor_count[j] += 1;
                        }
                    }
                }

            },
            LockdownIndicator::NotLockdown =>{
                let graph = &mut self.base_model.ensemble;
                for &i in new_infected.iter(){   
                    let contained_mut = graph.at_mut(i);
                    *contained_mut = InfectionState::Infected;
                    self.currently_infected_count += 1;
                    self.dangerous_neighbor_count[i] = 0;
                    let contained_index_iter = graph
                    .contained_iter_neighbors_with_index(i);

                    for (j, contained) in contained_index_iter
                    {
                        if contained.sus_check(){
                            self.dangerous_neighbor_count[j] += 1;
                        }
                    }
                }
            }

        }

        new_infected.clear();
    }

    pub fn ld_iterate(&mut self) -> u32{

        self.total_simulations_counter += 1;
        self.reset_ld_sir_simulation();
        //let mut vec = Vec::new();
        //for (index,count) in self.dangerous_neighbor_count.iter().enumerate(){
        //    if *count > 0{
        //        vec.push(index);
        //    }
        //}
        //println!("{:?}",vec);
        let lockdown_threshold = self.lockdownparams.lock_threshold;
        let release_threshold = self.lockdownparams.release_threshold;

        let mut max_infected = self.currently_infected_count;
        
        //println!("{max_infected}");
        debug_assert_eq!(max_infected,self.patient_zero_vec.len() as u32);
        for i in 0..self.time_steps.get(){

            self.offset.set_time(i);
            let inf = self.currently_infected_count as f64 /self.system_size.get() as f64;

            //println!("inf {}",inf);
            if inf > lockdown_threshold && self.lockdown_indicator.is_not_in_lockdown(){
                self.lockdown_indicator = LockdownIndicator::Lockdown;                
                //println!("lock");
                self.create_dangerous_neighbours_trans_sir();
            }
            if inf < release_threshold && self.lockdown_indicator.is_in_lockdown(){
                self.lockdown_indicator = LockdownIndicator::NotLockdown;                
                //println!("rel");
                self.create_dangerous_neighbours_trans_sir();

            }

            self.ld_iterate_once();
            max_infected = max_infected.max(self.currently_infected_count);

            if self.currently_infected_count == 0{
                self.last_extinction_index = i+1;

                
                if self.lockdown_indicator.is_in_lockdown(){
                    self.lock_graph.contained_iter()
                        .zip(
                            self.base_model.ensemble.contained_iter_mut()
                        ).for_each(
                            |(old,new)|
                            {
                                *new= *old
                            }
                        );
                }                
                return max_infected;
            }
        }
        self.last_extinction_index = usize::MAX;
        self.unfinished_simulations_counter += 1;
        max_infected
    }


    pub fn ld_iterate_printing(&mut self,writer: &mut SirWriter)
    {   
        self.reset_ld_sir_simulation();
        let lockdown_threshold = self.lockdownparams.lock_threshold;
        let release_threshold = self.lockdownparams.release_threshold;
        let x = self.base_model.ensemble.graph().contained_iter().filter(|state| state.inf_check()).count();

        debug_assert_eq!(x,self.patient_zero_vec.len());
        //let mut max_infected = self.currently_infected_count;

        writer.write_current(self.ensemble().graph())
            .unwrap();
        
        for i in 0..self.time_steps.get(){

            self.offset.set_time(i);
            let inf = self.currently_infected_count as f64 /self.system_size.get() as f64;
            //println!("inf {}",inf);
            if inf > lockdown_threshold && self.lockdown_indicator.is_not_in_lockdown(){
                self.lockdown_indicator = LockdownIndicator::Lockdown;                
                //println!("lock");
                self.create_dangerous_neighbours_trans_sir();
            }
            if inf < release_threshold && self.lockdown_indicator.is_in_lockdown(){
                self.lockdown_indicator = LockdownIndicator::NotLockdown;                
                //println!("rel");
                self.create_dangerous_neighbours_trans_sir();
            }


            self.ld_iterate_once();
            


            //println!("{i}");
            
            match self.lockdown_indicator{
                LockdownIndicator::Lockdown => {
                    writer.write_current(&self.lock_graph)
                    .unwrap();
                },
                LockdownIndicator::NotLockdown =>{
                    writer.write_current(self.ensemble().graph())
                    .unwrap();
                }
            }

            

            if self.currently_infected_count == 0 {
                if self.lockdown_indicator.is_in_lockdown(){
                    self.lock_graph.contained_iter()
                        .zip(
                            self.base_model.ensemble.contained_iter_mut()
                        ).for_each(
                            |(old,new)|
                            {
                                *new= *old
                            }
                        );
                }
                break;
            }
        }
        writer.write_line().unwrap();
    }
    pub fn reset_ld_sir_simulation(&mut self){

        self.lockdown_indicator = LockdownIndicator::NotLockdown;

        self.dangerous_neighbor_count
            .iter_mut()
            .for_each(|v| *v = 0);


        let pzerovec = &self.patient_zero_vec;

        self.base_model.infect_many_patients(pzerovec);

        //Need all infected to be complete so we can calc dangerous neighbour list correctly.

        for j in &mut self.patient_zero_vec{
            let neighbor_iter = self
                .base_model
                .ensemble
                .contained_iter_neighbors_with_index(*j)
                .filter(|(_, state)| state.sus_check());
        
            for (i, _) in neighbor_iter
            {
                self.dangerous_neighbor_count[i] += 1;
            }

        }
        self.currently_infected_count = self.patient_zero_vec.len() as u32;
    }
    pub fn unfinished_counter(&self) -> u64 
    {
        self.unfinished_simulations_counter
    }

    pub fn total_simulations_counter(&self) -> u64
    {
        self.total_simulations_counter
    }

    pub fn last_extinction_index(&self) -> usize
    {
        self.last_extinction_index
    }
}



impl MarkovChain<MarkovStep, ()> for SWLargeDeviation
{
    fn m_step(&mut self) -> MarkovStep
    {
        unimplemented!()
    }

    fn undo_step(&mut self, step: &MarkovStep) 
    {
        self.undo_step_quiet(step)
    }

    fn undo_step_quiet(&mut self, step: &MarkovStep)
    {
        match step
        {
            
            MarkovStep::Transmission(info) => {
                self.transmission_rand_vec[info.index] = info.old_val;
            },
            MarkovStep::Recovery(info) => {
                self.recovery_rand_vec[info.index] = info.old_val;
            },
            MarkovStep::SwapTrans((a, b)) => {
                self.transmission_rand_vec.swap(*a, *b);
            },
            MarkovStep::SwapRec((a, b)) => {
                self.recovery_rand_vec.swap(*a, *b);
            }
            MarkovStep::RotateLeft => {
                // rotate right
                self.offset.minus_1();
            },
            MarkovStep::RotateRight => {
                // rotate left
                self.offset.plus_1();
            },
            MarkovStep::MovePatientZero(old_patient,index,..) | MarkovStep::MovePatientZeroRandom(old_patient,index) => {
                self.patient_zero_vec[*index] = *old_patient;
            }
            
        }
    }

    /* 
    fn m_steps_retired(&mut self, count: usize, steps: &mut Vec<MarkovStep>)
    {   
        //self.patient_move_boolean = false;
        steps.clear();
        self.markov_changed = true;
        let uniform = Uniform::new_inclusive(0.0_f64, 1.0);
        let which = uniform.sample(&mut self.markov_rng);
        //change to vec.len
        

        if which <= ROTATE_LEFT {
            //println!("rotating left");
            self.offset.plus_1();
            steps.push(MarkovStep::RotateLeft);
        } else if which <= ROTATE_RIGHT
        {
            //println!("rotating r");
            self.offset.minus_1();
            steps.push(MarkovStep::RotateRight);
        } else if which <= PATIENT_MOVE//maybe make this high, big degrees in BA
        {
            let f = 1.0 / self.max_degree as f64;
            let mut old = 0.0;
            let p = uniform.sample(&mut self.markov_rng);
            //Choose the patient zero in question at random from the array
            let uniform_usize = Uniform::new(0,self.initial_infected);
            let index_of_patient_to_be_changed = uniform_usize.sample(&mut self.markov_rng);
            //let neighbours_of_patient_zero:Vec<usize> = self.base_model.ensemble
            //.contained_iter_neighbors_with_index(self.patient_zero_vec[index_of_patient_to_be_changed]).map(|(index,_)| index).collect()
            //let neighbours_of_patient_zero_2 = self.base_model.ensemble.graph().container(self.patient_zero).neighbors()
            //change this to the iterator.
            for n in self.base_model.ensemble
            .contained_iter_neighbors_with_index(self.patient_zero_vec[index_of_patient_to_be_changed]).map(|(index,_)| index)
            {
                //let mut old = 0.0;
                let new_prob = f  + old;// <- less efficient
                if (old..new_prob).contains(&p)
                {   
                    if self.patient_zero_vec.contains(&n){
                        //println!("rejecting a patient move");
                        break
                        //continue
                    }
                    else{
                        //println!("doing a patient move");
                        //check if new index is already a p0
                        //yes: break
                        //no: below
                        let old_patient = self.patient_zero_vec[index_of_patient_to_be_changed];
                        self.patient_zero_vec[index_of_patient_to_be_changed]  = n;
                        //let x = self.patient_zero_vec[index_of_patient_to_be_changed];
                        //println!("old {old_patient}, new {x}");

                        //self.hist_patient_zero.increment_quiet(self.patient_zero_vec[index_of_patient_to_be_changed]);
                        for i in self.patient_zero_vec.iter(){
                            self.hist_patient_zero.increment_quiet(i);
                        }
                        
                        steps.push(MarkovStep::MovePatientZero(old_patient,index_of_patient_to_be_changed,true));
                        return;

                    }
                    
                }
                old = new_prob;
            }   
            self.markov_changed = false;
            //self.hist_patient_zero.increment_quiet(self.patient_zero_vec[index_of_patient_to_be_changed]);
            for i in self.patient_zero_vec.iter(){
                self.hist_patient_zero.increment_quiet(i);
            }
            steps.push(MarkovStep::MovePatientZero(self.patient_zero_vec[index_of_patient_to_be_changed],index_of_patient_to_be_changed,false));
            


        }
        else if which < PATIENT_MOVE_RANDOM{
            let index_of_patient_to_be_changed = self.markov_rng.gen_range(0..self.initial_infected);

            loop{
                let node_index = self.markov_rng.gen_range(0..self.system_size.get());
                if !self.patient_zero_vec.contains(&node_index){
                    steps.push(MarkovStep::MovePatientZeroRandom(self.patient_zero_vec[index_of_patient_to_be_changed],index_of_patient_to_be_changed));

                    self.patient_zero_vec[index_of_patient_to_be_changed] = node_index;
                    return;

                }

            }
            
        } 
        else {
            let amount:usize = Binomial::new(count as u64, 0.5).unwrap().sample(&mut self.markov_rng) as usize;
            let index_uniform = Uniform::new(0, self.system_size.get());
            let w = uniform.sample(&mut self.markov_rng);

            let tmp = count /2;
            let amount = self.markov_rng.gen_range(1..=tmp);
            let count = amount *2;

            let time_step = self.markov_rng.gen_range(0..self.time_steps.get());


            if w < 0.5 {
                steps.extend(
                    (0..count)
                        .map(
                            |_|
                            {
                                let index = index_uniform.sample(&mut self.markov_rng);
                                // Transmission 
                                let index = index + time_step*self.system_size.get();
                                    let old_val = std::mem::replace(
                                        &mut self.transmission_rand_vec[index],
                                        uniform.sample(&mut self.markov_rng)
                                    ); 
                                        MarkovStep::Transmission(
                                        ExchangeInfo{
                                            index,
                                            old_val
                                            }
                                        )
                                        
                                
                                
                            }
                        )
                );
            } else {
                steps.extend(
                    (0..(count))
                        .map(
                            |_|
                            {
                                let index = index_uniform.sample(&mut self.markov_rng);
                                // Transmission 
                                let index = index + time_step*self.system_size.get();

                                let old_val = std::mem::replace(
                                    &mut self.recovery_rand_vec[index],
                                    uniform.sample(&mut self.markov_rng)
                                ); 
                                        MarkovStep::Recovery(
                                        ExchangeInfo{
                                            index,
                                            old_val
                                            }
                                        )
                                
                            }
                        )
                );
            }
            

            
        }
    }
    */
    fn m_steps(&mut self, count: usize, steps: &mut Vec<MarkovStep>)
    {   
        //self.patient_move_boolean = false;
        steps.clear();

        let uniform = Uniform::new_inclusive(0.0_f64, 1.0);
        let which = uniform.sample(&mut self.markov_rng);
        //change to vec.len
        

        if which <= ROTATE_LEFT {
            //println!("rotating left");
            self.offset.plus_1();
            steps.push(MarkovStep::RotateLeft);
        } else if which <= ROTATE_RIGHT
        {
            //println!("rotating r");
            self.offset.minus_1();
            steps.push(MarkovStep::RotateRight);
        } else if which <= PATIENT_MOVE//maybe make this high, big degrees in BA
        {
            let f = 1.0 / self.max_degree as f64;
            let mut old = 0.0;
            let p = uniform.sample(&mut self.markov_rng);
            //Choose the patient zero in question at random from the array
            let uniform_usize = Uniform::new(0,self.patient_zero_vec.len());
            let index_of_patient_to_be_changed = uniform_usize.sample(&mut self.markov_rng);
            //let neighbours_of_patient_zero:Vec<usize> = self.base_model.ensemble
            //.contained_iter_neighbors_with_index(self.patient_zero_vec[index_of_patient_to_be_changed]).map(|(index,_)| index).collect()
            //let neighbours_of_patient_zero_2 = self.base_model.ensemble.graph().container(self.patient_zero).neighbors()
            //change this to the iterator.
            for n in self.base_model.ensemble
            .contained_iter_neighbors_with_index(self.patient_zero_vec[index_of_patient_to_be_changed]).map(|(index,_)| index)
            {
                //let mut old = 0.0;
                let new_prob = f  + old;// <- less efficient
                if (old..new_prob).contains(&p)
                {   
                    if self.patient_zero_vec.contains(&n){
                        //println!("rejecting a patient move");
                        break
                        //continue
                    }
                    else{
                        //println!("doing a patient move");
                        //check if new index is already a p0
                        //yes: break
                        //no: below
                        let old_patient = self.patient_zero_vec[index_of_patient_to_be_changed];
                        self.patient_zero_vec[index_of_patient_to_be_changed]  = n;
                        //let x = self.patient_zero_vec[index_of_patient_to_be_changed];
                        //println!("old {old_patient}, new {x}");

                        //self.hist_patient_zero.increment_quiet(self.patient_zero_vec[index_of_patient_to_be_changed]);
                        for i in self.patient_zero_vec.iter(){
                            self.hist_patient_zero.increment_quiet(i);
                        }
                        steps.push(MarkovStep::MovePatientZero(old_patient,index_of_patient_to_be_changed,true));
                        return;

                    }
                    
                }
                old = new_prob;
            }   
            self.markov_changed = false;
            //self.hist_patient_zero.increment_quiet(self.patient_zero_vec[index_of_patient_to_be_changed]);
            for i in self.patient_zero_vec.iter(){
                self.hist_patient_zero.increment_quiet(i);
            }
            steps.push(MarkovStep::MovePatientZero(self.patient_zero_vec[index_of_patient_to_be_changed],index_of_patient_to_be_changed,false));
            


        }
        else if which < PATIENT_MOVE_RANDOM{
            let index_of_patient_to_be_changed = self.markov_rng.gen_range(0..self.patient_zero_vec.len());

            loop{
                let node_index = self.markov_rng.gen_range(0..self.system_size.get());
                if !self.patient_zero_vec.contains(&node_index){
                    steps.push(MarkovStep::MovePatientZeroRandom(self.patient_zero_vec[index_of_patient_to_be_changed],index_of_patient_to_be_changed));

                    self.patient_zero_vec[index_of_patient_to_be_changed] = node_index;
                    return;

                }

            }
            
        } 
        else {
            let amount = Binomial::new(count as u64, 0.5).unwrap().sample(&mut self.markov_rng);
            let index_uniform = Uniform::new(0, self.recovery_rand_vec.len());
            
            steps.extend(
                (0..amount)
                    .map(
                        |_|
                        {
                            let index = index_uniform.sample(&mut self.markov_rng);
                            // Transmission 
                            let old_val = std::mem::replace(
                                &mut self.transmission_rand_vec[index], 
                                uniform.sample(&mut self.markov_rng)
                            );
                            MarkovStep::Transmission(
                                ExchangeInfo{
                                    index,
                                    old_val
                                }
                            )
                        }
                    )
            );

            steps.extend(
                (0..(count as u64 - amount))
                    .map(
                        |_|
                        {
                            let index = index_uniform.sample(&mut self.markov_rng);
                            // Recovery 
                            let old_val = std::mem::replace(
                                &mut self.recovery_rand_vec[index], 
                                uniform.sample(&mut self.markov_rng)
                            );
                            MarkovStep::Recovery(
                                ExchangeInfo{
                                    index,
                                    old_val
                                }
                            )
                        }
                    )
            );
        }
   
}
}





#[derive(Clone, Serialize, Deserialize)]
pub struct SWLargeDeviationWithLocks
{
    pub ld_model: SWLargeDeviation,
    markov_workaround: Vec<MarkovStep>,
    lockdown: LockdownParameters,
    pub tracker: AcceptanceTracker
}

impl HasRng<Pcg64> for SWLargeDeviationWithLocks{
    fn rng(&mut self) -> &mut Pcg64
    {
        &mut self.ld_model.markov_rng
    }

    fn swap_rng(&mut self, other: &mut Pcg64)
    {
        std::mem::swap(&mut self.ld_model.markov_rng, other)
    } 
}

impl SWLargeDeviationWithLocks
{
    pub fn new(ld_model: SWLargeDeviation) -> Self
        {
            Self{
                lockdown: ld_model.lockdownparams,

                ld_model,
                markov_workaround: Vec::new(),
                tracker: AcceptanceTracker::new()
            
            }
        }
    pub fn infect_initial_patients(&mut self)
        {
            let pzerovec = &self.ld_model.patient_zero_vec;

            self.ld_model.base_model.infect_many_patients(pzerovec);
            
        }


    pub fn ld_energy_m(&mut self) -> u32{
        self.ld_model.ld_iterate() //this gives M
    }


    pub fn ld_energy_m_and_print(&mut self, writer: &mut SirWriter)
    {
        //list should already be calculated!
        // That means that the following function should short-circut
        
        self.ld_model.ld_iterate_printing(writer);

    }
    pub fn random_lockdown_monte_carlo(&mut self){
        //I guess to randomie the montecarlo for the random lockdown.. we just have to redraw the edges to be removed?
        let pairs_struct = create_lock_pairs_lists(self.lockdown, self.ld_model.base_model.ensemble.graph(),&mut self.ld_model.markov_rng);
        
        self.ld_model.lock_graph = create_locked_down_network_from_pair_list(&pairs_struct, self.ld_model.base_model.ensemble.graph());
        self.ld_model.lockdown_pairs = pairs_struct.to_be_kept;
        self.ld_model.not_lockdown_pairs = pairs_struct.to_be_removed;
    }

    pub fn randomise_monte_carlo(&mut self,rng: &mut Pcg64){

        match self.lockdown.lock_style{
            LockdownType::Random(_) => self.random_lockdown_monte_carlo(),
            //other lockdown types, such as neone and circuit breaker, do not need to be randomised.
            LockdownType::CircuitBreaker => {},
            LockdownType::None => {},
            _ => unimplemented!() 
        };

        let uniform = Uniform::new_inclusive(0.0_f64, 1.0);
        let mut iter = uniform.sample_iter(rng);

        //randomising order of transmission vector.
        self.transmission_rand_vec.iter_mut().for_each(|elem| *elem = iter.next().unwrap());


        let mut patient_zero_vec = Vec::new();
        
        let uniform = Uniform::new(0_usize, self.base_model.ensemble.vertex_count());
        while patient_zero_vec.len() < self.ld_model.patient_zero_vec.len(){
            let index = uniform.sample(&mut self.markov_rng);
            if !patient_zero_vec.iter().any(|i| *i == index){
                //patient_zero_vec.infect_patient(index);
                patient_zero_vec.push(index);
            }
        }
        self.patient_zero_vec = patient_zero_vec;
        self.recovery_rand_vec
            .iter_mut()
            .zip(iter)
            .for_each(
                |(elem, random)|
                *elem = random
            );
    }
    

    pub fn find_lockdown_markovmove(&mut self) -> LockdownMarkovMove{
        match self.ld_model.lockdownparams.lock_style{
            LockdownType::Random(_) => self.find_replace_edge_for_random_lock(),
            _ => unimplemented!()
        }
    }

    pub fn find_replace_edge_for_random_lock(&mut self) -> LockdownMarkovMove{

        let uniform = Uniform::new(0,self.not_lockdown_pairs.len());
        let not_lockdown_index = uniform.sample(&mut self.markov_rng);
        //let pair_to_be_restored = self.not_lockdown_pairs[not_lockdown_index];

        let uniform = Uniform::new(0,self.lockdown_pairs.len());
        let lockdown_index = uniform.sample(&mut self.markov_rng);
        //let pair_to_be_removed = self.lockdown_pairs[lockdown_index];

        
        LockdownMarkovMove{
            lockdown_index,
            not_lockdown_index
        }

    
    }
    pub fn find_replace_edge_for_limiting_lock(&mut self) -> LockdownMarkovMove{
        //not yet implemented.
        let lockdown_index = 0;
        let not_lockdown_index = 0;
        let uniform = Uniform::new(0,self.system_size.get());
        let _node_to_be_changed = uniform.sample(&mut self.markov_rng);
        

        LockdownMarkovMove{
            lockdown_index,
            not_lockdown_index
        }
    }



    pub fn change_edge(&mut self, lockmarkovmove:&LockdownMarkovMove){

        let b = &mut self.ld_model.lockdown_pairs[lockmarkovmove.lockdown_index];
        let a = &mut self.ld_model.not_lockdown_pairs[lockmarkovmove.not_lockdown_index];

        let lockgraph = &mut self.ld_model.lock_graph;
        lockgraph.remove_edge(b[0],b[1]).unwrap();
        lockgraph.add_edge(a[0],a[1]).unwrap();

       
        std::mem::swap(a, b);
    }


}
pub struct Measure{
    pub max_infected: u32,
    pub ever_infected: u32
}
impl Deref for SWLargeDeviationWithLocks
{
    type Target = SWLargeDeviation;

    fn deref(&self) -> &Self::Target
    {
        &self.ld_model
    }
}

impl DerefMut for SWLargeDeviationWithLocks
{
    fn deref_mut(&mut self) -> &mut Self::Target
    {
        &mut self.ld_model
    }
}
impl MarkovChain<MarkovStepWithLocks, ()> for SWLargeDeviationWithLocks
{
    fn m_step(&mut self) -> MarkovStepWithLocks
    {
        unimplemented!()
    }

    fn undo_step(&mut self, step: &MarkovStepWithLocks) 
    {
        self.undo_step_quiet(step)
    }

    fn undo_step_quiet(&mut self, step: &MarkovStepWithLocks)
    {
        match step{
            MarkovStepWithLocks::BaseMarkovStep(base_step) => {
                self.deref_mut().undo_step_quiet(base_step)
            },
            MarkovStepWithLocks::LockdownStep(to_step) => {
                //
                //possibly have a match in here to the appropriate undo of the lockdown step corresponding to the correct lockdown.
                //the method of doing swap in the 
                self.change_edge(to_step)
                //self.high_degree_helper.swap_order(to_step.0, to_step.1)
            }
        }
    }
    
    fn undo_steps_quiet(&mut self, steps: &[MarkovStepWithLocks]) {
        self.energy = self.old_energy;
        steps.iter()
            .rev()
            .for_each( |step| self.undo_step_quiet(step));
    }

    fn steps_accepted(&mut self, steps: &[MarkovStepWithLocks])
    {
        self.tracker.accept(&steps[0]);
    }

    fn steps_rejected(&mut self, steps: &[MarkovStepWithLocks])
    {
        self.tracker.reject(&steps[0])
    }

    fn m_steps(&mut self, count: usize, steps: &mut Vec<MarkovStepWithLocks>)
    {   
        self.ld_model.markov_changed = true;

        steps.clear();

        match self.lockdown.lock_style{
            LockdownType::Random(_) => 
            {
                let uniform = Uniform::new_inclusive(0.0_f64, 1.0);
                let which = uniform.sample(&mut self.markov_rng);


                if which <= LOCKDOWN_CHANGE
                {   
                    //println!("lockdown markov move");
                    //let rng  = &mut self.markov_rng;
                    let num_moves = self.markov_rng.gen_range(1..=15);
                    for _ in 0..num_moves{
                        let ld_step = self.find_lockdown_markovmove();
                        self.change_edge(&ld_step);
                        steps.push(MarkovStepWithLocks::LockdownStep(ld_step));
                    }
                    
                

                } else {
                    //println!("not lockdown markov move");
                    self.ld_model.m_steps(count, &mut self.markov_workaround);
                
                    steps.extend(
                        self.markov_workaround
                            .iter()
                            .copied()
                            .map(MarkovStepWithLocks::from)
                    );
                    self.markov_workaround.clear();
                }

            },
            _ =>{ 
                //other lockdowns have no markov moves as they are not random. hence we just do the 'else' of the other step.
                self.ld_model.m_steps(count, &mut self.markov_workaround);

                steps.extend(
                    self.markov_workaround
                        .iter()
                        .copied()
                        .map(MarkovStepWithLocks::from)
                );
                self.markov_workaround.clear();

                }
        }

        
    }

    
}


#[cfg(test)]
mod tests {
    //use crate::large_deviations::*;
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    use crate::sir_model::*;
    #[test]
    fn markov_test() {

        let param = LargeDeviationParam {
            time_steps: ONE,
            markov_seed: DEFAULT_MARKOV_SEED,
            initial_infected:DEFAULT_INITIAL_INFECTED,
            zero_starting_conditions:false
        };

        let lockparams = LockdownParameters{
            lock_style: LockdownType::Random(0.6),
            lock_threshold: 0.1,
            release_threshold: 0.05
        };

        let base_options = SWOptions{
            graph_seed: DEFAULT_GRAPH_SEED,
            lambda: DEFAULT_LAMBDA,
            gamma: DEFAULT_RECOVERY_PROB,
            system_size: DEFAULT_SYSTEM_SIZE,
            rewire_prob: 0.1}; 
    
        let ba:SWModel = base_options.into();
        let ld  = SWLargeDeviation::new(ba,param,lockparams);

        let mut test_model = SWLargeDeviationWithLocks::new(ld);
        
        //let markov = Randomize::default();

        let num_steps = 1000;
        for _i in 0..10000

        {   
            let energy = test_model.ld_energy_m();
            let mut markov_vec = Vec::new();

            
            test_model.m_steps(num_steps, &mut markov_vec);
            //assert_eq!(num_steps,markov_vec.len());
            
            println!("{}",test_model.ld_energy_m());    
            
            //println!("{}",markov_vec.len());

            test_model.undo_steps_quiet(&markov_vec);
            //println!("{energy}");
            /*markov_vec.reverse();
            for step in &mut markov_vec
                {
                    test_model.undo_step(&step)
                }
                
        */
            let energy_2 = test_model.ld_energy_m();
            println!("{energy},{energy_2}");
            assert_eq!(energy,energy_2)

    
        }

    }

    
}