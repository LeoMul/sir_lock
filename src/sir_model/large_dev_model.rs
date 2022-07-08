use {
    super::*,
    std::{num::*, ops::{Deref, DerefMut}},
    serde::{Serialize, Deserialize},
    net_ensembles::rand::{distributions::{Uniform, Distribution}, SeedableRng},
    net_ensembles::{GraphIterators, MarkovChain, MeasurableGraphQuantities, WithGraph},
    rand_pcg::Pcg64,
    net_ensembles::sampling::{HasRng, histogram::*},
    rand_distr::Binomial,
    crate::lockdown_methods::*
};

//might need but leaving it here.
//use net_ensembles::Contained,


const ROTATE_LEFT: f64 =  0.005;
const ROTATE_RIGHT: f64 =  0.01;
const PATIENT_MOVE: f64 = 0.03;

pub struct LockdownMarkovMove{
    pair_to_be_removed: [usize;2],
    removed_index: usize,
    pair_to_be_restored: [usize;2],
    restored_index: usize
}



use net_ensembles::GraphIteratorsMut;

const LOCKDOWN_CHANGE: f64 = 0.20;

#[derive(Clone, Serialize, Deserialize)]
pub struct BALargeDeviation
{
    base_model: BarabasiModel,
    lock_graph: GenGraphSIR,
    removed_pairs: Vec<[usize;2]>,
    remaining_pairs: Vec<[usize;2]>,
    lockdownparams:LockdownParameters,
    transmission_rand_vec: Vec<f64>,
    recovery_rand_vec: Vec<f64>,
    time_steps: NonZeroUsize,
    unfinished_simulations_counter: u64,
    total_simulations_counter: u64,
    markov_rng: Pcg64,
    offset: Offset,
    // temporary storage, allocated here so that 
    // I do not need to allocate it again and again
    new_infected: Vec<usize>,
    // Number of infected neighbors of a node
    // but only if the node is susceptible,
    // otherwise this is 0, as there is no danger
    dangerous_neighbor_count: Vec<i32>,
    one_minus_lambda: f64,
    system_size: NonZeroUsize,
    /// counting how many nodes are currently infected
    currently_infected_count: u32,
    last_extinction_index: usize,
    patient_zero: usize,
    max_degree: usize,
    hist_patient_zero: HistUsizeFast
}
impl Deref for BALargeDeviation
{
    type Target = BarabasiModel;
    fn deref(&self) -> &Self::Target
    {
        &self.base_model
    }
}

impl DerefMut for BALargeDeviation
{
    fn deref_mut(&mut self) -> &mut Self::Target
    {
        &mut self.base_model
    }
}
#[derive(Clone, Serialize, Deserialize, Copy)]
pub struct LargeDeviationParam
{
    pub time_steps: NonZeroUsize,
    pub markov_seed: u64
}

impl BALargeDeviation
{
    pub fn new(base_model: BarabasiModel, param: LargeDeviationParam,lockdown:LockdownParameters) -> Self
    {

        //right now, using two different dangerous neighbour lists. Perhaps it is better to use only one.

        //do I need the lockdown parameters here also? yes
        let pairs_struct = create_lock_pairs_lists(lockdown, base_model.ensemble.graph());

        


        let mut markov_rng = Pcg64::seed_from_u64(param.markov_seed);

        let lock_graph =create_locked_down_network_from_pair_list(&pairs_struct, base_model.ensemble.graph());

        let removed_pairs = pairs_struct.to_be_removed;
        let remaining_pairs = pairs_struct.to_be_kept;


        let system_size = NonZeroUsize::new(base_model.ensemble.graph().vertex_count())
            .unwrap();
        let rng_vec_len = system_size.get() * param.time_steps.get();

        let uniform = Uniform::new_inclusive(0.0_f64, 1.0);

        let transmission_rand_vec: Vec<_> = (0..rng_vec_len)
            .map(|_| uniform.sample(&mut markov_rng))
            .collect();

        let recovery_rand_vec: Vec<_> = (0..rng_vec_len)
            .map(|_| uniform.sample(&mut markov_rng))
            .collect();

        let offset = Offset::new(param.time_steps.get(), system_size.get());
        let dangerous_neighbor_count = vec![0; system_size.get()];
        //let dangerous_neighbor_count_lock = vec![0; system_size.get()];
        let patient_zero = Uniform::new(0, base_model.ensemble.vertex_count()).sample(&mut markov_rng);
        let mut hist = HistUsizeFast::new(0, base_model.ensemble.vertex_count()).unwrap();
        hist.increment_quiet(patient_zero);

        let max_degree = base_model.ensemble.graph().degree_iter().max().unwrap();

        Self{
            one_minus_lambda: 1.0 - base_model.lambda,
            base_model,
            lock_graph,
            removed_pairs,
            remaining_pairs,
            lockdownparams:lockdown,
            markov_rng,
            time_steps: param.time_steps,
            offset,
            transmission_rand_vec,
            recovery_rand_vec,
            unfinished_simulations_counter: 0,
            total_simulations_counter: 0,
            new_infected: Vec::new(),
            dangerous_neighbor_count,
            system_size,
            currently_infected_count: 0,
            last_extinction_index: usize::MAX,
            patient_zero,
            max_degree,
            hist_patient_zero: hist
        }
    }

    pub fn hist_patient_zero(&self) -> &HistUsizeFast
    {
        &self.hist_patient_zero
    }

    //Change this to use lockdown indicator. Only one implementation.
    pub fn create_dangerous_neighbours_trans_sir(&mut self,lockdown_indicator:bool){

        let n = self.system_size.get();
        self.dangerous_neighbor_count.iter_mut().for_each(|item| *item = 0);
        //let graph = self.base_model.ensemble.graph();

        if !lockdown_indicator{
            //this transfers the sir information to the base graph and creates the new dangerous neighbour count.
            self.lock_graph.contained_iter()
            .zip(
                self.base_model.ensemble.contained_iter_mut()
            ).for_each(
                |(old,new)|
                {
                    *new= *old
                }
            );


            for index in 0..n{
                let contained_mut = unsafe{self.base_model.ensemble.graph().get_contained_unchecked(index)};
                if !contained_mut.inf_check(){
                    continue;
                }
                else{
                    let contained_index_iter = self.base_model.ensemble.contained_iter_neighbors_with_index(index);
                    for (j, contained) in contained_index_iter
                    {
                        // only susceptible nodes are relevant later
                        if contained.sus_check(){
                            // keep track of neighbors of infected nodes
                            self.dangerous_neighbor_count[j] += 1;
                        }
                    }

                }

            }
        }
        else{
            self.base_model.ensemble.contained_iter()
            .zip(
                self.lock_graph.contained_iter_mut()
            ).for_each(
                |(old,new)|
                {
                    *new= *old
                }
            );


            for index in 0..n{
                let contained_mut = unsafe{self.lock_graph.get_contained_unchecked(index)};
                if !contained_mut.inf_check(){
                    continue;
                }
                else{
                    let contained_index_iter = self.lock_graph.contained_iter_neighbors_with_index(index);
                    for (j, contained) in contained_index_iter
                    {
                        // only susceptible nodes are relevant later
                        if contained.sus_check(){
                            // keep track of neighbors of infected nodes
                           self.dangerous_neighbor_count[j] += 1;
                        }
                    }
    
                }
    
            }

        }
    }
    pub fn create_dangerous_neighbours_from_lock(&mut self){
        unimplemented!()
    }



    pub fn ld_iterate_once(&mut self,lockdown_indicator:bool){
        let new_infected = &mut self.new_infected;
        debug_assert!(new_infected.is_empty());


        //current implementation: If the system goes into lockdown, the dangerous neighbour list is recalculated before this function is called.
        //Hence, this iter contains the correct SIR information.
        let iter_danger = self.dangerous_neighbor_count.iter().enumerate();

        //New Infections. The choice of topology is handled here.
        if !lockdown_indicator{
            for (index, &count) in iter_danger{
                if count > 0{
                    debug_assert!(self.base_model.ensemble.at(index).sus_check());
                    // probability of NOT infecting this node
                    let prob = self.one_minus_lambda.powi(count);
                    // infect if random number is bigger
                    if self.transmission_rand_vec[self.offset.lookup_index(index)] >= prob {
                        new_infected.push(index);
                
                }
            }
        }
        } else{
            for (index,&count) in iter_danger{
                if count > 0{
                    debug_assert!(self.lock_graph.at(index).sus_check());
                    let prob = self.one_minus_lambda.powi(count);
                    if self.transmission_rand_vec[self.offset.lookup_index(index)]>= prob{
                        new_infected.push(index);
                    }

                }
            }
        }
        
        let gamma = self.base_model.gamma;

        //if we run into problems, this is probably part of it
        //let graph = &mut self.base_model.ensemble;
        
        //doing the recoveries.
        if !lockdown_indicator{
            let graph = &mut self.base_model.ensemble;
            for index in 0..self.system_size.get() {
                let contained_mut = graph.at_mut(index);
                if !contained_mut.inf_check(){
                    continue;
                }
                if self.recovery_rand_vec[self.offset.lookup_index(index)] < gamma {
                    *contained_mut = InfectionState::Recovered;
                    self.currently_infected_count -= 1;
                }
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
        else{
            
            for index in 0..self.system_size.get() {
                let contained_mut = self.lock_graph.at_mut(index);
                if !contained_mut.inf_check(){
                    continue;
                }
                if self.recovery_rand_vec[self.offset.lookup_index(index)] < gamma {
                    *contained_mut = InfectionState::Recovered;
                    self.currently_infected_count -= 1;
                }
                let contained_index_iter = self.lock_graph
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




        //transfering SIR infor
        

        
        //we have already dealt with the choice of topology, so we need only update the graphs's sir information.

        
        if !lockdown_indicator{
            let graph = &mut self.base_model.ensemble;
            for &i in new_infected.iter()
            {   let contained_mut = graph.at_mut(i);
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
        else{
            //let graph = &mut self.base_model.ensemble;
            for &i in new_infected.iter()
            {   let contained_mut = self.lock_graph.at_mut(i);
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
        }

        new_infected.clear();
    }

    pub fn ld_iterate(&mut self) -> u32{
        self.total_simulations_counter += 1;
        self.reset_ld_sir_simulation();

        let lockdown_threshold = self.lockdownparams.lock_threshold;
        let release_threshold = self.lockdownparams.release_threshold;

        let mut max_infected = self.currently_infected_count;
        let mut lockdown_indicator = false;

        for i in 0..self.time_steps.get(){
            self.offset.set_time(i);
            let inf = self.currently_infected_count as f64 /self.system_size.get() as f64;
            if inf > lockdown_threshold && !lockdown_indicator{
                lockdown_indicator = true;
                self.create_dangerous_neighbours_trans_sir(lockdown_indicator);
            }
            if inf < release_threshold && lockdown_indicator{
                lockdown_indicator = false;
                self.create_dangerous_neighbours_trans_sir(lockdown_indicator);

            }

            self.ld_iterate_once(lockdown_indicator);
            if self.currently_infected_count == 0{
                self.last_extinction_index = i+1;
                return max_infected;
            }
            max_infected = max_infected.max(self.currently_infected_count);
        }
        self.last_extinction_index = usize::MAX;
        self.unfinished_simulations_counter += 1;
        max_infected
    }
    pub fn reset_ld_sir_simulation(&mut self){
        let p0 = self.patient_zero;
        self.infect_patient(p0);

        self.dangerous_neighbor_count
            .iter_mut()
            .for_each(|v| *v = 0);
        
        
        let neighbor_iter = self
            .base_model
            .ensemble
            .contained_iter_neighbors_with_index(self.patient_zero)
            .filter(|(_, state)| state.sus_check());
        
        for (i, _) in neighbor_iter
        {
            self.dangerous_neighbor_count[i] += 1;
        }
        self.currently_infected_count = 1;
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
    MovePatientZero(usize)
}
impl MarkovChain<MarkovStep, ()> for BALargeDeviation
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
            MarkovStep::MovePatientZero(old_patient) => {
                self.patient_zero = *old_patient;
            }
        }
    }

    fn m_steps(&mut self, count: usize, steps: &mut Vec<MarkovStep>)
    {
        steps.clear();

        let uniform = Uniform::new_inclusive(0.0_f64, 1.0);
        let which = uniform.sample(&mut self.markov_rng);

        if which <= ROTATE_LEFT {
            self.offset.plus_1();
            steps.push(MarkovStep::RotateLeft);
        } else if which <= ROTATE_RIGHT
        {
            self.offset.minus_1();
            steps.push(MarkovStep::RotateRight);
        } else if which <= PATIENT_MOVE//maybe make this high, big degrees in BA
        {
            let f = 1.0 / self.max_degree as f64;
            let mut old = 0.0;
            let p = uniform.sample(&mut self.markov_rng);

            let neighbours_of_patient_zero:Vec<usize> = self.base_model.ensemble
            .contained_iter_neighbors_with_index(self.patient_zero).map(|(index,_)| index).collect();

            //let neighbours_of_patient_zero_2 = self.base_model.ensemble.graph().container(self.patient_zero).neighbors();
            for n in neighbours_of_patient_zero
            {
                let new_prob = f  + old;// <- less efficient
                if (old..new_prob).contains(&p)
                {
                    let old_patient = self.patient_zero;
                    self.patient_zero  = n;
                    self.hist_patient_zero.increment_quiet(self.patient_zero);
                    steps.push(MarkovStep::MovePatientZero(old_patient));
                    return;
                }
                old = new_prob;
            }
            self.hist_patient_zero.increment_quiet(self.patient_zero);
            steps.push(MarkovStep::MovePatientZero(self.patient_zero));
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





pub enum MarkovStepWithLocks{
    BaseMarkovStep(MarkovStep),
    LockdownStep(LockdownMarkovMove) 
}
impl From<MarkovStep> for MarkovStepWithLocks{
    fn from(other:MarkovStep) -> Self{
        Self::BaseMarkovStep(other)
    }
}


#[derive(Clone, Serialize, Deserialize)]
pub struct BALargeDeviationWithLocks
{
    ld_model: BALargeDeviation,
    markov_workaround: Vec<MarkovStep>,
    lockdown: LockdownParameters,
}

impl HasRng<Pcg64> for BALargeDeviationWithLocks{
    fn rng(&mut self) -> &mut Pcg64
    {
        &mut self.markov_rng
    }

    fn swap_rng(&mut self, other: &mut Pcg64)
    {
        std::mem::swap(&mut self.markov_rng, other)
    } 
}

impl BALargeDeviationWithLocks
{
    pub fn new(ld_model: BALargeDeviation) -> Self
        {
            Self{
                lockdown: ld_model.lockdownparams,

                ld_model,
                markov_workaround: Vec::new(),
            
            }
        }
    pub fn infect_patient(&mut self)
        {
            let p0 = self.patient_zero;
            self.ld_model.infect_patient(p0);
        }
    pub fn ld_energy_m(&mut self) -> u32{
        self.ld_model.ld_iterate() //this gives M
    }
    
    

    pub fn randomise_monte_carlo(&mut self,rng: &mut Pcg64){

        //need to do something with lockdowns here
        match self.lockdown.lock_style{
            LockdownType::Random(_,_) => random_lockdown_monte_carlo(),
            _ => unimplemented!()
        };

        let uniform = Uniform::new_inclusive(0.0_f64, 1.0);
        let mut iter = uniform.sample_iter(rng);

        //randomising order of transmission vector.
        self.transmission_rand_vec.iter_mut().for_each(|elem| *elem = iter.next().unwrap());

        let uniform = Uniform::new(0_usize, self.base_model.ensemble.vertex_count());

        let patient_zero = uniform.sample(&mut self.markov_rng);
        self.patient_zero = patient_zero;
        self.recovery_rand_vec
            .iter_mut()
            .zip(iter)
            .for_each(
                |(elem, random)|
                *elem = random
            );
    }
    #[inline]
    pub fn patient_zero(&self) -> usize
    {
        self.ld_model.patient_zero
    }

    pub fn find_lockdown_markovmove(&mut self) -> LockdownMarkovMove{
        match self.lockdownparams.lock_style{
            LockdownType::Random(_,_) => self.find_replace_edge_for_random_lock(),
            _ => unimplemented!()
        }
    }

    pub fn find_replace_edge_for_random_lock(&mut self) -> LockdownMarkovMove{

        let uniform = Uniform::new(0,self.removed_pairs.len());
        let restored_index = uniform.sample(&mut self.markov_rng);
        let pair_to_be_restored = self.removed_pairs[restored_index];

        let uniform = Uniform::new(0,self.remaining_pairs.len());
        let removed_index = uniform.sample(&mut self.markov_rng);
        let pair_to_be_removed = self.remaining_pairs[removed_index];

        
        LockdownMarkovMove{
            pair_to_be_removed,
            removed_index,
            pair_to_be_restored,
            restored_index
        }

    
    }
    pub fn change_edge(&mut self, lockmarkovmove:&LockdownMarkovMove){

        self.removed_pairs.swap_remove(lockmarkovmove.restored_index);
        self.removed_pairs.push(lockmarkovmove.pair_to_be_removed);
        self.remaining_pairs.swap_remove(lockmarkovmove.removed_index);
        self.removed_pairs.push(lockmarkovmove.pair_to_be_restored);



        let mut lockgraph = self.lock_graph.clone();
        let rem = lockmarkovmove.pair_to_be_removed;
        let rep = lockmarkovmove.pair_to_be_restored;
        lockgraph.remove_edge(rem[0],rem[1]).unwrap();
        lockgraph.add_edge(rep[0],rep[1]).unwrap();
        self.lock_graph = lockgraph;

    }


}
pub struct Measure{
    pub max_infected: u32,
    pub ever_infected: u32
}
impl Deref for BALargeDeviationWithLocks
{
    type Target = BALargeDeviation;

    fn deref(&self) -> &Self::Target
    {
        &self.ld_model
    }
}

impl DerefMut for BALargeDeviationWithLocks
{
    fn deref_mut(&mut self) -> &mut Self::Target
    {
        &mut self.ld_model
    }
}
impl MarkovChain<MarkovStepWithLocks, ()> for BALargeDeviationWithLocks
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
                self.change_edge(to_step)
                //self.high_degree_helper.swap_order(to_step.0, to_step.1)
            }
        }
    }

    fn m_steps(&mut self, count: usize, steps: &mut Vec<MarkovStepWithLocks>)
    {
        steps.clear();
        let uniform = Uniform::new_inclusive(0.0_f64, 1.0);

        let which = uniform.sample(&mut self.markov_rng);


        //first we orient ourselves with the random lockdown.
        if which <= LOCKDOWN_CHANGE
        {
            //let rng  = &mut self.markov_rng;
            let ld_step = self.find_lockdown_markovmove();
            steps.push(MarkovStepWithLocks::LockdownStep(ld_step));
    
            
        } else {
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

pub fn random_lockdown_monte_carlo(){
    unimplemented!()
}