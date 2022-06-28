use net_ensembles::*;
use serde::{Serialize, Deserialize};
//use crate::sir_model::sir_states::*;
use rand::*;
use rand::distributions::*;
use rand_pcg::Pcg64;
use net_ensembles::rand::prelude::Distribution;
#[derive(Debug, Clone,Serialize, Deserialize,Copy)]

pub enum LockdownType{
    CircuitBreaker,
    None,
    Random,
    Invalid

}



pub fn lockdown<InfectionState,A>(lockdowntype:LockdownType,graph:GenericGraph<InfectionState,A>) -> GenericGraph<InfectionState,A> 
where A: net_ensembles::AdjContainer<InfectionState>, InfectionState: net_ensembles::Node, A: std::clone::Clone{
    match lockdowntype{
        LockdownType::CircuitBreaker => circuit_breaking_lockdown(graph),
        LockdownType::None => no_lockdown(graph),
        LockdownType::Random => random_lockdown(graph),
        _ => unimplemented!()
    }


}
pub fn lockdown_naming_string(lockdowntype:LockdownType) -> String{
    match lockdowntype{
        LockdownType::CircuitBreaker => "CircuitBreaker".to_owned(),
        LockdownType::Random => "Random".to_owned(),
        LockdownType::None => "None".to_owned(),
        _ => unimplemented!()
    }

}

pub fn circuit_breaking_lockdown<T,A>(mut graph:GenericGraph<T,A>)-> GenericGraph<T,A> where A: net_ensembles::AdjContainer<T>{
    graph.clear_edges();
    
    graph

}

pub fn no_lockdown<T,A>(graph:GenericGraph<T,A>)-> GenericGraph<T,A> where A: net_ensembles::AdjContainer<T>{graph}

pub fn random_lockdown<InfectionState,A>(graph:GenericGraph<InfectionState,A>)-> GenericGraph<InfectionState,A> 
where A: net_ensembles::AdjContainer<InfectionState>, InfectionState: net_ensembles::Node, A: std::clone::Clone{
    //note to self. add some more enums and structs, eg like running the execute files in main to better capture the parameters of each of the individual lockdowns.__rust_force_expr!;
    let mut clone = graph.clone();
    let pairs = pair_finder(graph);
    let prob_dist= Uniform::new_inclusive(0.0,1.0);
    let mut rng = Pcg64::seed_from_u64(12312313131231414);
    let prob = 0.9;
    for pair in pairs{

        let i1 = pair[0];
        let i2 = pair[1];
        //println!("{},{}",i1,i2);
        let p = prob_dist.sample(&mut rng);
        if p < prob{
            clone.remove_edge(i1,i2 ).unwrap();
        }

    }


    clone}


pub fn lockdown_target_infections(){
    //the goal of this will be to target those who are infected, and lockdown the surround clique of individuals. for instance, suppose 
    //an infected as three susceptible neighbours. those neighbours have their links removed to their neighbours (not the infected though)
    //this perhaps simulates the close-contact dynamics of the lockdown.
    //some credit is due to Peter Werner for discussing this with me.
}
//remember the remove_edge(index1,index2 ) command!

pub fn pair_finder<InfectionState,A>(mut graph:GenericGraph<InfectionState,A>)-> Vec<Vec<usize>> 
where A: net_ensembles::AdjContainer<InfectionState>, InfectionState: net_ensembles::Node, A: std::clone::Clone{
    //This function finds all the pairs of indicies connected in the graph. 

    let mut vec:Vec<Vec<_>> = Vec::new();
    for j in 0..graph.vertex_count(){
        let adj = graph.contained_iter_neighbors_mut_with_index(j);
        for (n_index,_neighbor) in adj{
            let mut new_vec = vec![j,n_index];
            new_vec.sort_unstable();
            vec.push(new_vec);
        } 
    }
    vec.sort();
    vec.dedup();
    vec
    
}

