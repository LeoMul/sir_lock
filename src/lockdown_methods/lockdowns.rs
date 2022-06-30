use net_ensembles::*;
use serde::{Serialize, Deserialize};
//use crate::sir_model::sir_states::*;
use rand::*;
use rand::distributions::*;
use rand_pcg::Pcg64;
use net_ensembles::rand::prelude::Distribution;
use crate::misc_types::*;
#[derive(Debug, Clone,Serialize, Deserialize,Copy)]

pub enum LockdownType{
    CircuitBreaker,
    None,
    Random(u64,f64),
    Targeted,
    Invalid

}

#[derive(Debug, Clone,Serialize, Deserialize,Copy)]
pub struct LockdownParameters{
    pub lock_style: LockdownType,
    pub release_threshold: f64,
    pub lock_threshold: f64,
    pub dynamic_bool: bool,
}




pub fn lockdown(lockdownparams: LockdownParameters,graph:GenGraphSIR,inflist:&Vec<usize>) -> GenGraphSIR
{
    
    match lockdownparams.lock_style{
        LockdownType::CircuitBreaker => circuit_breaking_lockdown(graph),
        LockdownType::None => no_lockdown(graph),
        LockdownType::Random(seed,prob) => random_lockdown(graph,seed,prob),
        LockdownType::Targeted => target_infection_clusters(graph, inflist),
        _ => unimplemented!()
    }


}
pub fn lockdown_naming_string(lockdowntype:LockdownType) -> String{
    match lockdowntype{
        LockdownType::CircuitBreaker => "CircuitBreaker".to_owned(),
        LockdownType::Random(_,_) => "Random".to_owned(),
        LockdownType::None => "None".to_owned(),
        LockdownType::Targeted => "Targeted".to_owned(),
        _ => unimplemented!()
    }

}

pub fn circuit_breaking_lockdown<T,A>(mut graph:GenericGraph<T,A>)-> GenericGraph<T,A> where A: net_ensembles::AdjContainer<T>{
    graph.clear_edges();
    
    graph

}

pub fn no_lockdown<T,A>(graph:GenericGraph<T,A>)-> GenericGraph<T,A> where A: net_ensembles::AdjContainer<T>{graph}

pub fn random_lockdown<InfectionState,A>(graph:GenericGraph<InfectionState,A>,seed:u64,prob:f64)-> GenericGraph<InfectionState,A> 
where A: net_ensembles::AdjContainer<InfectionState>, InfectionState: net_ensembles::Node, A: std::clone::Clone{
    //note to self. add some more enums and structs, eg like running the execute files in main to better capture the parameters of each of the individual lockdowns.__rust_force_expr!;
    
    
    let mut clone = graph.clone();
    let pairs = pair_finder(graph);
    let prob_dist= Uniform::new_inclusive(0.0,1.0);
    let mut rng = Pcg64::seed_from_u64(seed);

    //let prob = 0.6;

    for pair in pairs{

        let i1 = pair[0];
        let i2 = pair[1];
        //println!("{},{}",i1,i2);
        let p = prob_dist.sample(&mut rng);
        if p < prob{
            clone.remove_edge(i1,i2).unwrap();
        }

    }


    clone}



pub fn target_infection_clusters(graph:GenGraphSIR,inflist:&Vec<usize>) -> GenGraphSIR{

     //the goal of this will be to target those who are infected, and lockdown the surround clique of individuals. for instance, suppose 
    //an infected as three susceptible neighbours. those neighbours have their links removed to their neighbours (not the infected though)
    //this perhaps simulates the close-contact dynamics of the lockdown.
    //some credit is due to Peter Werner for discussing this with me.

    //can't quite get it to work yet.... get duplicates of certain infected when updating the graph in realtime
    let mut pairs:Vec<Vec<usize>> = Vec::new();
    let mut locked_down = graph.clone();

    for &ind in inflist{

        for (n_ind,_neigh) in graph.contained_iter_neighbors_with_index(ind).filter(|(_,neigh)|neigh.sus_check()){
            //This finds all the susceptible neighbours to the infected.

            //you made a mistake! leave edge between sus and infected.
            for (n_ind_2,_neigh_2) in graph.contained_iter_neighbors_with_index(n_ind).filter(|(_,neigh_2)|neigh_2.sus_check()){
                let mut v = vec![n_ind,n_ind_2];
                v.sort();
                pairs.push(v);

            }
        }

    }
    pairs.sort();
    pairs.dedup();
    for pair in pairs{
        locked_down.remove_edge(pair[0],pair[1]).unwrap();
    }



    locked_down
}























pub fn pair_finder<InfectionState,A>(mut graph:GenericGraph<InfectionState,A>)-> Vec<[usize;2]> 
where A: net_ensembles::AdjContainer<InfectionState>, InfectionState: net_ensembles::Node, A: std::clone::Clone{
    //This function finds all the pairs of indicies connected in the graph. 

    let mut vec:Vec<[usize;2]> = Vec::new();
    for j in 0..graph.vertex_count(){
        let adj = graph.contained_iter_neighbors_mut_with_index(j);
        for (n_index,_neighbor) in adj{
            //let mut new_vec = vec![j,n_index];
            if j < n_index{
                vec.push([j,n_index]);
            }

            //new_vec.sort_unstable();
            //vec.push(new_vec);
        } 
    }
    //vec.sort_unstable();
    //vec.dedup();
    vec
    //return iterator maybe, sort dedepu not needed if ordering is taken into account
    //size 2 slice instead
}

