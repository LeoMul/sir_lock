use net_ensembles::*;
use serde::{Serialize, Deserialize};
//use crate::sir_model::sir_states::*;
use rand::seq::SliceRandom;
use rand::*;
use rand::distributions::*;
use rand_pcg::Pcg64;
use net_ensembles::rand::prelude::Distribution;
//use crate::misc_types::*;
#[derive(Debug, Clone,Serialize, Deserialize,Copy)]

pub enum LockdownType{
    CircuitBreaker,
    None,
    Random(f64),
    Targeted,
    LimitContacts(usize),
    Invalid

}

pub type GenGraphSIR = net_ensembles::GenericGraph<crate::sir_model::sir_states::InfectionState, net_ensembles::graph::NodeContainer<crate::sir_model::sir_states::InfectionState>>;
pub type SwSIR = net_ensembles::GenericGraph<crate::sir_model::sir_states::InfectionState, net_ensembles::sw_graph::SwContainer<crate::sir_model::sir_states::InfectionState>>;


#[derive(Debug, Clone,Serialize, Deserialize,Copy)]
pub struct LockdownParameters{
    pub lock_style: LockdownType,
    pub release_threshold: f64,
    pub lock_threshold: f64
}

impl LockdownParameters{
    pub fn set_lock_thresh(&mut self,lock_thresh:f64){
        self.lock_threshold = lock_thresh;

    }
    pub fn set_rel_thresh(&mut self,rel_thresh:f64){
        self.release_threshold = rel_thresh;

    }
}

#[derive(Debug, Clone,Serialize, Deserialize)]
pub struct LockdownPairStorage{
    pub to_be_removed: Vec<[usize;2]>,
    pub to_be_kept: Vec<[usize;2]>
}
pub fn lockdown<InfectionState,A,R>(lockdownparams: LockdownParameters,graph:GenericGraph<InfectionState,A>,rng:R) 
-> GenericGraph<InfectionState,A> 
    where 
    A: net_ensembles::AdjContainer<InfectionState>, 
    InfectionState: net_ensembles::Node, 
    A: std::clone::Clone, 
    R:Rng
{
    let pairs_struct = create_lock_pairs_lists(lockdownparams, &graph,rng);
    create_locked_down_network_from_pair_list(&pairs_struct, &graph)

}
pub fn create_lock_pairs_lists<InfectionState,A,R>(lockdownparams: LockdownParameters,graph:&GenericGraph<InfectionState,A>,rng:R) 
-> LockdownPairStorage  where A: net_ensembles::AdjContainer<InfectionState>, InfectionState: net_ensembles::Node, A: std::clone::Clone, R: Rng {
    match lockdownparams.lock_style{
        LockdownType::Random(prob) => create_pair_lists_for_random(graph, prob,rng),
        LockdownType::CircuitBreaker => create_pair_lists_for_circuitbreaker(graph),
        LockdownType::None => create_pairs_lists_for_nolock(graph),
        _ => unimplemented!()
    }
}

pub fn create_locked_down_network_from_pair_list<InfectionState,A>(pairs:&LockdownPairStorage, graph:&GenericGraph<InfectionState,A>) -> GenericGraph<InfectionState,A>
where A: net_ensembles::AdjContainer<InfectionState>, InfectionState: net_ensembles::Node, A: std::clone::Clone{
    //I think the lockdowns will be reformulated to follow this pattern.
    let pairs_to_be_removed = &pairs.to_be_removed;
    let mut clone = graph.clone();
    for pair in pairs_to_be_removed{

        let i1 = pair[0];
        let i2 = pair[1];
        //println!("{},{}",i1,i2);
        clone.remove_edge(i1,i2).unwrap();
    

    }
    clone
    
}

pub fn create_pairs_lists_for_nolock<T,A>(graph:&GenericGraph<T,A>) -> LockdownPairStorage 
where A: net_ensembles::AdjContainer<T>, T: net_ensembles::Node, A: std::clone::Clone {
    let pairs_to_be_kept = pair_finder(graph);
    LockdownPairStorage{
        to_be_removed: Vec::new(),
        to_be_kept: pairs_to_be_kept
    }
}


pub fn create_pair_lists_for_circuitbreaker<T,A>(graph:&GenericGraph<T,A>) ->LockdownPairStorage
where A: net_ensembles::AdjContainer<T>, T: net_ensembles::Node, A: std::clone::Clone {
    let pairs = pair_finder(graph);
    LockdownPairStorage{
        to_be_removed:pairs,
        to_be_kept: Vec::new()
    }
}

pub fn create_pair_lists_for_random<T,A,R>(graph:&GenericGraph<T,A>,prob:f64,mut rng:R)-> LockdownPairStorage
where A: net_ensembles::AdjContainer<T>, T: net_ensembles::Node, A: std::clone::Clone, 
    R: Rng{

    let mut pairs = pair_finder(graph);
    let amount = (pairs.len()as f64 *prob).ceil() as usize;
    //let mut rng = Pcg64::seed_from_u64(seed);
    //let (shuffled,not_shuffled) = pairs.partial_shuffle(&mut rng, amount);
    //unimplemented!();
    pairs.shuffle(&mut rng);
    let to_be_kept:Vec<_> = pairs.drain(amount..).collect();
    //This computation only happens once so happy to have it be done this way,
    LockdownPairStorage{
        to_be_removed: pairs,
        to_be_kept
    }

}

pub fn pair_finder<InfectionState,A>(graph:&GenericGraph<InfectionState,A>)-> Vec<[usize;2]> 
where A: net_ensembles::AdjContainer<InfectionState>, InfectionState: net_ensembles::Node, A: std::clone::Clone{
    //This function finds all the pairs of indicies connected in the graph. 
    let clone = graph;
    let mut vec:Vec<[usize;2]> = Vec::new();
    for j in 0..graph.vertex_count(){
        let adj = clone.contained_iter_neighbors_with_index(j);
        for (n_index,_neighbor) in adj{
            //let mut new_vec = vec![j,n_index];
            if j < n_index{
                vec.push([j,n_index]);
            }

            
        } 
    }

    vec
    //return iterator maybe
}
pub fn lockdown_naming_string(lockdowntype:LockdownType) -> String{

    
    match lockdowntype{
        LockdownType::CircuitBreaker => "CircuitBreaker".to_owned(),
        LockdownType::Random(_) => "Random".to_owned(),
        LockdownType::None => "None".to_owned(),
        LockdownType::Targeted => "Targeted".to_owned(), //retired
        LockdownType::LimitContacts(_) => "LimitContacts".to_owned(), //retired
        _ => unimplemented!()
    }

}
//retired
pub fn limit_contacts<InfectionState,A>(mut graph:GenericGraph<InfectionState,A>,max_new_edges:usize, rng:&mut Pcg64)-> GenericGraph<InfectionState,A> 
where A: net_ensembles::AdjContainer<InfectionState>, InfectionState: net_ensembles::Node, A: std::clone::Clone{
    
    let mut index_list:Vec<usize> = (0..graph.vertex_count()).collect();
    index_list.shuffle(rng);

    for j in index_list{

        if graph.degree(j).unwrap() > max_new_edges{
            let adj = graph.contained_iter_neighbors_mut_with_index(j);
            let mut vec = Vec::new();
            for (ind,_neighbour) in adj{
                vec.push(ind);
            }
            
            while graph.degree(j).unwrap() > max_new_edges{
                let un = Uniform::new(0,graph.degree(j).unwrap());  
                let rand_index = un.sample(rng);
                graph.remove_edge(j,vec[rand_index]).unwrap();
                vec.remove(rand_index);
            }
            

        }
        
    }


    graph

}
