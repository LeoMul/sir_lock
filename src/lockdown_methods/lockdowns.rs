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
    Random(u64,f64),
    Targeted,
    LimitContacts(usize),
    Invalid

}

#[derive(Debug, Clone,Serialize, Deserialize,Copy)]

pub enum GraphToBeLockedDownType{
    SmallWorldToBeLockedDown,
    OtherToBeLockedDown,
    Invalid}


pub type GenGraphSIR = net_ensembles::GenericGraph<crate::sir_model::sir_states::InfectionState, net_ensembles::graph::NodeContainer<crate::sir_model::sir_states::InfectionState>>;
pub type SwSIR = net_ensembles::GenericGraph<crate::sir_model::sir_states::InfectionState, net_ensembles::sw_graph::SwContainer<crate::sir_model::sir_states::InfectionState>>;




#[derive(Debug, Clone,Serialize, Deserialize,Copy)]
pub struct LockdownParameters{
    pub lock_style: LockdownType,
    pub release_threshold: f64,
    pub lock_threshold: f64,
    pub dynamic_bool: bool,
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



pub fn lockdown<InfectionState,A>(lockdownparams: LockdownParameters,graph:GenericGraph<InfectionState,A>, rng:&mut Pcg64) 
-> GenericGraph<InfectionState,A> where A: net_ensembles::AdjContainer<InfectionState>, InfectionState: net_ensembles::Node, A: std::clone::Clone{
    match lockdownparams.lock_style{
        LockdownType::CircuitBreaker => circuit_breaking_lockdown(graph),
        LockdownType::None => no_lockdown(graph),
        LockdownType::Random(seed,prob) => random_lockdown(graph,seed,prob),
        //LockdownType::Targeted => target_infection_clusters(&graph, inflist),
        LockdownType::LimitContacts(m) => limit_contacts(graph, m, rng),
        _ => unimplemented!()
    }

}



//reformulating how lockdowns are calculated.
pub fn create_lock_pairs_lists(lockdownparams:LockdownParameters,graph:&GenGraphSIR) -> LockdownPairStorage{
    match lockdownparams.lock_style{
        LockdownType::Random(seed,prob) => create_pair_lists_for_random(graph,seed,prob),
        _ => unimplemented!()
    }
}


pub fn lockdown_naming_string(lockdowntype:LockdownType) -> String{

    
    match lockdowntype{
        LockdownType::CircuitBreaker => "CircuitBreaker".to_owned(),
        LockdownType::Random(_,_) => "Random".to_owned(),
        LockdownType::None => "None".to_owned(),
        LockdownType::Targeted => "Targeted".to_owned(),
        LockdownType::LimitContacts(_) => "LimitContacts".to_owned(),
        _ => unimplemented!()
    }

}

pub fn create_locked_down_network_from_pair_list(pairs:&LockdownPairStorage, graph:&GenGraphSIR) -> GenGraphSIR{
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



pub fn circuit_breaking_lockdown<T,A>(mut graph:GenericGraph<T,A>)-> GenericGraph<T,A> where A: net_ensembles::AdjContainer<T>{
    graph.clear_edges();
    
    graph

}


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

pub fn no_lockdown<T,A>(graph:GenericGraph<T,A>)-> GenericGraph<T,A> where A: net_ensembles::AdjContainer<T>{graph}





pub fn create_pair_lists_for_random(graph:&GenGraphSIR,seed:u64,prob:f64)-> LockdownPairStorage{
    let mut pairs = pair_finder(graph);
    let amount = (pairs.len()as f64 *prob).ceil() as usize;
    let mut rng = Pcg64::seed_from_u64(seed);
    let (shuffled,not_shuffled) = pairs.partial_shuffle(&mut rng, amount);

    //This computation only happens once so happy to have it be done this way,
    LockdownPairStorage{
        to_be_removed: shuffled.to_vec(),
        to_be_kept: not_shuffled.to_vec()
    }

}


pub fn random_lockdown<T,A>(graph:GenericGraph<T,A>,seed:u64,prob:f64)-> GenericGraph<T,A> 
where A: net_ensembles::AdjContainer<T>, T: net_ensembles::Node, A: std::clone::Clone{
    //note to self. add some more enums and structs, eg like running the execute files in main to better capture the parameters of each of the individual lockdowns.__rust_force_expr!;
    
    
    let mut clone = graph.clone();
    let mut pairs = pair_finder(&graph);
    
    let amount = (pairs.len()as f64 *prob).ceil() as usize;
    //let prob_dist= Uniform::new_inclusive(0.0,1.0);
    let mut rng = Pcg64::seed_from_u64(seed);
    let (rando,_) = pairs.partial_shuffle(&mut rng,amount);

    for pair in rando{

        let i1 = pair[0];
        let i2 = pair[1];
        //println!("{},{}",i1,i2);
        clone.remove_edge(i1,i2).unwrap();
    

    }


    clone
    }

//pub type Temp = crate::sir_model::sir_states::InfectionState;
//
//pub fn target_infection_clusters<Temp,A>(graph:&GenericGraph<Temp,A>,inflist:&Vec<usize>) -> GenericGraph<Temp,A>
//where A: net_ensembles::AdjContainer<Temp>, Temp: net_ensembles::Node, A: std::clone::Clone{
//
//     //the goal of this will be to target those who are infected, and lockdown the surround clique of individuals. for instance, suppose 
//    //an infected as three susceptible neighbours. those neighbours have their links removed to their neighbours (not the infected though)
//    //this perhaps simulates the close-contact dynamics of the lockdown.
//    //some credit is due to Peter Werner for discussing this with me.
//
//    //can't quite get it to work yet.... get duplicates of certain infected when updating the graph in realtime
//    let mut pairs:Vec<(usize,usize)> = Vec::new();
//    let mut locked_down = graph.clone();
//
//    for &ind in inflist{
//
//        for (n_ind,_neigh) in graph.contained_iter_neighbors_with_index(ind).filter(|(_,neigh)|neigh.sus_check()){
//            //This finds all the susceptible neighbours to the infected.
//
//            //you made a mistake! leave edge between sus and infected.
//            for (n_ind_2,_neigh_2) in graph.contained_iter_neighbors_with_index(n_ind).filter(|(_,neigh_2)|neigh_2.sus_check()){
//                
//                if n_ind < n_ind_2{
//                    pairs.push((n_ind,n_ind_2));
//                }
//
//            }
//        }
//
//    }
//    pairs.sort_unstable();
//    pairs.dedup();
//    for pair in pairs{
//        locked_down.remove_edge(pair.0,pair.1).unwrap();
//    }
//
//
//
//    locked_down
//}
//
//
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

