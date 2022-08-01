use{
    super::*,
    std::num::*,
    serde_json::Value,

};
use std::{ io::Write};
use rand::seq::SliceRandom;
use rand::Rng;
//use serde_json::from_slice;

use{
    crate::{GraphType, sir_model::*},
    std::{ fs::File, io::BufWriter},
    rand_pcg::Pcg64,
    net_ensembles::rand::SeedableRng,
    crate::lockdown_methods::pair_finder,
};
use net_ensembles::{rand};
use net_ensembles::WithGraph;

pub fn run_simulation(param:ConnectedComponentParams, json: Value, num_threads: Option<NonZeroUsize>){
    match param.graph_type{
        GraphType::SmallWorld(_) => sim_small_world(param, json, num_threads),
        GraphType::Barabasi(_,_) => sim_ba(param, json, num_threads),
        _ => unimplemented!()
    }

}

fn sim_small_world(param: ConnectedComponentParams,json:Value,num_threads:Option<NonZeroUsize>){
    let num_chunks = param.num_chunks;  
    //let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    let mut var_vec:Vec<f64> = vec![0.;num_chunks as usize];
    let mut avg_vec:Vec<f64> = var_vec.clone();
    let mut percent_vec = Vec::new();
    let mut graph_rng = Pcg64::seed_from_u64(param.graph_seed);

    //let iter = (0..param.num_networks).into_iter();

    for k in 0..param.num_networks{
        let new_graph_seed = graph_rng.gen::<u64>();
        let opt = SWOptions::from_connectedcomponent_param(&param,new_graph_seed);
        let world:SWModel = opt.into(); 
        let mut graph = world.ensemble.graph().clone();

        let mut pairs = pair_finder(&graph);
        pairs.shuffle(&mut graph_rng);

        

        //let chunk_size = pairs.len()/num_chunks as usize;

        let chunk_len = (pairs.len() as f64/num_chunks as f64).ceil() as usize;

        let chunked_vectors:Vec<&[_]> = pairs.chunks(chunk_len).collect();
        if k == 0{
            let percent_removed_per_chunk = (chunk_len as f64)/(pairs.len() as f64);
            for p in 0..(chunked_vectors.len()-1){
                percent_vec.push(p as f64 *percent_removed_per_chunk);
            }
            percent_vec.push(1.);
        }
        //let vec = graph.connected_components();
        assert_eq!(chunked_vectors.len(),num_chunks as usize);

        for i in 0..chunked_vectors.len(){

            let pairs_to_be_removed = chunked_vectors[i];
            for pair in pairs_to_be_removed{
                graph.remove_edge(pair[0],pair[1]).unwrap();
            }
            //let num_edges = graph.edge_count();
            //println!("{num_edges}");
            let largest_connected_component = graph.connected_components()[0];

            avg_vec[i] += largest_connected_component as f64;
            var_vec[i] += largest_connected_component as f64 *largest_connected_component as f64;

        }

    };

    for i in 0..var_vec.len(){
        avg_vec[i] /= param.num_networks as f64;
        var_vec[i] /= param.num_networks as f64;
        var_vec[i] -= avg_vec[i] *avg_vec[i]
    }
    writing(&param,&json,num_threads,avg_vec,var_vec,percent_vec);




}
fn writing(param:&ConnectedComponentParams,json:&Value,num_threads: Option<NonZeroUsize>,avg:Vec<f64>,var:Vec<f64>,percentages:Vec<f64>){
    let name = param.name("dat", num_threads);
    println!("creating: {name}");
    let file = File::create(name).expect("unable to create file");
    let mut buf = BufWriter::new(file);
    write!(buf, "#").unwrap();
    serde_json::to_writer(&mut buf, &json).unwrap();
    writeln!(buf).unwrap();
    // let data_vec = &data[n];
    
    //let data = &data_master[j];
    writeln!(buf,"percent_removed average_largest var_largest").unwrap();

    for k in 0..avg.len(){
        
        writeln!(buf,"{} {} {}",percentages[k],avg[k],var[k]).unwrap();
    }
    
}
fn sim_ba(param: ConnectedComponentParams,json:Value,num_threads:Option<NonZeroUsize>){
    let num_chunks = param.num_chunks;  
    //let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    let mut var_vec:Vec<f64> = vec![0.;num_chunks as usize];
    let mut avg_vec:Vec<f64> = var_vec.clone();
    let mut percent_vec = Vec::new();
    let mut graph_rng = Pcg64::seed_from_u64(param.graph_seed);

    //let iter = (0..param.num_networks).into_iter();

    for k in 0..param.num_networks{
        let new_graph_seed = graph_rng.gen::<u64>();
        let opt = SWOptions::from_connectedcomponent_param(&param,new_graph_seed);
        let world:SWModel = opt.into(); 
        let mut graph = world.ensemble.graph().clone();

        let mut pairs = pair_finder(&graph);
        pairs.shuffle(&mut graph_rng);

        

        //let chunk_size = pairs.len()/num_chunks as usize;

        let chunk_len = (pairs.len() as f64/num_chunks as f64).ceil() as usize;

        let chunked_vectors:Vec<&[_]> = pairs.chunks(chunk_len).collect();
        if k == 0{
            let percent_removed_per_chunk = (chunk_len as f64)/(pairs.len() as f64);
            for p in 0..(chunked_vectors.len()-1){
                percent_vec.push(p as f64 *percent_removed_per_chunk);
            }
            percent_vec.push(1.);
        }
        //let vec = graph.connected_components();
        assert_eq!(chunked_vectors.len(),num_chunks as usize);

        for i in 0..chunked_vectors.len(){

            let pairs_to_be_removed = chunked_vectors[i];
            for pair in pairs_to_be_removed{
                graph.remove_edge(pair[0],pair[1]).unwrap();
            }
            //let num_edges = graph.edge_count();
            //println!("{num_edges}");
            let largest_connected_component = graph.connected_components()[0];

            avg_vec[i] += largest_connected_component as f64;
            var_vec[i] += largest_connected_component as f64 *largest_connected_component as f64;

        }

    };

    for i in 0..var_vec.len(){
        avg_vec[i] /= param.num_networks as f64;
        var_vec[i] /= param.num_networks as f64;
        var_vec[i] -= avg_vec[i] *avg_vec[i]
    }
    writing(&param,&json,num_threads,avg_vec,var_vec,percent_vec);




}