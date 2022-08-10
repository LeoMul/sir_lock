use{
    super::*,
    std::num::*,
    serde_json::Value,

};

use rayon::iter::{ParallelIterator};
use rayon::iter::IntoParallelRefMutIterator;
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
        GraphType::SmallWorld(_) => sim_small_world_percent_instead_of_chunk_size(param, json, num_threads),
        GraphType::Barabasi(_,_) => unimplemented!(),//sim_ba(param, json, num_threads),
        _ => unimplemented!()
    }

}
fn sim_small_world_percent_instead_of_chunk_size(param: ConnectedComponentParams,json:Value,num_threads:Option<NonZeroUsize>){
    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());

    let desired_percent = param.desired_percent_step;

    //let chunk_len = param.chunk_len;



    let mut graph_rng = Pcg64::seed_from_u64(param.graph_seed);

    //let iter = (0..param.num_networks).into_iter();
    let new_graph_seed = graph_rng.gen::<u64>();
    let opt = SWOptions::from_connectedcomponent_param(&param,new_graph_seed);
    let world:SWModel = opt.into(); 
    let graph = world.ensemble.graph().clone();
    let pairs = pair_finder(&graph);
    
    let res = (desired_percent*pairs.len() as f64).floor() as usize;
    let chunk_len = if res == 0{
            println!("Desired percent is too low. Setting chunk size to one.");
            1
        }
        else{
            res
    };
    let percent_removed_per_chunk = (chunk_len as f64)/(pairs.len() as f64);
    let start_index= (param.percent_start/percent_removed_per_chunk).floor() as usize;
    let end_index = (param.percent_end/percent_removed_per_chunk).ceil() as usize;
    
    //println!("start {start_index}, end {end_index}");

    let vec_length = end_index - start_index;
    println!("Start percent {}, end percent {}, desired increment {}",param.percent_start, param.percent_end,param.desired_percent_step);
    println!("Actual increment {} with {} data points sampled over {} networks.",percent_removed_per_chunk,vec_length,param.num_networks);
    println!("Note: data points are normalised by systemsize");

    let mut rngs: Vec<_> = (0..k.get())
        .map(
            |_| 
                {      
                        Pcg64::from_rng(&mut graph_rng).unwrap()
                }
            )
        .collect();

    
    

    let mut var_vec:Vec<f64> = vec![0.;vec_length];
    let mut avg_vec:Vec<f64> = var_vec.clone();
    let mut percent_vec = Vec::new();


    for p in start_index..end_index{
        percent_vec.push(p as f64 *percent_removed_per_chunk);
    }
    //percent_vec.push(1.);

    let per_thread = param.num_networks/k.get() as u64;
    let bar = crate::indication_bar(param.num_networks);

    let vector_of_core_vectors:Vec<_> = rngs.par_iter_mut().map(|rng|{

        let mut core_av_vec:Vec<f64> = vec![0.;vec_length as usize];
        let mut core_var_vec = core_av_vec.clone();
        //println!("test");
        let iter = 0..per_thread;
        iter.for_each(|_|{
            let new_graph_seed = rng.gen::<u64>();
            let opt = SWOptions::from_connectedcomponent_param(&param,new_graph_seed);
            let world:SWModel = opt.into(); 
            let mut graph = world.ensemble.graph().clone();

            let mut pairs = pair_finder(&graph);
            pairs.shuffle(rng);
            let chunked_vectors:Vec<&[_]> = pairs.chunks(chunk_len as usize).collect();
            //println!("j");

            for i in start_index..end_index{

                if i ==start_index{
                    //for j in 0..i{
                    //    let pairs_to_be_removed = chunked_vectors[j];
                    //    for pair in pairs_to_be_removed{
                    //        graph.remove_edge(pair[0],pair[1]).unwrap();
                    //    }
                //
//
                    //}
                    for item in chunked_vectors.iter().take(i){
                        for pair in *item{
                            graph.remove_edge(pair[0], pair[1]).unwrap();
                        }
                    }

                    let largest_connected_component = graph.connected_components()[0] as f64/param.system_size.get() as f64;
                    //println!("{}",largest_connected_component);
                    core_av_vec[i-start_index] += largest_connected_component ;
                    core_var_vec[i-start_index] += largest_connected_component*largest_connected_component;
                }
                else
                {


                let pairs_to_be_removed = chunked_vectors[i];
                for pair in pairs_to_be_removed{
                    graph.remove_edge(pair[0],pair[1]).unwrap();
                }
                //let num_edges = graph.edge_count();
                //println!("{num_edges}");
                let largest_connected_component = graph.connected_components()[0] as f64/param.system_size.get() as f64;
                    //println!("{}",largest_connected_component);
                    core_av_vec[i-start_index] += largest_connected_component ;
                    core_var_vec[i-start_index] += largest_connected_component*largest_connected_component;
                }
    
            }
            bar.inc(1);
        });
        (core_av_vec,core_var_vec)
    }).collect();
    bar.finish();
    println!("");
    for (averages,variances) in vector_of_core_vectors{
        for i in 0..averages.len(){
            avg_vec[i] += averages[i];
            var_vec[i] += variances[i];
        }
    }


    let actual_samples = per_thread*k.get() as u64;
    for i in 0..var_vec.len(){
        avg_vec[i] /= actual_samples as f64;
        var_vec[i] /= actual_samples as f64;
        var_vec[i] -= avg_vec[i] *avg_vec[i]
    }
    
    writing(&param,&json,num_threads,avg_vec,var_vec,percent_vec);
}
/* 
fn sim_ba(param: ConnectedComponentParams,json:Value,num_threads:Option<NonZeroUsize>){

    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    let chunk_len = param.chunk_len;
    let mut graph_rng = Pcg64::seed_from_u64(param.graph_seed);

    //let iter = (0..param.num_networks).into_iter();
    let new_graph_seed = graph_rng.gen::<u64>();
    let opt = BarabasiOptions::from_connectedcomponent_param(&param,new_graph_seed);
    let world:BarabasiModel = opt.into(); 
    let graph = world.ensemble.graph().clone();
    let pairs = pair_finder(&graph);



    

    let mut rngs: Vec<_> = (0..k.get())
        .map(
            |_| 
                {      
                        Pcg64::from_rng(&mut graph_rng).unwrap()
                }
            )
        .collect();

    let percent_removed_per_chunk = (chunk_len as f64)/(pairs.len() as f64);
    
    let start_index= (param.percent_start/percent_removed_per_chunk).floor() as usize;
    let end_index = (param.percent_end/percent_removed_per_chunk).ceil() as usize;
    
    println!("start {start_index}, end {end_index}");

    let vec_length = end_index - start_index;

    let mut var_vec:Vec<f64> = vec![0.;vec_length];
    let mut avg_vec:Vec<f64> = var_vec.clone();
    let mut percent_vec = Vec::new();


    for p in start_index..end_index{
        percent_vec.push(p as f64 *percent_removed_per_chunk);
    }
    //percent_vec.push(1.);

    let per_thread = param.num_networks/k.get() as u64;
    let bar = crate::indication_bar(param.num_networks);

    let vector_of_core_vectors:Vec<_> = rngs.par_iter_mut().map(|rng|{

        let mut core_av_vec:Vec<f64> = vec![0.;vec_length as usize];
        let mut core_var_vec = core_av_vec.clone();
        //println!("test");
        let iter = 0..per_thread;
        iter.for_each(|_|{
            let new_graph_seed = rng.gen::<u64>();
            let opt = BarabasiOptions::from_connectedcomponent_param(&param,new_graph_seed);
            let world:BarabasiModel = opt.into(); 
            let mut graph = world.ensemble.graph().clone();

            let mut pairs = pair_finder(&graph);
            pairs.shuffle(rng);
            let chunked_vectors:Vec<&[_]> = pairs.chunks(chunk_len as usize).collect();
            //println!("j");

            for i in start_index..end_index{

                if i ==start_index{

                    //for j in 0..i{
                    //    let pairs_to_be_removed = chunked_vectors[j];
                    //    for pair in pairs_to_be_removed{
                    //        graph.remove_edge(pair[0],pair[1]).unwrap();
                    //    }
                //
//
                    //}
                    for item in chunked_vectors.iter().take(i){
                        for pair in *item{
                            graph.remove_edge(pair[0], pair[1]).unwrap();
                        }
                    }

                    let largest_connected_component = graph.connected_components()[0];
                    //println!("{}",largest_connected_component);
                    core_av_vec[i-start_index] += largest_connected_component as f64;
                    core_var_vec[i-start_index] += largest_connected_component as f64 *largest_connected_component as f64;
                }
                else{
                    let pairs_to_be_removed = chunked_vectors[i];
                    for pair in pairs_to_be_removed{
                        graph.remove_edge(pair[0],pair[1]).unwrap();
                    }
                    //let num_edges = graph.edge_count();
                    //println!("{num_edges}");
                    let largest_connected_component = graph.connected_components()[0];
                    //println!("{}",largest_connected_component);
                    core_av_vec[i-start_index] += largest_connected_component as f64;
                    core_var_vec[i-start_index] += largest_connected_component as f64 *largest_connected_component as f64;
                }
    
            }
            bar.inc(1);
        });
        (core_av_vec,core_var_vec)
    }).collect();
    bar.finish();
    println!("");
    for (averages,variances) in vector_of_core_vectors{
        for i in 0..averages.len(){
            avg_vec[i] += averages[i];
            var_vec[i] += variances[i];
        }
    }


    let actual_samples = per_thread*k.get() as u64;
    for i in 0..var_vec.len(){
        avg_vec[i] /= actual_samples as f64;
        var_vec[i] /= actual_samples as f64;
        var_vec[i] -= avg_vec[i] *avg_vec[i]
    }
    
    writing(&param,&json,num_threads,avg_vec,var_vec,percent_vec);




}
fn sim_small_world_new(param: ConnectedComponentParams,json:Value,num_threads:Option<NonZeroUsize>){
    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    let chunk_len = param.chunk_len;
    let mut graph_rng = Pcg64::seed_from_u64(param.graph_seed);

    //let iter = (0..param.num_networks).into_iter();
    let new_graph_seed = graph_rng.gen::<u64>();
    let opt = SWOptions::from_connectedcomponent_param(&param,new_graph_seed);
    let world:SWModel = opt.into(); 
    let graph = world.ensemble.graph().clone();
    let pairs = pair_finder(&graph);



    

    let mut rngs: Vec<_> = (0..k.get())
        .map(
            |_| 
                {      
                        Pcg64::from_rng(&mut graph_rng).unwrap()
                }
            )
        .collect();

    let percent_removed_per_chunk = (chunk_len as f64)/(pairs.len() as f64);
    
    let start_index= (param.percent_start/percent_removed_per_chunk).floor() as usize;
    let end_index = (param.percent_end/percent_removed_per_chunk).ceil() as usize;
    
    println!("start {start_index}, end {end_index}");

    let vec_length = end_index - start_index;

    let mut var_vec:Vec<f64> = vec![0.;vec_length];
    let mut avg_vec:Vec<f64> = var_vec.clone();
    let mut percent_vec = Vec::new();


    for p in start_index..end_index{
        percent_vec.push(p as f64 *percent_removed_per_chunk);
    }
    //percent_vec.push(1.);

    let per_thread = param.num_networks/k.get() as u64;
    let bar = crate::indication_bar(param.num_networks);

    let vector_of_core_vectors:Vec<_> = rngs.par_iter_mut().map(|rng|{

        let mut core_av_vec:Vec<f64> = vec![0.;vec_length as usize];
        let mut core_var_vec = core_av_vec.clone();
        //println!("test");
        let iter = 0..per_thread;
        iter.for_each(|_|{
            let new_graph_seed = rng.gen::<u64>();
            let opt = SWOptions::from_connectedcomponent_param(&param,new_graph_seed);
            let world:SWModel = opt.into(); 
            let mut graph = world.ensemble.graph().clone();

            let mut pairs = pair_finder(&graph);
            pairs.shuffle(rng);
            let chunked_vectors:Vec<&[_]> = pairs.chunks(chunk_len as usize).collect();
            //println!("j");

            for i in start_index..end_index{

                if i ==start_index{
                    //for j in 0..i{
                    //    let pairs_to_be_removed = chunked_vectors[j];
                    //    for pair in pairs_to_be_removed{
                    //        graph.remove_edge(pair[0],pair[1]).unwrap();
                    //    }
                //
//
                    //}
                    for item in chunked_vectors.iter().take(i){
                        for pair in *item{
                            graph.remove_edge(pair[0], pair[1]).unwrap();
                        }
                    }

                    let largest_connected_component = graph.connected_components()[0];
                    //println!("{}",largest_connected_component);
                    core_av_vec[i-start_index] += largest_connected_component as f64;
                    core_var_vec[i-start_index] += largest_connected_component as f64 *largest_connected_component as f64;
                }
                else
                {


                let pairs_to_be_removed = chunked_vectors[i];
                for pair in pairs_to_be_removed{
                    graph.remove_edge(pair[0],pair[1]).unwrap();
                }
                //let num_edges = graph.edge_count();
                //println!("{num_edges}");
                let largest_connected_component = graph.connected_components()[0];
                //println!("{}",largest_connected_component);
                core_av_vec[i-start_index] += largest_connected_component as f64;
                core_var_vec[i-start_index] += largest_connected_component as f64 *largest_connected_component as f64;
                }
    
            }
            bar.inc(1);
        });
        (core_av_vec,core_var_vec)
    }).collect();
    bar.finish();
    println!("");
    for (averages,variances) in vector_of_core_vectors{
        for i in 0..averages.len(){
            avg_vec[i] += averages[i];
            var_vec[i] += variances[i];
        }
    }


    let actual_samples = per_thread*k.get() as u64;
    for i in 0..var_vec.len(){
        avg_vec[i] /= actual_samples as f64;
        var_vec[i] /= actual_samples as f64;
        var_vec[i] -= avg_vec[i] *avg_vec[i]
    }
    
    writing(&param,&json,num_threads,avg_vec,var_vec,percent_vec);




}
*/


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
