

use rayon::iter::{ ParallelIterator};
use{
    crate::{ sir_model::*},
    std::{num::*},
    rand_pcg::Pcg64,
    net_ensembles::rand::SeedableRng,
    rayon::prelude::*,
    crate::lockdown_methods::*,
    crate::lifespanhist::parser::LifeSpanParams
    
};
use indicatif::*;

use net_ensembles::sampling::histogram::*;

pub fn acquire_sorted_data(param: &LifeSpanParams,k:NonZeroUsize,sir_seed:u64,samples:u64,lockdown:LockdownParameters)-> Vec<u32> {

    

    let mut rng = Pcg64::seed_from_u64(sir_seed);

    let mut container: Vec<_> = (0..k.get()).map(
        |_|
        {
            Pcg64::from_rng(&mut rng).unwrap()
        }
    ).collect();

    let per_threads = samples/k.get() as u64;

    let bar = crate::indication_bar(samples);

    let mut data_chunked:Vec<Vec<u32>> = container.par_iter_mut().map(|r|{

        let vec_per_thread:Vec<_> = (0..per_threads).into_iter().map(|_|{
            let opt = SWOptions::from_lifespan_param(&param);
            let small_world = opt.into();
            let mut model = SimpleSampleSW::from_base(small_world, param.sir_seed,param.initial_infected);
            model.reseed_sir_rng(r);
            //let mut lockgraph = model.create_locked_down_network(lockdown);
           model.propagate_until_completion_time_with_locks_new_lockgraph_for_each_lockdown(lockdown)
        }).progress_with(bar.clone()).collect();
        vec_per_thread

    }
    ).collect();
    //data is a vector of four vectors. we now decompose this into a single vector.
    let mut data:Vec<u32> = Vec::new();
    
    for j in &mut data_chunked{
        data.append(j);
    }
    //let max = data.iter().max();

    data.sort_unstable();
    data
    


}

pub fn convert_sorted_to_hist(data:&Vec<u32>) -> Vec<(u32, usize)>{
    //Converts already sorted

    //data.sort_unstable();
    //println!("{:?}",data);
    let mut hist = HistU32Fast::new_inclusive(0,data[data.len()-1]).unwrap();

    for item in data{
        hist.count_val(item).unwrap();
    }
    let vector_data: Vec<(u32, usize)> =  hist.bin_hits_iter().collect();
    vector_data

}


pub fn acquire_percent_life_span(sorted_data:&Vec<u32>,fraction:f64)-> u32{
    
    
    let index = (sorted_data.len() as f64 * fraction).round() as usize;
    sorted_data[index]


}