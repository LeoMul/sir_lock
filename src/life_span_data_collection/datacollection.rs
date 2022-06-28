

use rayon::iter::{ ParallelIterator};
use{
    crate::{ sir_model::*},
    std::{num::*},
    rand_pcg::Pcg64,
    net_ensembles::rand::SeedableRng,
    rayon::prelude::*
};

use net_ensembles::sampling::histogram::*;

pub fn acquire_sorted_data(model:SimpleSample,k:NonZeroUsize,sir_seed:u64,samples:u64)-> Vec<u32> {

    

    let mut rng = Pcg64::seed_from_u64(sir_seed);

    let mut container: Vec<_> = (0..k.get()).map(
        |_|
        {
            let mut model = model.clone();
            model.reseed_sir_rng(&mut rng);
            model
        }
    ).collect();

    let per_threads = samples/k.get() as u64;

    
    let mut data_chunked:Vec<Vec<u32>> = container.par_iter_mut().map(|model|{
        let mut vec_per_thread:Vec<u32> = Vec::new();
        for _i in 0..per_threads{
            let time = model.propagate_until_completion_time();
            vec_per_thread.push(time)
        }
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