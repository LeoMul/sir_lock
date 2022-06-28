use std::{ io::Write};

//***************READ ME***************** */
//I had intended to get some kind of infections vs time graph, after discussing this with yannick maybe we need to take more care in our sampling. perhaps we can find an average functional??

use rayon::iter::{ ParallelIterator};



use{
    super::*,
    crate::{GraphType, sir_model::*},
    serde_json::Value,
    std::{num::*, fs::File, io::BufWriter},
    rand_pcg::Pcg64,
    net_ensembles::rand::SeedableRng,
    rayon::prelude::*,
    crate::stats_methods::MyVariance,
    crate::grid::*,
   
};

pub fn run_simulation(param:TimeGraphParams, json: Value, num_threads: Option<NonZeroUsize>){
    match param.graph_type{
        GraphType::SmallWorld(_) => sim_small_world(param, json, num_threads),
        _ => unimplemented!()
    }

}

fn sim_small_world(param: TimeGraphParams, json: Value, num_threads:Option<NonZeroUsize>){
    let opt = BaseSwOptions::from_time_graph_param(&param);
    let small_world = opt.into();
    let model = SimpleSample::from_base(small_world, param.sir_seed);

    let system_size_fraction = if param.fraction{
        param.system_size.get() as f64
    }else{
        1.
    };
    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    rayon::ThreadPoolBuilder::new().num_threads(k.get()).build_global().unwrap();

    let mut rng = Pcg64::seed_from_u64(param.sir_seed);

    let mut container: Vec<_> = (0..k.get()).map(
        |_|
        {
            let mut model = model.clone();
            model.reseed_sir_rng(&mut rng);
            model
        }
    ).collect();

    let per_threads = param.samples/k.get() as u64;

}
