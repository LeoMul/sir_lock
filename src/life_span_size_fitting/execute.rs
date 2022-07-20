use std::{ io::Write};

use indicatif::ProgressIterator;

use{
    super::*,
    crate::{GraphType, sir_model::*},
    serde_json::Value,
    std::{num::*, fs::File, io::BufWriter},
    crate::stats_methods::MyVariance,
    crate::life_span_data_collection::*,
    //crate::grid::*
    rayon::iter::{ParallelIterator},

};
use rand::Rng;
use{
    
    rand_pcg::Pcg64,
    net_ensembles::rand::SeedableRng,
    rayon::prelude::*,
};

use crate:: sir_model::base_model_options::BaseSwOptions;

pub fn run_simulation(param:LifespanSizeFittingParams, json: Value, num_threads: Option<NonZeroUsize>){
    match param.graph_type{
        GraphType::SmallWorld(_) => sim_small_world(param, json, num_threads),
        GraphType::Barabasi(_,_) => sim_barabasi(param, json, num_threads),
        _ => unimplemented!()
    }

}


fn sim_barabasi(param: LifespanSizeFittingParams, json: Value, num_threads:Option<NonZeroUsize>){
    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    rayon::ThreadPoolBuilder::new().num_threads(k.get()).build_global().unwrap();

    let mut n_size:Vec<_> = param.system_size_range.to_vec();
    println!("{:?}",n_size);
    let lambda_range = param.trans_prob_range.get_range();
    let lambda_vec:Vec<_> = lambda_range.iter().collect();
    let lockparams = param.lockdown;
    let mut sir_rng_2 = Pcg64::seed_from_u64(param.sir_seed);
    let mut graph_rng = Pcg64::seed_from_u64(param.graph_seed);


    let _data_vec:Vec<_> = n_size.iter_mut().map(|n|{
        println!("{}",n);
        
        
        

        //let mut container: Vec<_> = (0..k.get()).map(|_|{}).collect();

        let per_thread = param.num_networks/k.get() as u64;
        let bar = crate::indication_bar(param.num_networks);
        let mut rngs: Vec<_> = (0..k.get())
        .map(
            |_| 
                {
                    
                        (Pcg64::from_rng(&mut sir_rng_2).unwrap(),
                        Pcg64::from_rng(&mut graph_rng).unwrap())
                    
                }
            )
        .collect();
        //let bar = crate::indication_bar(per_thread);

        //this is a vector of matrices. there is a matrix for each thread, with each row being the lifespan for a particular network spanned over lambda.
        // 
        let intermediate_data_vecs: Vec<_> = rngs.par_iter_mut().map(|(rng,graph_rng)|{
            //let mut new_vec:Vec<_> = Vec::with_capacity(per_thread as usize);
            //let mut rng = Pcg64::seed_from_u64(param.sir_seed);

            let iter = (0..per_thread).into_iter().map(|_|{});
            let mut rng = rng;

            //new_vec has num_networks 
            let new_vec:Vec<Vec<u32>> = iter.map(|_|{
                let new_graph_seed = graph_rng.gen::<u64>();
                let opt = BarabasiOptions::from_life_span_fiting_params(&param,NonZeroUsize::new(*n).unwrap(),new_graph_seed);
                let small_world = opt.into();
                let mut model = SimpleSampleBarabasi::from_base(small_world, param.sir_seed);
                let lockgraph = model.create_locked_down_network(lockparams);

                let data_point_for_each_lambda:Vec<_> = lambda_vec.iter().map(|lambda|{
                    model.set_lambda(*lambda);
                    model.reseed_sir_rng(&mut rng);
                    model.propagate_until_completion_time_with_locks(lockgraph.clone(),lockparams)
                    }
                ).collect(); 
                data_point_for_each_lambda
            }).progress_with(bar.clone()).collect();

            new_vec
            }
            
        ).collect();

       let x = acquire_life_span_data_from_vec_of_matrices(intermediate_data_vecs, param.lifespanpercent);
       alternate_writing(&param, &json, num_threads, x, *n, &lambda_vec);

    }).collect();

}
fn sim_small_world(param: LifespanSizeFittingParams, json: Value, num_threads:Option<NonZeroUsize>){
    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    rayon::ThreadPoolBuilder::new().num_threads(k.get()).build_global().unwrap();

    let mut n_size:Vec<_> = param.system_size_range.to_vec();
    println!("{:?}",n_size);
    let lambda_range = param.trans_prob_range.get_range();
    let lambda_vec:Vec<_> = lambda_range.iter().collect();
    let lockparams = param.lockdown;
    let mut sir_rng_2 = Pcg64::seed_from_u64(param.sir_seed);
    let mut graph_rng = Pcg64::seed_from_u64(param.graph_seed);


    let _data_vec:Vec<_> = n_size.iter_mut().map(|n|{
        println!("{}",n);
        
        
        

        //let mut container: Vec<_> = (0..k.get()).map(|_|{}).collect();

        let per_thread = param.num_networks/k.get() as u64;
        let bar = crate::indication_bar(param.num_networks);
        let mut rngs: Vec<_> = (0..k.get())
        .map(
            |_| 
                {
                    
                        (Pcg64::from_rng(&mut sir_rng_2).unwrap(),
                        Pcg64::from_rng(&mut graph_rng).unwrap())
                    
                }
            )
        .collect();
        //let bar = crate::indication_bar(per_thread);

        //this is a vector of matrices. there is a matrix for each thread, with each row being the lifespan for a particular network spanned over lambda.
        // 
        let intermediate_data_vecs: Vec<_> = rngs.par_iter_mut().map(|(rng,graph_rng)|{
            //let mut new_vec:Vec<_> = Vec::with_capacity(per_thread as usize);
            //let mut rng = Pcg64::seed_from_u64(param.sir_seed);

            let iter = (0..per_thread).into_iter().map(|_|{});
            let mut rng = rng;

            //new_vec has num_networks 
            let new_vec:Vec<Vec<u32>> = iter.map(|_|{
                let new_graph_seed = graph_rng.gen::<u64>();
                let opt = BaseSwOptions::from_lifespan_size_fitting_param(&param,NonZeroUsize::new(*n).unwrap(),new_graph_seed);
                let small_world = opt.into();
                let mut model = SimpleSample::from_base(small_world, param.sir_seed);
                let lockgraph = model.create_locked_down_network(lockparams);

                let data_point_for_each_lambda:Vec<_> = lambda_vec.iter().map(|lambda|{
                    model.set_lambda(*lambda);
                    model.reseed_sir_rng(&mut rng);
                    model.propagate_until_completion_time_with_locks(lockgraph.clone(),lockparams)
                    }
                ).collect(); 
                data_point_for_each_lambda
            }).progress_with(bar.clone()).collect();

            new_vec
            }
            
        ).collect();

       let x = acquire_life_span_data_from_vec_of_matrices(intermediate_data_vecs, param.lifespanpercent);
       alternate_writing(&param, &json, num_threads, x, *n, &lambda_vec);

    }).collect();

}






fn transpose(matrix:&Vec<Vec<u32>>) -> Vec<Vec<u32>>{
    
    let num_rows = matrix.len();
    let num_cols = matrix[0].len();
    let mut new_matrix = Vec::with_capacity(num_cols);
    for j in 0..num_cols{
        let mut new_row = Vec::with_capacity(num_rows);
        for i in matrix.iter().take(num_rows){
            let data = &i;
            new_row.push(data[j]);
        }
        new_matrix.push(new_row);
    }
    new_matrix
    

}
fn decompose_thread_split_data(mut vecofmatrices:Vec<Vec<Vec<u32>>>) -> Vec<Vec<u32>>{
    let length = vecofmatrices.len();
    
    for i in 0..length{
        vecofmatrices[i] = transpose(&vecofmatrices[i])
    }

    let num_lams = vecofmatrices[0].len();
    let mut matrix_new = Vec::with_capacity(num_lams);
    for i in 0..num_lams{
        let mut new_row:Vec<u32> = Vec::new();
        for item in &vecofmatrices{
            let matrix = item;
            new_row.extend(&matrix[i]);
        }
        matrix_new.push(new_row);
    }
    matrix_new
}

fn acquire_life_span_data_from_vec_of_matrices(vecofmatrices:Vec<Vec<Vec<u32>>>,fraction:f64) -> Vec<u32>{
    let combined_data = decompose_thread_split_data(vecofmatrices);

    let mut lifespans = Vec::new();
    for item in &combined_data{
        let row = item;
        let mut new_row = row.clone();
        new_row.sort_unstable();
        let percent_life_span = acquire_percent_life_span(&new_row,fraction);
        lifespans.push(percent_life_span)
    }

    lifespans
}


















pub struct Measured
{
    pub var_t: MyVariance,
}


fn alternate_writing(param:&LifespanSizeFittingParams,json:&Value,num_threads: Option<NonZeroUsize>,data_master:Vec<u32>,n:usize,lambda:&Vec<f64>){
    
    
        let name = param.name("dat", num_threads,n);
        println!("creating: {name}");
        let file = File::create(name).expect("unable to create file");
        let mut buf = BufWriter::new(file);
        write!(buf, "#").unwrap();
        serde_json::to_writer(&mut buf, &json).unwrap();
        writeln!(buf).unwrap();
        // let data_vec = &data[n];
        writeln!(buf, "#lambda percentlifespan").unwrap();
        let data = &data_master;
        for k in 0..data.len(){
            
            writeln!(buf,"{} {}",lambda[k],data[k]).unwrap();
        }

    
    
    

}





