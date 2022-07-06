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

use{
    
    rand_pcg::Pcg64,
    net_ensembles::rand::SeedableRng,
    rayon::prelude::*,
};

use crate:: sir_model::base_model_options::BaseSwOptions;

pub fn run_simulation(param:LifespanSizeFittingParams, json: Value, num_threads: Option<NonZeroUsize>){
    match param.graph_type{
        GraphType::SmallWorld(_) => sim_small_world(param, json, num_threads),
        GraphType::Barabasi(_,_) => sim_barabasi_new(param, json, num_threads),
        _ => unimplemented!()
    }

}


fn sim_barabasi_new(param: LifespanSizeFittingParams, json: Value, num_threads:Option<NonZeroUsize>){
    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    rayon::ThreadPoolBuilder::new().num_threads(k.get()).build_global().unwrap();

    let mut n_size:Vec<_> = param.system_size_range.iter().map(|n|{*n}).collect();
    println!("{:?}",n_size);
    let lambda_range = param.trans_prob_range.get_range();
    let lambda_vec:Vec<_> = lambda_range.iter().collect();
    let lockparams = param.lockdown;



    let data_vec:Vec<_> = n_size.iter_mut().map(|n|{
        println!("{}",n);
        
        
        

        let mut container: Vec<_> = (0..k.get()).map(|_|{}).collect();

        let per_thread = param.num_networks/k.get() as u64;
        
        //let bar = crate::indication_bar(per_thread);

        //this is a vector of matrices. there is a matrix for each thread, with each row being the lifespan for a particular network spanned over lambda.
        // 
        let intermediate_data_vecs: Vec<_> = container.par_iter_mut().map(|_|{
            //let mut new_vec:Vec<_> = Vec::with_capacity(per_thread as usize);
            let mut rng = Pcg64::seed_from_u64(param.sir_seed);

            let iter = (0..per_thread).into_iter().map(|_|{});
            

            //new_vec has num_networks 
            let new_vec:Vec<Vec<u32>> = iter.map(|_|{
                let opt = BarabasiOptions::from_life_span_fiting_params(&param,NonZeroUsize::new(*n).unwrap());
                let small_world = opt.into();
                let mut model = SimpleSampleBarabasi::from_base(small_world, param.sir_seed);
                let lockgraph = model.create_locked_down_network(lockparams);

                let data_point_for_each_lambda:Vec<_> = lambda_vec.iter().map(|lambda|{
                    model.set_lambda(*lambda);
                    model.reseed_sir_rng(&mut rng);
                    let t = model.propagate_until_completion_time_with_locks(lockgraph.clone(),lockparams);
                    t }
                ).collect(); 
                data_point_for_each_lambda
            }).collect();

            new_vec
            }
            
        ).collect();

       let lifespans = acquire_life_span_data_from_vec_of_matrices(intermediate_data_vecs, param.lifespanpercent);
        lifespans
    }).collect();
    alternate_writing(param, json, num_threads, data_vec, n_size, lambda_vec);

}

fn transpose(matrix:&Vec<Vec<u32>>) -> Vec<Vec<u32>>{
    
    let num_rows = matrix.len();
    let num_cols = matrix[0].len();
    let mut new_matrix = Vec::with_capacity(num_cols);
    for j in 0..num_cols{
        let mut new_row = Vec::with_capacity(num_rows);
        for i in 0..num_rows{
            let data = &matrix[i];
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
        for j in 0..vecofmatrices.len(){
            let matrix = &vecofmatrices[j];
            new_row.extend(&matrix[i]);
        }
        matrix_new.push(new_row);
    }
    matrix_new
}

fn acquire_life_span_data_from_vec_of_matrices(vecofmatrices:Vec<Vec<Vec<u32>>>,fraction:f64) -> Vec<u32>{
    let combined_data = decompose_thread_split_data(vecofmatrices);

    let mut lifespans = Vec::new();
    for i in 0..combined_data.len(){
        let row = &combined_data[i];
        let mut new_row = row.clone();
        new_row.sort_unstable();
        let percent_life_span = acquire_percent_life_span(&new_row,fraction);
        lifespans.push(percent_life_span)
    }

    lifespans
}





fn _sim_barabasi(param: LifespanSizeFittingParams, json: Value, num_threads:Option<NonZeroUsize>){

    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    rayon::ThreadPoolBuilder::new().num_threads(k.get()).build_global().unwrap();

    let mut n_size:Vec<_> = param.system_size_range.iter().map(|n|{*n}).collect();
    println!("{:?}",n_size);
    let lambda_range = param.trans_prob_range.get_range();
    let lambda_vec:Vec<_> = lambda_range.iter().collect();
    let lockparams = param.lockdown;
    let rng = Pcg64::seed_from_u64(param.sir_seed);


    let data_vec:Vec<_> = n_size.iter_mut().map(|n|{
        println!("{}",n);
        let bar = crate::indication_bar(lambda_vec.len() as u64);

        let n_lambdaspanned:Vec<_> = lambda_vec.iter().map(|lambda|{

            let per_thread = param.num_networks/k.get() as u64;

            let mut container: Vec<_> = (0..k.get()).map(
                |_|
                {
                    //need this to avoid recreating rng several thousand times. this way ony num_thread times.
                    rng.clone()
                    
                }
            ).collect();

            //this is the vector of lifespans for the given (N,lambda)
            let mut time_vec:Vec<_> = container.par_iter_mut().flat_map(|mut r|{
                
                let iterator = 0..per_thread;

                let vec_temp:Vec<_> = iterator.map(|_|{
                    let opt = BarabasiOptions::from_life_span_fiting_params(&param,NonZeroUsize::new(*n).unwrap());
                    let small_world = opt.into();
                    let mut model = SimpleSampleBarabasi::from_base(small_world, param.sir_seed);
                    let lockgraph = model.create_locked_down_network(lockparams);

                    model.set_lambda(*lambda);
                    model.reseed_sir_rng(&mut r);
                    let t = model.propagate_until_completion_time_with_locks(lockgraph,lockparams) as u32;
                    t

                }).collect();
                vec_temp

                }).collect();
            //sort time_vec, so we can easily histogram it.
            time_vec.sort_unstable();
            //let hist = convert_sorted_to_hist(&time_vec);
            let percent_life_span = acquire_percent_life_span(&time_vec,param.lifespanpercent);
            percent_life_span
        }).progress_with(bar).collect();
        n_lambdaspanned



    }).collect();


    alternate_writing(param, json, num_threads, data_vec, n_size, lambda_vec);}




    fn sim_small_world(param: LifespanSizeFittingParams, json: Value, num_threads:Option<NonZeroUsize>){

    //this is retired. see barabasi for proper function.
    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());

    rayon::ThreadPoolBuilder::new().num_threads(k.get()).build_global().unwrap();

    let mut n_size:Vec<_> = param.system_size_range.iter().map(|n|{*n}).collect();
    println!("{:?}",n_size);
    let lambda_range = param.trans_prob_range.get_range();
    let mut lambda_vec:Vec<_> = lambda_range.iter().collect();

    //let bar = crate::indication_bar(n_size.len() as u64);

    let list_of_data_vec:Vec<_> = n_size.iter_mut().map(|n|{
        println!("{}",n);
        let bar = crate::indication_bar(lambda_vec.len() as u64);

        let lambda_t_percent_avg_and_var_vec:Vec<_> = lambda_vec.iter_mut().map(|lambda|{
            println!("{}",lambda);
            let mut avg = 0.;
            let mut var = 0.;
            for _b in 0..param.num_networks{
                let opt = BaseSwOptions::from_lifespan_size_fitting_param_and_size(&param,NonZeroUsize::new(*n).unwrap());
                let small_world = opt.into();
                let mut model = SimpleSample::from_base(small_world, param.sir_seed);
                
                
                model.set_lambda(*lambda);

        
                let sorted_data = acquire_sorted_data(model,k,param.sir_seed,10);
                let t_90 = acquire_percent_life_span(&sorted_data,0.9);
                avg += t_90 as f64;
                var += t_90 as f64 *t_90 as f64;
            }
            avg /=param.num_networks as f64;
            var /= param.num_networks as f64;
            var -= avg*avg;

            Measured{
                var_t: MyVariance{
                    mean: avg,
                    var
                }
            }
        }).progress_with(bar).collect();
        lambda_t_percent_avg_and_var_vec

    }).collect();

    write_data_files(&list_of_data_vec, &param, &json, num_threads, n_size, lambda_vec);





}










pub struct Measured
{
    pub var_t: MyVariance,
}

fn alternate_writing(param:LifespanSizeFittingParams,json:Value,num_threads: Option<NonZeroUsize>,data_master:Vec<Vec<u32>>,n_list:Vec<usize>,lambda:Vec<f64>){
    
    for j in 0..n_list.len(){
        let name = param.name("dat", num_threads,n_list[j]);
        println!("creating: {name}");
        let file = File::create(name).expect("unable to create file");
        let mut buf = BufWriter::new(file);
        write!(buf, "#").unwrap();
        serde_json::to_writer(&mut buf, &json).unwrap();
        writeln!(buf).unwrap();
        // let data_vec = &data[n];
        writeln!(buf, "#lambda percentlifespan").unwrap();
        let data = &data_master[j];
        for k in 0..data.len(){
            
            writeln!(buf,"{} {}",lambda[k],data[k]).unwrap();
        }

    }
    
    

}





fn write_data_files(data:&[Vec<Measured>],param: &LifespanSizeFittingParams, json: &Value, num_threads:Option<NonZeroUsize>,nvec:Vec<usize>,lambda_vec:Vec<f64>){
    

    for n in 0..nvec.len(){
        let name = param.name("dat", num_threads,nvec[n]);
        println!("creating: {name}");



        let file = File::create(name).expect("unable to create file");
        let mut buf = BufWriter::new(file);
        write!(buf, "#").unwrap();
        serde_json::to_writer(&mut buf, &json).unwrap();
        writeln!(buf).unwrap();
        let data_vec = &data[n];
        writeln!(buf, "#lambda time_steps_avg time_steps_var").unwrap();
        for j in 0..data_vec.len(){
            let var = &data_vec[j];
            writeln!(buf,"{} {} {}",lambda_vec[j],var.var_t.mean,var.var_t.var).unwrap();
        }

    }


    
    
}