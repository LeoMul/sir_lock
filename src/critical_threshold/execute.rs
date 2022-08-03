use{
    super::*,
    std::num::*,
    serde_json::Value,

};
use std::{ io::Write};
use indicatif::*;
use rayon::iter::{ParallelIterator};
use rand::Rng;
//use serde_json::from_slice;
use std::sync::Mutex;
use std::ops::DerefMut;
use{
    crate::{GraphType, sir_model::*},
    std::{ fs::File, io::BufWriter},
    rand_pcg::Pcg64,
    net_ensembles::rand::SeedableRng,
    rayon::prelude::*,
    crate::stats_methods::MyVariance,
    crate::misc_types::MeasureType,
    crate::lockdown_methods::*,
};

use crate::stats_methods::*;
pub struct Measured
{
    pub var_m: MyVariance,
    pub var_c: MyVariance,
}

pub type MyVarianceInBoots = MyVarianceBootstrap;
pub struct MeasuredWithErrors{
    pub m: MyVarianceInBoots,
    pub c: MyVarianceInBoots
}




pub fn run_simulation(param:CriticalThreshParams, json: Value, num_threads: Option<NonZeroUsize>){
    match param.graph_type{
        GraphType::SmallWorld(_) => sim_small_world_both_m_and_c(param, json, num_threads),
        GraphType::Barabasi(_,_) => sim_barabasi_both_m_and_c(param, json, num_threads),
        _ => unimplemented!()
    }

}
fn sim_small_world_both_m_and_c(param: CriticalThreshParams,json:Value,num_threads:Option<NonZeroUsize>){
    
    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    rayon::ThreadPoolBuilder::new().num_threads(k.get()).build_global().unwrap();
    let n_size:Vec<_> = param.system_size_range.to_vec();
    println!("Sampling energies C and M for system sizes{:?} over thresh-values {}-{} with {} points.\n Each data point is averaged over {} networks. \n Bootstrapping is {}."
    ,n_size,param.thresh_range.start,param.thresh_range.end,param.thresh_range.steps,param.num_networks,param.bootbool);
    let thresh_range = param.thresh_range.get_range();
    let thresh_vec:Vec<_> = thresh_range.iter().collect();

    let mut sir_rng_2 = Pcg64::seed_from_u64(param.sir_seed);
    let mut graph_rng = Pcg64::seed_from_u64(param.graph_seed);
    for n in n_size{
        println!("Commencing sampling for {}",n);
        let system_size_frac = if param.fraction{
            Some(n as f64)
        } else{
            None
        };
        let mut rngs: Vec<_> = (0..k.get())
        .map(
            |_| 
                {      (Pcg64::from_rng(&mut sir_rng_2).unwrap(),
                        Pcg64::from_rng(&mut graph_rng).unwrap())
                    
                }
            )
        .collect();
        let per_thread = param.num_networks/k.get() as u64;
        
        let bar = crate::indication_bar(param.num_networks);
        //This vector contains 4 vectors of length per_thread. These constituent vectors contain (c.m) spanned over lambda.
        let intermediate_data_vecs: Vec<Vec<Vec<(u32,u32)>>> = rngs.par_iter_mut().map(
            |(r,graph_rng)|
            {
                //let mut new_vec:Vec<_> = Vec::with_capacity(per_thread as usize);
                let iter = (0..per_thread).into_iter().map(|_|{});
                iter.map(
                    |_|
                    {
                    let mut lockparams = LockdownParameters{
                        lock_style: param.lockdowntype,
                        release_threshold: 0.0,
                        lock_threshold: 0.0
                    };

                    let new_graph_seed = graph_rng.gen::<u64>();
                    let opt = SWOptions::from_critical_thresh_params(&param,NonZeroUsize::new(n).unwrap(),new_graph_seed);
                    let small_world = opt.into();
                    let mut model = SimpleSampleSW::from_base(small_world, param.sir_seed,param.initial_infected);
                    let mut lockgraph = model.create_locked_down_network(lockparams);

                    //(c,m) for each lambda. 
                    let data_point_for_each_lambda:Vec<_> = thresh_vec.iter().map(|thresh|{
                        //model.set_lambda(*lambda);
                        let rel = get_release_threshold(param.releasetype,*thresh);
                        lockparams.set_lock_thresh(*thresh);
                        lockparams.set_rel_thresh(rel);
                        //println!("{rel},{thresh}");
                        model.reseed_sir_rng(r);
                        let m = model.propagate_until_completion_max_with_lockdown(&mut lockgraph,lockparams);
                        //let c = 
                        //println!("{},{},{},{},{}",c,m,lambda,model.ensemble.graph().edge_count(),new_graph_seed);
                        let c = model.calculate_ever_infected() as u32;
                        (m as u32,c)
                        }).collect(); 

                    data_point_for_each_lambda
                    }).progress_with(bar.clone()).collect()//this produces per_thread vectors of (c,m) for each lambda
            }
        ).collect();
        //goal: to convert this to two vectors 

        let mut merged_vec = Vec::new();
        for vec in intermediate_data_vecs{
            merged_vec.extend(vec);
        }
        let (m_matrix,c_matrix) = transpose_into_two_vecs(&merged_vec);
        //This is now all of the data decomposed for this n.
        
        if param.bootbool{
            
            let m_measure_data = produce_measure_data_with_boots(m_matrix, &mut sir_rng_2, &mut graph_rng, k, system_size_frac);
            let c_measure_data = produce_measure_data_with_boots(c_matrix, &mut sir_rng_2, &mut graph_rng, k, system_size_frac);

            writing_with_boot_single_energy(&param, &json, num_threads, &m_measure_data, n, &thresh_vec,MeasureType::M);
            writing_with_boot_single_energy(&param, &json, num_threads, &c_measure_data, n, &thresh_vec,MeasureType::C)

            }
        else{
            //do not resample and produce data file.print
            let m_measure_data:Vec<MyVariance> = m_matrix.par_iter().map(|item|{
                crate::stats_methods::MyVariance::from_slice(item,system_size_frac)
            }).collect();
            let c_measure_data:Vec<MyVariance> = c_matrix.par_iter().map(|item|{
                crate::stats_methods::MyVariance::from_slice(item,system_size_frac)
            }).collect();
            writing_with_no_boot_single_energy(&param, &json, num_threads, &m_measure_data, n, &thresh_vec,MeasureType::M);
            writing_with_no_boot_single_energy(&param, &json, num_threads, &c_measure_data, n, &thresh_vec,MeasureType::C)

        }

    }
}


fn sim_barabasi_both_m_and_c(param: CriticalThreshParams,json:Value,num_threads:Option<NonZeroUsize>){
    
    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    rayon::ThreadPoolBuilder::new().num_threads(k.get()).build_global().unwrap();
    let n_size:Vec<_> = param.system_size_range.to_vec();
    println!("Sampling energies C and M for system sizes{:?} over thresh-values {}-{} with {} points.\n Each data point is averaged over {} networks. \n Bootstrapping is {}."
    ,n_size,param.thresh_range.start,param.thresh_range.end,param.thresh_range.steps,param.num_networks,param.bootbool);
    let thresh_range = param.thresh_range.get_range();
    let thresh_vec:Vec<_> = thresh_range.iter().collect();

    let mut sir_rng_2 = Pcg64::seed_from_u64(param.sir_seed);
    let mut graph_rng = Pcg64::seed_from_u64(param.graph_seed);
    for n in n_size{
        println!("Commencing sampling for {}",n);
        let system_size_frac = if param.fraction{
            Some(n as f64)
        } else{
            None
        };
        let mut rngs: Vec<_> = (0..k.get())
        .map(
            |_| 
                {      (Pcg64::from_rng(&mut sir_rng_2).unwrap(),
                        Pcg64::from_rng(&mut graph_rng).unwrap())
                    
                }
            )
        .collect();
        let per_thread = param.num_networks/k.get() as u64;
        
        let bar = crate::indication_bar(param.num_networks);
        //This vector contains 4 vectors of length per_thread. These constituent vectors contain (c.m) spanned over lambda.
        let intermediate_data_vecs: Vec<Vec<Vec<(u32,u32)>>> = rngs.par_iter_mut().map(
            |(r,graph_rng)|
            {
                //let mut new_vec:Vec<_> = Vec::with_capacity(per_thread as usize);
                let iter = (0..per_thread).into_iter().map(|_|{});
                iter.map(
                    |_|
                    {
                    let mut lockparams = LockdownParameters{
                        lock_style: param.lockdowntype,
                        release_threshold: 0.0,
                        lock_threshold: 0.0
                    };

                    let new_graph_seed = graph_rng.gen::<u64>();
                    let opt = BarabasiOptions::from_critical_thresh_params(&param,NonZeroUsize::new(n).unwrap(),new_graph_seed);
                    let world = opt.into();
                    let mut model = SimpleSampleBarabasi::from_base(world, param.sir_seed,param.initial_infected);
                    let mut lockgraph = model.create_locked_down_network(lockparams);

                    //(c,m) for each lambda. 
                    let data_point_for_each_lambda:Vec<_> = thresh_vec.iter().map(|thresh|{
                        //model.set_lambda(*lambda);
                        let rel = get_release_threshold(param.releasetype,*thresh);
                        lockparams.set_lock_thresh(*thresh);
                        lockparams.set_rel_thresh(rel);
                        //println!("{rel},{thresh}");
                        model.reseed_sir_rng(r);
                        let m = model.propagate_until_completion_max_with_lockdown(&mut lockgraph,lockparams);
                        //let c = 
                        //println!("{},{},{},{},{}",c,m,lambda,model.ensemble.graph().edge_count(),new_graph_seed);
                        let c = model.calculate_ever_infected() as u32;
                        (m as u32,c)
                        }).collect(); 

                    data_point_for_each_lambda
                    }).progress_with(bar.clone()).collect()//this produces per_thread vectors of (c,m) for each lambda
            }
        ).collect();
        //goal: to convert this to two vectors 

        let mut merged_vec = Vec::new();
        for vec in intermediate_data_vecs{
            merged_vec.extend(vec);
        }
        let (m_matrix,c_matrix) = transpose_into_two_vecs(&merged_vec);
        //This is now all of the data decomposed for this n.
        
        if param.bootbool{
            
            let m_measure_data = produce_measure_data_with_boots(m_matrix, &mut sir_rng_2, &mut graph_rng, k, system_size_frac);
            let c_measure_data = produce_measure_data_with_boots(c_matrix, &mut sir_rng_2, &mut graph_rng, k, system_size_frac);

            writing_with_boot_single_energy(&param, &json, num_threads, &m_measure_data, n, &thresh_vec,MeasureType::M);
            writing_with_boot_single_energy(&param, &json, num_threads, &c_measure_data, n, &thresh_vec,MeasureType::C)

            }
        else{
            //do not resample and produce data file.print
            let m_measure_data:Vec<MyVariance> = m_matrix.par_iter().map(|item|{
                crate::stats_methods::MyVariance::from_slice(item,system_size_frac)
            }).collect();
            let c_measure_data:Vec<MyVariance> = c_matrix.par_iter().map(|item|{
                crate::stats_methods::MyVariance::from_slice(item,system_size_frac)
            }).collect();
            writing_with_no_boot_single_energy(&param, &json, num_threads, &m_measure_data, n, &thresh_vec,MeasureType::M);
            writing_with_no_boot_single_energy(&param, &json, num_threads, &c_measure_data, n, &thresh_vec,MeasureType::C)

        }

    }
}




fn produce_measure_data_with_boots<R>(merged_vec:Vec<Vec<u32>>,mut rng1:&mut R,mut rng2:&mut  R,num_threads:NonZeroUsize,system_size_frac:Option<f64>) 
-> Vec<MyVarianceBootstrap> where R: Rng + rand::RngCore{
    println!("Commencing bootstrapping");
    let bar = crate::indication_bar(merged_vec.len() as u64);
    
    //let measure_data = rngs.par_iter()
    let end = merged_vec.len()/num_threads.get();
    let rngs: Vec<_> = (0..num_threads.get()).map(
        |_|
        {
            
            Mutex::new(
                (Pcg64::from_rng(&mut rng1).unwrap(),
                Pcg64::from_rng(&mut rng2).unwrap())
            )

        }
    ).collect();
    let iterator = (0..=end)
        .into_par_iter()
        .flat_map(
        |i|
        {
            let start = i*rngs.len();
            rngs.par_iter()
                .zip(merged_vec[ start..].par_iter())
        }
        );
    let vec = iterator.progress_with(bar).map(
        |(container,slice)|{
            let mut lock = container.lock()
            .expect("unable to lock");
        let inner = lock.deref_mut();
        //let vaccine_rng = &mut inner.0;
        let rng1 = &mut inner.0;
        //let vaccine_list_helper = &mut inner.2;
        let rng2 = &mut inner.1;
        crate::stats_methods::MyVarianceBootstrap::from_slice(slice,system_size_frac,rng1 ,rng2)
        }
    ).collect();
    println!("Bootstrapping complete.");
    vec
    

}

fn _transpose_into_a_vec(merged_vec:&Vec<Vec<u32>>) -> Vec<Vec<u32>>{
    let num_lam = merged_vec[0].len();
    let num_cols = merged_vec.len();
    let mut c_matrix:Vec<Vec<u32>> = Vec::with_capacity(num_lam);
    
    for i in 0..num_lam{
        let mut c_vec = Vec::with_capacity(num_cols);
        for j in 0..num_cols{
            let col = &merged_vec[j];
            let point = col[i];
            c_vec.push(point as u32);
        }

        c_matrix.push(c_vec);

    }   
    c_matrix        
}
fn transpose_into_two_vecs(merged_vec:&Vec<Vec<(u32,u32)>>) -> (Vec<Vec<u32>>,Vec<Vec<u32>>){
    let num_lam = merged_vec[0].len();
    let num_cols = merged_vec.len();
    let mut c_matrix:Vec<Vec<u32>> = Vec::with_capacity(num_lam);
    let mut m_matrix:Vec<Vec<u32>> = Vec::with_capacity(num_lam);

    for i in 0..num_lam{
        let mut c_vec = Vec::with_capacity(num_cols);
        let mut m_vec = Vec::with_capacity(num_cols);

        for j in 0..num_cols{
            let col = &merged_vec[j];
            let point = col[i];
            c_vec.push(point.1 as u32);
            m_vec.push(point.0 as u32);
        }

        c_matrix.push(c_vec);
        m_matrix.push(m_vec);

    }   
    (m_matrix,c_matrix)        
}
fn writing_with_boot_single_energy(param:&CriticalThreshParams,json:&Value,num_threads: Option<NonZeroUsize>,data:&Vec<MyVarianceInBoots>,n:usize,lambda:&Vec<f64>,energy:MeasureType){
    let name = param.name("dat", num_threads,n,energy);
    println!("creating: {name}");
    let file = File::create(name).expect("unable to create file");
    let mut buf = BufWriter::new(file);
    write!(buf, "#").unwrap();
    serde_json::to_writer(&mut buf, &json).unwrap();
    writeln!(buf).unwrap();
    // let data_vec = &data[n];
    if energy.is_c(){
        writeln!(buf, "#thresh c c_var c_err c_var_err").unwrap();
    }
    else{
        writeln!(buf, "#thresh m m_var m_err m_var_err").unwrap();

    }
    //let data = &data_master[j];
    for k in 0..data.len(){
        
        writeln!(buf,"{} {} {} {} {}",lambda[k],data[k].mean,data[k].var,data[k].mean_err,data[k].var_err).unwrap();
    }
    
}
fn writing_with_no_boot_single_energy(param:&CriticalThreshParams,json:&Value,num_threads: Option<NonZeroUsize>,data:&Vec<MyVariance>,n:usize,lambda:&Vec<f64>,energy:MeasureType){
    let name = param.name("dat", num_threads,n,energy);
    println!("creating: {name}");
    let file = File::create(name).expect("unable to create file");
    let mut buf = BufWriter::new(file);
    write!(buf, "#").unwrap();
    serde_json::to_writer(&mut buf, &json).unwrap();
    writeln!(buf).unwrap();
    // let data_vec = &data[n];
    if energy.is_c(){
        writeln!(buf, "#thresh c c_var ").unwrap();
    }
    else{
        writeln!(buf, "#thresh m m_var ").unwrap();

    }
    //let data = &data_master[j];
    for k in 0..data.len(){
        
        writeln!(buf,"{} {} {}",lambda[k],data[k].mean,data[k].var).unwrap();
    }
    
}





