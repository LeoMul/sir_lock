use{
    super::*,
    std::num::*,
    serde_json::Value,

};
use crate::{ sir_model::*};

use std::{io::Write};
use rayon::iter::{ParallelIterator};



use{
    crate::{GraphType},
    
    std::{fs::File, io::BufWriter},
    rand_pcg::Pcg64,
    net_ensembles::rand::SeedableRng,

    rayon::prelude::*,
    crate::stats_methods::MyVariance
};

//use crate:: sir_model::base_model_options::BaseSwOptions;

pub fn run_simulation(param:CriticalLambdaParams, json: Value, num_threads: Option<NonZeroUsize>){
    match param.graph_type{
        GraphType::SmallWorld(_) => sim_small_world(param, json, num_threads),
        _ => unimplemented!()
    }

}

//returns scan-lambda-esque data about a PARTICULAR GRAPH for a particular lambda range.


fn sim_small_world(param:CriticalLambdaParams, json: Value, num_threads: Option<NonZeroUsize>){
    
    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    rayon::ThreadPoolBuilder::new().num_threads(k.get()).build_global().unwrap();
    //let n_range = param.system_size_range.get_range();
    let mut n_size:Vec<_> = param.system_size_range.iter().map(|n|{*n}).collect();
    println!("{:?}",n_size);
    let lambda_range = param.lambda_range.get_range();
    let lambda_vec:Vec<_> = lambda_range.iter().collect();

    //let bar = crate::indication_bar(lambda_vec.len() as u64);
    //contains a vector for each value of n. 
    let data_vec:Vec<_> = n_size.iter_mut().map(|n|{
        //let mut netvec = Vec::new();
        println!("{}",n);
        let n_lambdaspanned:Vec<_> = lambda_vec.par_iter().map(|lambda|{
            let mut rng = Pcg64::seed_from_u64(param.sir_seed);
            let mut avg_c = 0.;
            let mut var_c = 0.;
            //might be better to have the parallel here? good enough for now tho
            

            for _j in 0..param.num_networks{
                let opt = BaseSwOptions::from_critical_lambda_params(&param,NonZeroUsize::new(*n).unwrap());
                let small_world = opt.into();
                let mut model = SimpleSample::from_base(small_world, param.sir_seed);

                model.set_lambda(*lambda);
                model.reseed_sir_rng(&mut rng);
                model.propagate_until_completion_max();

                let c = model.calculate_ever_infected() as f64/ *n as f64;
                avg_c += c as f64;
                //avg_c /= *n as f64;
                var_c += (c*c) as f64;
            }



            avg_c /= (param.num_networks) as f64;
            var_c /= param.num_networks as f64;
            var_c = var_c - avg_c*avg_c;
            MyVariance{
                mean: avg_c,
                var: var_c
            }
        }).collect();
        n_lambdaspanned



    }).collect();

    alternate_writing(param, json, num_threads, data_vec,n_size,lambda_vec);

}



fn _find_critical_lambda_one_data_set(vec:Vec<MyVariance>,lambdavals:&[f64])->f64{
    let mut var_vec:Vec<f64> = vec![0.;vec.len()];
    for i in 0..vec.len(){
        var_vec[i] = vec[i].var;
    }
    let abs_max = var_vec
        .iter()
        .max_by(|x, y| x.abs().partial_cmp(&y.abs()).unwrap())
        .unwrap();
    let index = var_vec.iter().position(|&x| x == *abs_max).unwrap();

    lambdavals[index]

}



fn alternate_writing(param:CriticalLambdaParams,json:Value,num_threads: Option<NonZeroUsize>,data_master:Vec<Vec<MyVariance>>,n_list:Vec<usize>,lambda:Vec<f64>){
    
    for j in 0..n_list.len(){
        let name = param.name("dat", num_threads,n_list[j]);
        println!("creating: {name}");
        let file = File::create(name).expect("unable to create file");
        let mut buf = BufWriter::new(file);
        write!(buf, "#").unwrap();
        serde_json::to_writer(&mut buf, &json).unwrap();
        writeln!(buf).unwrap();
        // let data_vec = &data[n];
        writeln!(buf, "#lambda c var_c").unwrap();
        let data = &data_master[j];
        for k in 0..data.len(){
            
            writeln!(buf,"{} {} {}",lambda[k],data[k].mean,data[k].var).unwrap();
        }

    }
    
    

}


