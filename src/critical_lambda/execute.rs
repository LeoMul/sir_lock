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
fn scan_lambda_function(param:&CriticalLambdaParams, _json: &Value, num_threads: Option<NonZeroUsize>,n:usize) -> Vec<MyVariance>{

    let opt = BaseSwOptions::from_critical_lambda_params(param,NonZeroUsize::new(n).unwrap());
    let small_world = opt.into();
    let model = SimpleSample::from_base(small_world, param.sir_seed);
    let lambda_range = param.lambda_range.get_range();
    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    let lambda_vec:Vec<_> = lambda_range.iter().collect();

    let system_size_fracion_indicator = if param.fraction{
        n as f64
    } else{
        1.
    };
    let lambda_data_vec:Vec<_> = lambda_vec.iter().map(|lambda|{
        let mut rng = Pcg64::seed_from_u64(param.sir_seed);

        let mut container: Vec<_> = (0..k.get()).map(
        |_|
        {
            let mut model = model.clone();
            model.reseed_sir_rng(&mut rng);
            model
        }
        ).collect();
        let per_threads = param.samples_per_lambda/k.get() as u64;
        let vec:Vec<_> = container.par_iter_mut().map(|model|{
            model.set_lambda(*lambda);
            //println!("model lambda {}",model.lambda);
            let mut sum_c = 0.;
            let mut var_c = 0.;
            for _j in 0..per_threads{
                model.propagate_until_completion_max();
                let c= model.calculate_ever_infected() as f64/system_size_fracion_indicator;
                sum_c += c;
                var_c += c*c;

            }
            (sum_c,var_c)

        }).collect();
        let mut avg_c = 0.;
        let mut var_c = 0.;
        for (i,j) in vec{
            avg_c += i;
            var_c += j;

        }
        avg_c /= k.get() as f64 *per_threads as f64;
        var_c /=  k.get() as f64 *per_threads as f64;
        var_c -= avg_c*avg_c;
        MyVariance{
            mean:avg_c,
            var:var_c,
        }

    }).collect();
    lambda_data_vec
}


fn _sim_small_world_testing(param:CriticalLambdaParams, json: Value, num_threads: Option<NonZeroUsize>){
    let n = 3200;
    let opt = BaseSwOptions::from_critical_lambda_params(&param,NonZeroUsize::new(n).unwrap());
    let small_world = opt.into();
    let _model = SimpleSample::from_base(small_world, param.sir_seed);
    let lambda_range = param.lambda_range.get_range();
    let _k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    let lambda_vec:Vec<_> = lambda_range.iter().collect();
    let lambda_data = scan_lambda_function(&param, &json, num_threads, n);

    let name = "testingfile.dat";
        println!("creating: {name}");
        let file = File::create(name).expect("unable to create file");
        let mut buf = BufWriter::new(file);
        write!(buf, "#").unwrap();
        serde_json::to_writer(&mut buf, &json).unwrap();
        writeln!(buf).unwrap();
       // let data_vec = &data[n];
        writeln!(buf, "#lambda time_steps_avg time_steps_var").unwrap();
        for j in 0..lambda_vec.len(){
    
            writeln!(buf,"{} {} {}",lambda_vec[j],lambda_data[j].mean,lambda_data[j].var).unwrap();
        }

}

fn sim_small_world(param:CriticalLambdaParams, json: Value, num_threads: Option<NonZeroUsize>){
    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());

    rayon::ThreadPoolBuilder::new().num_threads(k.get()).build_global().unwrap();

    let n_range = param.system_size_range.get_range();
    let mut n_size:Vec<_> = n_range.iter().collect();
    println!("{:?}",n_size);
    let lambda_range = param.lambda_range.get_range();
    let lambda_vec:Vec<_> = lambda_range.iter().collect();
    //let bar = crate::indication_bar(lambda_vec.len() as u64);
    //.progress_with(bar)

    let data_vec:Vec<_> = n_size.iter_mut().map(|n|{
        println!("{}",n);
        let mut avg_crit_lambda = 0.;
        let mut var_crit_lambda = 0.;
        for _j in 0..param.num_networks{
            let lambda_data = scan_lambda_function(&param, &json, num_threads, *n);
            let new_crit = find_critical_lambda_one_data_set(lambda_data, &lambda_vec);
            //println!("{}",new_crit);
            avg_crit_lambda += new_crit;
            var_crit_lambda += new_crit*new_crit;
            
        }
        avg_crit_lambda /= param.num_networks as f64;
        var_crit_lambda /=param.num_networks as f64;
        var_crit_lambda -= avg_crit_lambda*avg_crit_lambda;
        (avg_crit_lambda,var_crit_lambda)

    }).collect();
    write_data_files(param, json, num_threads, n_size, data_vec);

}


fn find_critical_lambda_one_data_set(vec:Vec<MyVariance>,lambdavals:&[f64])->f64{
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

fn write_data_files(param:CriticalLambdaParams, json: Value, num_threads: Option<NonZeroUsize>,n_vector:Vec<usize>,data:Vec<(f64,f64)>){
    
        let name = param.name("dat", num_threads);
        println!("creating: {name}");
        let file = File::create(name).expect("unable to create file");
        let mut buf = BufWriter::new(file);
        write!(buf, "#").unwrap();
        serde_json::to_writer(&mut buf, &json).unwrap();
        writeln!(buf).unwrap();
       // let data_vec = &data[n];
        writeln!(buf, "#lambda time_steps_avg time_steps_var").unwrap();
        for j in 0..data.len(){
    
            writeln!(buf,"{} {} {}",n_vector[j],data[j].0,data[j].1).unwrap();
        }

    

}



