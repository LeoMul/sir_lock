use std::{ io::Write};


use{
    super::*,
    crate::{GraphType, sir_model::*},
    serde_json::Value,
    std::{num::*, fs::File, io::BufWriter},
    crate::stats_methods::MyVariance,
    crate::life_span_data_collection::*,
    //crate::grid::*
};


use crate:: sir_model::base_model_options::BaseSwOptions;

pub fn run_simulation(param:LifespanSizeFittingParams, json: Value, num_threads: Option<NonZeroUsize>){
    match param.graph_type{
        GraphType::SmallWorld(_) => sim_small_world(param, json, num_threads),
        _ => unimplemented!()
    }

}

fn sim_small_world(param: LifespanSizeFittingParams, json: Value, num_threads:Option<NonZeroUsize>){
    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());

    rayon::ThreadPoolBuilder::new().num_threads(k.get()).build_global().unwrap();

    let n_range = param.system_size_range.get_range();
    let mut n_size:Vec<_> = n_range.iter().collect();
    println!("{:?}",n_size);
    let lambda_range = param.trans_prob_range.get_range();
    let mut lambda_vec:Vec<_> = lambda_range.iter().collect();

    //let bar = crate::indication_bar(n_size.len() as u64);

    let list_of_data_vec:Vec<_> = n_size.iter_mut().map(|n|{
        println!("{}",n);

        let lambda_t_percent_avg_and_var_vec:Vec<_> = lambda_vec.iter_mut().map(|lambda|{
            println!("{}",lambda);
            let mut avg = 0.;
            let mut var = 0.;
            for _b in 0..param.num_networks{
                let opt = BaseSwOptions::from_lifespan_size_fitting_param_and_size(&param,NonZeroUsize::new(*n).unwrap());
                let small_world = opt.into();
                let mut model = SimpleSample::from_base(small_world, param.sir_seed);
                
                
                model.set_lambda(*lambda);

        
                let sorted_data = acquire_sorted_data(model,k,param.sir_seed,param.samples_per_step);
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
        }).collect();
        lambda_t_percent_avg_and_var_vec

    }).collect();

    write_data_files(&list_of_data_vec, &param, &json, num_threads, n_size, lambda_vec);





}

pub struct Measured
{
    pub var_t: MyVariance,
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