use{
    super::*,
    std::num::*,
    serde_json::Value,

};
use std::{ io::Write};
use indicatif::*;
use rayon::iter::{ParallelIterator};
use rand::Rng;
use{
    crate::{GraphType, sir_model::*},
    std::{ fs::File, io::BufWriter},
    rand_pcg::Pcg64,
    net_ensembles::rand::SeedableRng,
    rayon::prelude::*,
    crate::stats_methods::MyVariance
};








//use crate:: sir_model::base_model_options::BaseSwOptions;

pub fn run_simulation(param:CriticalLambdaParams, json: Value, num_threads: Option<NonZeroUsize>){
    match param.graph_type{
        GraphType::SmallWorld(_) => sim_small_world(param, json, num_threads),
        GraphType::Barabasi(_,_) => sim_ba(param, json, num_threads),
        _ => unimplemented!()
    }

}

fn sim_ba(param:CriticalLambdaParams,json:Value,num_threads:Option<NonZeroUsize>){

    //this is the fastest program I have for this yet. network size 3000 with 50,000 networks only takes 8 minutes.
    // This same computation would have taken 
    //2,5 hrs previously

    //some optimisation could be done with the random numbers.. 
    

    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    rayon::ThreadPoolBuilder::new().num_threads(k.get()).build_global().unwrap();
    let mut n_size:Vec<_> = param.system_size_range.to_vec();
    let lockparams = param.lockdown;
    println!("{:?}",n_size);
    let lambda_range = param.lambda_range.get_range();
    let lambda_vec:Vec<_> = lambda_range.iter().collect();
    //println!("Progress will take approx 3-4 times what the bar says");
    let mut sir_rng_2 = Pcg64::seed_from_u64(param.sir_seed);
    let mut graph_rng = Pcg64::seed_from_u64(param.graph_seed);
    let _data_vec:Vec<_> = n_size.iter_mut().map(|n|{
        println!("{}",n);
        
        let divisor = if param.fraction{
            *n as f64
        } else{ 1.};
        

        //let mut container: Vec<_> = (0..k.get()).map(|_|{}).collect();

        let mut rngs: Vec<_> = (0..k.get())
        .map(
            |_| 
                {
                    
                        (Pcg64::from_rng(&mut sir_rng_2).unwrap(),
                        Pcg64::from_rng(&mut graph_rng).unwrap())
                    
                }
            )
        .collect();


        let per_thread = param.num_networks/k.get() as u64;
        
        let bar = crate::indication_bar(param.num_networks);

        //this is a vector of vectors. let let the elements of this vector be rows. for the measurements we want, we require summing over the columns.
        let intermediate_data_vecs: Vec<_> = rngs.par_iter_mut().map(|(r,graph_rng)|{
            //let mut new_vec:Vec<_> = Vec::with_capacity(per_thread as usize);
            let mut rng = r;

            let iter = (0..per_thread).into_iter().map(|_|{});
            


            let new_vec = iter.map(|_|{
                let new_graph_seed = graph_rng.gen::<u64>();
                let opt = BarabasiOptions::from_critical_lambda_params(&param,NonZeroUsize::new(*n).unwrap(),new_graph_seed);
                let small_world = opt.into();
                let mut model = SimpleSampleBarabasi::from_base(small_world, param.sir_seed);
                let lockgraph = model.create_locked_down_network(lockparams);

                let data_point_for_each_lambda:Vec<_> = lambda_vec.iter().map(|lambda|{
                    model.set_lambda(*lambda);
                    model.reseed_sir_rng(&mut rng);
                    let m = model.propagate_until_completion_max_with_lockdown(lockgraph.clone(),lockparams) as f64/divisor;
                    let c = model.calculate_ever_infected() as f64/divisor;
                    //println!("{},{},{},{},{}",c,m,lambda,model.ensemble.graph().edge_count(),new_graph_seed);
                    (c,m)

                }).collect(); data_point_for_each_lambda}).progress_with(bar.clone()).collect();

                new_vec
            }
            
        ).collect();

        let mut measure_vec:Vec<Measured> = Vec::with_capacity(lambda_vec.len());
        let num = (per_thread*k.get() as u64) as f64;
        let mut copy_vec = Vec::with_capacity(k.get());
        for item in intermediate_data_vecs.iter().take(copy_vec.capacity()){
            copy_vec.push(collect_columns_sums(item));
        }
        for j in 0..measure_vec.capacity(){
            let mut new_tuple = (0.,0.,0.,0.);
            for x in &copy_vec{
                //let x  = &copy_vec[i];
                new_tuple = add_four_tuples(new_tuple, x[j]);}

            let av_c = new_tuple.0/num;
            let var_c = new_tuple.1/num- av_c*av_c;
            let av_m = new_tuple.2/num ;
            let var_m = new_tuple.3/num - av_m*av_m;
            
            measure_vec.push(Measured{
                var_c: MyVariance{
                    mean: av_c,
                    var: var_c,
                },
                var_m: MyVariance{
                    mean: av_m,
                    var: var_m
                }
            });
            //println!("{}",av_c);
            }
    
    writing(param.clone(),json.clone(),num_threads,&measure_vec,*n,lambda_vec.clone());
    measure_vec
    }
    
    ).collect();
    //alternate_writing(param, json, num_threads, data_vec,n_size,lambda_vec);

}
fn sim_small_world(param:CriticalLambdaParams,json:Value,num_threads:Option<NonZeroUsize>){

    //this is the fastest program I have for this yet. network size 3000 with 50,000 networks only takes 8 minutes.
    // This same computation would have taken 
    //2,5 hrs previously

    //some optimisation could be done with the random numbers.. 
    

    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    rayon::ThreadPoolBuilder::new().num_threads(k.get()).build_global().unwrap();
    let mut n_size:Vec<_> = param.system_size_range.to_vec();
    let lockparams = param.lockdown;
    println!("{:?}",n_size);
    let lambda_range = param.lambda_range.get_range();
    let lambda_vec:Vec<_> = lambda_range.iter().collect();
    //println!("Progress will take approx 3-4 times what the bar says");
    let mut sir_rng_2 = Pcg64::seed_from_u64(param.sir_seed);
    let mut graph_rng = Pcg64::seed_from_u64(param.graph_seed);
    let _data_vec:Vec<_> = n_size.iter_mut().map(|n|{
        println!("{}",n);
        
        let divisor = if param.fraction{
            *n as f64
        } else{ 1.};
        

        //let mut container: Vec<_> = (0..k.get()).map(|_|{}).collect();

        let mut rngs: Vec<_> = (0..k.get())
        .map(
            |_| 
                {
                    
                        (Pcg64::from_rng(&mut sir_rng_2).unwrap(),
                        Pcg64::from_rng(&mut graph_rng).unwrap())
                    
                }
            )
        .collect();


        let per_thread = param.num_networks/k.get() as u64;
        
        let bar = crate::indication_bar(param.num_networks);

        //this is a vector of vectors. let let the elements of this vector be rows. for the measurements we want, we require summing over the columns.
        let intermediate_data_vecs: Vec<_> = rngs.par_iter_mut().map(|(r,graph_rng)|{
            //let mut new_vec:Vec<_> = Vec::with_capacity(per_thread as usize);
            let mut rng = r;

            let iter = (0..per_thread).into_iter().map(|_|{});
            


            let new_vec = iter.map(|_|{
                let new_graph_seed = graph_rng.gen::<u64>();
                let opt = BaseSwOptions::from_critical_lambda_params(&param,NonZeroUsize::new(*n).unwrap(),new_graph_seed);
                let small_world = opt.into();
                let mut model = SimpleSample::from_base(small_world, param.sir_seed);
                let lockgraph = model.create_locked_down_network(lockparams);

                let data_point_for_each_lambda:Vec<_> = lambda_vec.iter().map(|lambda|{
                    model.set_lambda(*lambda);
                    model.reseed_sir_rng(&mut rng);
                    let m = model.propagate_until_completion_max_with_lockdown(lockgraph.clone(),lockparams) as f64/divisor;
                    let c = model.calculate_ever_infected() as f64/divisor;
                    //println!("{},{},{},{},{}",c,m,lambda,model.ensemble.graph().edge_count(),new_graph_seed);
                    (c,m)

                }).collect(); data_point_for_each_lambda}).progress_with(bar.clone()).collect();

                new_vec
            }
            
        ).collect();

        let mut measure_vec:Vec<Measured> = Vec::with_capacity(lambda_vec.len());
        let num = (per_thread*k.get() as u64) as f64;
        let mut copy_vec = Vec::with_capacity(k.get());
        for item in intermediate_data_vecs.iter().take(copy_vec.capacity()){
            copy_vec.push(collect_columns_sums(item));
        }
        for j in 0..measure_vec.capacity(){
            let mut new_tuple = (0.,0.,0.,0.);
            for x in &copy_vec{
                //let x  = &copy_vec[i];
                new_tuple = add_four_tuples(new_tuple, x[j]);}

            let av_c = new_tuple.0/num;
            let var_c = new_tuple.1/num- av_c*av_c;
            let av_m = new_tuple.2/num ;
            let var_m = new_tuple.3/num - av_m*av_m;
            
            measure_vec.push(Measured{
                var_c: MyVariance{
                    mean: av_c,
                    var: var_c,
                },
                var_m: MyVariance{
                    mean: av_m,
                    var: var_m
                }
            });
            //println!("{}",av_c);
            }
    
    writing(param.clone(),json.clone(),num_threads,&measure_vec,*n,lambda_vec.clone());
    measure_vec
    }
    
    ).collect();
    //alternate_writing(param, json, num_threads, data_vec,n_size,lambda_vec);

}
fn collect_columns_sums(data:&Vec<Vec<(f64,f64)>>)-> Vec<(f64,f64,f64,f64)>{
    
    let num_columns = data[0].len();
    let num_rows = data.len();
    let mut return_vec:Vec<(f64,f64,f64,f64)> = Vec::with_capacity(num_columns);
    for j in 0..num_columns{
        let mut avg_c = 0.;
        let mut var_c = 0.;
        let mut avg_m = 0.;
        let mut var_m = 0.;
        for item in data.iter().take(num_rows){
            let row_i = item;
            let c = row_i[j].0;
            let m = row_i[j].1;
            avg_c += c;
            var_c += c*c;
            avg_m += m;
            var_m += m*m;
            
        }
        return_vec.push((avg_c,var_c,avg_m,var_m));
    }
    return_vec
}

fn add_four_tuples(a:(f64,f64,f64,f64),b:(f64,f64,f64,f64))-> (f64,f64,f64,f64){
    (a.0+b.0,a.1+b.1,a.2+b.2,a.3+b.3)
}










//
//
//fn sim_small_world(param:CriticalLambdaParams, json: Value, num_threads: Option<NonZeroUsize>){
//    
//    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
//    rayon::ThreadPoolBuilder::new().num_threads(k.get()).build_global().unwrap();
//    //let n_range = param.system_size_range.get_range();
//    let mut n_size:Vec<_> = param.system_size_range.to_vec();
//    println!("{:?}",n_size);
//    let lambda_range = param.lambda_range.get_range();
//    let lambda_vec:Vec<_> = lambda_range.iter().collect();
//
//    //let bar = crate::indication_bar(lambda_vec.len() as u64);
//    //contains a vector for each value of n. 
//    let data_vec:Vec<_> = n_size.iter_mut().map(|n|{
//        println!("{}",n);
//        let n_lambdaspanned:Vec<_> = lambda_vec.par_iter().map(|lambda|{
//            let mut rng = Pcg64::seed_from_u64(param.sir_seed);
//            let mut avg_c = 0.;
//            let mut var_c = 0.;
//            //might be better to have the parallel here? good enough for now tho
//            
//
//            for _j in 0..param.num_networks{
//                let opt = BaseSwOptions::from_critical_lambda_params(&param,NonZeroUsize::new(*n).unwrap());
//                let small_world = opt.into();
//                let mut model = SimpleSample::from_base(small_world, param.sir_seed);
//
//                model.set_lambda(*lambda);
//                model.reseed_sir_rng(&mut rng);
//                model.propagate_until_completion_max();
//
//                let c = model.calculate_ever_infected() as f64/ *n as f64;
//                avg_c += c as f64;
//                //avg_c /= *n as f64;
//                var_c += (c*c) as f64;
//            }
//
//
//
//            avg_c /= (param.num_networks) as f64;
//            var_c /= param.num_networks as f64;
//            var_c -= avg_c*avg_c;
//            MyVariance{
//                mean: avg_c,
//                var: var_c
//            }
//        }).collect();
//        n_lambdaspanned
//
//
//
//    }).collect();
//
//    alternate_writing_old(param, json, num_threads, data_vec,n_size,lambda_vec);
//
//}
//fn _find_critical_lambda_one_data_set(vec:Vec<MyVariance>,lambdavals:&[f64])->f64{
//    let mut var_vec:Vec<f64> = vec![0.;vec.len()];
//    for i in 0..vec.len(){
//        var_vec[i] = vec[i].var;
//    }
//    let abs_max = var_vec
//        .iter()
//        .max_by(|x, y| x.abs().partial_cmp(&y.abs()).unwrap())
//        .unwrap();
//    let index = var_vec.iter().position(|&x| x == *abs_max).unwrap();
//
//    lambdavals[index]
//
//}


fn writing(param:CriticalLambdaParams,json:Value,num_threads: Option<NonZeroUsize>,data:&Vec<Measured>,n:usize,lambda:Vec<f64>){
    let name = param.name("dat", num_threads,n);
    println!("creating: {name}");
    let file = File::create(name).expect("unable to create file");
    let mut buf = BufWriter::new(file);
    write!(buf, "#").unwrap();
    serde_json::to_writer(&mut buf, &json).unwrap();
    writeln!(buf).unwrap();
    // let data_vec = &data[n];
    writeln!(buf, "#lambda c var_c  m var_m").unwrap();
    //let data = &data_master[j];
    for k in 0..data.len(){
        
        writeln!(buf,"{} {} {} {} {}",lambda[k],data[k].var_c.mean,data[k].var_c.var,data[k].var_m.mean,data[k].var_m.var).unwrap();
    }
    
}
fn _alternate_writing(param:CriticalLambdaParams,json:Value,num_threads: Option<NonZeroUsize>,data_master:Vec<Vec<Measured>>,n_list:Vec<usize>,lambda:Vec<f64>){
    
    for j in 0..n_list.len(){
        let name = param.name("dat", num_threads,n_list[j]);
        println!("creating: {name}");
        let file = File::create(name).expect("unable to create file");
        let mut buf = BufWriter::new(file);
        write!(buf, "#").unwrap();
        serde_json::to_writer(&mut buf, &json).unwrap();
        writeln!(buf).unwrap();
        // let data_vec = &data[n];
        writeln!(buf, "#lambda c var_c  m var_m").unwrap();
        let data = &data_master[j];
        for k in 0..data.len(){
            
            writeln!(buf,"{} {} {} {} {}",lambda[k],data[k].var_c.mean,data[k].var_c.var,data[k].var_m.mean,data[k].var_m.var).unwrap();
        }

    }
    
    

}





fn _alternate_writing_old(param:CriticalLambdaParams,json:Value,num_threads: Option<NonZeroUsize>,data_master:Vec<Vec<MyVariance>>,n_list:Vec<usize>,lambda:Vec<f64>){
    
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
pub struct Measured
{
    pub var_m: MyVariance,
    pub var_c: MyVariance,
}