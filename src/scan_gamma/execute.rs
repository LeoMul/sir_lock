use std::{ops::DerefMut, io::Write};
use crate::json_parsing::*;
use indicatif::ParallelProgressIterator;
use rayon::iter::{IntoParallelIterator, ParallelIterator, IntoParallelRefIterator};

use{
    super::*,
    crate::{GraphType, sir_model::*},
    serde_json::Value,
    std::{num::*, sync::Mutex, fs::File, io::BufWriter},
    rand_pcg::Pcg64,
    net_ensembles::rand::SeedableRng,
    crate::lockdown_methods::*,
    rayon::prelude::*,
    crate::stats_methods::MyVariance
};

//use crate:: sir_model::*;

pub fn run_simulation(param:ScanGammaParams, json: Value, num_threads: Option<NonZeroUsize>){
    match param.graph_type{
        GraphType::SmallWorld(_) => sim_small_world(param, json, num_threads),
        GraphType::Barabasi(_,_) => sim_barabasi(param, json, num_threads),
        _ => unimplemented!()
    }

}

fn sim_barabasi(param:ScanGammaParams,json:Value,num_threads:Option<NonZeroUsize>){
    let opt = BarabasiOptions::from_gamma_scan_param(&param);
    let barabasi_world = opt.into();
    let mut model = SimpleSampleBarabasi::from_base(barabasi_world, param.sir_seed,param.initial_infected);
    let range = param.gamma_range.get_range();
    let gamma_range:Vec<_> = range.iter().collect();
    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    rayon::ThreadPoolBuilder::new().num_threads(k.get()).build_global().unwrap();

    let mut lock_graph  = model.create_locked_down_network(param.lockdown);

    //let y = scanning_lambda_function_dynamic(&param, &json, num_threads, &model);
    let y = scanning_gamma_function_static(&param, &json, num_threads, &model,param.lockdown,&mut lock_graph);
    let y_clone = y.clone();
    let samples = Samples{
        gamma: gamma_range,
        var: y
    };
    if !param.compare{
        let name = param.name("dat", num_threads);
        println!("Creating: {}", &name);
        let file = File::create(name)
            .expect("unable to create file");
        let mut buf = BufWriter::new(file);
        write_json(&mut buf, &json);
        samples.write(buf).unwrap()

    }
    else{
        println!("comparing");

        let gamma_range:Vec<_> = range.iter().collect();
        let mut lock_graph = model.create_locked_down_network(param.complock);
        println!("{}",lock_graph.edge_count());
        //let y2 = scanning_lambda_function_dynamic(&param, &json, num_threads, &model);
        let y2 = scanning_gamma_function_static(&param, &json, num_threads, &model,param.complock,&mut lock_graph);

        let samples2 = Samples2{
            gamma:gamma_range,
            var: y_clone,
            var2: y2
        };
        let name = param.name("dat", num_threads);
        println!("Creating: {}", &name);
        let file = File::create(name)
            .expect("unable to create file");
        let mut buf = BufWriter::new(file);
        write_json(&mut buf, &json);
        samples2.write(buf).unwrap()

    }

    


}

fn sim_small_world(param: ScanGammaParams, json: Value, num_threads:Option<NonZeroUsize>){
    let opt = SWOptions::from_gamma_scan_param(&param);
    let small_world = opt.into();
    let mut model = SimpleSampleSW::from_base(small_world, param.sir_seed,param.initial_infected);
    let locked_down_graph = model.create_locked_down_network(param.lockdown);
    let range = param.gamma_range.get_range();
    let gamma_range:Vec<_> = range.iter().collect();
    //println!("locked down {}",locked_down_graph.edge_count());


    //The following allows for the multi-core functionality.

    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    rayon::ThreadPoolBuilder::new().num_threads(k.get()).build_global().unwrap();

    let mut rng = Pcg64::seed_from_u64(param.sir_seed);

    let container: Vec<_> = (0..k.get()).map(
        |_|
        {
            let mut model = model.clone();
            let locked_down_graph = locked_down_graph.clone();
            model.reseed_sir_rng(&mut rng);
            Mutex::new(
                (model,locked_down_graph )
            )

        }
    ).collect();

    let end = gamma_range.len()/k.get();
    let bar = crate::indication_bar(gamma_range.len() as u64);

    let iterator = (0..=end)
        .into_par_iter()
        .flat_map(
                |i|
                {
                    let start = i*container.len();
                    container.par_iter()
                        .zip(gamma_range[ start..].par_iter())
                }
            );
    let system_size_fraction_bool = if param.fraction{
        Some(param.system_size.get() as f64)
    }else{
        None
    };
    let factor = if param.fraction{
        param.system_size.get() as f64
    }
    else{ 
        1.
    };
    let y: Vec<MyVariance> = iterator
        .progress_with(bar)
        .map(
            |(container, &gamma)|
            {
                let mut lock = container.lock()
                    .expect("unable to lock");
                let inner = lock.deref_mut();
                //let vaccine_rng = &mut inner.0;
                let model = &mut inner.0;
                let mut locked_down_graph = &mut inner.1;
                //let vaccine_list_helper = &mut inner.2;

                model.set_gamma(gamma);
                //println!("{lambda}");
                let vals: Vec<_> = (0..param.samples_per_step)
                    .map(
                        |_|
                        {
                            //vaccine_list_helper.randomize(vaccine_rng);
                            if gamma == 0.{
                                (param.system_size.get() as f64/factor) as u32
                            }
                            else{
                            //let vaccine_list = vaccine_list_helper.get_vaccine_list(param.vaccine_doses, model.ensemble().graph());
                            let res = model.propagate_until_completion_max_with_lockdown(&mut locked_down_graph,param.lockdown) as u32;
                            //println!("locked down {}",locked_down_graph.edge_count());
                            if param.measure.is_c() 
                            {
                                model.calculate_ever_infected() as u32
                            } else {
                                res
                            }}
                        }
                    ).collect();

                MyVariance::from_slice(&vals, system_size_fraction_bool)
            }
        ).collect();    
    //let x = find_critical_lambda_one_data_set(&y, &lambda_range);
    //println!("{}",x);
    let y_clone = y.clone();
    let samples = Samples{
        gamma: gamma_range,
        var: y
    };
    //let y_clone = y.clone();
    if !param.compare{
    let name = param.name("dat", num_threads);
    println!("Creating: {}", &name);
    let file = File::create(name)
        .expect("unable to create file");
    let mut buf = BufWriter::new(file);
    write_json(&mut buf, &json);
    samples.write(buf).unwrap()}
    else{
        println!("comparing");
        

        let gamma_range:Vec<_> = range.iter().collect();
        let new_locked_down_graph = model.create_locked_down_network(param.complock);
        let container: Vec<_> = (0..k.get()).map(
            |_|
            {
                let mut model = model.clone();
                let new_locked_down_graph = new_locked_down_graph.clone();
                model.reseed_sir_rng(&mut rng);
                Mutex::new(
                    (model,new_locked_down_graph )
                )
    
            }
        ).collect();
        //println!("{}",locked_down_graph.edge_count());
        //let y2 = scanning_lambda_function_dynamic(&param, &json, num_threads, &model);
        let iterator = (0..=end)
        .into_par_iter()
        .flat_map(
                |i|
                {
                    let start = i*container.len();
                    container.par_iter()
                        .zip(gamma_range[ start..].par_iter())
                }
            );
        let system_size_fraction_bool = if param.fraction{
            Some(param.system_size.get() as f64)
        }else{
            None
        };
        let bar = crate::indication_bar(gamma_range.len() as u64);
        let factor = if param.fraction{
            param.system_size.get() as f64
        }
        else{ 
            1.
        };
        let y2: Vec<MyVariance> = iterator
        .progress_with(bar)
        .map(
            |(container, &gamma)|
            {
                let mut lock = container.lock()
                    .expect("unable to lock");
                let inner = lock.deref_mut();
                //let vaccine_rng = &mut inner.0;
                let model = &mut inner.0;
                let mut locked_down_graph = &mut inner.1;
                //let vaccine_list_helper = &mut inner.2;

                model.set_lambda(gamma);

                let vals: Vec<_> = (0..param.samples_per_step)
                    .map(
                        |_|
                        {
                            //vaccine_list_helper.randomize(vaccine_rng);
                            
                            //let vaccine_list = vaccine_list_helper.get_vaccine_list(param.vaccine_doses, model.ensemble().graph());
                            //println!("sim");
                            //println!("locked down {}",locked_down_graph.edge_count());
                            if gamma == 0.{
                                (param.system_size.get() as f64/factor) as u32
                            }
                            else{
                            let res = model.propagate_until_completion_max_with_lockdown(&mut locked_down_graph,param.lockdown) as u32;
                            if param.measure.is_c() 
                            {
                                model.calculate_ever_infected() as u32
                            } else {
                                res
                            }}
                        }
                    ).collect();

                MyVariance::from_slice(&vals, system_size_fraction_bool)
            }
        ).collect();

        let samples2 = Samples2{
            gamma:gamma_range,
            var: y_clone,
            var2: y2
        };
        let name = param.name("dat", num_threads);
        println!("Creating: {}", &name);
        let file = File::create(name)
            .expect("unable to create file");
        let mut buf = BufWriter::new(file);
        write_json(&mut buf, &json);
        samples2.write(buf).unwrap()

    }

}
fn scanning_gamma_function_static(param:&ScanGammaParams,_json:&Value,num_threads:Option<NonZeroUsize>,model:&SimpleSampleBarabasi,lockparams:LockdownParameters,lockgraph:&mut GenGraphSIR) -> Vec<MyVariance>{
    //,lockgraph:&GenGraphSIR
    let mut rng = Pcg64::seed_from_u64(param.sir_seed);

    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    let range = param.gamma_range.get_range();

    let gamma_range:Vec<_> = range.iter().collect();
    let container: Vec<_> = (0..k.get()).map(
        |_|
        {
            let mut model = model.clone();
            model.reseed_sir_rng(&mut rng);
            let lock_down_graph = lockgraph.clone();
            Mutex::new(
                (model,lock_down_graph)
            )

        }
    ).collect();

    let end = gamma_range.len()/k.get();
    let bar = crate::indication_bar(gamma_range.len() as u64);

    let iterator = (0..=end)
        .into_par_iter()
        .flat_map(
                |i|
                {
                    let start = i*container.len();
                    container.par_iter()
                        .zip(gamma_range[ start..].par_iter())
                }
            );
    
    let system_size_fraction_bool = if param.fraction{
        Some(param.system_size.get() as f64)
    }else{
        None
    };

    let factor = if param.fraction{
        param.system_size.get() as f64
    }
    else{ 
        1.
    };



    let y: Vec<MyVariance> = iterator
        .progress_with(bar)
        .map(
            |(container, &gamma)|
            {
                let mut lock = container.lock()
                    .expect("unable to lock");
                let inner = lock.deref_mut();
                //let vaccine_rng = &mut inner.0;
                let model = &mut inner.0;
                //let vaccine_list_helper = &mut inner.2;
                //let locked_down = &mut inner.1;
                let lockgraph = &mut inner.1;
                model.set_gamma(gamma);

                let vals: Vec<_> = (0..param.samples_per_step)
                    .map(
                        |_|
                        {
                            //vaccine_list_helper.randomize(vaccine_rng);
                            if gamma == 0.{
                                (param.system_size.get() as f64/factor) as u32
                            }
                            else{ 
                            //let vaccine_list = vaccine_list_helper.get_vaccine_list(param.vaccine_doses, model.ensemble().graph());
                            let res = model.propagate_until_completion_max_with_lockdown(lockgraph,lockparams) as u32;
                            if param.measure.is_c() 
                            {
                                model.calculate_ever_infected() as u32
                            } else {
                                res
                            }}
                        }
                    ).collect();

                MyVariance::from_slice(&vals, system_size_fraction_bool)
            }
        ).collect();
        y
}










fn _find_critical_lambda_one_data_set(vec:&Vec<MyVariance>,lambdavals:&[f64])->f64{
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


pub struct Samples{
    gamma: Vec<f64>,
    var: Vec<MyVariance>
}

impl Samples{
    fn write<W>(&self, mut writer: W) -> std::io::Result<()>
    where W: Write
    {
        writeln!(writer, "#gammaa mean variance")?;

        for (lambda, var) in self.gamma.iter().zip(self.var.iter())
        {
            writeln!(writer, "{:e} {:e} {:e}", lambda, var.mean(), var.variance_of_mean())?
        }
        Ok(())
    }
}
pub struct Samples2{
    gamma: Vec<f64>,
    var: Vec<MyVariance>,
    var2: Vec<MyVariance>
}

impl Samples2{
    fn write<W>(&self, mut writer: W) -> std::io::Result<()>
    where W: Write
    {
        writeln!(writer, "#gamma meanWL varianceWL meanNL varianceNL")?;

        for j in 0..self.gamma.len(){
            writeln!(writer, "{:e} {:e} {:e} {:e} {:e}", self.gamma[j], self.var[j].mean(), self.var[j].variance_of_mean(),self.var2[j].mean(),self.var2[j].variance_of_mean())?

        }
        println!("{}",self.gamma.len()-self.var.len());
        
        Ok(())
    }
}