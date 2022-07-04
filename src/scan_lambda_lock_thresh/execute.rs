use std::{ io::Write};


use rayon::iter::{ ParallelIterator};
use indicatif::ProgressIterator;



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
use crate::misc_types::*;
pub fn run_simulation(param:ScanLambdaThreshParams, json: Value, num_threads: Option<NonZeroUsize>){
    match param.graph_type{
        //GraphType::SmallWorld(_) => sim_small_world(param, json, num_threads),
        GraphType::Barabasi(_,_) => sim_barabasi(param, json, num_threads),
        _ => unimplemented!()
    }

}

pub fn sim_barabasi(param:ScanLambdaThreshParams, json: Value, num_threads: Option<NonZeroUsize>){
    let opt = BarabasiOptions::from_lambda_thresh_param(&param);
    let barabasi_world = opt.into();
    let mut model = SimpleSampleBarabasi::from_base(barabasi_world, param.sir_seed);

    //let range_lam = param.lambda_range.get_range();
    //let range_gam = param.gamma_range.get_range();
    //let lambda_range:Vec<_> = range_lam.iter().collect();
    //let gamma_range: Vec<_> = range_gam.iter().collect();

    let grid_2 = GridF64::new(
        param.lambda_range.get_range(),
        param.lockt_range.get_range(),   
    );

    //let grid_size:usize = param.lambda_range.get_range().iter().collect().len();

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

    let lockparams = param.lockdown;

    let lock_graph  = model.create_locked_down_network(lockparams);

    let per_threads = param.samples_per_step/k.get() as u64;
    let bar = crate::indication_bar(grid_2.len() as u64);


    let data:Vec<_> = grid_2.grid_point2d_iter().map(|pair|{

        let vec:Vec<_> = container.par_iter_mut().map(|model|{
            
            
            let mut lockparams = param.lockdown;

            //model.set_gamma(pair.y);
            model.set_lambda(pair.x);
            lockparams.set_lock_thresh(pair.y);
            //let t = lockparams.lock_threshold;
            //println!("{}",t);
            let mut sum_m_pt = 0.;
            let mut var_m_pt = 0.;

            let mut sum_c_pt = 0.;
            let mut var_c_pt = 0.;

            

            for _i in 0..per_threads{
                //let  = if pair.y == 0. && pair.x == 0.{
                  //  (1.,1.)
                    //let c = 1./system_size_fraction;
                
                //}
               // else if pair.y == 0.{
                    //    (param.system_size.get() as f64,param.system_size.get() as f64)
                  //  }
                //else{
               // (model.propagate_until_completion_max_with_lockdown(lock_graph.clone(),lockparams) as f64,model.calculate_ever_infected() as f64)};
               let (mut m,mut c) = (model.propagate_until_completion_max_with_lockdown(lock_graph.clone(),lockparams) as f64,model.calculate_ever_infected() as f64);
                
               
                m/= system_size_fraction;
                c/= system_size_fraction;
                
                sum_m_pt += m;
                sum_c_pt += c;
                
                var_m_pt += m*m;
                var_c_pt += c*c;
                


            }
            
            
                

            (sum_c_pt,sum_m_pt,var_c_pt,var_m_pt)
        }).collect(); //change goes here**;
        

        
        let mut sum_c = 0.;
        let mut var_c = 0.;
        let mut sum_m = 0.;
        let mut var_m = 0.;
        //println!("{}",vec.len());
        for (sum_c_i,sum_m_i,var_c_i,var_m_i) in vec{
            //println!("{}",var_c_i);
            sum_c += sum_c_i as f64;
            var_c += var_c_i as f64;
            sum_m += sum_m_i as f64;
            var_m += var_m_i as f64;
        }
        let p = (k.get() as u64 *per_threads) as f64 ;
        //println!("{}",p);
        sum_c /= p;
        var_c /= p;
        sum_m /= p;
        var_m /= p;

        //println!("sum c{}",sum_c);
        //println!("var c b{}",var_c);
        let var_c = var_c - sum_c*sum_c;
        let var_m = var_m - sum_m*sum_m;
        //println!("var c{}",var_c);
        //mean varian
        //let (c_vec, m_vec): (Vec<_>, Vec<_>) = vec.into_iter().map(|(a,b)| (a, b)).unzip();
        Measured{
            var_m: MyVariance{
                mean: sum_m,
                var: var_m
            },
            var_c: MyVariance{
                mean: sum_c,
                var: var_c}
            
        }

    }).progress_with(bar).collect();

    let map = GridMapF64Generic::from_vec_unchecked(grid_2,data);
    let measure = MeasureType::C;
    let name = param.name(measure,"mean.dat",num_threads);
    println!("Creating: {}", &name);
    let file = File::create(name).expect("unable to create file");  
    let mut buf = BufWriter::new(file);
    write!(buf,"#").unwrap();
    serde_json::to_writer(&mut buf, &json).expect("unable to write");
    writeln!(buf).unwrap();

    //let _ = writelin!(buf)

    map.write(
        &mut buf, 
        |measured| measured.var_c.mean()
    ).expect("unable to write");

    let name = param.name(measure, "var.dat", num_threads);
    println!("Creating: {}", &name);
    let file =  File::create(name)
        .expect("unable to create file");
    let mut buf = BufWriter::new(file);
    write!(buf, "#").unwrap();
    serde_json::to_writer(&mut buf, &json)
        .expect("unable to write");
    writeln!(buf).unwrap();

    map.write(
        &mut buf, 
        |measured| measured.var_c.variance_of_mean()
    ).expect("unable to write");

    // now write M
    let measure = MeasureType::M;

    let name = param.name(measure, "mean.dat", num_threads);
    println!("Creating: {}", &name);
    let file =  File::create(name)
        .expect("unable to create file");
    let mut buf = BufWriter::new(file);
    write!(buf, "#").unwrap();
    serde_json::to_writer(&mut buf, &json)
        .expect("unable to write");
    writeln!(buf).unwrap();


    map.write(
        &mut buf, 
        |measured| measured.var_m.mean()
    ).expect("unable to write");

    let name = param.name(measure, "var.dat", num_threads);
    println!("Creating: {}", &name);
    let file =  File::create(name)
        .expect("unable to create file");
    let mut buf = BufWriter::new(file);
    write!(buf, "#").unwrap();
    serde_json::to_writer(&mut buf, &json)
        .expect("unable to write");
    writeln!(buf).unwrap();


    map.write(
        &mut buf, 
        |measured| measured.var_m.variance_of_mean()
    ).expect("unable to write");

}
