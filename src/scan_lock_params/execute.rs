use{
    super::*,
    std::num::*,
    serde_json::Value,

};
use std::{ io::Write};
use rayon::prelude::*;
//use rayon::iter::{ParallelIterator, IntoParallelRefIterator};
use indicatif::ProgressIterator;

use{
    crate::{GraphType, sir_model::*},
    std::{fs::File, io::BufWriter},
    rand_pcg::Pcg64,
    net_ensembles::rand::SeedableRng,
    crate::stats_methods::MyVariance,
    crate::lockdown_methods::*,
};

//use serde_json::from_slice;


pub struct Measured
{
    pub var_m: MyVariance,
    pub var_c: MyVariance,
}

pub fn run_simulation(param:ScanLockParams, json: Value, num_threads: Option<NonZeroUsize>){
    match param.graph_type{
        GraphType::SmallWorld(_) => sim_small_world(param, json, num_threads),
        GraphType::Barabasi(_,_) => sim_ba(param, json, num_threads),
        _ => unimplemented!()
    }

}

pub fn sim_small_world(param:ScanLockParams,json:Value,num_threads:Option<NonZeroUsize>){


    let opt = SWOptions::from_lock_scan_param(&param);
    let world = opt.into();
    let model = SimpleSampleSW::from_base(world, param.sir_seed,param.initial_infected);


    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    rayon::ThreadPoolBuilder::new().num_threads(k.get()).build_global().unwrap();
    let severities:Vec<_> = param.lockdown_severities.to_vec();

    //let mut lockparams = param.lockdown;
    let thresh_range = param.thresholdrange.get_range();
    let thresh_vec:Vec<_> = thresh_range.iter().collect();

    let mut sir_rng_2 = Pcg64::seed_from_u64(param.sir_seed);
    //let mut graph_rng = Pcg64::seed_from_u64(param.graph_seed);

    let mut container: Vec<_> = (0..k.get()).map(
        |_|
        {
            let mut model = model.clone();
            model.reseed_sir_rng(&mut sir_rng_2);
            
            model
            

        }
    ).collect();


    for severity in severities{
        println!("Commencing sampling for {}",severity);
        let system_size_frac = if param.fraction{
            param.system_size.get() as f64
        } else{
            1.
        };
        
        let per_thread = param.num_samples/k.get() as u64;
        
        let bar = crate::indication_bar(thresh_vec.len() as u64);
        
        let data:Vec<_> = thresh_vec.iter().map(|t|{
            let rel = get_release_threshold(param.releaseparams,*t);
            //println!("{rel},{t}");
            let lockparams = LockdownParameters{
                lock_style:LockdownType::Random(severity),
                release_threshold: rel,
                lock_threshold:*t,
                };

            let thread_vecs:Vec<_> = container.par_iter_mut().map(|cont|{
                //let mut lock = cont.lock()
                //    .expect("unable to lock");
                //let inner = lock.deref_mut();
                //let vaccine_rng = &mut inner.0;
                let model = cont;//&mut inner.0;
                //let iter = (0..per_thread).into_iter();
                //individual_thread data
                let mut av_m = 0.;
                let mut var_m = 0.;
                let mut av_c = 0.;
                let mut var_c =0.;
                for _ in 0..per_thread{

                    let mut lockgraph = model.create_locked_down_network(lockparams);
                    let m = model.propagate_until_completion_max_with_lockdown(&mut lockgraph,lockparams) as f64 /system_size_frac;
                    let c = model.calculate_ever_infected() as f64 /system_size_frac;

                    av_m += m;
                    av_c += c;
                    var_m += m*m;
                    var_c += c*c;
                }
                (av_m,var_m,av_c,var_c)
            }
            ).collect();
            let mut average_m = 0.;
            let mut variance_m = 0.;
            let mut average_c = 0.;
            let mut variance_c = 0.;
            for (i,j,ke,l) in thread_vecs{
                average_m += i;
                variance_m += j;
                average_c += ke;
                variance_c += l;
            }
            let num_samples = (per_thread*k.get() as u64) as f64;
            average_m /= num_samples;
            average_c /= num_samples;
            variance_c /= num_samples;
            variance_m /= num_samples;

            variance_c -= average_c*average_c;
            variance_m -= average_m*average_m;
            Measured{
                var_c:MyVariance{
                    var:variance_c,
                    mean:average_c,
                },
                var_m:MyVariance{
                    var:variance_m,
                    mean:average_m
                }
            }

        }).progress_with(bar).collect();

        writing(&param,&json,num_threads,&data,severity,&thresh_vec);
        

    }


}

pub fn sim_ba(param:ScanLockParams,json:Value,num_threads:Option<NonZeroUsize>){


    let opt = BarabasiOptions::from_lock_scan_param(&param);
    let world = opt.into();
    let model = SimpleSampleBarabasi::from_base(world, param.sir_seed,param.initial_infected);


    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    rayon::ThreadPoolBuilder::new().num_threads(k.get()).build_global().unwrap();
    let severities:Vec<_> = param.lockdown_severities.to_vec();

    //let mut lockparams = param.lockdown;
    let thresh_range = param.thresholdrange.get_range();
    let thresh_vec:Vec<_> = thresh_range.iter().collect();

    let mut sir_rng_2 = Pcg64::seed_from_u64(param.sir_seed);
    //let mut graph_rng = Pcg64::seed_from_u64(param.graph_seed);

    let mut container: Vec<_> = (0..k.get()).map(
        |_|
        {
            let mut model = model.clone();
            model.reseed_sir_rng(&mut sir_rng_2);
            
            model
            

        }
    ).collect();


    for severity in severities{
        println!("Commencing sampling for {}",severity);
        let system_size_frac = if param.fraction{
            param.system_size.get() as f64
        } else{
            1.
        };
        
        let per_thread = param.num_samples/k.get() as u64;
        
        let bar = crate::indication_bar(thresh_vec.len() as u64);
        
        let data:Vec<_> = thresh_vec.iter().map(|t|{

            let lockparams = LockdownParameters{
                lock_style:LockdownType::Random(severity),
                release_threshold: get_release_threshold(param.releaseparams,*t),
                lock_threshold:*t,
                };

            let thread_vecs:Vec<_> = container.par_iter_mut().map(|cont|{
                //let mut lock = cont.lock()
                //    .expect("unable to lock");
                //let inner = lock.deref_mut();
                //let vaccine_rng = &mut inner.0;
                let model = cont;//&mut inner.0;
                //let iter = (0..per_thread).into_iter();
                //individual_thread data
                let mut av_m = 0.;
                let mut var_m = 0.;
                let mut av_c = 0.;
                let mut var_c =0.;
                for _ in 0..per_thread{

                    let mut lockgraph = model.create_locked_down_network(lockparams);
                    let m = model.propagate_until_completion_max_with_lockdown(&mut lockgraph,lockparams) as f64 /system_size_frac;
                    let c = model.calculate_ever_infected() as f64 /system_size_frac;

                    av_m += m;
                    av_c += c;
                    var_m += m*m;
                    var_c += c*c;
                }
                (av_m,var_m,av_c,var_c)
            }
            ).collect();
            let mut average_m = 0.;
            let mut variance_m = 0.;
            let mut average_c = 0.;
            let mut variance_c = 0.;
            for (i,j,ke,l) in thread_vecs{
                average_m += i;
                variance_m += j;
                average_c += ke;
                variance_c += l;
            }
            let num_samples = (per_thread*k.get() as u64) as f64;
            average_m /= num_samples;
            average_c /= num_samples;
            variance_c /= num_samples;
            variance_m /= num_samples;

            variance_c -= average_c*average_c;
            variance_m -= average_m*average_m;
            Measured{
                var_c:MyVariance{
                    var:variance_c,
                    mean:average_c,
                },
                var_m:MyVariance{
                    var:variance_m,
                    mean:average_m
                }
            }

        }).progress_with(bar).collect();

        writing(&param,&json,num_threads,&data,severity,&thresh_vec);
        

    }


}

fn writing(param:&ScanLockParams,json:&Value,num_threads: Option<NonZeroUsize>,data:&Vec<Measured>,severity:f64,thrsh:&[f64]){
    let name = param.name("dat", num_threads,severity);
    println!("creating: {name}");
    let file = File::create(name).expect("unable to create file");
    let mut buf = BufWriter::new(file);
    write!(buf, "#").unwrap();
    serde_json::to_writer(&mut buf, &json).unwrap();
    writeln!(buf).unwrap();
    // let data_vec = &data[n];
    writeln!(buf, "#thresh c c_var m m_var ").unwrap();
    
    for k in 0..data.len(){
        
        writeln!(buf,"{} {} {} {} {}",thrsh[k],data[k].var_c.mean,data[k].var_c.var,data[k].var_m.mean,data[k].var_m.var).unwrap();
    }}


