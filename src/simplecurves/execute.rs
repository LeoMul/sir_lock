use {
    super::parser::*,
    serde_json::Value,
    std::{num::*},
    crate::*,
    crate::sir_model::*,
    //rand_pcg::Pcg64,
    //net_ensembles::rand::SeedableRng,
};

pub fn execute_sir(
    param: SimpleCurvesParam,
    json: Value,
    num_threads: Option<NonZeroUsize>
)
{
    match param.graph_type{
        GraphType::Barabasi(_,_) => execute_ba(param, json, num_threads),
        GraphType::SmallWorld(_) => execute_sw(param, json, num_threads),
        _ => unimplemented!()
    }
}

pub fn execute_sw(param: SimpleCurvesParam,json: Value,num_threads: Option<NonZeroUsize>){
    let opt = SWOptions::from_simplecurves(&param);
    let ba = opt.into();
    let mut model = SimpleSampleSW::from_base(ba,param.sir_seed,param.initial_infected);

    let j = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    // limit number of threads to j
    rayon::ThreadPoolBuilder::new().num_threads(j.get()).build_global().unwrap();
    //let mut sir_rng_2 = Pcg64::seed_from_u64(param.sir_seed);


    
    let lockdownparams = param.lockdown;
    //use this in other programs.
    //let mut rngs: Vec<_> = (0..j.get())
    //    .map(
    //        |_| 
    //            {
    //                
    //                    Pcg64::from_rng(&mut sir_rng_2).unwrap()
    //                
    //            }
    //        )
    //    .collect();
    //let bar = indication_bar(param.samples as u64);
    let name = param.quick_name();
    let mut writer = SirWriter::new(&name, 1);
                writer.write_header(&[json]).unwrap();
    for _ in 0..param.samples{
        let mut post_locked_down_graph = model.create_locked_down_network(lockdownparams);
        writer.write_energy(5, 10).unwrap();
        model.propagate_until_completion_max_with_lockdown_printing(&mut post_locked_down_graph,lockdownparams,&mut writer)
    }
}
pub fn execute_ba(param: SimpleCurvesParam,json: Value,num_threads: Option<NonZeroUsize>){
    let opt = BarabasiOptions::from_simplecurves(&param);
    let ba = opt.into();
    let mut model = SimpleSampleBarabasi::from_base(ba,param.sir_seed,param.initial_infected);

    let j = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    // limit number of threads to j
    rayon::ThreadPoolBuilder::new().num_threads(j.get()).build_global().unwrap();
    //let mut sir_rng_2 = Pcg64::seed_from_u64(param.sir_seed);


    
    let lockdownparams = param.lockdown;
    //use this in other programs.
    //let mut rngs: Vec<_> = (0..j.get())
    //    .map(
    //        |_| 
    //            {
    //                
    //                    Pcg64::from_rng(&mut sir_rng_2).unwrap()
    //                
    //            }
    //        )
    //    .collect();
    //let bar = indication_bar(param.samples as u64);
    let name = param.quick_name();
    let mut writer = SirWriter::new(&name, 1);
                writer.write_header(&[json]).unwrap();
    for _ in 0..param.samples{
        let mut post_locked_down_graph = model.create_locked_down_network(lockdownparams);
        writer.write_energy(5, 10).unwrap();
        model.propagate_until_completion_max_with_lockdown_printing(&mut post_locked_down_graph,lockdownparams,&mut writer)
    }
}