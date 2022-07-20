use {
    super::parser::*,
    serde_json::Value,
    std::{num::*},
    crate::*,
    crate::sir_model::*,
    rayon::prelude::*,
    rand_pcg::Pcg64,
    net_ensembles::rand::SeedableRng,
    net_ensembles::sampling::{HistU32Fast},
    crate::large_deviations::*,
};

pub fn execute_sir(
    param: SimpleSampleParam,
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

pub fn execute_ba(param: SimpleSampleParam,json: Value,num_threads: Option<NonZeroUsize>){
    let opt = BarabasiOptions::from_simple_sample(&param);
    let ba = opt.into();
    let model = SimpleSampleBarabasi::from_base(ba,param.sir_seed);

    let j = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    // limit number of threads to j
    rayon::ThreadPoolBuilder::new().num_threads(j.get()).build_global().unwrap();
    let mut sir_rng_2 = Pcg64::seed_from_u64(param.sir_seed);

    let samples_per_thread = param.samples / j.get();

    

    //use this in other programs.
    let mut rngs: Vec<_> = (0..j.get())
        .map(
            |_| 
                {
                    
                        Pcg64::from_rng(&mut sir_rng_2).unwrap()
                    
                }
            )
        .collect();
    let bar = indication_bar(param.samples as u64);

    let mut histograms: Vec<_> = rngs.par_iter_mut()
        .map(
            |sir_rng|
            {
                let mut model = model.clone();
                model.reseed_sir_rng(sir_rng);
                let mut hist_c = HistU32Fast::new_inclusive(1, param.system_size.get() as u32)
                    .unwrap();
                let mut hist_m = hist_c.clone();
                
                for i in 0..samples_per_thread
                {
                    let lockgraph = model.create_locked_down_network(param.lockdown);
                    //let lockgraph = create_locked_down_network_from_pair_list(&pairs_struct, model.ensemble.graph());

                    let m = model.propagate_until_completion_max_with_lockdown(lockgraph,param.lockdown) as u32;
                    hist_m.increment_quiet(m);
                    let c = model.calculate_ever_infected() as u32;
                    hist_c.increment_quiet(c);
                    if i % 1000 == 0 {
                        bar.inc(1000)
                    }
                }
                (hist_c, hist_m)
            }
        ).collect();
    bar.finish_with_message("Done");

    let mut combined_hist = histograms.pop()
        .unwrap();
    
    for hist in histograms{
        combined_hist.0.try_add(&hist.0)
            .unwrap();
        combined_hist.1.try_add(&hist.1)
            .unwrap();
    }

    let (hist_c, hist_m) = combined_hist;

    let name = param.quick_name(
        num_threads,
        MeasureType::C
    );

    hist_to_file(&hist_c, name, &json);

    let name = param.quick_name(
        num_threads,
        MeasureType::M
    );

    hist_to_file(&hist_m, name, &json);
}

pub fn execute_sw(param: SimpleSampleParam,json: Value,num_threads: Option<NonZeroUsize>){
    let opt = SWOptions::from_simple_sample(&param);
    let ba = opt.into();
    let model = SimpleSampleSW::from_base(ba,param.sir_seed);

    let j = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    // limit number of threads to j
    rayon::ThreadPoolBuilder::new().num_threads(j.get()).build_global().unwrap();
    let mut sir_rng_2 = Pcg64::seed_from_u64(param.sir_seed);

    let samples_per_thread = param.samples / j.get();

    

    //use this in other programs.
    let mut rngs: Vec<_> = (0..j.get())
        .map(
            |_| 
                {
                    
                        Pcg64::from_rng(&mut sir_rng_2).unwrap()
                    
                }
            )
        .collect();
    let bar = indication_bar(param.samples as u64);

    let mut histograms: Vec<_> = rngs.par_iter_mut()
        .map(
            |sir_rng|
            {
                let mut model = model.clone();
                model.reseed_sir_rng(sir_rng);
                let mut hist_c = HistU32Fast::new_inclusive(1, param.system_size.get() as u32)
                    .unwrap();
                let mut hist_m = hist_c.clone();
                
                for i in 0..samples_per_thread
                {
                    let lockgraph = model.create_locked_down_network(param.lockdown);

                    let m = model.propagate_until_completion_max_with_lockdown(lockgraph,param.lockdown) as u32;
                    hist_m.increment_quiet(m);
                    let c = model.calculate_ever_infected() as u32;
                    hist_c.increment_quiet(c);
                    if i % 1000 == 0 {
                        bar.inc(1000)
                    }
                }
                (hist_c, hist_m)
            }
        ).collect();
    bar.finish_with_message("Done");

    let mut combined_hist = histograms.pop()
        .unwrap();
    
    for hist in histograms{
        combined_hist.0.try_add(&hist.0)
            .unwrap();
        combined_hist.1.try_add(&hist.1)
            .unwrap();
    }

    let (hist_c, hist_m) = combined_hist;

    let name = param.quick_name(
        num_threads,
        MeasureType::C
    );

    hist_to_file(&hist_c, name, &json);

    let name = param.quick_name(
        num_threads,
        MeasureType::M
    );

    hist_to_file(&hist_m, name, &json);

}