use net_ensembles::{sampling::{RewlBuilder, HistogramFast}, rand::{distributions::Uniform, distributions::Distribution}};


use {
    super::*,
    crate::json_parsing::*,
    crate::sir_model::*,
    rand_pcg::Pcg64,
    net_ensembles::rand::SeedableRng,
    humantime::format_duration,
    std::fs::File,
    std::io::{BufWriter, Write, BufReader},
    net_ensembles::sampling::Rewl,
    std::num::*,
    serde::{Serialize, Deserialize},
    crate::misc_types::*
};


pub fn execute_large_dev(opt: BALDOptsLD, instant: std::time::Instant){
    //Initialising number of threads
    if let Some(num) = opt.num_threads
    {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num.get())
            .build_global()
            .unwrap();
    }

    let (param, value): (BALDLDparam, _) = parse(opt.json.as_ref());

    let (m,n) = match &param.graph_type {
        GraphType::Barabasi(m,n) => (*m,*n),
        _ => unimplemented!()
    };

    let base_options = BarabasiOptions{
        graph_seed: param.graph_seed,
        lambda: param.lambda,
        gamma: param.recovery_prob,
        system_size: param.system_size,
        m,
        source_n:n}; 

    let ba:BarabasiModel = base_options.into();
    let model = BALargeDeviation::new(ba.clone(),param.large_deviation_param,param.lockdownparams);
    let mut ld_model = BALargeDeviationWithLocks::new(model);
    ld_model.infect_patient();
    let cluster_size = 1;
    //let mut cluster_size = ld_model.cluster_size_patient(ld_model.patient_zero());


        //ask Yannick if we need this, the function involves vaccinating, but it might have some use. Ask exactly what it is supposed to be doing.


    //if let Some(steps) = param.histograms.greedy_search_steps()
    //{
      //  println!("Start greedy cluster size search. Current estimate: {cluster_size}, steps: {}", steps);
        //cluster_size = ld_model.greedy_search_cluster_size(steps.get(), false);
    //}
    //let now = instant.elapsed();
    //println!(
      //  "{} have passed. At least {cluster_size} nodes can still be infected at once",
        //humantime::format_duration(now)
    //);

    let histograms = param.histograms.create(
        1, 
        cluster_size as u32, 
        param.walkers_per_interval
    );
    print_min_max_intervalsize(&histograms);

    let mut markov_seeder = Pcg64::seed_from_u64(param.large_deviation_param.markov_seed);
    let uniform = Uniform::new(0, 90357094723984_u64);

    let ensembles = histograms
        .iter()
        .map(
            |_|
            {
                let mut l_param = param.large_deviation_param;
                l_param.markov_seed = uniform.sample(&mut markov_seeder);
                let ld = BALargeDeviation::new(ba.clone(), l_param,param.lockdownparams);
                BALargeDeviationWithLocks::new(ld)
            }
        ).collect();

    let rewl_builder = RewlBuilder::from_ensemble_vec(
        ensembles,
        histograms,
        param.step_size,
        param.sweep_size,
        param.walkers_per_interval,
        param.f_threshold
    ).expect("unable to create rewl builder");

    println!("Start greedy build");

    let energy = |model: &mut BALargeDeviationWithLocks|
    {
        let measure = model.ld_energy_m();
        Some(
            match param.energy{
                MeasureType::M => measure,
                MeasureType::C => model.calculate_ever_infected() as u32
            }
        )
    };

    let rewl = rewl_builder.greedy_build(energy);
    let duration = format_duration(instant.elapsed());
    println!("Finished greedy build after {duration}");

}

fn print_min_max_intervalsize(hists: &[HistogramFast<u32>])
{
    let min_interval_size = hists.iter()
        .map(|h| h.range_inclusive().size_hint().0)
        .min()
        .unwrap();

    let max_interval_size = hists.iter()
        .map(|h| h.range_inclusive().size_hint().0)
        .max()
        .unwrap();

    
    println!("min_interval len: {min_interval_size} max interval size {max_interval_size}");
}