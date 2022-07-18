use net_ensembles::{sampling::{RewlBuilder, HistogramFast}, rand::{distributions::Uniform, distributions::Distribution}};


use {
    super::*,
    crate::json_parsing::*,
    crate::sir_model::*,
    rand_pcg::Pcg64,
    net_ensembles::rand::SeedableRng,
    humantime::format_duration,
    std::fs::File,
    std::io::{BufWriter, Write},
    net_ensembles::sampling::Rewl,
    crate::misc_types::*,
    //bincode::*
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
    

    let histograms = param.histograms.create(
        1, 
        param.system_size.get() as u32, 
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

    let energy = energy_function_returner(param.energy);

    let rewl = rewl_builder.greedy_build(energy);
    let duration = format_duration(instant.elapsed());
    println!("Finished greedy build after {duration}");
    
    //println!("{}",rewl.walkers().len());
    let allowed = param.allowed_seconds();
    let e = param.energy;
    let name_fn = move |interval: Option<usize>, end: &str, mode: LargeDeviationMode|
    {
        param.quick_name(interval, end, mode)
    };

    let opts = BALDLdOpts{
        quick_name: Box::new(name_fn),
        allowed_seconds: allowed,
        energy: e,
        value: vec![value],
        no_save: opt.no_save
    };

    execute_high_degree_helper(rewl, instant, opts)

}

fn print_min_max_intervalsize(hists: &[HistogramFast<u32>]){
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

pub fn execute_high_degree_helper(mut rewl: Rewl<BALargeDeviationWithLocks, Pcg64, HistogramFast<u32>, u32, MarkovStepWithLocks, ()>, 
    instant: std::time::Instant,
    opts: BALDLdOpts){

    let energy = energy_function_returner(opts.energy);

    rewl.simulate_while(energy, |_| instant.elapsed().as_secs() < opts.allowed_seconds);

    let unfinished_count = unsafe{
        rewl.ensemble_iter_mut()
            .fold(
                0, 
                |acc, m| m.unfinished_counter() + acc
            )};

    let total_simulations_count = unsafe{
        rewl.ensemble_iter_mut()
            .fold(
                0, 
                |acc, m| m.total_simulations_counter() + acc
            )};

    let unfinished_frac = unfinished_count as f64 / total_simulations_count as f64;

    rewl.walkers()
        .iter()
        .enumerate()
        .for_each(
            |(index, walker)|
            {
                println!("walker {index}, steps: {}, log_f {}", walker.step_count(), walker.log_f());
                let density = walker.log10_density();
                let name = (opts.quick_name)(Some(index), "dat", LargeDeviationMode::Rewl);
                println!("Created: {}", &name);
                let file = File::create(name)
                    .expect("unable to create file");
                let mut buf = BufWriter::new(file);
                for value in opts.value.iter()
                {
                    write!(buf, "#").unwrap();
                    serde_json::to_writer(&mut buf, value)
                        .unwrap();
                    writeln!(buf).unwrap();
                }

                writeln!(buf, "#hist log10").unwrap();
                writeln!(buf, "#walker {index}, steps: {}, log_f {}", walker.step_count(), walker.log_f())
                    .unwrap();

                writeln!(
                    buf, 
                    "# Replica_exchanges {}, proposed_replica_exchanges {} acceptance_rate {}",
                    walker.replica_exchanges(),
                    walker.proposed_replica_exchanges(),
                    walker.replica_exchange_frac()
                ).unwrap();

                writeln!(
                    buf,
                    "# Acceptance_rate_markov: {}",
                    walker.acceptance_rate_markov()
                ).unwrap();

                write!(buf, "#roundtrips of all walkers:").unwrap();

                for roundtrip in rewl.roundtrip_iter()
                {
                    write!(buf, " {roundtrip}").unwrap();
                }
                writeln!(buf).unwrap();
                writeln!(buf, "#All walker: Unfinished Sim: {unfinished_count} total: {total_simulations_count}, unfinished_frac {unfinished_frac}")
                    .unwrap();

                let hist = walker.hist();
                for (bin, density) in hist.bin_iter().zip(density)
                {
                    writeln!(buf, "{bin} {:e}", density).unwrap();
                }

                let name = (opts.quick_name)(Some(index), "hist", LargeDeviationMode::Rewl);
                let file = File::create(name)
                    .expect("unable to create file");
                let mut buf = BufWriter::new(file);
                let ensemble = rewl.get_ensemble(walker.id()).unwrap();
                for (bin, hits) in ensemble.hist_patient_zero().bin_hits_iter()
                {
                    let _ = writeln!(buf, "{} {}", bin, hits);
                }
            }
        );

    for roundtrip in rewl.roundtrip_iter() 
    {
        println!("rountrips: {roundtrip}")
    }

    if !opts.no_save{
        let save_name = (opts.quick_name)(None, "bincode", LargeDeviationMode::Rewl);
        println!("Create save: {save_name}");
        let string_vec: Vec<_> = opts.value.iter().map(serde_json::Value::to_string).collect();
        
        let save = (rewl, string_vec);

        let file = File::create(save_name)
            .expect("Unable to create save file");
        let buf = BufWriter::new(file);

        bincode::serialize_into(buf, &save)
            .unwrap();
    }
}

pub fn load_high_degree_rewl(opts: BALDContinueOpts, instant: std::time::Instant, no_save: bool){
    let allowed_seconds = opts.seconds();

    let (mut rewl, mut json_string) = deserialize_from_file(&opts.file_name);




    if let Some(step_size) = opts.change_step_size
    {
        for i in 0..
        {
            if rewl.change_step_size_of_interval(i, step_size.get()).is_err()
            {
                break;
            }
        }
    }

    // add the json from this call
    let json_val =  serde_json::to_value(opts).unwrap();
    json_string.push(json_val.to_string());
    

    let mut old_opts: BALDLDparam = serde_json::from_str(&json_string[0])
        .expect("Unable to deserialize old params");
    old_opts.times_plus_1();

    let mut jsons = vec![serde_json::to_value(old_opts.clone()).expect("unable to create json")];

    jsons.extend(
        json_string.iter()
            .skip(1)
            .filter_map(
                |s|
                {
                    let res: Result<BALDContinueOpts, _> = serde_json::from_str(s);
                    
                    match res
                    {
                        Ok(val) => Some(serde_json::to_value(val).expect("Unable to create json in continue")),
                        Err(_) => {
                            eprintln!("Cannot read old continue opts. Ignoring the error.");
                            None
                        }
                    }
                }
            )
    );

    let e = old_opts.energy;
    let name_fn = move |interval: Option<usize>, end: &str, mode: LargeDeviationMode|
    {
        old_opts.quick_name(interval, end, mode)
    };

    let opts = BALDLdOpts{
        quick_name: Box::new(name_fn),
        allowed_seconds,
        energy: e,
        value: jsons,
        no_save
    };

    execute_high_degree_helper(rewl, instant, opts)

}