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

pub fn execute_large_dev(opt: LDOptsLD, instant: std::time::Instant){
    let (param, value): (LDLDparam, _) = parse(opt.json.as_ref());

    if param.change_energy_res && param.large_deviation_param.initial_infected%2 == param.system_size.get()%2{
        println!("Error! The changing of resolution requires two energies per bin. The initial infected and system size must then have opposite parity");
        panic!()

    }
    else{
        match param.graph_type{
            GraphType::SmallWorld(r) => execute_sw(opt, instant,r,param,value),
            GraphType::Barabasi(m,n) => execute_ba(opt, instant,m,n,param,value),
            _ => unimplemented!()
        }
    }
}



pub fn execute_ba(opt: LDOptsLD, instant: std::time::Instant,m:usize,n:usize,param:LDLDparam,value:serde_json::Value){
    //Initialising number of threads
    if let Some(num) = opt.num_threads
    {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num.get())
            .build_global()
            .unwrap();
    }
    
    
    let boolean_res:bool = (param.change_energy_res).clone();

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
    ld_model.infect_initial_patients();
    

    let histograms =  if !boolean_res{
        param.histograms.create(
            ld_model.initial_infected as u32, 
            param.system_size.get() as u32, 
            param.walkers_per_interval
        )}
    else{
        param.histograms.create(
            (ld_model.initial_infected as u32-1)/2, 
            (param.system_size.get() as u32-1)/2, 
            param.walkers_per_interval
            )
    };
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

    let energy = energy_function_returner_ba(param.energy,param.change_energy_res);

    let rewl = rewl_builder.greedy_build(energy);
    let duration = format_duration(instant.elapsed());
    println!("Finished greedy build after {duration}");
    //panic!();
    //println!("{}",rewl.walkers().len());
    let allowed = param.allowed_seconds();
    let e = param.energy;
    let name_fn = move |interval: Option<usize>, end: &str, mode: LargeDeviationMode|
    {
        param.quick_name(interval, end, mode)
    };

    let opts = LDLdOpts{
        quick_name: Box::new(name_fn),
        allowed_seconds: allowed,
        energy: e,
        change_energy_res:boolean_res,
        value: vec![value],
        no_save: opt.no_save
    };

    execute_high_degree_helper_ba(rewl, instant, opts)

}

pub fn execute_sw(opt: LDOptsLD, instant: std::time::Instant,r:f64,param:LDLDparam,value:serde_json::Value){
    if let Some(num) = opt.num_threads
    {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num.get())
            .build_global()
            .unwrap();
    }

    let boolean_res:bool = (param.change_energy_res).clone();

    let base_options = SWOptions{
        graph_seed: param.graph_seed,
        lambda: param.lambda,
        gamma: param.recovery_prob,
        system_size: param.system_size,
        rewire_prob:r}; 

    let sw:SWModel = base_options.into();
    let model = SWLargeDeviation::new(sw.clone(),param.large_deviation_param,param.lockdownparams);
    let mut ld_model = SWLargeDeviationWithLocks::new(model);
    ld_model.infect_initial_patients();
    

    let histograms =  if !boolean_res{
            param.histograms.create(
                ld_model.patient_zero_vec.len() as u32, 
                param.system_size.get() as u32, 
                param.walkers_per_interval
            )}
        else{
            println!("yes");
            param.histograms.create(
                (ld_model.patient_zero_vec.len() as u32+1)/2, 
                (param.system_size.get() as u32+1)/2, 
                param.walkers_per_interval
                )
    };
    print_min_max_intervalsize(&histograms);
    print_intervals_hists(&histograms);

    let mut markov_seeder = Pcg64::seed_from_u64(param.large_deviation_param.markov_seed);
    let uniform = Uniform::new(0, 90357094723984_u64);

    let ensembles = histograms
        .iter()
        .map(
            |_|
            {
                let mut l_param = param.large_deviation_param;
                l_param.markov_seed = uniform.sample(&mut markov_seeder);
                let ld = SWLargeDeviation::new(sw.clone(), l_param,param.lockdownparams);
                SWLargeDeviationWithLocks::new(ld)
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

    let name = "curve";

    let mut writer = SirWriter::new(&name, 1);
    writer.write_header(&[value.clone()]).unwrap();

    let energy = energy_function_returner_sw(param.energy,param.change_energy_res);

    let rewl = rewl_builder.greedy_build(energy);
    let duration = format_duration(instant.elapsed());
    println!("Finished greedy build after {duration}");
    //panic!();
    
    rewl.ensemble_iter().for_each(
        |ensemble|
            {
                ensemble.tracker.write_stats(std::io::stdout());
            }
        );
    //println!("{}",rewl.walkers().len());
    let allowed = param.allowed_seconds();
    let e = param.energy;
    let name_fn = move |interval: Option<usize>, end: &str, mode: LargeDeviationMode|
    {
        param.quick_name(interval, end, mode)
    };

    let opts = LDLdOpts{
        quick_name: Box::new(name_fn),
        allowed_seconds: allowed,
        energy: e,
        change_energy_res:boolean_res,
        value: vec![value],
        no_save: opt.no_save
    };

    execute_high_degree_helper_sw(rewl, instant, opts)}


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

fn print_intervals_hists(hists: &[HistogramFast<u32>]){
    hists.iter().for_each(|his|{
        let range = his.range_inclusive();
        println!("beg {}, end {}",range.start(),range.end() );
    }
    )
}

pub fn execute_high_degree_helper_sw(mut rewl: Rewl<SWLargeDeviationWithLocks, Pcg64, HistogramFast<u32>, u32, MarkovStepWithLocks, ()>, 
    instant: std::time::Instant,
    opts: LDLdOpts){

    let energy = energy_function_returner_sw(opts.energy,opts.change_energy_res);

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
                rewl.ensemble_iter()
                    .enumerate()
                    .for_each(
                        |(index, ensemble)|
                        {
                            let _ = writeln!(buf, "#Ensemble {index}");
                            ensemble.tracker.write_stats(&mut buf);
                        }
                    );
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
pub fn execute_high_degree_helper_ba(mut rewl: Rewl<BALargeDeviationWithLocks, Pcg64, HistogramFast<u32>, u32, MarkovStepWithLocks, ()>, 
    instant: std::time::Instant,
    opts: LDLdOpts){

    let energy = energy_function_returner_ba(opts.energy,opts.change_energy_res);

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
                rewl.ensemble_iter()
                    .enumerate()
                    .for_each(
                        |(index, ensemble)|
                        {
                            let _ = writeln!(buf, "#Ensemble {index}");
                            ensemble.tracker.write_stats(&mut buf);
                        }
                    );

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

pub type BALDRewl = Rewl<BALargeDeviationWithLocks, Pcg64, HistogramFast<u32>, u32, MarkovStepWithLocks, ()>;
pub type SWLDRewl = Rewl<SWLargeDeviationWithLocks, Pcg64, HistogramFast<u32>, u32, MarkovStepWithLocks, ()>;






pub fn load_high_degree_rewl(opts: LDContinueOpts, instant: std::time::Instant, no_save: bool){

    match opts.rewltype{
        RewlType::SmallWorld => load_small_world_rewl(opts, instant, no_save),
        RewlType::Barabasi => load_barabasi_rewl(opts,instant,no_save),
        _ => println!("Enter the graph type under the RewlType field. SmallWorld or Barabasi")
    }



    
    

}

pub fn load_small_world_rewl(opts: LDContinueOpts, instant: std::time::Instant, no_save: bool){
    let allowed_seconds = opts.seconds();

    let (mut rewl, mut json_string):(Rewl<SWLargeDeviationWithLocks, _, _, _, _, ()>, Vec<_>) = deserialize_from_file(&opts.file_name);




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
    

    let mut old_opts: LDLDparam = serde_json::from_str(&json_string[0])
        .expect("Unable to deserialize old params");
    old_opts.times_plus_1();
    let boolean_res = old_opts.change_energy_res;

    let mut jsons = vec![serde_json::to_value(old_opts.clone()).expect("unable to create json")];

    jsons.extend(
        json_string.iter()
            .skip(1)
            .filter_map(
                |s|
                {
                    let res: Result<LDContinueOpts, _> = serde_json::from_str(s);
                    
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
    let opts = LDLdOpts{
        quick_name: Box::new(name_fn),
        allowed_seconds,
        energy: e,
        change_energy_res: boolean_res,
        value: jsons,
        no_save
    };
    
    execute_high_degree_helper_sw(rewl, instant, opts)
        

    

}
pub fn load_barabasi_rewl(opts: LDContinueOpts, instant: std::time::Instant, no_save: bool){
    let allowed_seconds = opts.seconds();

    let (mut rewl, mut json_string):(Rewl<BALargeDeviationWithLocks, _, _, _, _, ()>, Vec<_>) = deserialize_from_file(&opts.file_name);

    


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
    

    let mut old_opts: LDLDparam = serde_json::from_str(&json_string[0])
        .expect("Unable to deserialize old params");
    old_opts.times_plus_1();

    let mut jsons = vec![serde_json::to_value(old_opts.clone()).expect("unable to create json")];

    jsons.extend(
        json_string.iter()
            .skip(1)
            .filter_map(
                |s|
                {
                    let res: Result<LDContinueOpts, _> = serde_json::from_str(s);
                    
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
    let boolean_res:bool = old_opts.change_energy_res;

    let e = old_opts.energy;
    let name_fn = move |interval: Option<usize>, end: &str, mode: LargeDeviationMode|
    {
        old_opts.quick_name(interval, end, mode)
    };

    let opts = LDLdOpts{
        quick_name: Box::new(name_fn),
        allowed_seconds,
        energy: e,
        change_energy_res: boolean_res,
        value: jsons,
        no_save
    };
    
    execute_high_degree_helper_ba(rewl, instant, opts)
        

    

}


