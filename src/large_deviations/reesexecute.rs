use core::panic;

use{
    std::{
        fs::File,
        io::{BufWriter, Write},
        time::Instant
    },
    crate::sir_model::*,
    serde_json::Value,
    super::*,

};

pub fn execute_entropic_sampling(opt: ReesOpts, start_time: Instant)
{

    let mut log_file = LogfilePrinter::new("rees.log");
    let allowed = opt.allowed_seconds();
    writeln!(log_file, "Oppening {}", &opt.save_file).unwrap();

    // read in save state
    let (rewl, old_json_string) = deserialize_from_file(&opt.save_file);

    let old_jsons: Vec<Value> = 
    old_json_string.iter()
        .filter_map(
            |st|
            {
                match serde_json::from_str(st)
                 {
                    Ok(o) => Some(o),
                    Err(e) => {
                        eprintln!("Cannot read old json - skipping due to error {e}");
                        None
                    }
                 }
            }
        ).collect();

    let old_param: BALDLDparam = serde_json::from_value(old_jsons[0].clone())
        .expect("Unable to parse old json");

    let energy = energy_function_returner(old_param.energy);

    let writer: Vec<_> = rewl.hists()
        .iter()
        .enumerate()
        .zip(rewl.walkers())
        .map(
            |((index, hist), walker)|
            {
                let left = hist.left();
                let right = hist.right();
                let name = old_param.quick_name(None, "", LargeDeviationMode::Rees);
                let name = format!("{name}_{left}_{right}");

                let ever = walker.step_count() / opt.print_count;

                if ever == 0 {
                    panic!("Ever == 0! Abbort");
                }

                let mut writer = SirWriter::new(&name, index);
                writer.write_header(&old_jsons).unwrap();

                (ever, writer)
            }
        ).collect();

    // convert to rees
    let mut rees = match rewl.into_rees_with_extra(writer)
    {
        Ok(rewl) => rewl,
        Err(_) => panic!("Unable to create rees with extra!")
    };

    rees
        .simulate_while(energy, 
            |_| start_time.elapsed().as_secs() < allowed,
            |rees_walker, ensemble, extra| {
                if rees_walker.step_count() % extra.0 == 0{
                    let extinction_index = ensemble.last_extinction_index();
                    let energy = rees_walker.energy_copy();
                    extra.1.write_energy(energy, extinction_index).unwrap();
                    ensemble.ld_energy_m_and_print(&mut extra.1)
                }
            }
        );

    let unfinished_count = unsafe{
        rees.ensemble_iter_mut()
            .fold(
                0, 
                |acc, m| m.unfinished_counter() + acc
            )};

    let total_simulations_count = unsafe{
        rees.ensemble_iter_mut()
            .fold(
                0, 
                |acc, m| m.total_simulations_counter() + acc
            )};

    let unfinished_frac = unfinished_count as f64 / total_simulations_count as f64;

    rees.walkers()
        .iter()
        .enumerate()
        .for_each(
            |(index, walker)|
            {
                writeln!(log_file, "walker {index}, steps: {}", walker.step_count())
                    .unwrap();
                let density = walker.log10_density();
                let name = old_param.quick_name(Some(index), "dat", LargeDeviationMode::Rees);
                writeln!(log_file, "Created: {}", &name).unwrap();
                let file = File::create(name)
                    .expect("unable to create file");
                let mut buf = BufWriter::new(file);
                write!(buf, "#old json: ").unwrap();
                write_jsons(&old_jsons, &mut buf).unwrap();
                writeln!(buf, "#hist log10").unwrap();
                writeln!(buf, "#walker {index}, steps: {}", walker.step_count())
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

                write!(buf, "#rewl roundtrips of all walkers:").unwrap();

                for roundtrip in rees.rewl_roundtrip_iter()
                {
                    write!(buf, " {roundtrip}").unwrap();
                }
                writeln!(buf).unwrap();
                write!(buf, "#rees roundtrips of all walkers:").unwrap();

                for roundtrip in rees.rees_roundtrip_iter()
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
            }
        );

    writeln!(log_file, "OLD rewl roundtrips")
        .unwrap();
    for roundtrip in rees.rewl_roundtrip_iter()
    {
        write!(log_file, " {roundtrip}").unwrap();
    }
    write!(log_file, "\nNew rees roundtrips").unwrap();
    for roundtrip in rees.rees_roundtrip_iter()
    {
        write!(log_file, " {roundtrip}").unwrap();
    }
    writeln!(log_file).unwrap();

    if !opt.no_save{
        let (rees, _) = rees.unpack_extra();
        let save_name = old_param.quick_name(None, "bincode", LargeDeviationMode::Rees);
        writeln!(log_file, "Create save: {save_name}").unwrap();
        let save = (rees, old_json_string);

        let file = File::create(save_name)
            .expect("Unable to create save file");
        let buf = BufWriter::new(file);

        bincode::serialize_into(buf, &save)
            .unwrap();
    }

}