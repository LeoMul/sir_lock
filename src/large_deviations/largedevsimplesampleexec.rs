use{
    net_ensembles::sampling::{HistU32Fast, MarkovChain},
    super::*,
    serde_json::Value,

    crate::{json_parsing::*, GraphType},
    crate::sir_model::*,
    rand_pcg::Pcg64,
    net_ensembles::rand::{SeedableRng},
    net_ensembles::sampling::{Histogram},

    //crate::simple_sample::execute::hist_to_file
};
use std::{fs::File, io::BufWriter};
use std::{ io::Write};


pub fn typical_event_sampling(opts: SimpleSampleCMDopts)
{
    let (param, value): (SimpleSampleldParam, _) = parse(opts.json.as_ref());


    let (m,n) = match &param.graph_type {
        GraphType::Barabasi(p,m) => (*p,*m),
       _ =>unimplemented!()
    };

    let base_options = BarabasiOptions{
        graph_seed: param.graph_seed,
        lambda: param.lambda,
        gamma: param.recovery_prob,
        system_size: param.system_size,
        m,
        source_n:n
    };

    let base_sw: BarabasiModel = base_options.into();

    let model = BALargeDeviation::new(base_sw, param.large_deviation_param, param.lockdown);

    let mut ld_model = BALargeDeviationWithLocks::new(model);

    let system_size = param.system_size.get() as u32; 
    let mut hist = HistU32Fast::new_inclusive(1, system_size)
        .unwrap();

    let mut rng = Pcg64::seed_from_u64(param.large_deviation_param.markov_seed);

    //let mut vaccine_rng = Pcg64::seed_from_u64(4896709264107025);
    let mut markov_vec = Vec::new();

    let bar = crate::indication_bar(param.samples as u64);
    for i in 0..param.samples
    {
        match param.randomize{
            Randomize::Random => ld_model.randomise_monte_carlo(&mut rng),
            Randomize::Markov(markov) => {
                for _ in 0..markov.every.get()
                {
                    ld_model.m_steps(markov.step_size.get(), &mut markov_vec);
                }
            }
        }
        

        let mut energy = ld_model.ld_energy_m();

        if matches!(param.energy, crate::misc_types::MeasureType::C)
        {
            energy = ld_model.calculate_ever_infected() as u32;
        }

        hist.increment(energy)
            .unwrap();

        if i % 1000 == 0 
        {
            bar.inc(1000);
        }
    }
    bar.finish_and_clear();

    let name = param.quick_name();

    hist_to_file(&hist, name, &value);

    println!("unfinished: {}", ld_model.unfinished_counter())
}

pub fn hist_to_file(hist: &HistU32Fast, file_name: String, json: &Value)
{
    let normed = norm_hist(hist);

    println!("Creating {}", &file_name);
    let file = File::create(file_name)
        .unwrap();
    let mut buf = BufWriter::new(file);
    write!(buf, "#").unwrap();
    serde_json::to_writer(&mut buf, json)
        .unwrap();
    writeln!(buf).unwrap();
    writeln!(buf, "#bin log10_prob hits").unwrap();

    hist.bin_hits_iter()
        .zip(normed)
        .for_each(
            |((bin, hits), log_prob)|
            {
                writeln!(buf, "{} {} {}", bin, log_prob, hits).unwrap()
            }
        );
}

pub fn norm_hist(hist: &HistU32Fast) -> Vec<f64>
{
    let mut density: Vec<_> = hist.hist()
        .iter()
        .map(|&hits| (hits as f64).log10())
        .collect();

    subtract_max(density.as_mut_slice());
    let int = integrate_log(density.as_slice(), hist.hist().len());
    let sub = int.log10();
    density.iter_mut()
        .for_each(|v| *v -= sub);
    density
}
pub fn subtract_max(slice: &mut[f64])
{
    let mut max = std::f64::NEG_INFINITY;
    slice.iter()
        .for_each(
            |v|
            {
                if *v > max {
                    max = *v;
                }
            }
        );
    slice.iter_mut()
        .for_each(|v| *v -= max);
}

pub fn integrate_log(curve: &[f64], n: usize) -> f64
{
    let delta = 1.0 / n as f64;
    curve.iter()
        .map(|&val| delta * 10_f64.powf(val))
        .sum()
}