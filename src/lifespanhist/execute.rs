use std::{ io::Write};

use crate::life_span_data_collection::*;

use{
    super::*,
    crate::{GraphType, sir_model::*},
    serde_json::Value,
    std::{num::*, fs::File, io::BufWriter},
};



pub fn run_simulation(param:LifeSpanParams, json: Value, num_threads: Option<NonZeroUsize>){
    match param.graph_type{
        GraphType::SmallWorld(_) => sim_small_world(param, json, num_threads),
        _ => unimplemented!()
    }

}

fn sim_small_world(param: LifeSpanParams, json: Value, num_threads:Option<NonZeroUsize>){
    let opt = BaseSwOptions::from_lifespan_param(&param);
    let small_world = opt.into();
    let model = SimpleSample::from_base(small_world, param.sir_seed);
    let k = num_threads.unwrap_or_else(|| NonZeroUsize::new(1).unwrap());

    rayon::ThreadPoolBuilder::new().num_threads(k.get()).build_global().unwrap();
    
    let data = acquire_sorted_data(model, k, param.sir_seed, param.samples);
    let vector_data = convert_sorted_to_hist(&data);
    
    

   
    
    write_cumulative(&data, &param, &json, num_threads);
    write_histogram(&vector_data,&param,&json,num_threads);
    let life95 = acquire_percent_life_span(&data,0.95);
    println!{"{}",life95};

}



fn write_histogram(data:&Vec<(u32,usize)>,param: &LifeSpanParams, json: &Value, num_threads:Option<NonZeroUsize>){
    //takes properly binned data vector and outs a data file for histogram plotting in gnuplot using with boxes
    let name = param.name("hist","dat", num_threads);
    println!("creating: {name}");
    let norm_fraction = if param.fraction{
        param.samples as f64
    }else{
        1.
    };
    let file = File::create(name).expect("unable to create file");
    let mut buf = BufWriter::new(file);
    write!(buf, "#").unwrap();
    serde_json::to_writer(&mut buf, &json)
        .unwrap();
    writeln!(buf).unwrap();

    writeln!(buf, "#frac time_steps").unwrap();
    for j in data{
        writeln!(buf,"{} {}",j.0,j.1 as f64/norm_fraction).unwrap();
    }
    
}

fn write_cumulative(data:&[u32],param: &LifeSpanParams, json: &Value, num_threads:Option<NonZeroUsize>){
    let name = param.name("cum","dat", num_threads);
    println!("creating: {name}");

    let file = File::create(name)
        .expect("unable to create file");
    let mut buf = BufWriter::new(file);

    write!(buf, "#").unwrap();
    serde_json::to_writer(&mut buf, &json)
        .unwrap();
    writeln!(buf).unwrap();

    writeln!(buf, "#frac time_steps").unwrap();
    let total = param.samples as f64;
    let mut last_val = 0;
    let mut last_index = 0;
    for (&val, i) in data.iter().zip(1_usize..)
    {   
        // skip if value does not change
        if val == last_val{
            continue;
        }

        // value changed. If we did not write the 
        // datapoint right before yet, then write it now
        if i-1 != last_index{
            let percent = ((i-1) as f64) / total;
            writeln!(buf, "{:e} {last_val}", percent).unwrap()  
        }
        last_val = val;
        last_index = i;
        // write new value
        let percent = (i as f64) / total;
        writeln!(buf, "{:e} {val}", percent).unwrap()
    }

}

