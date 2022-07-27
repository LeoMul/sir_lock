
use{
    super::*,
    crate::{GraphType, sir_model::*},
    serde_json::Value,
    std::{num::*},
    crate::lockdown_methods::*,
    crate::misc_types::*,

};


pub fn run_simulation(param:PropTestParams, json: Value, num_threads: Option<NonZeroUsize>){
    match param.graph_type{
        GraphType::SmallWorld(_) => sim_small_world(param, json, num_threads),
        GraphType::Barabasi(_,_) => sim_barabasi(param, json, num_threads),
        _ => unimplemented!()
    }

}

fn sim_barabasi(_param:PropTestParams, _json: Value, _num_threads: Option<NonZeroUsize>){}

fn sim_small_world(param:PropTestParams, _json: Value, _num_threads: Option<NonZeroUsize>){
    let opt = SWOptions::from_prop_test_param(&param);
    let barabasi_world = opt.into();
    let mut model = SimpleSampleSW::from_base(barabasi_world, param.sir_seed,param.initial_infected);

    let lock = LockdownParameters{
        lock_style: LockdownType::Random(DEFAULT_RANDOM_LOCKDOWN_FRAC),
        lock_threshold: 0.1,
        release_threshold: 0.05
    };
    let mut lock_graph = model.create_locked_down_network(lock);
    let x = model.propagate_until_completion_max_with_lockdown(&mut lock_graph,lock);
    println!("{}",x);

}