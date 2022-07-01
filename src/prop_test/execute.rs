
use{
    super::*,
    crate::{GraphType, sir_model::*},
    serde_json::Value,
    std::{num::*},
    crate::lockdown_methods::*,

};


pub fn run_simulation(param:PropTestParams, json: Value, num_threads: Option<NonZeroUsize>){
    match param.graph_type{
        GraphType::SmallWorld(_) => sim_small_world(param, json, num_threads),
        GraphType::Barabasi(_,_) => sim_barabasi(param, json, num_threads),
        _ => unimplemented!()
    }

}

fn sim_small_world(_param:PropTestParams, _json: Value, _num_threads: Option<NonZeroUsize>){}

fn sim_barabasi(param:PropTestParams, _json: Value, _num_threads: Option<NonZeroUsize>){
    let opt = BarabasiOptions::from_prop_test_param(&param);
    let barabasi_world = opt.into();
    let mut model = SimpleSampleBarabasi::from_base(barabasi_world, param.sir_seed);

    let lock = LockdownParameters{
        lock_style: LockdownType::LimitContacts(2),
        dynamic_bool: true,
        lock_threshold: 0.1,
        release_threshold: 0.05
    };
    let lock_graph = model.create_locked_down_network(lock);
    model.propagate_until_completion_max_with_lockdown(lock_graph,lock);

}