use std::fmt::Display;

use{
    super::*,
    structopt::StructOpt,
    std::num::*,
    crate::json_parsing::*,
    serde::{Serialize, Deserialize},
    serde_json::Value,
   
    crate::misc_types::*,
};


#[derive(Debug, StructOpt, Clone)]
/// Largest connected component measurement.
pub struct ConnectedComponent{
    #[structopt(long)]
    json: Option<String>,

    #[structopt(long)]
    num_threads:Option<NonZeroUsize>
}

impl ConnectedComponent{
    pub fn parse(&self) -> (ConnectedComponentParams, Value){
        parse(self.json.as_ref())
    }
    pub fn execute(&self){
        let (opt, json) = self.parse();
        run_simulation(opt,json,self.num_threads)
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ConnectedComponentParams{

    pub system_size: NonZeroUsize,
    pub graph_type: GraphType,
    pub graph_seed: u64,
    pub num_networks: u64,
    pub chunk_len: u64,
    pub percent_start:f64,
    pub percent_end:f64,

}
impl Default for ConnectedComponentParams{
    fn default() -> Self{
    Self{
        system_size: NonZeroUsize::new(2048).unwrap(),
        graph_type: GraphType::SmallWorld(0.1),
        graph_seed:DEFAULT_GRAPH_SEED,
        num_networks: 100000,
        chunk_len: 10,
        percent_start:0.4,
        percent_end:0.6,
    }}
}
impl ConnectedComponentParams{
    pub fn name<E>(&self, file_ending:E , num_threads:Option<NonZeroUsize>) -> String where E:Display{
        let k = match num_threads{
            None => "".to_owned(),
            Some(v) => format!("k{}",v)
        };
        
        format!(
            "ver{}ConnectedComponentSize{}NumNet{}chunk_len{}_GT{}_GS{}_THR{}.{}",
            crate::VERSION,
            self.system_size.get(),
            self.num_networks,
            self.chunk_len,
            self.graph_type.name(),
            self.graph_seed,
            k,
            file_ending


        )
    }
}