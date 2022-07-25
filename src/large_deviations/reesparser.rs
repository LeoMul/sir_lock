use core::panic;
use std::num::NonZeroU64;

use{
    structopt::StructOpt,
    std::{
        fs::File,
        io::{BufWriter, Write, BufReader},
    },
    //crate::sir_model::sir_states::Deserialize,
    super::*,
    std::{num::*},
    serde::{Serialize, Deserialize},
};
use net_ensembles::{sampling::{HistogramFast}};
use crate::misc_types::*;

use {
    rand_pcg::Pcg64,
    net_ensembles::sampling::Rewl,
};

#[derive(Clone, Serialize, Deserialize)]
pub struct REESJsonOpts
{
    pub file_name: String,
    pub seconds: Option<NonZeroU64>,
    pub minutes: Option<NonZeroU64>,
    pub days: Option<NonZeroU64>,
    pub hours: Option<NonZeroU64>,
    pub change_step_size: Option<NonZeroUsize>,
    pub rewltype: RewlType,
    pub no_save: bool,
    pub print_count:Option<NonZeroU64>
}
impl REESJsonOpts{
    pub fn allowed_seconds(&self) -> u64
    {
        let mut secs = 0;
        if let Some(sec) = self.seconds{
            secs += sec.get();
        }
        if let Some(min) = self.minutes
        {
            secs += min.get() * 60;
        }

        if let Some(hours) = self.hours
        {
            secs += hours.get() * 60 * 60;
        }

        if let Some(days) = self.days
        {
            secs += days.get() * 60 * 60 * 24;
        }

        if secs == 0 {
            panic!("Specified time is 0. Unable to process anything without time")
        }

        secs
    }
}
impl Default for REESJsonOpts
{
    fn default() -> Self
    {
        Self{
            file_name: "insert file with binary save".to_owned(),
            seconds: None,
            minutes: None,
            hours: None,
            days: None,
            change_step_size: None,
            rewltype: RewlType::None,
            no_save: false,
            print_count: NonZeroU64::new(2000)
        }
    }
}


#[derive(StructOpt, Debug, Clone)]
/// Replica exchange entropic sampling for lockdowns
pub struct ReesOpts
{   
    #[structopt(long)]
    pub json: Option<String>,
}

impl ReesOpts
{
    
    pub fn execute_old(self, _instant: std::time::Instant)
    {
        //execute_entropic_sampling(self,instant)
    }
    pub fn execute(&self, instant: std::time::Instant)
    {
    
        let json_opts: REESJsonOpts = match self.json.as_ref()
        {
            None => {
                let def = REESJsonOpts::default();
                serde_json::to_writer_pretty(
                    std::io::stdout(), 
                    &def
                ).unwrap();
                std::process::exit(0)
            },
            Some(filename) => {
                let file = File::open(filename).expect("Unable to open json file");
                let reader = BufReader::new(file);
                serde_json::from_reader(reader).expect("Parsing error for json")
            }
        };
        //let rewlt = json_opts.rewltype;
        execute_entropic_sampling(&json_opts,instant,json_opts.rewltype.clone())
        
    

    }


}

pub struct LogfilePrinter
{
    printer: BufWriter<File>
}

impl LogfilePrinter
{
    pub fn new(name: &str) -> Self
    {
        let file = File::create(name)
            .unwrap();
        let buf = BufWriter::new(file);
        
        Self{
            printer: buf
        }
    }
}

impl Write for LogfilePrinter
{
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize>
    {
        std::io::stdout().write_all(buf)?;
        self.printer.write(buf)
    }

    fn flush(&mut self) -> std::io::Result<()>
    {
        self.printer.flush()
    }
}







//use {serde::{Serialize, Deserialize}};

type RewlCompact<T> = (Rewl<T,Pcg64,HistogramFast<u32>,u32,MarkovStepWithLocks,()>, Vec<String>);

pub fn deserialize_from_file<T>(filename: &str) -> RewlCompact<T>
where T: serde::de::DeserializeOwned
{
    let file = File::open(filename).expect("unable to open save file");
    let reader = BufReader::new(file);

    let res: Result<RewlCompact<T>, _> =  bincode::deserialize_from(reader);

        match res
        {
            Ok((rewl, json_string)) => {
                (
                    rewl,
                    json_string
                )
            },
            _ => 
            {
                let file = File::open(filename).expect("unable to open save file");
                let reader = BufReader::new(file);
                bincode::deserialize_from(reader).expect("unable to parse bincode file")
            }
        }
}